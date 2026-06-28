# Peaky — Chemistry core (masses, DBE/Senior gates, the formula grid)

This document explains **how a region of m/z space becomes a finite list of
plausible neutral formulas** — the exact-mass arithmetic, the structural gates
(integer DBE, Senior's rule, the oxygen cap), and the enumerated grid every
candidate flows from. It is a module deep-dive companion to
[`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole pipeline, §5 *Chemistry gates*),
[`ISOTOPES.md`](ISOTOPES.md) (which *shrinks* this grid from spectral evidence),
and [`ASSIGNMENT.md`](ASSIGNMENT.md) (which scores + arbitrates the candidates).

**Code:** `peaky/chem/chemistry.py` — the lowest layer. It knows nothing about
Mascope, contexts, or peaks; it is pure arithmetic + combinatorics.

> Keep this in sync with the code. Every threshold below is a named constant or a
> literal in `chemistry.py`; if you change one there, change it here.

---

## 1. What this stage does

This module turns *element ranges* (a box like `C0-40 H0-80 O0-18 …`) into the
set of **neutral** formulas that (a) are valence-legal and (b) match an observed
peak m/z under some adduct, within ppm tolerance. It is **structural, not
statistical** — the gates are valence facts applied identically every run, which
is why a "mass fit" alone never wins a peak. It also supplies the per-element
**complexity penalty** arbitration uses to break score ties toward the simpler
interpretation.

```
element ranges (C0-40 H0-80 O0-18 N0-3 …)   +   mass window [30, 900] Da
        │  enumerate_grid:  derive H from integer DBE, not iterate it
        ▼
   for each (C,Si,N,P,F,Cl,Br,I):  DBE 0..Senior-cap → H ;  then S, O (O-cap)
        │   keep iff DBE≥0 ∧ integer ∧ ≤ Senior cap ∧ O ≤ 2(C+N+S+P)+4
        ▼
   sorted [(neutral_mass, formula)]  (memoised grid cache)
        │  candidates_for_peaks: m_neutral = peak_mz − ADDUCT_SHIFTS[adduct]
        ▼   bisect the grid masses within ±ppm of each peak
   set of neutral formulas matching at least one peak under at least one adduct
```

---

## 2. Inputs

- **Element ranges** — `parse_ranges("C0-30 H0-60 O0-15 N0-5 Br0-2")` →
  `{'C':(0,30), …}`. The box usually comes from a reagent profile
  ([`REAGENTS.md`](REAGENTS.md)) after [`ISOTOPES.md`](ISOTOPES.md) tightens it.
- **Adducts** — labels keyed in `ADDUCT_SHIFTS` (`ion_mz = neutral_mass + shift`).
- **Peak m/z** + a **ppm tolerance** for `candidates_for_peaks`.

---

## 3. The transformation, stage by stage

1. **Parse / format** (`parse_formula`, `format_formula`). Strings ↔ count dicts.
   Output is **Hill notation** (C, then H, then alphabetical) to match exactly how
   Mascope formats formula strings — string-level comparison against scorer output
   depends on it.

2. **Exact mass** (`neutral_mass`). Σ of monoisotopic element masses from `M`
   (C = 12.0 by definition; H = 1.0078250319; etc.). The **electron mass**
   `M_E = 0.0005485799` enters only at the ion layer.

3. **Ion m/z** (`ion_mz`). `neutral_mass + ADDUCT_SHIFTS[adduct]`. Each shift is
   built from the same exact masses ± `M_E` (an anion *gains* an electron, a cation
   *loses* one). Unknown adduct → `KeyError`.

4. **DBE on the NEUTRAL** (`dbe`). `DBE = 1 + (C+Si) + (N+P)/2 − (H+F+Cl+Br+I)/2`.
   O and S are divalent (zero contribution); Si is tetravalent (carbon-equivalent),
   P trivalent (like N), **halogens count as hydrogens**. A real neutral has a
   **non-negative integer** DBE; a half-integer DBE is an *ion-only* artifact of
   (de)protonation, so it is rejected here and handled at the ion layer (organic
   nitrates pass as neutrals).

5. **Senior cap** (`seniors_cap`). Max DBE = `C + Si + N/2 + 1`.

6. **Oxygen cap** (`oxygen_ok`). `O ≤ 2·(C+N+S+P) + 4`. Pure valence (every O
   needs a bonding site on the skeleton; +4 headroom for peroxide chains /
   terminal acids), **not** a Van Krevelen prior. DBE can't constrain O (divalent
   → zero), so without this the O-rich corner hits *any* mass within tolerance
   (`C3H5ClO17`); real HOMs (`C10H18O7`) and inorganic acids (`H2SO4`, `HNO3`) pass.

7. **The hard gate** (`dbe_ok`, tol `1e-9`). A formula survives iff DBE ≥ 0,
   DBE is integer, **and** DBE ≤ Senior cap. Returns `(ok, reason_if_not)`.

8. **Enumerate the grid** (`enumerate_grid`, mass `[30, 900]` Da). The clever part:
   **H is derived, not iterated.** For each `(C,Si,N,P,F,Cl,Br,I)` it walks integer
   DBE `0..⌊C+Si+N/2+1⌋` and solves `H = 2(1+C+Si) + (N+P) − 2·DBE − halo`, so only
   valid integer-DBE neutrals are ever visited (correct *and* fast). The O loop is
   capped inline at `min(O_max, 2(C+N+S+P)+4)`. `_grid_cached` memoises the sorted
   grid (LRU-ish: cleared at **`_GRID_CACHE_MAX` 16** entries).

9. **Complexity penalty** (`complexity_penalty`). Per-atom weights
   `{N:3, S:8, P:25, Cl:50, Br:50, Si:80, I:80, F:30}`, scaled by `0.01`, capped at
   `0.20`. Arbitration subtracts this so a heteroatom-rich neutral must out-score a
   CHO competitor by a real margin (a neutral halogen never wins on a mass tie).

10. **Pre-filter to peaks** (`candidates_for_peaks`). For each peak m/z and adduct,
    `m_neutral = mz − shift`; `bisect` the sorted grid masses within
    `±(m_neutral·ppm·1e-6)`; union the matching formulas. The result is the
    candidate universe handed to scoring.

---

## 4. Constants reference

All in `peaky/chem/chemistry.py`.

| constant | value | role |
| --- | --- | --- |
| `M["C"]` … `M["I"]` | exact monoisotopic masses | C=12.0, H=1.0078250319, O=15.9949146221, N=14.0030740052, … |
| `M_E` | 0.0005485799 Da | electron mass (anion +e, cation −e in shifts) |
| `_M_15N` | 15.0001088989 | ¹⁵N mass for the labelled-nitrate adduct shift |
| `_DBE_PLUS` | `(C, Si)` | +1 each in DBE (valence 4) |
| `_DBE_HALF_PLUS` | `(N, P)` | +½ each (valence 3) |
| `_DBE_HALF_MINUS` | `(H, F, Cl, Br, I)` | −½ each (valence 1; halogens = H) |
| Senior cap | `C + Si + N/2 + 1` | max permitted DBE |
| oxygen cap | `2·(C+N+S+P) + 4` | structural O ceiling |
| `dbe_ok` `tol` | 1e-9 | integer/Senior comparison tolerance |
| `enumerate_grid` `mass_min` / `mass_max` | 30.0 / 900.0 Da | grid mass window |
| `_COMPLEXITY_WEIGHT` | N3 S8 P25 Cl50 Br50 Si80 I80 F30 | per-atom heteroatom cost |
| `complexity_penalty` `scale` / `cap` | 0.01 / 0.20 | penalty = min(Σweight·scale, cap) |
| `_GRID_CACHE_MAX` | 16 | memoised-grid entries before the cache is cleared |
| `ADDUCT_SHIFTS` | per-adduct Δm | `ion_mz = neutral_mass + shift` (see §5) |

---

## 5. Metrics, defined

- **DBE** — `1 + (C+Si) + (N+P)/2 − (H+F+Cl+Br+I)/2`, on the *neutral*; must be a
  non-negative integer.
- **Senior cap** — `C + Si + N/2 + 1`; the largest DBE the atom count can support.
- **O cap** — `2·(C+N+S+P) + 4`; the largest oxygen count valence allows.
- **complexity penalty** — `min(Σ weightₑ·nₑ · 0.01, 0.20)` ∈ [0, 0.20]; the
  score handicap that biases ties to CHO < CHON < heteroatom-rich.
- **`ADDUCT_SHIFTS`** — e.g. `[M-H]⁻ = −H + e`, `[M+Br]⁻ = Br + e`,
  `[M+NO3]⁻ = N+3O + e`, `[M+^NO3]⁻ = ¹⁵N+3O + e` (+62.9855, not +61.9885),
  `[M+H]⁺ = H − e`, `[M+NH4]⁺ = N+4H − e`, `[M+(CH4N2O)H]⁺ = C+5H+2N+O − e`
  (+61.0396, the uronium analyte channel).

---

## 6. Outputs

| artifact | content |
| --- | --- |
| `enumerate_grid` / `_grid_cached` | sorted `[(neutral_mass, formula)]` for a box — the valence-legal universe |
| `candidates_for_peaks` | `set[str]` of neutral formulas matching ≥1 peak under ≥1 adduct within ppm |
| `neutral_mass` / `ion_mz` | exact mass of a neutral / theoretical m/z of its singly-charged ion |
| `dbe` / `dbe_ok` / `oxygen_ok` | DBE value; `(ok, reason)` gate verdicts |
| `complexity_penalty` | a 0–0.20 score handicap for arbitration |

---

## 7. Properties, invariants & gotchas

- **DBE is enforced on the neutral, as an integer.** Half-integer DBE is *expected*
  on the observed ion (`C8H14NO12⁻`) and handled at the ion layer — never here.
- **The O cap is valence, not a VK prior.** It is the one gate that constrains the
  oxygen axis (DBE can't), and it is what kills O-rich mass-coincidences while
  letting genuine HOMs through.
- **Halogens count as hydrogens** in DBE; Si is a carbon-equivalent backbone atom
  (so the VK H/C ratio downstream uses `(H+F+Cl+Br+I)/(C+Si)`).
- **Hill notation is load-bearing.** `format_formula` must match Mascope's string
  formatting or string-level joins against scorer output silently miss.
- **H is derived from DBE**, not iterated — both a correctness guarantee (only
  valid DBE visited) and the reason the grid is fast.
- **The complexity penalty is capped at 0.20** and is a *tie-breaker*, not a veto:
  strong isotope/series evidence can still carry a heteroatom formula.
- **Electron mass is in the shift, not the neutral.** `neutral_mass` is electron-
  neutral; `ion_mz` adds/removes `M_E` per the adduct's charge.
- **Grid cache is unbounded-then-flushed**: at 16 distinct boxes it clears wholesale
  (not true LRU) — fine for a run that reuses a handful of boxes.

---

## 8. Code map

| function | role |
| --- | --- |
| `parse_formula` / `format_formula` | string ↔ count dict (Hill notation) |
| `neutral_mass` / `ion_mz` | exact neutral mass; theoretical ion m/z via `ADDUCT_SHIFTS` |
| `dbe` | Double-Bond-Equivalents of the neutral |
| `seniors_cap` / `dbe_ok` | Senior max DBE; the integer-DBE + Senior hard gate |
| `oxygen_ok` | structural oxygen cap `O ≤ 2(C+N+S+P)+4` |
| `parse_ranges` | `"C0-30 …"` → element box dict |
| `enumerate_grid` | valence-legal `(mass, formula)` enumeration (H derived from DBE) |
| `_grid_cached` | memoised sorted grid (clears at 16 entries) |
| `complexity_penalty` | per-atom heteroatom score handicap (≤ 0.20) |
| `candidates_for_peaks` | grid pre-filtered to observed peak m/z within ppm |
