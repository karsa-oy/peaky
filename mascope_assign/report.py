"""Reporting: turn the finished ledger into human-facing outputs.

  * build_sheets(ledger)  -> dict of DataFrames (pure; unit-tested offline)
  * write_excel(...)      -> multi-sheet .xlsx (needs openpyxl)
  * write_markdown(...)   -> narrative summary with commentary + alternatives

The commentary and close-alternatives the user asked for are generated
mechanically from ledger columns, so they are reproducible.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import contexts as X
from . import ledger as L

__version__ = "0.1.0"


def _alts_to_text(cell) -> str:
    try:
        alts = json.loads(cell) if isinstance(cell, str) else (cell or [])
    except Exception:
        return ""
    parts = []
    for a in alts[:3]:
        s = a.get("ion_score") or a.get("raw_score")
        ppm = a.get("ppm")
        seg = a.get("formula", "?")
        if s is not None:
            seg += f" (score {s:.2f}"
            if ppm is not None:
                seg += f", {ppm:.1f} ppm"
            seg += ")"
        parts.append(seg)
    return "; ".join(parts)


def _iso_to_text(cell) -> str:
    try:
        iso = json.loads(cell) if isinstance(cell, str) else (cell or [])
    except Exception:
        return ""
    return "; ".join(f"{i.get('label')}={i.get('score'):.2f}"
                     for i in iso if i.get("score") is not None)


def build_sheets(ledger: pd.DataFrame, context: str = "ambient-air") -> dict[str, pd.DataFrame]:
    """Return the report sheets as DataFrames."""
    led = ledger.copy()
    m0 = led[led["role"] == L.ROLE_M0].copy()

    # enrich assignments with class + readable alternatives/isotopologues
    if len(m0):
        cls = m0["neutral_formula"].map(lambda f: X.classify_compound(f))
        m0["compound_class"] = [c[0] for c in cls]
        m0["oxidation"] = [c[1] for c in cls]
        m0["heteroatoms"] = [c[2] for c in cls]
        m0["alternatives_text"] = m0["alternatives"].map(_alts_to_text)
        m0["isotopologues_text"] = m0["isotopologues"].map(_iso_to_text)

    assignments = m0[[
        "peak_id", "mz", "height", "neutral_formula", "ion_formula", "adduct",
        "dbe", "ion_score", "compound_score", "ppm_error", "confidence",
        "compound_class", "oxidation", "heteroatoms", "pass_no", "method",
        "isotopologues_text", "alternatives_text", "commentary",
    ]].sort_values(["confidence", "ion_score"], ascending=[True, False]) if len(m0) else m0

    # by class
    by_class = (m0.groupby(["compound_class", "heteroatoms"])
                .agg(n=("peak_id", "count"), signal=("height", "sum"))
                .reset_index()) if len(m0) else pd.DataFrame()

    # unique formulas
    uniq = (m0.groupby("neutral_formula")
            .agg(n_peaks=("peak_id", "count"), best_score=("ion_score", "max"),
                 signal=("height", "sum"))
            .reset_index().sort_values("best_score", ascending=False)) if len(m0) else pd.DataFrame()

    # target list (formula + adduct + best ppm)
    target = (m0[["neutral_formula", "adduct", "ion_formula", "mz", "ppm_error",
                  "ion_score", "confidence"]]
              .sort_values("mz")) if len(m0) else pd.DataFrame()

    # isotopologues (one row per attributed child)
    iso = led[led["role"] == L.ROLE_ISO][[
        "peak_id", "mz", "height", "parent_peak_id", "iso_label", "iso_match_score"]].copy()

    # unassigned
    un = led[led["role"] == L.ROLE_UNEXPLAINED][["peak_id", "mz", "height"]].copy()
    un = un.sort_values("height", ascending=False)

    # reagent
    reag = led[led["role"] == L.ROLE_REAGENT][["peak_id", "mz", "height", "commentary"]].copy()

    # ownership audit (one row per physical peak)
    ownership = led[["peak_id", "mz", "height", "role", "neutral_formula",
                     "adduct", "ion_score", "ppm_error", "confidence",
                     "parent_peak_id", "iso_label", "pass_no", "method",
                     "commentary"]].sort_values("height", ascending=False)

    return {
        "Assignments": assignments,
        "By class": by_class,
        "Unique formulas": uniq,
        "Target list": target,
        "Isotopologues": iso,
        "Peak ownership": ownership,
        "Unassigned": un,
        "Reagent ions": reag,
    }


def summary_stats(ledger: pd.DataFrame) -> pd.DataFrame:
    st = L.stats(ledger)
    rows = [{"metric": "n_peaks", "value": st["n_peaks"]}]
    for role, n in st["by_role"].items():
        rows.append({"metric": f"peaks[{role}]", "value": n})
        rows.append({"metric": f"peaks%[{role}]",
                     "value": round(100 * st["count_frac_by_role"][role], 1)})
        rows.append({"metric": f"signal%[{role}]",
                     "value": round(100 * st["signal_by_role"][role], 1)})
    for conf, n in st.get("by_confidence", {}).items():
        rows.append({"metric": f"confidence[{conf}]", "value": n})
    return pd.DataFrame(rows)


def write_excel(ledger: pd.DataFrame, path: str | Path, context: str = "ambient-air"):
    sheets = build_sheets(ledger, context)
    sheets["Summary"] = summary_stats(ledger)
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            (df if len(df) else pd.DataFrame({"(empty)": []})).to_excel(
                xl, sheet_name=name[:31], index=False)
    return path


def write_markdown(result: dict, path: str | Path) -> Path:
    led = result["ledger"]
    st = result["stats"]
    m0 = led[led["role"] == L.ROLE_M0]
    lines = [
        f"# Peak assignment — sample {result['sample_id']}",
        "",
        f"- Context: **{result['context']}**",
        f"- Peaks: {st['n_peaks']}  |  assigned (M0): {st['by_role']['M0']}  "
        f"|  isotopologues: {st['by_role']['iso_child']}  "
        f"|  unexplained: {st['by_role']['unexplained']}",
        (f"- Peaks explained: "
         f"{100*(1 - st['count_frac_by_role']['unexplained']):.1f}% "
         f"({st['by_role']['unexplained']}/{st['n_peaks']} unexplained)  |  "
         if "count_frac_by_role" in st else "- ")
        + f"Signal explained: "
        f"{100*(st['signal_by_role']['M0']+st['signal_by_role']['iso_child']):.1f}%",
        f"- Confidence: {st.get('by_confidence', {})}",
        f"- Prescan: {result['prescan']}",
        "",
        "## Top assignments",
        "",
    ]
    top = m0.sort_values("ion_score", ascending=False).head(20)
    for _, r in top.iterrows():
        lines.append(f"- **{r['neutral_formula']}** {r['adduct']} "
                     f"(m/z {r['mz']:.4f}, {r['confidence']}) — {r['commentary']}")
    if result.get("problems"):
        lines += ["", "## ⚠ Ledger validation problems", ""]
        lines += [f"- {p}" for p in result["problems"]]
    Path(path).write_text("\n".join(lines))
    return Path(path)
