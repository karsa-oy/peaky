"""Standalone interactive rotating-GKA widget generator.

Reads an assignment ledger CSV (the pipeline's `<id>_ledger.csv`) and writes a
self-contained HTML file: a Generalized Kendrick mass-defect plot with a slider
that rotates the scaling factor X so homologous series flatten into horizontal
rows. Peaks are colored by status (backbone / low-confidence / unassigned) so
you can hunt for structure in the questionable tail.

The band detector uses the mass-accuracy-derived tolerance (Alton et al. trick):

    delta_GKA ~= (X / mass(R)) * delta_m,    delta_m = ppm * m/z * 1e-6

so the horizontal-row width tracks the instrument, not a magic constant.

Usage:
    python3 scripts/gka_widget.py LEDGER.csv [-o out.html] [--ppm 2]
Open the resulting HTML in any browser. No server, no dependencies.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# repeat-unit bases offered in the UI (label -> exact mass, traditional X)
BASES = {
    "CH2": (14.015650, 14), "O": (15.994915, 16), "H2O": (18.010565, 18),
    "CO": (27.994915, 28), "C2H2O": (42.010565, 42), "CO2": (43.989830, 44),
    "CF2": (49.996806, 50), "C2H4O": (44.026215, 44),
    "C3H6O": (58.041865, 58), "C2H4O2": (60.021130, 60),
    "C2H6OSi": (74.018792, 74), "HBr": (79.926160, 80),
}


def _category(row) -> int:
    role = row.get("role")
    if role == "M0":
        conf = str(row.get("confidence", ""))
        return 0 if conf.startswith("High") or conf.startswith("Good") else 1
    if role == "unexplained":
        return 2
    return -1   # iso_child / reagent -> not plotted


def build_points(ledger: pd.DataFrame) -> list[list]:
    rows = []
    for r in ledger.itertuples(index=False):
        d = r._asdict() if hasattr(r, "_asdict") else dict(zip(ledger.columns, r))
        c = _category(d)
        if c < 0:
            continue
        h = d.get("height")
        h = float(h) if h and not pd.isna(h) else 1.0
        rows.append([round(float(d["mz"]), 4),
                     int(round(np.log10(max(h, 1)) * 10)), c])
    return rows


def render_html(points: list[list], title: str, ppm: float) -> str:
    pts = json.dumps(points, separators=(",", ":"))
    bases = json.dumps({k: {"m": v[0], "a": v[1]} for k, v in BASES.items()})
    n0 = sum(1 for p in points if p[2] == 0)
    n1 = sum(1 for p in points if p[2] == 1)
    n2 = sum(1 for p in points if p[2] == 2)
    return _TEMPLATE.format(title=title, pts=pts, bases=bases, ppm=ppm,
                            n0=n0, n1=n1, n2=n2)


_TEMPLATE = r"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>GKA rotating — {title}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:18px;color:#1a1a1a;background:#fff;}}
 h1{{font-size:18px;font-weight:500;}}
 .ctl{{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:10px 0;}}
 .leg{{display:flex;flex-wrap:wrap;gap:14px;font-size:13px;color:#555;margin:6px 0;}}
 label{{font-size:13px;color:#555;}} select,input,button{{font-size:13px;}}
 #wrap{{position:relative;width:100%;max-width:1000px;height:440px;}}
 canvas{{width:100%;height:100%;}}
 #tip{{position:absolute;display:none;pointer-events:none;background:#fff;border:1px solid #ccc;border-radius:6px;padding:6px 9px;font-size:12px;line-height:1.5;}}
 #rows span{{display:inline-block;padding:3px 9px;margin:3px;border:1px solid #ddd;border-radius:6px;color:#555;font-size:12.5px;}}
 .sw{{width:10px;height:10px;border-radius:2px;display:inline-block;}}
</style></head><body>
<h1>Rotating GKA — {title}</h1>
<div class="ctl">
 <label>Base</label><select id="base"></select>
 <label>X</label><input type="range" id="xs" min="4" max="80" step="1" value="14" style="width:240px">
 <span id="xv">14</span>
 <button id="trad">A(R)</button>
 <button id="play">▶ rotate</button>
 <label>ppm</label><input type="number" id="ppm" value="{ppm}" step="0.5" style="width:56px">
</div>
<div class="leg">
 <label><input type="checkbox" id="c0" checked> <span class="sw" style="background:#888780"></span> backbone ({n0})</label>
 <label><input type="checkbox" id="c1" checked> <span class="sw" style="background:#e09b18"></span> low / suspect ({n1})</label>
 <label><input type="checkbox" id="c2" checked> <span class="sw" style="background:#d8453f"></span> unassigned ({n2})</label>
 <label style="margin-left:auto"><input type="checkbox" id="bands" checked> highlight rows</label>
</div>
<div id="wrap"><canvas id="cv"></canvas><div id="tip"></div></div>
<p style="font-size:13px;color:#555">Detected rows (low/unassigned within the mass-accuracy band):</p>
<div id="rows"></div>
<script>
const P={pts}, BASES={bases};
const sel=document.getElementById('base'),xs=document.getElementById('xs'),xv=document.getElementById('xv'),
 cv=document.getElementById('cv'),wrap=document.getElementById('wrap'),tip=document.getElementById('tip'),rowsEl=document.getElementById('rows'),ppmEl=document.getElementById('ppm');
for(const k in BASES){{const o=document.createElement('option');o.value=k;o.textContent=k+' ('+BASES[k].m.toFixed(4)+')';sel.appendChild(o);}}
sel.value='CH2';
const COLS=['#888780','#e09b18','#d8453f'];
const ML=46,MR=12,MT=10,MB=30,MZ0=80,MZ1=760;
let pts=[],curBands=[];
function gkd(mz,bm,X){{const g=mz*X/bm;return g-Math.round(g);}}
function vis(){{return[c0.checked,c1.checked,c2.checked];}}
function draw(){{
 const X=+xs.value,bm=BASES[sel.value].m,ppm=+ppmEl.value;xv.textContent=X;
 const W=wrap.clientWidth,H=440,dpr=devicePixelRatio||1;
 cv.width=W*dpr;cv.height=H*dpr;const ctx=cv.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,W,H);
 const px=mz=>ML+(mz-MZ0)/(MZ1-MZ0)*(W-ML-MR), py=d=>MT+(0.5-d)*(H-MT-MB);
 ctx.strokeStyle='#eee';ctx.fillStyle='#888';ctx.font='11px sans-serif';
 for(let d=-0.4;d<=0.41;d+=0.2){{const y=py(d);ctx.beginPath();ctx.moveTo(ML,y);ctx.lineTo(W-MR,y);ctx.stroke();ctx.fillText(d.toFixed(1),8,y+4);}}
 for(let m=100;m<=700;m+=100){{const x=px(m);ctx.beginPath();ctx.moveTo(x,MT);ctx.lineTo(x,H-MB);ctx.stroke();ctx.fillText(m,x-10,H-12);}}
 const v=vis();pts=[];const qd=[];
 for(const[mz,lh,c]of P){{if(!v[c])continue;const d=gkd(mz,bm,X);pts.push([px(mz),py(d),mz,lh,c,d]);if(c>=1)qd.push([d,mz]);}}
 curBands=[];
 if(document.getElementById('bands').checked&&qd.length){{
  // mass-accuracy-derived band half-width: dGKA = (X/mass(R))*dm, dm=ppm*m/z*1e-6
  const mref=440, bw=(X/bm)*ppm*mref*1e-6;
  const bins={{}};for(const[d,mz]of qd){{const b=Math.round(d/(bw*2));(bins[b]=bins[b]||[]).push(d);}}
  curBands=Object.entries(bins).map(([b,a])=>[a.reduce((x,y)=>x+y,0)/a.length,a.length]).filter(e=>e[1]>=5).sort((a,b)=>b[1]-a[1]).slice(0,10);
  ctx.fillStyle='rgba(224,155,24,0.13)';for(const[d0,n]of curBands){{const y=py(d0);ctx.fillRect(ML,y-3,W-ML-MR,6);}}
 }}
 for(const p of pts){{if(p[4]!==0)continue;ctx.globalAlpha=.45;ctx.fillStyle=COLS[0];ctx.beginPath();ctx.arc(p[0],p[1],1.2+(p[3]-21)*.1,0,6.28);ctx.fill();}}
 ctx.globalAlpha=.85;for(const p of pts){{if(p[4]===0)continue;ctx.fillStyle=COLS[p[4]];ctx.beginPath();ctx.arc(p[0],p[1],1.4+(p[3]-21)*.12,0,6.28);ctx.fill();}}
 ctx.globalAlpha=1;
 rowsEl.innerHTML=curBands.length?curBands.map(([d,n])=>'<span>GKD '+d.toFixed(3)+' · '+n+'</span>').join(''):'<span>none at this X</span>';
}}
xs.oninput=draw;ppmEl.oninput=draw;sel.onchange=()=>{{xs.value=BASES[sel.value].a;draw();}};
for(const id of['c0','c1','c2','bands'])document.getElementById(id).onchange=draw;
document.getElementById('trad').onclick=()=>{{xs.value=BASES[sel.value].a;draw();}};
let t=null;const pb=document.getElementById('play');
pb.onclick=()=>{{if(t){{clearInterval(t);t=null;pb.textContent='▶ rotate';return;}}pb.textContent='⏸ pause';t=setInterval(()=>{{let v=+xs.value+1;if(v>+xs.max)v=+xs.min;xs.value=v;draw();}},300);}};
cv.onmousemove=e=>{{const r=cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;let b=null,bd=100;
 for(const p of pts){{const dd=(p[0]-mx)**2+(p[1]-my)**2;if(dd<bd){{bd=dd;b=p;}}}}
 if(b&&bd<=81){{tip.style.display='block';tip.style.left=Math.min(b[0]+10,wrap.clientWidth-160)+'px';tip.style.top=(b[1]-40)+'px';
  tip.innerHTML='m/z '+b[2].toFixed(4)+'<br>'+['backbone','low/suspect','unassigned'][b[4]]+'<br>GKD '+b[5].toFixed(4);}}else tip.style.display='none';}};
cv.onmouseleave=()=>tip.style.display='none';
new ResizeObserver(draw).observe(wrap);draw();
</script></body></html>"""


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate an interactive rotating-GKA HTML from a ledger CSV")
    ap.add_argument("ledger_csv")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--ppm", type=float, default=2.0, help="mass accuracy for the band width")
    a = ap.parse_args(argv)
    led = pd.read_csv(a.ledger_csv)
    pts = build_points(led)
    out = a.out or (Path(a.ledger_csv).with_suffix("").as_posix() + "_gka.html")
    Path(out).write_text(render_html(pts, Path(a.ledger_csv).stem, a.ppm))
    print(f"wrote {out}  ({len(pts)} points: "
          f"{sum(1 for p in pts if p[2]==0)} backbone / "
          f"{sum(1 for p in pts if p[2]==1)} low / "
          f"{sum(1 for p in pts if p[2]==2)} unassigned)")


if __name__ == "__main__":
    main()
