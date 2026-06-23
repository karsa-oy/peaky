"""Standalone INTERACTIVE analyte widgets: a Van Krevelen scatter and a raw
log time-series, each a self-contained HTML file (canvas + hover tooltip, no
server, no external library) where hovering a point/trace shows the neutral
FORMULA and the ION channel(s) it was detected as.

Data comes from peaky.analyte_viz so the figures are computed
identically for any instrument (one row per neutral, Si excluded, RAW intensity,
changing = cv>=0.30).

    python3 scripts/analyte_widgets.py \
        --ledger <LEDGER.csv> --ts-parquet <BATCH_peaks.parquet> \
        --adducts '[M+Br]-,[M-H]-,[M+HBr+Br]-' \
        --out-prefix ~/mascope-output/<name>/<name> --label 'Br⁻ CIMS' --batch <SAMPLE_ID>

--label is the reagent ('Br⁻ CIMS' / 'Ur⁺ CIMS'); --batch is the batch/sample id;
both appear in the widget header as "<label> · <batch>".
Writes <prefix>_vankrevelen.html and <prefix>_timeseries.html.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from peaky import analyte_viz as V  # noqa: E402


def html_van_krevelen(payload: dict, title: str) -> str:
    vk = payload["vk"]
    nf = sum(1 for p in vk if not p[3]); nch = len(vk) - nf
    return (_VK.replace("__DATA__", json.dumps(vk, separators=(",", ":")))
            .replace("__TITLE__", title).replace("__NF__", str(nf)).replace("__NCH__", str(nch)))


def html_timeseries(payload: dict, title: str) -> str:
    ts = payload["ts"] or {"grid": [], "series": []}
    return (_TS.replace("__DATA__", json.dumps(ts, separators=(",", ":")))
            .replace("__TITLE__", title).replace("__N__", str(len(ts["series"]))))


_HEAD = """<!DOCTYPE html><html><head><meta charset="utf-8"><title>__TITLE__</title>
<style>
 body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:18px;color:#1a1a1a;background:#fff;}
 h1{font-size:18px;font-weight:500;}
 .leg{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:#555;margin:8px 0;}
 .sw{width:11px;height:11px;border-radius:50%;display:inline-block;vertical-align:-1px;margin-right:4px;}
 .ln{width:18px;height:0;border-top:3px solid;display:inline-block;vertical-align:3px;margin-right:4px;}
 #wrap{position:relative;width:100%;max-width:980px;height:480px;}
 canvas{width:100%;height:100%;}
 #tip{position:absolute;display:none;pointer-events:none;background:#fff;border:1px solid #bbb;border-radius:6px;
   padding:6px 9px;font-size:12px;line-height:1.55;box-shadow:0 1px 4px rgba(0,0,0,.12);}
 #tip b{font-weight:600;} #search{font-size:13px;padding:4px 8px;width:200px;margin:4px 0;}
</style></head><body>"""


_VK = _HEAD + """
<h1>Van Krevelen — __TITLE__</h1>
<div class="leg">
 <span><span class="sw" style="background:#B4B2A9"></span>flat background (__NF__)</span>
 <span><span class="sw" style="background:#1D9E75"></span>changing CHO</span>
 <span><span class="sw" style="background:#7F77DD"></span>changing CHON</span>
 <span><span class="sw" style="background:#D85A30"></span>changing CHOS</span>
 <span style="margin-left:auto;color:#888">size &prop; log intensity &middot; hover for formula/ion (__NCH__ changing)</span>
</div>
<div id="wrap"><canvas id="cv"></canvas><div id="tip"></div></div>
<script>
const D=__DATA__;
const cv=document.getElementById('cv'),wrap=document.getElementById('wrap'),tip=document.getElementById('tip');
const KCOL={CHO:'#1D9E75',CHON:'#7F77DD',CHOS:'#D85A30'};
const ML=48,MR=14,MT=12,MB=34,OCX=1.3,HCY0=0.3,HCY1=3.0;
let P=[];
function draw(){
 const W=wrap.clientWidth,H=480,dpr=devicePixelRatio||1;
 cv.width=W*dpr;cv.height=H*dpr;const x=cv.getContext('2d');x.setTransform(dpr,0,0,dpr,0,0);x.clearRect(0,0,W,H);
 const px=oc=>ML+Math.min(oc,OCX)/OCX*(W-ML-MR), py=hc=>MT+(HCY1-Math.max(Math.min(hc,HCY1),HCY0))/(HCY1-HCY0)*(H-MT-MB);
 x.strokeStyle='#eee';x.fillStyle='#999';x.font='11px sans-serif';
 for(let h=0.5;h<=3.0;h+=0.5){const y=py(h);x.beginPath();x.moveTo(ML,y);x.lineTo(W-MR,y);x.stroke();x.fillText(h.toFixed(1),10,y+4);}
 for(let o=0;o<=1.2;o+=0.2){const xx=px(o);x.beginPath();x.moveTo(xx,MT);x.lineTo(xx,H-MB);x.stroke();x.fillText(o.toFixed(1),xx-8,H-14);}
 x.fillStyle='#555';x.fillText('O/C',(W)/2,H-2);x.save();x.translate(12,H/2);x.rotate(-Math.PI/2);x.fillText('H/C',0,0);x.restore();
 P=[];
 const order=[...D.keys()].sort((a,b)=>D[a][3]-D[b][3]); // flat first, changing on top
 for(const i of order){const d=D[i];const oc=d[0],hc=d[1],ch=d[3],logI=d[4];
  const X=px(oc),Y=py(hc),r=Math.max(2,(logI-2.5)*1.6);
  const kl=/S\\d|S$/.test(d[5])?'CHOS':(d[2]?'CHON':'CHO'); const col=ch?KCOL[kl]:'#B4B2A9';
  x.globalAlpha=ch?0.85:0.4;x.fillStyle=col;x.beginPath();x.arc(X,Y,r,0,6.283);x.fill();
  P.push([X,Y,r,d[5],d[6],oc,hc,logI,d[7]]);}
 x.globalAlpha=1;
}
cv.onmousemove=e=>{const r=cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;let b=null,bd=1e9;
 for(const p of P){const dd=(p[0]-mx)**2+(p[1]-my)**2;if(dd<bd&&dd<=(p[2]+4)**2){bd=dd;b=p;}}
 if(b){tip.style.display='block';tip.style.left=Math.min(b[0]+12,wrap.clientWidth-180)+'px';tip.style.top=(b[1]-52)+'px';
  tip.innerHTML='<b>'+b[3]+'</b><br>'+(b[4]||'')+'<br>O/C '+b[5].toFixed(2)+' &middot; H/C '+b[6].toFixed(2)+
   '<br>~'+Math.round(10**b[7]).toLocaleString()+' cps'+(b[8]?' &middot; '+b[8]:'');}else tip.style.display='none';};
cv.onmouseleave=()=>tip.style.display='none';
new ResizeObserver(draw).observe(wrap);draw();
</script></body></html>"""


_TS = _HEAD + """
<h1>Raw time series — __TITLE__</h1>
<div class="leg">
 <span><span class="ln" style="border-color:#1D9E75"></span>CHO</span>
 <span><span class="ln" style="border-color:#7F77DD"></span>CHON</span>
 <span style="margin-left:auto;color:#888">__N__ changing analytes &middot; raw cps, log y &middot; hover a trace for formula/ion</span>
</div>
<div id="wrap"><canvas id="cv"></canvas><div id="tip"></div></div>
<script>
const D=__DATA__;
const cv=document.getElementById('cv'),wrap=document.getElementById('wrap'),tip=document.getElementById('tip');
const ML=58,MR=14,MT=12,MB=32;
let SR=[],hov=-1,GMAX=24;
function bounds(){let mn=1e18,mx=1;for(const s of D.series)for(const v of s.y)if(v){mn=Math.min(mn,v);mx=Math.max(mx,v);}
 return [Math.max(10,10**Math.floor(Math.log10(mn))),10**Math.ceil(Math.log10(mx))];}
function draw(){
 const W=wrap.clientWidth,H=480,dpr=devicePixelRatio||1;
 cv.width=W*dpr;cv.height=H*dpr;const x=cv.getContext('2d');x.setTransform(dpr,0,0,dpr,0,0);x.clearRect(0,0,W,H);
 GMAX=D.grid.length?D.grid[D.grid.length-1]:24;const [Y0,Y1]=bounds();const L0=Math.log10(Y0),L1=Math.log10(Y1);
 const px=h=>ML+h/GMAX*(W-ML-MR), py=v=>MT+(L1-Math.log10(Math.max(v,Y0)))/(L1-L0)*(H-MT-MB);
 x.strokeStyle='#eee';x.fillStyle='#999';x.font='11px sans-serif';
 for(let e=L0;e<=L1+0.01;e+=1){const v=10**e,y=py(v);x.beginPath();x.moveTo(ML,y);x.lineTo(W-MR,y);x.stroke();
  x.fillText(v>=1e6?(v/1e6)+'M':v>=1e3?(v/1e3)+'k':v,8,y+4);}
 for(let h=0;h<=GMAX;h+=3){const xx=px(h);x.beginPath();x.moveTo(xx,MT);x.lineTo(xx,H-MB);x.stroke();x.fillText(h,xx-5,H-12);}
 x.fillStyle='#555';x.fillText('hour of experiment (UTC)',W/2-70,H-1);
 SR=[];
 for(let si=0;si<D.series.length;si++){const s=D.series[si];const col=s.n?'#7F77DD':'#1D9E75';
  x.strokeStyle=col;x.globalAlpha=si===hov?1:0.6;x.lineWidth=si===hov?2.4:1.1;
  x.beginPath();let pen=false;const pts=[];
  for(let i=0;i<s.y.length;i++){const v=s.y[i];if(!v){pen=false;continue;}const X=px(D.grid[i]),Y=py(v);
   if(!pen){x.moveTo(X,Y);pen=true;}else x.lineTo(X,Y);pts.push([X,Y,D.grid[i],v]);}
  x.stroke();SR.push(pts);}
 x.globalAlpha=1;x.lineWidth=1;
}
cv.onmousemove=e=>{const r=cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;
 let bs=-1,bi=-1,bd=1e9;
 for(let si=0;si<SR.length;si++)for(const p of SR[si]){const dd=(p[0]-mx)**2+(p[1]-my)**2;if(dd<bd){bd=dd;bs=si;bi=p;}}
 if(bs>=0&&bd<=144){if(hov!==bs){hov=bs;draw();}const s=D.series[bs];
  tip.style.display='block';tip.style.left=Math.min(mx+12,wrap.clientWidth-180)+'px';tip.style.top=(my-48)+'px';
  tip.innerHTML='<b>'+s.f+'</b><br>'+(s.ch||'')+'<br>hour '+bi[2].toFixed(1)+' &middot; '+Math.round(bi[3]).toLocaleString()+' cps';
 }else{tip.style.display='none';if(hov!==-1){hov=-1;draw();}}};
cv.onmouseleave=()=>{tip.style.display='none';if(hov!==-1){hov=-1;draw();}};
new ResizeObserver(draw).observe(wrap);draw();
</script></body></html>"""


def main(argv=None):
    ap = argparse.ArgumentParser(description="interactive analyte VK + time-series HTML widgets")
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--ts-parquet", required=True)
    ap.add_argument("--adducts", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--batch", default="", help="batch / sample id, shown in widget header")
    ap.add_argument("--mode", default="raw", choices=["raw", "reagent", "tic"])
    ap.add_argument("--reagent-mzs", default=None)
    ap.add_argument("--top-ts", type=int, default=28)
    args = ap.parse_args(argv)

    adducts = [a.strip() for a in args.adducts.split(",") if a.strip()]
    reagent_mzs = ([float(x) for x in args.reagent_mzs.split(",")] if args.reagent_mzs else None)
    an = V.analyte_table(pd.read_csv(args.ledger))
    ts = pd.read_parquet(args.ts_parquet)
    an = V.attach_dynamics(an, ts, adducts, mode=args.mode, reagent_mzs=reagent_mzs)
    grid, traces = V.time_traces(ts, an["neutral_formula"].tolist(), adducts,
                                 mode=args.mode, reagent_mzs=reagent_mzs)
    payload = V.widget_payload(an, grid, traces, top_ts=args.top_ts)

    pre = Path(args.out_prefix).expanduser(); pre.parent.mkdir(parents=True, exist_ok=True)
    lab = " · ".join(x for x in (args.label, args.batch) if x) or "analytes"
    Path(f"{pre}_vankrevelen.html").write_text(html_van_krevelen(payload, lab))
    Path(f"{pre}_timeseries.html").write_text(html_timeseries(payload, lab))
    print(f"[{lab}] wrote {pre}_{{vankrevelen,timeseries}}.html "
          f"({len(payload['vk'])} analytes, {len(payload['ts']['series'])} changing traces)")


if __name__ == "__main__":
    main()
