"""Build a self-contained rotatable cube from the R-generated mesh data.
Only the Python standard library is required. No model calculations in Python.
Run after 06B_many_tournament_figures.R. The R figures do not need Python.
"""
import json
import csv
from pathlib import Path

root = Path(__file__).resolve().parent
data = json.loads((root / '06B_many_results/allocation_cube_meshes.json').read_text())
rows = list(csv.DictReader((root / '06B_many_results/three_tournament_population.csv').open()))
checks = list(csv.DictReader((root / '06B_many_results/population_convergence.csv').open()))
beta = ', '.join(f"{float(r['beta']):.2f}" for r in rows)
phi = ', '.join(f"{float(r['phi']):.2f}" for r in rows)
n = max(int(r['population_nodes']) for r in checks)
description = f'β = ({beta}), φ = ({phi}), budget B = 50; {n:,} population nodes.'

# Reuse vertices instead of repeating all three coordinates in every triangle.
for parts in data:
    for part in parts:
        lookup, vertices, faces = {}, [], []
        for triangle in part.pop('triangles'):
            face = []
            for point in triangle:
                key = tuple(round(float(v), 6) for v in point)
                if key not in lookup:
                    lookup[key] = len(vertices)
                    vertices.append(key)
                face.append(lookup[key])
            faces.append(face)
        part.update(vertices=vertices, faces=faces)

html = '<!doctype html><html lang="en"><meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n<title>Three asymmetric tournaments — allocation cube</title>\n<style>body{font:16px system-ui,sans-serif;color:#222;max-width:1250px;margin:24px auto;padding:0 18px}h1{font-size:24px;font-weight:600}p{max-width:980px;line-height:1.5}.controls{display:flex;gap:24px;align-items:center;flex-wrap:wrap}main{display:grid;grid-template-columns:1fr 1fr;gap:16px}.panel h2{font-size:19px;text-align:center}canvas{width:100%;height:470px;touch-action:none;cursor:grab;background:#fff}label{white-space:nowrap}.key{display:flex;gap:24px;flex-wrap:wrap;margin:12px 0}.sw{display:inline-block;width:13px;height:13px;margin-right:7px}@media(max-width:760px){main{grid-template-columns:1fr}canvas{height:390px}}</style>\n<h1>Three asymmetric tournaments</h1>\n<p>__CUBE_DESCRIPTION__ Drag either cube to rotate both views. All talent axes run from 0 to 2.5. Gold surfaces show equal positive investment in two jointly leading tournaments; they do not require half of the total budget in either one.</p>\n<div class="controls"><label><input id="specialists" type="checkbox" checked> Specialist boundaries</label><label><input id="equal" type="checkbox" checked> Equal leading investments</label><label>Surface opacity <input id="opacity" type="range" min=".05" max=".65" step=".01" value=".28"></label><button id="reset">Reset view</button></div>\n<div class="key"><span><i class="sw" style="background:#EE0066"></i>Tournament 1</span><span><i class="sw" style="background:#0088FF"></i>Tournament 2</span><span><i class="sw" style="background:#40A56A"></i>Tournament 3</span><span><i class="sw" style="background:#D9AE00"></i>Equal leading investments</span></div>\n<main><section class="panel"><h2>Fixed equal Eₜ reference</h2><canvas aria-label="First-response allocation cube"></canvas></section><section class="panel"><h2>Responsive Eₜ</h2><canvas aria-label="Responsive allocation cube"></canvas></section></main>\n<p>The population reference uses independent lognormal talents with median 1 and SD 0.5. Initial allocations are uniform over the budget simplex. Population integration uses deterministic talent vectors, with the count recorded above. The visible cube is a bounded view of the full talent space.</p>\n<script>\nconst data=PAYLOAD;\nconst canvases=[...document.querySelectorAll(\'canvas\')];let az=-.65,el=.40,pending=false;\nconst vertices=[];for(let x of [0,2.5])for(let y of [0,2.5])for(let z of [0,2.5])vertices.push([x,y,z]);\nconst edges=[];for(let i=0;i<8;i++)for(let j=i+1;j<8;j++)if(vertices[i].filter((v,k)=>v!==vertices[j][k]).length===1)edges.push([i,j]);\nfunction project(p){let x=p[0]-1.25,y=p[1]-1.25,z=p[2]-1.25;return [Math.cos(az)*x-Math.sin(az)*y,Math.sin(az)*Math.sin(el)*x+Math.cos(az)*Math.sin(el)*y+Math.cos(el)*z,Math.sin(az)*Math.cos(el)*x+Math.cos(az)*Math.cos(el)*y-Math.sin(el)*z]}\nfunction draw(){pending=false;const op=+document.querySelector(\'#opacity\').value;\ncanvases.forEach((canvas,row)=>{let box=canvas.getBoundingClientRect(),ratio=devicePixelRatio||1;canvas.width=box.width*ratio;canvas.height=box.height*ratio;let ctx=canvas.getContext(\'2d\');ctx.scale(ratio,ratio);let w=box.width,h=box.height,scale=Math.min(w,h)/4.7;ctx.clearRect(0,0,w,h);\nlet map=p=>[w/2+p[0]*scale,h/2-p[1]*scale];let shapes=[];\ndata[row].forEach((part,index)=>{if(index<3&&!document.querySelector(\'#specialists\').checked)return;if(index>=3&&!document.querySelector(\'#equal\').checked)return;let triangles=part.vertices?part.faces.map(f=>f.map(i=>part.vertices[i])):(part.triangles||[]);triangles.forEach(tri=>{let p=tri.map(project);shapes.push({p,depth:p.reduce((s,v)=>s+v[2],0)/3,col:part.colour,alpha:index<3?Math.min(.9,op*1.8):op})})});\nshapes.sort((a,b)=>a.depth-b.depth);for(let s of shapes){ctx.globalAlpha=s.alpha;ctx.fillStyle=s.col;ctx.beginPath();s.p.forEach((p,i)=>{let q=map(p);if(i)ctx.lineTo(...q);else ctx.moveTo(...q)});ctx.closePath();ctx.fill()}ctx.globalAlpha=1;ctx.strokeStyle=\'#555\';ctx.lineWidth=.8;\nlet V=vertices.map(project);for(let [a,b] of edges){ctx.beginPath();ctx.moveTo(...map(V[a]));ctx.lineTo(...map(V[b]));ctx.stroke()}\nctx.fillStyle=\'#222\';ctx.font=\'17px system-ui\';ctx.textAlign=\'center\';for(let t=0;t<3;t++){let p=[0,0,0];p[t]=2.78;let q=map(project(p));ctx.fillText([\'α₁\',\'α₂\',\'α₃\'][t],q[0],q[1]);ctx.font=\'12px system-ui\';for(let a of [0,1,2]){p=[0,0,0];p[t]=a;if(a===0&&t)continue;q=map(project(p));ctx.fillText(a,q[0],q[1]+14)}ctx.font=\'17px system-ui\'}\n})}\nfunction schedule(){if(!pending){pending=true;requestAnimationFrame(draw)}}\nfor(let c of canvases){let last=null;c.addEventListener(\'pointerdown\',e=>{last=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId)});c.addEventListener(\'pointermove\',e=>{if(!last)return;az+=(e.clientX-last[0])*.009;el=Math.max(-1.35,Math.min(1.35,el-(e.clientY-last[1])*.009));last=[e.clientX,e.clientY];schedule()});c.addEventListener(\'pointerup\',()=>last=null);c.addEventListener(\'pointercancel\',()=>last=null)}\ndocument.querySelectorAll(\'input\').forEach(e=>e.addEventListener(\'input\',schedule));document.querySelector(\'#reset\').onclick=()=>{az=-.65;el=.4;schedule()};window.addEventListener(\'resize\',schedule);draw();\n</script></html>'

(root / 'Figure_many_asymmetric_allocation_cube.html').write_text(
    html.replace('PAYLOAD', json.dumps(data, separators=(',', ':'))).replace('__CUBE_DESCRIPTION__', description), encoding='utf-8')
print('Standalone HTML cube written.')
