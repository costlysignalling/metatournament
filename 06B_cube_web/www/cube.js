(()=>{function init(){

let data=[],appliedBeta=null,appliedPhi=null;
const canvases=[...document.querySelectorAll('canvas')];let az=-.65,el=.40,pending=false;
const vertices=[];for(let x of [0,2.5])for(let y of [0,2.5])for(let z of [0,2.5])vertices.push([x,y,z]);
const edges=[];for(let i=0;i<8;i++)for(let j=i+1;j<8;j++)if(vertices[i].filter((v,k)=>v!==vertices[j][k]).length===1)edges.push([i,j]);
function project(p){let x=p[0]-1.25,y=p[1]-1.25,z=p[2]-1.25;return [Math.cos(az)*x-Math.sin(az)*y,Math.sin(az)*Math.sin(el)*x+Math.cos(az)*Math.sin(el)*y+Math.cos(el)*z,Math.sin(az)*Math.cos(el)*x+Math.cos(az)*Math.cos(el)*y-Math.sin(el)*z]}
function draw(){pending=false;const op=+document.querySelector('#opacity').value;
canvases.forEach((canvas,row)=>{let box=canvas.getBoundingClientRect(),ratio=devicePixelRatio||1;canvas.width=box.width*ratio;canvas.height=box.height*ratio;let ctx=canvas.getContext('2d');ctx.scale(ratio,ratio);let w=box.width,h=box.height,scale=Math.min(w,h)/4.7;ctx.clearRect(0,0,w,h);
let map=p=>[w/2+p[0]*scale,h/2-p[1]*scale];let shapes=[];
(data[row]||[]).forEach((part,index)=>{if(index<3&&!document.querySelector('#specialists').checked)return;if(index>=3&&!document.querySelector('#equal').checked)return;let triangles=part.vertices?part.faces.map(f=>f.map(i=>part.vertices[i])):(part.triangles||[]);triangles.forEach(tri=>{let p=tri.map(project);shapes.push({p,depth:p.reduce((s,v)=>s+v[2],0)/3,col:part.colour,alpha:index<3?Math.min(.9,op*1.8):op})})});
shapes.sort((a,b)=>a.depth-b.depth);for(let s of shapes){ctx.globalAlpha=s.alpha;ctx.fillStyle=s.col;ctx.beginPath();s.p.forEach((p,i)=>{let q=map(p);if(i)ctx.lineTo(...q);else ctx.moveTo(...q)});ctx.closePath();ctx.fill()}ctx.globalAlpha=1;ctx.strokeStyle='#555';ctx.lineWidth=.8;
let V=vertices.map(project);for(let [a,b] of edges){ctx.beginPath();ctx.moveTo(...map(V[a]));ctx.lineTo(...map(V[b]));ctx.stroke()}
ctx.fillStyle='#222';ctx.font='17px system-ui';ctx.textAlign='center';for(let t=0;t<3;t++){let p=[0,0,0];p[t]=2.78;let q=map(project(p));ctx.fillText(['α₁','α₂','α₃'][t],q[0],q[1]);ctx.font='12px system-ui';for(let a of [0,1,2]){p=[0,0,0];p[t]=a;if(a===0&&t)continue;q=map(project(p));ctx.fillText(a,q[0],q[1]+14)}ctx.font='17px system-ui'}
})}
function schedule(){if(!pending){pending=true;requestAnimationFrame(draw)}}
for(let c of canvases){let last=null;c.addEventListener('pointerdown',e=>{last=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId)});c.addEventListener('pointermove',e=>{if(!last)return;az+=(e.clientX-last[0])*.009;el=Math.max(-1.35,Math.min(1.35,el-(e.clientY-last[1])*.009));last=[e.clientX,e.clientY];schedule()});c.addEventListener('pointerup',()=>last=null);c.addEventListener('pointercancel',()=>last=null)}
document.querySelectorAll('input').forEach(e=>e.addEventListener('input',schedule));document.querySelector('#reset').onclick=()=>{az=-.65;el=.4;schedule()};window.addEventListener('resize',schedule);draw();

// Read the currently displayed slider values at click time. Shiny's usual
// slider debounce must not let a rapid click calculate the previous values.
document.querySelector('#compute').addEventListener('click',()=>{
 document.querySelector('#compute').disabled=true;document.querySelector('#restore').disabled=true;
 const value=id=>$('#'+id).data('ionRangeSlider').result.from;
 Shiny.setInputValue('compute_request',{beta:[1,2,3].map(t=>value('b'+t)),
   phi:[1,2,3].map(t=>value('p'+t)),quality:document.querySelector('#quality').value,
   nonce:Date.now()},{priority:'event'});
});
function updateDirtyState(){
 if(!appliedBeta)return;
 const value=id=>{const slider=$('#'+id).data('ionRangeSlider');return slider?slider.result.from:Number(document.getElementById(id).value)};
 const dirty=[1,2,3].some((t,i)=>Math.abs(value('b'+t)-appliedBeta[i])>1e-9||Math.abs(value('p'+t)-appliedPhi[i])>1e-9);
 document.querySelector('#cube-status').textContent=dirty?'Parameters changed. Press Recalculate to update both cubes.':'';
}
Shiny.addCustomMessageHandler('cube-data',message=>{
 data=message.data;appliedBeta=message.beta;appliedPhi=message.phi;
 document.querySelector('#compute').disabled=false;document.querySelector('#restore').disabled=false;
 document.querySelector('#applied-parameters').textContent=message.parameters;
 updateDirtyState();schedule();
});
Shiny.addCustomMessageHandler('cube-status',message=>{
 document.querySelector('#cube-status').textContent=message.message;
 if(message.message.startsWith('Calculation stopped:')){document.querySelector('#compute').disabled=false;document.querySelector('#restore').disabled=false;}
});
$(document).on('change','#b1,#b2,#b3,#p1,#p2,#p3',updateDirtyState);

}if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();})();

