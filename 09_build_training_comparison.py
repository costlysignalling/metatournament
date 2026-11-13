"""Build the offline, self-contained figure toggle document. Run after the R script."""
from pathlib import Path
import base64, struct, re

root=Path(__file__).parent
panels=[]
def figure_label(document,label,fallback):
    path=root/document
    if not path.exists():
        return fallback
    source=path.read_text(encoding='utf-8-sig')
    labels=re.findall(r'^```\{r ([^, ]+).*fig\.cap=',source,re.M)
    prefix='Figure S' if document=='Supplement.Rmd' else 'Figure '
    return prefix+str(labels.index(label)+1)
for document,label_id,name,fallback in [
    ('Paper.Rmd','paper-fixed-performance-examples','fixed','Fixed-performance examples'),
    ('Supplement.Rmd','supp-harvesting-beta04','beta04','Mixed tournaments with training β = 0.4'),
    ('Supplement.Rmd','supp-harvesting-beta08','beta08','Mixed tournaments with training β = 0.8'),
]:
    title=figure_label(document,label_id,fallback)
    files=[root/f'Figure_harvesting_examples_{name}.png',root/f'Figure_harvesting_comparison_training_{name}.png']
    sizes=[struct.unpack('>II',p.read_bytes()[16:24]) for p in files]
    assert sizes[0]==sizes[1], (name,sizes)
    images=[]
    for k,p in enumerate(files):
        data=base64.b64encode(p.read_bytes()).decode('ascii')
        label='P = α' if k==0 else 'β = 1'
        images.append(f'<img data-model="{k}" alt="{title}: {label}" src="data:image/png;base64,{data}" '+('' if k==0 else 'hidden')+'>')
    panels.append(f'''<section aria-labelledby="title-{name}">
<div class="toolbar"><h2 id="title-{name}">{title}</h2>
<div class="switch" role="group" aria-label="Model for {title}">
<button type="button" data-choice="0" aria-pressed="true">P = α</button>
<button type="button" data-choice="1" aria-pressed="false">β = 1</button>
</div></div>
<div class="views" style="aspect-ratio:{sizes[0][0]}/{sizes[0][1]}">{''.join(images)}</div>
</section>''')
html='''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Fixed performance and β = 1 comparison</title>
<style>
*{box-sizing:border-box}body{margin:0;background:white;color:#181818;font:16px Arial,sans-serif}
main{max-width:1600px;margin:0 auto;padding:18px 20px 42px}
section{margin:0 0 32px}.toolbar{display:flex;align-items:center;gap:24px;margin:0 0 6px;padding:0 12px}
h2{font-size:19px;font-weight:500;margin:0}.switch{display:flex;gap:0}
button{font:17px Georgia,serif;cursor:pointer;padding:8px 20px;background:white;border:1px solid #aaa;color:#444}
button:first-child{border-radius:5px 0 0 5px}button:last-child{border-left:0;border-radius:0 5px 5px 0}
button[aria-pressed="true"]{background:#222;color:white;border-color:#222}button:focus-visible{outline:3px solid #489fff;outline-offset:3px}
.views{position:relative;width:100%}.views img{position:absolute;inset:0;width:100%;height:100%;object-fit:contain}
[hidden]{display:none!important}@media(max-width:600px){main{padding:12px 4px}.toolbar{gap:18px}button{padding:7px 14px}}
</style><main>'''+''.join(panels)+'''</main><script>
document.querySelectorAll('section').forEach(section=>{
 const buttons=Array.from(section.querySelectorAll('button'));
 const images=Array.from(section.querySelectorAll('img'));
 buttons.forEach(button=>button.addEventListener('click',()=>{
  const model=button.dataset.choice;
  buttons.forEach(b=>b.setAttribute('aria-pressed',String(b===button)));
  images.forEach(img=>{img.hidden=img.dataset.model!==model;});
 }));
});
</script></html>'''
file=root/'Compare_harvesting_P_alpha_vs_beta1.html'
file.write_text(html,encoding='utf-8')
print(f'Built {file.name}; {file.stat().st_size:,} bytes; six aligned images embedded.')
