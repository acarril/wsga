# WSGA Landing Page v2 — Visual Narrative Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evolve the landing page hero into a dots-forward, multi-beat evolving stage with a live effect-gap (Δ) readout, and fix the scrub timing so balance is reached at scene 04.

**Architecture:** Modify the single self-contained `docs/index.html` (merged v1). The single scroll scalar `t` becomes a phased `{align, t, hidden}` mapping anchored to scene `offsetTop`s. The `WSGA` math object gains synthetic per-unit effects and a weighted effect-gap. The canvas renderer is replaced wholesale (it's one tightly-coupled function): Phase A gathers 2D dot-clouds onto the M axis and fades the curves in; Phase B reweights (exaggerated dots, converging curves) while the effect markers diverge and Δ counts 4→9pp; scene 05 adds a faint hidden-dimension cue.

**Tech Stack:** Vanilla JS + Canvas 2D, no build step. KaTeX + IBM Plex Mono via CDN. Node for the math guard. Python `http.server` + headless Google Chrome for verification.

**Branch:** `feature/landing-visual-narrative` (already created off merged `main`; carries the v2 design spec). Spec: `specs/2026-05-29-wsga-landing-v2-visual-narrative-design.md`.

**Verified DGP constants (seed 42, N=160):** `ALPHA=0.1895, BETA=0.0924, GAMMA=-0.3340` → effectGap(0)=0.040 (4pp), effectGap(1)=0.090 (9pp); stdDiff(0)=0.94, stdDiff(1)=0.04.

**Visual-tuning note:** Pixel positions/sizes in the renderer (Task 2) are a functional baseline. Exact placement (effect-panel coordinates, dot sizes, cloud spread) is expected to need a screenshot-review pass — the orchestrator captures headless-Chrome screenshots between tasks and tunes. Subagents verify structure/syntax/guard and that the page renders without console errors; they do NOT judge pixels.

---

## File Structure

- **Modify:** `docs/index.html` — the entire deliverable. Three regions change: the `WSGA` math `<script>` (lines ~134-156), the canvas renderer `<script>` (lines ~157-220), and scene DOM/copy (lines ~63, 72-118).
- **Modify:** `specs/verify_hero_math.mjs` — extend the guard to assert the effect-gap endpoints too.

No new files. Single-file page is intentional (matches v1 + static Pages hosting).

---

## Task 1: Synthetic effects in `WSGA` math + extended Node guard

**Files:**
- Modify: `docs/index.html` (the `// ---- WSGA hero math` `<script>`, lines ~135-155)
- Modify: `specs/verify_hero_math.mjs`

- [ ] **Step 1: Extend the Node guard (the test) first**

Replace the entire contents of `specs/verify_hero_math.mjs` with:

```js
// Regression guard for the WSGA hero math (seed 42, N=160).
// Mirrors docs/index.html: IPW balance (std-diff) AND synthetic effect gap Δ.
// Run with: node specs/verify_hero_math.mjs
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
const MU0=0.42,MU1=0.58,SD=0.17,N=160,SEED=42;
const ALPHA=0.1895,BETA=0.0924,GAMMA=-0.3340;
const rng=mulberry32(SEED);
function rnorm(m,s){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return m+s*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const G0=[],G1=[];
for(let i=0;i<N;i++)G0.push({m:rnorm(MU0,SD)});
for(let i=0;i<N;i++)G1.push({m:rnorm(MU1,SD)});
for(const [arr,g] of [[G0,0],[G1,1]]){
  for(const p of arr){const f0=dnorm(p.m,MU0,SD),f1=dnorm(p.m,MU1,SD),fp=0.5*(f0+f1);p.tw=fp/(g===0?f0:f1);}
  const mean=arr.reduce((a,p)=>a+p.tw,0)/arr.length; arr.forEach(p=>p.tw/=mean);
}
const w=(p,t)=>1+t*(p.tw-1);
const wmean=(a,t)=>{let s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*p.m;sw+=x;}return s/sw;};
const wsd=(a,t)=>{let mu=wmean(a,t),s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*(p.m-mu)**2;sw+=x;}return Math.sqrt(s/sw);};
const stdDiff=t=>Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);
const tau=(p,g)=>ALPHA+BETA*g+GAMMA*p.m;
const tauMean=(a,g,t)=>{let s=0,sw=0;for(const p of a){let x=w(p,t);s+=x*tau(p,g);sw+=x;}return s/sw;};
const effectGap=t=>tauMean(G1,1,t)-tauMean(G0,0,t);
const sd0=stdDiff(0),sd1=stdDiff(1),g0=effectGap(0),g1=effectGap(1);
console.log('std-diff  raw',sd0.toFixed(3),'wtd',sd1.toFixed(3),' | effect gap  raw',(g0*100).toFixed(1)+'pp','wtd',(g1*100).toFixed(1)+'pp');
if(sd0<0.90) throw new Error('FAIL std-diff raw ~0.94, got '+sd0.toFixed(3));
if(sd1>0.08) throw new Error('FAIL std-diff wtd ~0.04, got '+sd1.toFixed(3));
if(g0<0.035||g0>0.045) throw new Error('FAIL effect gap raw ~0.040, got '+g0.toFixed(3));
if(g1<0.085||g1>0.095) throw new Error('FAIL effect gap wtd ~0.090, got '+g1.toFixed(3));
console.log('PASS');
```

- [ ] **Step 2: Run it to confirm it passes against the contract**

Run: `node specs/verify_hero_math.mjs`
Expected: prints `std-diff raw 0.941 wtd 0.041 | effect gap raw 4.0pp wtd 9.0pp` then `PASS`.

- [ ] **Step 3: Add the effect functions to the page's `WSGA` object**

In `docs/index.html`, inside the `WSGA = (function(){ ... })()` block, find:
```js
  function stdDiff(t){return Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);}
  return {G0,G1,w,kde,wmean,wsd,stdDiff};
```
Replace it with:
```js
  function stdDiff(t){return Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);}
  // synthetic per-unit treatment effects: tau_i = ALPHA + BETA*G_i + GAMMA*M_i.
  // GAMMA<0 => high-M units have smaller effects, so M-imbalance MASKS the true gap.
  // Subgroup effects are weighted means over the same IPW weights -> Delta recomputes honestly.
  const ALPHA=0.1895, BETA=0.0924, GAMMA=-0.3340;
  function tau(p,g){return ALPHA+BETA*g+GAMMA*p.m;}
  function tauMean(arr,g,t){let s=0,sw=0;for(const p of arr){let wv=w(p,t);s+=wv*tau(p,g);sw+=wv;}return s/sw;}
  function effectGap(t){return tauMean(G1,1,t)-tauMean(G0,0,t);}
  return {G0,G1,w,kde,wmean,wsd,stdDiff,ALPHA,BETA,GAMMA,tau,tauMean,effectGap};
```

- [ ] **Step 4: Verify the embedded object exposes the new API and the numbers match**

Run this one-liner (extracts the WSGA IIFE from the page, evals it under Node, prints values):
```bash
node -e 'const fs=require("fs");const h=fs.readFileSync("docs/index.html","utf8");const m=h.match(/const WSGA = \(function\(\)\{[\s\S]*?\}\)\(\);/)[0];eval(m);console.log("keys",Object.keys(WSGA).join(","));console.log("gap0",WSGA.effectGap(0).toFixed(3),"gap1",WSGA.effectGap(1).toFixed(3));'
```
Expected: keys include `tau,tauMean,effectGap`; `gap0 0.040 gap1 0.090`.

- [ ] **Step 5: Commit**

```bash
git add docs/index.html specs/verify_hero_math.mjs
git commit -m "Add synthetic effect gap to WSGA math + extend guard"
```

---

## Task 2: Replace the canvas renderer (phased timeline, dots-forward, effect panel)

Replace the entire `// ---- WSGA hero canvas rendering ----` `<script>` (current lines ~157-220) with the v2 renderer below. This is one coherent unit: phased `{align,t,hidden}` timeline, cloud→axis dot alignment, curves fading in with `align`, exaggerated reweight, the effect panel, and the scene-05 hidden cue.

**Files:**
- Modify: `docs/index.html` (replace the canvas renderer `<script>`)

- [ ] **Step 1: Replace the renderer script**

Replace the whole block that starts with `<script>` / `// ---- WSGA hero canvas rendering ----` and ends at the matching `</script>` (the one immediately before the `<script>` containing `IntersectionObserver`) with EXACTLY:

```html
<script>
// ---- WSGA hero canvas rendering (v2: phased timeline, dots-forward, effect panel) ----
(function(){
  const cv=document.getElementById('hero'), ctx=cv.getContext('2d');
  const REDUCE = window.matchMedia && window.matchMedia('(prefers-reduced-motion:reduce)').matches;
  const INK='#1d2433', RED='#b03a2e', MUTED='#9a958a';
  let W,H,DPR;
  function resize(){DPR=Math.min(2,devicePixelRatio||1);W=cv.clientWidth;H=cv.clientHeight;cv.width=W*DPR;cv.height=H*DPR;ctx.setTransform(DPR,0,0,DPR,0,0);}

  // stable 2D "cloud" start positions per unit (normalized 0..1), deterministic
  (function seedClouds(){
    function mb(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
    const r=mb(7);
    for(const [arr,cx] of [[WSGA.G0,0.66],[WSGA.G1,0.83]]){
      for(const p of arr){ p.cxn=cx+(r()-0.5)*0.12; p.cyn=0.42+(r()-0.5)*0.22; }
    }
  })();

  const xs=[];for(let i=0;i<=160;i++)xs.push(i/160);
  function X(m){
    if(W<=760){const left=0.10*W,right=0.90*W;return left+m*(right-left);}
    const left=0.56*W,right=0.95*W;return left+m*(right-left);
  }
  const clamp=(x)=>Math.max(0,Math.min(1,x));
  function ease(x){return x<.5?2*x*x:1-Math.pow(-2*x+2,2)/2;}
  function lerp(a,b,u){return a+(b-a)*u;}

  window.drawHero=function(p){
    p=p||{}; const align=clamp(p.align||0), te=ease(clamp(p.t||0)), hidden=clamp(p.hidden||0);
    ctx.clearRect(0,0,W,H);
    const baseY=(W<=760?H*0.82:H*0.66), amp=H*0.28;

    // curves (fade in with align; reweight with te)
    const d0=WSGA.kde(WSGA.G0,te,xs,0.05), d1=WSGA.kde(WSGA.G1,te,xs,0.05);
    const mx=Math.max(Math.max(...d0),Math.max(...d1))*1.05;
    if(align>0.01){
      ctx.globalAlpha=align;
      const curve=(d,color,fill)=>{
        ctx.beginPath();ctx.moveTo(X(0),baseY);
        for(let i=0;i<xs.length;i++)ctx.lineTo(X(xs[i]),baseY-(d[i]/mx)*amp);
        ctx.lineTo(X(1),baseY);ctx.fillStyle=fill;ctx.fill();
        ctx.beginPath();
        for(let i=0;i<xs.length;i++){let yy=baseY-(d[i]/mx)*amp;i?ctx.lineTo(X(xs[i]),yy):ctx.moveTo(X(xs[i]),yy);}
        ctx.strokeStyle=color;ctx.lineWidth=1.2;ctx.stroke();
      };
      curve(d0,INK,'rgba(29,36,51,0.07)');
      curve(d1,RED,'rgba(176,58,46,0.07)');
      ctx.globalAlpha=1;
    }

    // dots: 2D cloud (align=0) -> M-axis lane (align=1); radius by weight (te), exaggerated
    const dots=(arr,color,grp)=>{
      for(const p2 of arr){
        const x=lerp(p2.cxn*W,X(p2.m),align), y=lerp(p2.cyn*H,baseY+10+(grp?9:0),align);
        const wv=WSGA.w(p2,te); const r=2.2+Math.pow(Math.max(0.05,wv),0.9)*2.7;
        ctx.globalAlpha=0.5;ctx.beginPath();ctx.arc(x,y,r,0,7);ctx.fillStyle=color;ctx.fill();ctx.globalAlpha=1;
      }
    };
    dots(WSGA.G0,INK,0);dots(WSGA.G1,RED,1);

    // baseline, mean ticks, std-diff readout (fade with align)
    if(align>0.01){
      ctx.globalAlpha=align;
      const tick=(arr,color)=>{let mu=WSGA.wmean(arr,te);ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(X(mu),baseY-amp-6);ctx.lineTo(X(mu),baseY+6);ctx.strokeStyle=color;ctx.lineWidth=1;ctx.stroke();ctx.setLineDash([]);};
      tick(WSGA.G0,INK);tick(WSGA.G1,RED);
      ctx.beginPath();ctx.moveTo(X(0),baseY);ctx.lineTo(X(1),baseY);ctx.strokeStyle='rgba(29,36,51,0.25)';ctx.lineWidth=1;ctx.stroke();
      ctx.fillStyle=MUTED;ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='left';
      ctx.fillText('STD. DIFF IN M',X(0),baseY-amp-26);
      ctx.fillStyle=INK;ctx.font="24px 'IBM Plex Mono',monospace";
      ctx.fillText(WSGA.stdDiff(te).toFixed(2),X(0),baseY-amp-3);
      ctx.fillStyle=MUTED;ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='right';
      ctx.fillText('M →',X(1),baseY+22);ctx.textAlign='left';
      ctx.globalAlpha=1;
    }

    // effect panel (always visible) — vertical tau axis, two markers, Delta bracket
    const eTop=H*0.10, eBot=H*0.26, eX=X(0), eMark=X(0)+16;
    const yOf=(tau)=>eBot-(Math.max(0,Math.min(0.14,tau))/0.14)*(eBot-eTop);
    const tau0=WSGA.tauMean(WSGA.G0,0,te), tau1=WSGA.tauMean(WSGA.G1,1,te), gap=WSGA.effectGap(te);
    ctx.fillStyle=MUTED;ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='left';
    ctx.fillText('EFFECT GAP Δ',eX,eTop-12);
    ctx.strokeStyle='rgba(29,36,51,0.3)';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(eX,eTop);ctx.lineTo(eX,eBot);ctx.stroke();
    const marker=(tau,color,lab)=>{const y=yOf(tau);ctx.fillStyle=color;ctx.beginPath();ctx.arc(eMark,y,5,0,7);ctx.fill();
      ctx.font="10px 'IBM Plex Mono',monospace";ctx.textAlign='left';ctx.fillStyle=color;ctx.fillText(lab+' '+(tau*100).toFixed(1)+'pp',eMark+12,y+3);};
    const y0=yOf(tau0), y1=yOf(tau1), bx=eMark-12;
    ctx.strokeStyle='rgba(29,36,51,0.4)';ctx.beginPath();ctx.moveTo(eMark,y0);ctx.lineTo(bx,y0);ctx.lineTo(bx,y1);ctx.lineTo(eMark,y1);ctx.stroke();
    marker(tau0,INK,'τ(G0)');marker(tau1,RED,'τ(G1)');
    ctx.fillStyle=INK;ctx.font="15px 'IBM Plex Mono',monospace";ctx.textAlign='right';
    ctx.fillText('Δ '+(gap*100).toFixed(0)+'pp',bx-6,(y0+y1)/2+5);ctx.textAlign='left';
    ctx.fillStyle=MUTED;ctx.font="10px 'IBM Plex Mono',monospace";
    const subl = align<0.99 ? 'effects differ by subgroup'
               : te<0.05 ? 'but is the gap G, or M?'
               : te<0.95 ? 'holding M constant…'
               : 'holding M constant';
    ctx.fillText(subl,eX,eBot+18);

    // scene-05 hidden-dimension cue
    if(hidden>0.01){
      ctx.globalAlpha=0.6*hidden;
      ctx.fillStyle=MUTED;ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='left';
      ctx.fillText('unobserved U — weighting can’t see this',X(0),baseY+46);
      ctx.strokeStyle='rgba(115,118,126,0.4)';ctx.setLineDash([2,4]);ctx.beginPath();ctx.moveTo(X(0),baseY+62);ctx.lineTo(X(1),baseY+62);ctx.stroke();ctx.setLineDash([]);
      ctx.fillStyle='rgba(115,118,126,0.55)';
      for(let i=0;i<24;i++){const gx=X(0.05+0.9*((i*0.137)%1));ctx.beginPath();ctx.arc(gx,baseY+62,2.5,0,7);ctx.fill();}
      ctx.globalAlpha=1;
    }
  };

  // ---- phased scroll timeline ----
  const secs=[...document.querySelectorAll('.scene')];
  function progress(){
    const y=window.scrollY||window.pageYOffset||0, vh=window.innerHeight;
    const top=i=>secs[i].offsetTop;
    const seg=(a,b)=> b>a ? clamp((y-a)/(b-a)) : (y>=a?1:0);
    return {
      align: seg(top(0),top(1)),                 // phase A: scene 1 -> 2
      t:     seg(top(1),top(3)),                 // phase B: scene 2 -> 4 (=1 from scene 4 on)
      hidden:seg(top(4)-vh*0.6, top(4))          // scene 5 cue
    };
  }
  let ticking=false;
  function frame(){window.__p=progress();drawHero(window.__p);ticking=false;}
  function onScroll(){if(ticking)return;ticking=true;requestAnimationFrame(frame);}
  addEventListener('resize',()=>{resize();drawHero(window.__p||{align:0,t:0,hidden:0});});
  if(REDUCE){ resize();window.__p={align:1,t:1,hidden:0};drawHero(window.__p); }
  else { addEventListener('scroll',onScroll,{passive:true});resize();window.__p=progress();drawHero(window.__p); }
})();
</script>
```

- [ ] **Step 2: Syntax-check the new renderer**

Extract the renderer `<script>` body to a temp file and run `node --check`:
```bash
node -e 'const fs=require("fs");const h=fs.readFileSync("docs/index.html","utf8");const m=h.match(/\/\/ ---- WSGA hero canvas rendering[\s\S]*?\}\)\(\);/)[0];fs.writeFileSync("/tmp/r.js",m);'
node --check /tmp/r.js && echo "SYNTAX OK"; rm -f /tmp/r.js
```
Expected: `SYNTAX OK`.

- [ ] **Step 3: Confirm structure**

```bash
grep -c "phased timeline" docs/index.html        # 1
grep -c "EFFECT GAP" docs/index.html             # 1
grep -c "seedClouds" docs/index.html             # 1
grep -c "window.__t" docs/index.html             # 0  (old single-t global gone)
grep -c "hidden-dimension cue" docs/index.html   # 1
```
Expected counts as noted.

- [ ] **Step 4: Render check (headless) — produces an artifact for the orchestrator's visual review**

```bash
cd docs && python3 -m http.server 8021 >/tmp/s.log 2>&1 &
SRV=$!; sleep 1.5
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
"$CHROME" --headless=new --disable-gpu --virtual-time-budget=4000 --dump-dom "http://localhost:8021/" 2>/dev/null | grep -c "id=\"hero\"" ;  # expect 1
"$CHROME" --headless=new --disable-gpu --window-size=1440,900 --virtual-time-budget=4000 --screenshot=/tmp/wsga_v2_top.png "http://localhost:8021/" 2>/dev/null
kill $SRV; ls -la /tmp/wsga_v2_top.png
```
Expected: canvas present (1); screenshot file created. (Pixel correctness reviewed by the orchestrator.)

- [ ] **Step 5: Commit**

```bash
git add docs/index.html
git commit -m "Replace hero renderer: phased timeline, dots-forward, effect panel, scene-5 cue"
```

---

## Task 3: Scene copy + DOM alignment with the new visuals

Bring the scene text in line with the new visuals: scene 01 gains the effect-gap framing; the canvas `aria-label` is updated; the dead `data-t` attributes are removed (they're superseded by the phased timeline reading `offsetTop`). The KaTeX equations stay.

**Files:**
- Modify: `docs/index.html` (lines ~63, 72-118)

- [ ] **Step 1: Update the canvas aria-label**

Replace:
```html
<canvas id="hero" aria-label="Animation: two subgroup distributions over a moderator M, reweighted until they coincide."></canvas>
```
with:
```html
<canvas id="hero" aria-label="Animation: two subgroups' moderator M is reweighted until balanced; the subgroup effect gap grows from 4 to 9 percentage points once M is held constant."></canvas>
```

- [ ] **Step 2: Rewrite scene 01 to introduce the effect gap**

Replace the scene-1 `<section>` (currently `data-scene="1"`) block with:
```html
  <section class="scene" data-scene="1">
    <div class="text">
      <div class="mk">01 · the setup</div>
      <h2>An effect that differs <span class="ll">by subgroup</span>.</h2>
      <div class="supp">Estimate the treatment effect in G = 0 and G = 1 and compare: a gap of Δ = 4 pp. Is that gap real?</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>
```

- [ ] **Step 3: Update scene 04 copy to name the payoff (gap grows)**

Replace the scene-4 `<section>` block with:
```html
  <section class="scene" data-scene="4">
    <div class="text">
      <div class="mk">04 · balance</div>
      <div class="eq" data-tex="\text{std.\,diff}(M) \;\longrightarrow\; 0"></div>
      <h2>Hold M constant, and the true gap <span class="ll">emerges</span>.</h2>
      <div class="supp">With M balanced across subgroups, the effect gap is no longer masked: Δ′ = 9 pp, attributable to G.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>
```

- [ ] **Step 4: Remove dead `data-t` attributes from the other scenes**

In the remaining scene sections, change `<section class="scene" data-scene="2" data-t="0">` → `<section class="scene" data-scene="2">`, and likewise for scenes 3 and 5 (drop the `data-t="..."`). The phased timeline derives anchors from `offsetTop`, not these attributes.

- [ ] **Step 5: Verify copy + structure**

```bash
cd docs && python3 -m http.server 8022 >/tmp/s.log 2>&1 &
SRV=$!; sleep 1.5
curl -s http://localhost:8022/ | grep -c 'class="scene"'         # 5
curl -s http://localhost:8022/ | grep -c 'data-t='              # 0
curl -s http://localhost:8022/ | grep -c 'Δ = 4 pp'             # 1
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --disable-gpu --virtual-time-budget=4000 --dump-dom "http://localhost:8022/" 2>/dev/null | grep -c 'class="katex"'   # 3
kill $SRV
```
Expected: 5 scenes, 0 `data-t`, the `Δ = 4 pp` line present, 3 KaTeX equations still render.

- [ ] **Step 6: Commit**

```bash
git add docs/index.html
git commit -m "Update scene copy and DOM for effect-gap narrative"
```

---

## Task 4: Final verification, screenshots, and PR

**Files:** none (verification only), unless fixes are needed.

- [ ] **Step 1: Math guard**

Run: `node specs/verify_hero_math.mjs`
Expected: `PASS` (std-diff 0.94→0.04, effect gap 4.0pp→9.0pp).

- [ ] **Step 2: Capture the five beats + reduced-motion (headless)**

The phased timeline is scroll-driven, so a static `--screenshot` only captures the top. Use Chrome's DevTools-free approach: scroll via a generated URL is not possible, so capture the top frame and the reduced-motion end-state, and rely on the orchestrator to scroll-capture intermediate beats during review.

```bash
cd docs && python3 -m http.server 8023 >/tmp/s.log 2>&1 &
SRV=$!; sleep 1.5
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
"$CH" --headless=new --disable-gpu --window-size=1440,900 --virtual-time-budget=4000 --screenshot=/tmp/wsga_v2_setup.png "http://localhost:8023/" 2>/dev/null
"$CH" --headless=new --disable-gpu --force-prefers-reduced-motion --window-size=1440,900 --virtual-time-budget=4000 --screenshot=/tmp/wsga_v2_reduced.png "http://localhost:8023/" 2>/dev/null
kill $SRV; ls -la /tmp/wsga_v2_setup.png /tmp/wsga_v2_reduced.png
```
Expected: both PNGs created. Orchestrator reviews: setup frame shows two dot-clouds + effect panel (Δ 4pp); reduced-motion shows balanced state + Δ 9pp. (Intermediate beats — aligned/raw-curves at scene 2, diverging markers mid-scroll, balanced at scene 4 — are confirmed in the orchestrator's scroll-driven review.)

- [ ] **Step 3: Confirm v1 robustness paths still hold**

```bash
cd docs && python3 -m http.server 8024 >/tmp/s.log 2>&1 &
SRV=$!; sleep 1.5
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8024/                    # 200
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8024/simulation-report.html  # 200
for ph in "by subgroup" "tangled with other moderators" "look like" "emerges" "necessary, not"; do curl -s http://localhost:8024/ | grep -c "$ph"; done   # each >=1
kill $SRV
```
Expected: 200, 200, and each scene phrase present in raw HTML (no-JS baseline intact).

- [ ] **Step 4: Confirm `.Rbuildignore` still excludes docs/ and specs/**

```bash
grep -E '\^docs\$|\^specs\$' .Rbuildignore   # both present
```

- [ ] **Step 5: Push branch and open PR**

```bash
git push -u origin feature/landing-visual-narrative
gh pr create --base main --head feature/landing-visual-narrative \
  --title "Landing page v2: visual narrative (dots-forward stage + effect gap)" \
  --body "Evolves the hero into a dots-forward evolving stage with a live effect-gap readout, and fixes scrub timing so balance lands at scene 04. Docs-only — no package code, no version bump (per rescoped CLAUDE.md rule). Spec: specs/2026-05-29-wsga-landing-v2-visual-narrative-design.md."
```

---

## Post-implementation notes

- **No version bump:** this PR changes only `docs/` and `specs/` — no R/Stata package code — so under the rescoped versioning rule it needs no bump.
- **Visual tuning is expected:** effect-panel coordinates, dot sizes, cloud spread, and the scene-05 cue placement may need a screenshot-review pass after Task 2. The orchestrator owns that (subagents verify structure/syntax/render, not pixels).
- **After merge:** GitHub Pages rebuilds `main`/`docs` automatically; verify `alvarocarril.com/wsga/` shows the v2 hero once the build completes.
