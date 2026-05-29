# WSGA Landing Page Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single self-contained `docs/index.html` landing page that explains Weighted Subgroup Analysis (WSGA) to a broad academic audience via a 5-scene, scroll-scrubbed generative-ink visualization, published at `alvarocarril.com/wsga/`.

**Architecture:** One static HTML file. Real DOM text for all 5 scenes (works without JS). A sticky full-viewport `<canvas>` renders the hero "balancing" animation; the reader's scroll position over a 5×100vh track drives the animation clock `t ∈ [0,1]`. Per-scene layout is split (text left, contained plot right). Math (synthetic data + true IPW weights + standardized-difference readout) lives in an inline script; its correctness is guarded by a throwaway Node script. Equations rendered with KaTeX (CDN), display font IBM Plex Mono (CDN), both with graceful fallbacks.

**Tech Stack:** Vanilla JS, Canvas 2D, CSS (no framework, no build step). KaTeX + IBM Plex Mono via CDN. Python `http.server` for local serving. Node for the math sanity check.

**Branch:** `feature/landing-page` (already created off `main`). Spec: `specs/2026-05-29-wsga-landing-page-design.md`.

**Palette/type constants (use everywhere):**
- paper `#f1ede3`, ink `#1d2433`, red `#b03a2e`, muted `#73767e`, panel-border `#ddd7c8`
- mono: `'IBM Plex Mono', ui-monospace, 'Courier New', monospace`
- serif: `Georgia, 'Times New Roman', serif`

**Hero DGP/weights contract (validated, do not change):** two Gaussians `μ₀=0.42, μ₁=0.58, σ=0.17`, `N=160` per group; weights `wᵢ = f_pool(Mᵢ)/f_{G(i)}(Mᵢ)` using **analytic** Gaussian densities, normalized to mean 1 within group. Endpoints: unweighted std-diff ≈ 0.97 → fully weighted ≈ 0.05. **Do NOT compute weights from KDE densities** — KDE self-inclusion floors the weighted std-diff at ~0.25.

---

## File Structure

- **Modify → rename:** `docs/index.html` (currently the Quarto simulation report) → `docs/simulation-report.html`
- **Create:** `docs/index.html` (the new landing page — the entire deliverable)
- **Create (throwaway, committed for provenance):** `specs/verify_hero_math.mjs` (Node assertion of the IPW endpoints)

The landing page is intentionally one file. Sections within it, in order: `<head>` (meta, CDN links, `<style>`) → header → 5 `<section class="scene">` blocks inside a sticky-canvas scroll track → footer → inline `<script>` (math + draw + scroll).

---

## Task 1: Relocate the simulation report

**Files:**
- Rename: `docs/index.html` → `docs/simulation-report.html`

- [ ] **Step 1: Confirm current state**

Run: `git -C . ls-files docs/`
Expected: shows `docs/index.html` (the 2.7 MB report).

- [ ] **Step 2: Rename the report**

```bash
git mv docs/index.html docs/simulation-report.html
```

- [ ] **Step 3: Verify**

Run: `ls docs/`
Expected: `simulation-report.html` present, no `index.html`.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "Move simulation report to docs/simulation-report.html"
```

---

## Task 2: Page skeleton with no-JS baseline

Build the full DOM with all 5 scenes as real, readable text and the generative-ink styling. No canvas, no JS behavior yet — this is the accessible/no-JS fallback baseline.

**Files:**
- Create: `docs/index.html`

- [ ] **Step 1: Write the skeleton**

Create `docs/index.html` with this exact content:

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Weighted Subgroup Analysis</title>
<meta name="description" content="Weighted Subgroup Analysis: a reweighting method to hold moderators constant when comparing treatment effects across subgroups.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:ital,wght@0,400;0,500;1,400&display=swap" rel="stylesheet">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css" integrity="sha384-n8MVd4RsNIU0tAv4ct0nTaAbDJwPJzDEaqSD1odI+WdtXRGWt2kTvGFasHpSy3SV" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js" integrity="sha384-XjKyOOlGwcjNTAIQHIpgOno0Hl1YQqzUOEleOLALmuqehneUG+vnGctmUb0ZY0l8" crossorigin="anonymous"></script>
<style>
  :root{
    --paper:#f1ede3; --ink:#1d2433; --red:#b03a2e; --muted:#73767e; --panel-border:#ddd7c8;
    --mono:'IBM Plex Mono',ui-monospace,'Courier New',monospace;
    --serif:Georgia,'Times New Roman',serif;
  }
  *{box-sizing:border-box;}
  html{background:var(--paper);}              /* paper shows behind the fixed canvas */
  html,body{margin:0;color:var(--ink);-webkit-font-smoothing:antialiased;}
  body{font-family:var(--mono);background:transparent;}  /* transparent so canvas (z-index:-1) is visible */

  header.site{padding:6vh 6vw 2vh;background:var(--paper);position:relative;}
  header.site .kicker{font-size:13px;letter-spacing:.34em;text-transform:uppercase;color:var(--red);}
  header.site h1{font-family:var(--serif);font-weight:400;font-size:clamp(34px,6vw,68px);line-height:1.02;letter-spacing:-.02em;margin:18px 0 14px;}
  header.site h1 em{font-style:italic;color:var(--red);}
  header.site .authors{color:var(--muted);font-size:14px;letter-spacing:.04em;}

  /* scene template */
  .scene{min-height:100vh;display:flex;align-items:center;padding:8vh 6vw;gap:5vw;}
  .scene .text{flex:1.15;}
  .scene .panel{flex:1;}
  .scene .mk{font-size:13px;letter-spacing:.3em;text-transform:uppercase;color:var(--red);}
  .scene .eq{margin:18px 0 22px;border-left:2px solid var(--red);padding:5px 0 5px 15px;color:#3a3f4a;font-size:19px;}
  .scene h2{font-family:var(--serif);font-weight:400;font-size:clamp(26px,3.6vw,40px);line-height:1.08;letter-spacing:-.015em;margin:0;}
  .scene h2 .ll{color:var(--red);font-style:italic;}
  .scene .supp{color:var(--muted);font-size:14px;line-height:1.6;margin-top:20px;max-width:46ch;}

  /* contained plot panel (static placeholder until canvas added) */
  .plotframe{width:100%;aspect-ratio:16/11;background:var(--paper);border:1px solid var(--panel-border);border-radius:2px;}

  footer.site{padding:10vh 6vw 12vh;border-top:1px solid var(--panel-border);background:var(--paper);position:relative;}
  footer.site h2{font-family:var(--serif);font-weight:400;font-size:clamp(22px,3vw,32px);margin:0 0 28px;}
  footer.site .links{display:flex;flex-wrap:wrap;gap:32px;margin-bottom:40px;}
  footer.site .links a{color:var(--ink);text-decoration:none;border-bottom:1px solid var(--red);padding-bottom:2px;font-size:15px;}
  footer.site .install{font-size:13px;color:var(--muted);line-height:1.9;white-space:pre-wrap;}
  footer.site .cite{margin-top:36px;font-size:13px;color:var(--muted);max-width:60ch;line-height:1.6;}

  @media(max-width:760px){
    .scene{flex-direction:column;align-items:flex-start;gap:6vh;min-height:auto;padding:10vh 7vw;}
    .scene .panel{width:100%;}
  }
</style>
</head>
<body>
<header class="site">
  <div class="kicker">Weighted Subgroup Analysis</div>
  <h1>Two subgroups,<br>pulled into <em>balance</em>.</h1>
  <div class="authors">Carril · Cazor · Gerardino · Litschig · Pomeranz</div>
</header>

<main id="track">

  <section class="scene" data-scene="1" data-t="0">
    <div class="text">
      <div class="mk">01 · the setup</div>
      <h2>Many studies ask how an effect differs <span class="ll">by subgroup</span>.</h2>
      <div class="supp">Split the sample into G = 0 and G = 1, estimate the treatment effect in each, and compare.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>

  <section class="scene" data-scene="2" data-t="0">
    <div class="text">
      <div class="mk">02 · the problem</div>
      <div class="eq" data-tex="\Delta \;=\; \tau(G{=}1) \;-\; \tau(G{=}0)"></div>
      <h2>But G is tangled with other moderators M.</h2>
      <div class="supp">Is the effect different by G — or by what G is <em>correlated with</em>? The subgroups don't share the same distribution of M.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>

  <section class="scene" data-scene="3" data-t="0.6">
    <div class="text">
      <div class="mk">03 · the idea</div>
      <div class="eq" data-tex="w_i \;=\; \hat f_{\text{pool}}(M_i) \,/\, \hat f_{G(i)}(M_i)"></div>
      <h2>Up-weight the units that <span class="ll">look like</span> the other group.</h2>
      <div class="supp">— and down-weight the ones that look like their own. Estimate it from the propensity to belong to G = 1.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>

  <section class="scene" data-scene="4" data-t="1">
    <div class="text">
      <div class="mk">04 · balance</div>
      <div class="eq" data-tex="\text{std.\,diff}(M) \;\longrightarrow\; 0"></div>
      <h2>Now M is distributed <span class="ll">identically</span> across subgroups.</h2>
      <div class="supp">The remaining difference in effects is attributable to G, holding M constant.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>

  <section class="scene" data-scene="5" data-t="1">
    <div class="text">
      <div class="mk">05 · the honest caveat</div>
      <h2>Observable balance is necessary, not <span class="ll">sufficient</span>.</h2>
      <div class="supp">A causal reading still needs two assumptions: that unobservables are balanced too, and that the design's identifying assumptions hold within each subgroup.</div>
    </div>
    <div class="panel"><div class="plotframe"></div></div>
  </section>

</main>

<footer class="site">
  <h2>Use it.</h2>
  <div class="links">
    <a href="https://github.com/acarril/wsga">GitHub repository</a>
    <a href="simulation-report.html">Monte Carlo simulation report</a>
  </div>
  <div class="install">R:     devtools::install_github("acarril/wsga")
Stata: net install wsga, from("https://raw.githubusercontent.com/acarril/wsga/main/stata/")</div>
  <div class="cite">Carril, Alvaro, Andre Cazor, Maria Paula Gerardino, Stephan Litschig, and Dina Pomeranz. “Weighted Subgroup Analysis.” Working paper.</div>
</footer>

<script>
  // KaTeX render of inline equations (defer-loaded; run on load)
  window.addEventListener('load',function(){
    if(!window.katex)return;
    document.querySelectorAll('.eq[data-tex]').forEach(function(el){
      try{katex.render(el.getAttribute('data-tex'),el,{throwOnError:false,displayMode:false});}catch(e){el.textContent=el.getAttribute('data-tex');}
    });
  });
</script>
</body>
</html>
```

- [ ] **Step 2: Serve and verify the no-JS baseline**

Run: `cd docs && python3 -m http.server 8000`
Open: `http://localhost:8000`
Expected: header with serif title + red-italic "balance" + authors; five full-height scenes each with mono kicker, (scenes 2–4) a red-ruled equation rendered by KaTeX, a large serif headline with red-italic emphasis, and a supporting line; an empty contained plot frame on the right of each; footer with two working-looking links, install one-liners, and the citation. Page is fully readable top-to-bottom by scrolling. Stop server with Ctrl-C.

- [ ] **Step 3: Verify footer link target exists**

Open: `http://localhost:8000/simulation-report.html`
Expected: the Quarto simulation report loads (from Task 1).

- [ ] **Step 4: Commit**

```bash
git add docs/index.html
git commit -m "Add landing page skeleton (no-JS baseline, 5 scenes)"
```

---

## Task 3: Hero math + Node verification guard

Add the pure math (synthetic data, analytic IPW weights, weighted moments, standardized difference) as an inline script, and a throwaway Node script that asserts the validated endpoints so the bug-prone weighting can't silently regress.

**Files:**
- Modify: `docs/index.html` (add a `<script>` block defining `WSGA` math object — place it immediately before the existing KaTeX render script)
- Create: `specs/verify_hero_math.mjs`

- [ ] **Step 1: Write the Node verification script (the test)**

Create `specs/verify_hero_math.mjs`:

```js
// Regression guard for the WSGA hero IPW math.
// Mirrors the algorithm contract in docs/index.html. Run with: node specs/verify_hero_math.mjs
function rnorm(m,s){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return m+s*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
const MU0=0.42,MU1=0.58,SD=0.17,N=160;
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
const sd=t=>Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);
const raw=sd(0), wtd=sd(1);
console.log('unweighted std-diff:',raw.toFixed(3),' weighted std-diff:',wtd.toFixed(3));
if(raw < 0.7) throw new Error('FAIL: unweighted std-diff should be large (~0.97), got '+raw.toFixed(3));
if(wtd > 0.20) throw new Error('FAIL: weighted std-diff should collapse (<0.20, ~0.05), got '+wtd.toFixed(3));
console.log('PASS');
```

- [ ] **Step 2: Run it to verify the contract holds**

Run: `node specs/verify_hero_math.mjs`
Expected: prints `unweighted std-diff: 0.9xx  weighted std-diff: 0.0xx` then `PASS`. (Values vary slightly run-to-run; thresholds are loose.)

- [ ] **Step 3: Add the math to the page**

In `docs/index.html`, insert this `<script>` block immediately BEFORE the existing `<script>` that renders KaTeX:

```html
<script>
// ---- WSGA hero math (pure; no DOM) ----
const WSGA = (function(){
  const MU0=0.42, MU1=0.58, SD=0.17, N=160;
  function rnorm(m,s){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return m+s*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
  function dnorm(x,m,s){return Math.exp(-0.5*((x-m)/s)**2)/(s*Math.sqrt(2*Math.PI));}
  const G0=[],G1=[];
  for(let i=0;i<N;i++)G0.push({m:rnorm(MU0,SD)});
  for(let i=0;i<N;i++)G1.push({m:rnorm(MU1,SD)});
  // true IPW weights toward pooled (mode 2), ANALYTIC densities (NOT KDE)
  for(const [arr,g] of [[G0,0],[G1,1]]){
    for(const p of arr){const f0=dnorm(p.m,MU0,SD),f1=dnorm(p.m,MU1,SD),fp=0.5*(f0+f1);p.tw=fp/(g===0?f0:f1);}
    const mean=arr.reduce((a,p)=>a+p.tw,0)/arr.length; arr.forEach(p=>p.tw/=mean);
  }
  function w(p,t){return 1+t*(p.tw-1);}
  function kde(arr,t,xs,h){return xs.map(x=>{let s=0;for(const p of arr){let wv=w(p,t);let d=(x-p.m)/h;s+=wv*Math.exp(-0.5*d*d);}return s/arr.length;});}
  function wmean(arr,t){let s=0,sw=0;for(const p of arr){let wv=w(p,t);s+=wv*p.m;sw+=wv;}return s/sw;}
  function wsd(arr,t){let mu=wmean(arr,t),s=0,sw=0;for(const p of arr){let wv=w(p,t);s+=wv*(p.m-mu)**2;sw+=wv;}return Math.sqrt(s/sw);}
  function stdDiff(t){return Math.abs(wmean(G1,t)-wmean(G0,t))/Math.sqrt((wsd(G0,t)**2+wsd(G1,t)**2)/2);}
  return {G0,G1,w,kde,wmean,wsd,stdDiff};
})();
</script>
```

- [ ] **Step 4: Verify the page still loads (no behavior change yet)**

Run: `cd docs && python3 -m http.server 8000`
Open `http://localhost:8000`, open the browser console.
Expected: no errors; typing `WSGA.stdDiff(0).toFixed(2)` returns ~`"0.97"` and `WSGA.stdDiff(1).toFixed(2)` returns ~`"0.05"`. Stop server.

- [ ] **Step 5: Commit**

```bash
git add docs/index.html specs/verify_hero_math.mjs
git commit -m "Add hero IPW math with Node verification guard"
```

---

## Task 4: Canvas rendering at a fixed t

Add the sticky canvas and a `drawHero(t)` function. Wire a TEMPORARY range input to set `t` so rendering can be verified before scroll wiring.

**Files:**
- Modify: `docs/index.html` (add `<canvas>`, CSS, draw function, temporary slider)

- [ ] **Step 1: Add the sticky canvas element and CSS**

In `docs/index.html`, immediately after `<main id="track">`'s opening tag is NOT correct — instead add the canvas as the FIRST child of `<body>` (so it underlays content). Add right after `<body>`:

```html
<canvas id="hero" aria-label="Animation: two subgroup distributions over a moderator M, reweighted until they coincide."></canvas>
```

Add to `<style>`:

```css
  #hero{position:fixed;inset:0;width:100vw;height:100vh;z-index:-1;display:block;}
  /* scene panels now host nothing visible; the fixed canvas shows through.
     Keep .plotframe as a sizing spacer but make it invisible. */
  .plotframe{border-color:transparent;background:transparent;}
```

- [ ] **Step 2: Add the draw function and a temporary slider**

Add this `<script>` AFTER the `WSGA` math script and BEFORE the KaTeX render script:

```html
<script>
// ---- WSGA hero canvas rendering ----
(function(){
  const cv=document.getElementById('hero'), ctx=cv.getContext('2d');
  const INK='#1d2433', RED='#b03a2e';
  let W,H,DPR;
  function resize(){DPR=Math.min(2,devicePixelRatio||1);W=cv.clientWidth;H=cv.clientHeight;cv.width=W*DPR;cv.height=H*DPR;ctx.setTransform(DPR,0,0,DPR,0,0);}
  addEventListener('resize',()=>{resize();drawHero(window.__t||0);});

  // plot occupies the right ~42% of the viewport, vertically centered band
  const xs=[];for(let i=0;i<=160;i++)xs.push(i/160);
  function X(m){const left=0.56*W,right=0.95*W;return left+m*(right-left);}
  function ease(x){return x<.5?2*x*x:1-Math.pow(-2*x+2,2)/2;}

  window.drawHero=function(raw){
    const t=ease(Math.max(0,Math.min(1,raw)));
    ctx.clearRect(0,0,W,H);
    const baseY=H*0.66, amp=H*0.30;
    const d0=WSGA.kde(WSGA.G0,t,xs,0.05), d1=WSGA.kde(WSGA.G1,t,xs,0.05);
    const mx=Math.max(Math.max(...d0),Math.max(...d1))*1.05;
    function curve(d,color,fill){
      ctx.beginPath();ctx.moveTo(X(0),baseY);
      for(let i=0;i<xs.length;i++)ctx.lineTo(X(xs[i]),baseY-(d[i]/mx)*amp);
      ctx.lineTo(X(1),baseY);ctx.fillStyle=fill;ctx.fill();
      ctx.beginPath();
      for(let i=0;i<xs.length;i++){let yy=baseY-(d[i]/mx)*amp;i?ctx.lineTo(X(xs[i]),yy):ctx.moveTo(X(xs[i]),yy);}
      ctx.strokeStyle=color;ctx.lineWidth=1.2;ctx.stroke();
    }
    curve(d0,INK,'rgba(29,36,51,0.07)');
    curve(d1,RED,'rgba(176,58,46,0.07)');
    function dots(arr,color,grp){for(const p of arr){let wv=WSGA.w(p,t);let r=1.4+Math.sqrt(Math.max(0.05,wv))*1.7;ctx.globalAlpha=0.5;ctx.beginPath();ctx.arc(X(p.m),baseY+10+(grp?9:0),r,0,7);ctx.fillStyle=color;ctx.fill();ctx.globalAlpha=1;}}
    dots(WSGA.G0,INK,0);dots(WSGA.G1,RED,1);
    function tick(arr,color){let mu=WSGA.wmean(arr,t);ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(X(mu),baseY-amp-6);ctx.lineTo(X(mu),baseY+6);ctx.strokeStyle=color;ctx.lineWidth=1;ctx.stroke();ctx.setLineDash([]);}
    tick(WSGA.G0,INK);tick(WSGA.G1,RED);
    ctx.beginPath();ctx.moveTo(X(0),baseY);ctx.lineTo(X(1),baseY);ctx.strokeStyle='rgba(29,36,51,0.25)';ctx.lineWidth=1;ctx.stroke();
    // readout
    ctx.fillStyle='#9a958a';ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='left';
    ctx.fillText('STD. DIFF IN M',X(0),baseY-amp-26);
    ctx.fillStyle=INK;ctx.font="26px 'IBM Plex Mono',monospace";
    ctx.fillText(WSGA.stdDiff(t).toFixed(2),X(0),baseY-amp-2);
    ctx.fillStyle='#9a958a';ctx.font="11px 'IBM Plex Mono',monospace";ctx.textAlign='right';
    ctx.fillText('M →',X(1),baseY+22);ctx.textAlign='left';
  };

  // TEMPORARY: slider to scrub t (removed in Task 5)
  const sl=document.createElement('input');
  sl.type='range';sl.min=0;sl.max=1;sl.step=0.01;sl.value=0;
  sl.style.cssText='position:fixed;left:20px;bottom:20px;z-index:10;width:260px;';
  sl.id='__tmpslider';
  sl.addEventListener('input',()=>{window.__t=+sl.value;drawHero(+sl.value);});
  document.body.appendChild(sl);

  resize();window.__t=0;drawHero(0);
})();
</script>
```

- [ ] **Step 3: Verify rendering at both extremes**

Run: `cd docs && python3 -m http.server 8000`; open `http://localhost:8000`.
Expected: a contained plot on the right side of the viewport showing two density curves (ink + red) with weight-sized dots, dashed mean ticks, baseline, and a "STD. DIFF IN M" readout. Drag the bottom-left slider: at 0 the curves are apart and readout ≈ 0.97; at 1 the curves coincide and readout ≈ 0.05; dot sizes visibly change. Stop server.

- [ ] **Step 4: Commit**

```bash
git add docs/index.html
git commit -m "Add sticky canvas hero rendering with temporary t slider"
```

---

## Task 5: Scroll-scrubbing

Replace the temporary slider with scroll-driven `t`. The 5-scene `<main>` becomes the scroll track; vertical progress through it maps to `t ∈ [0,1]`.

**Files:**
- Modify: `docs/index.html` (remove slider, add scroll handler; the canvas is already `position:fixed`)

- [ ] **Step 1: Remove the temporary slider**

In the canvas script from Task 4, delete the entire `// TEMPORARY: slider...` block (the lines creating `sl` through `document.body.appendChild(sl);`).

- [ ] **Step 2: Add scroll → t mapping**

Replace the final line `resize();window.__t=0;drawHero(0);` of the canvas IIFE with:

```js
  const track=document.getElementById('track');
  let ticking=false;
  function progress(){
    const r=track.getBoundingClientRect();
    const total=r.height-window.innerHeight;          // scrollable distance within the track
    const scrolled=Math.min(Math.max(-r.top,0),total);
    return total>0 ? scrolled/total : 0;
  }
  function onScroll(){
    if(ticking)return;ticking=true;
    requestAnimationFrame(()=>{window.__t=progress();drawHero(window.__t);ticking=false;});
  }
  addEventListener('scroll',onScroll,{passive:true});
  resize();window.__t=progress();drawHero(window.__t);
```

- [ ] **Step 3: Verify scrubbing**

Run: `cd docs && python3 -m http.server 8000`; open `http://localhost:8000`.
Expected: no slider. Scrolling from top of scene 1 to bottom of scene 5 drives the readout smoothly from ≈0.97 down to ≈0.05; curves slide together; dots resize. Scrolling back up reverses it. Stop server.

- [ ] **Step 4: Commit**

```bash
git add docs/index.html
git commit -m "Drive hero animation by scroll position (scrub)"
```

---

## Task 6: Scene text cross-fade

Fade each scene's text in as it enters the viewport center and out as it leaves, so only the active beat reads strongly over the canvas.

**Files:**
- Modify: `docs/index.html` (CSS for `.scene .text` opacity + IntersectionObserver)

- [ ] **Step 1: Add fade CSS**

Add to `<style>`:

```css
  .scene .text{transition:opacity .5s ease, transform .5s ease;opacity:.12;transform:translateY(12px);}
  .scene .text.active{opacity:1;transform:none;}
  @media(prefers-reduced-motion:reduce){.scene .text{opacity:1;transform:none;transition:none;}}
```

- [ ] **Step 2: Add the observer**

Add this `<script>` after the canvas script, before the KaTeX render script:

```html
<script>
(function(){
  const texts=document.querySelectorAll('.scene .text');
  const io=new IntersectionObserver((entries)=>{
    entries.forEach(e=>{ if(e.isIntersecting) e.target.classList.add('active'); });
  },{rootMargin:'-35% 0px -35% 0px',threshold:0});
  texts.forEach(t=>io.observe(t));
})();
</script>
```

- [ ] **Step 3: Verify**

Run: `cd docs && python3 -m http.server 8000`; open `http://localhost:8000`.
Expected: as you scroll, each scene's text is faint until it reaches the vertical center band, then fades to full opacity and rises slightly; previous scenes fade back. The plot remains continuous throughout. Stop server.

- [ ] **Step 4: Commit**

```bash
git add docs/index.html
git commit -m "Cross-fade scene text on scroll"
```

---

## Task 7: prefers-reduced-motion + no-JS robustness

Ensure the page is sound when motion is disabled or JS is unavailable.

**Files:**
- Modify: `docs/index.html` (reduced-motion: render balanced end-state, no scrub; verify no-JS)

- [ ] **Step 1: Honor reduced-motion in the canvas script**

At the very top of the canvas IIFE body (right after `const cv=...,ctx=...;`), add:

```js
  const REDUCE = window.matchMedia && window.matchMedia('(prefers-reduced-motion:reduce)').matches;
```

Then change the scroll wiring added in Task 5 so that when `REDUCE` is true we do NOT bind scroll and we render the balanced end-state. Replace the Task 5 block's last lines:

```js
  addEventListener('scroll',onScroll,{passive:true});
  resize();window.__t=progress();drawHero(window.__t);
```

with:

```js
  if(REDUCE){ resize();window.__t=1;drawHero(1); }
  else { addEventListener('scroll',onScroll,{passive:true});resize();window.__t=progress();drawHero(window.__t); }
```

- [ ] **Step 2: Verify reduced-motion**

Run: `cd docs && python3 -m http.server 8000`.
On macOS enable System Settings → Accessibility → Display → Reduce motion, reload `http://localhost:8000`.
Expected: canvas shows the balanced end-state (readout ≈0.05, curves coincident) and does NOT change on scroll; all five scenes' text is fully visible (from Task 6 reduced-motion CSS). Disable Reduce motion afterward. Stop server.

- [ ] **Step 3: Verify no-JS**

In the browser, open devtools → command palette → "Disable JavaScript", reload.
Expected: header, all five scenes' text (kickers, headlines, supporting lines), and footer are all present and readable as a stacked document. Equations show their raw TeX text (acceptable fallback). The canvas is empty (fine — it's `z-index:-1` background). Re-enable JS. Stop server.

- [ ] **Step 4: Commit**

```bash
git add docs/index.html
git commit -m "Add reduced-motion end-state and confirm no-JS fallback"
```

---

## Task 8: Responsive / mobile pass

**Files:**
- Modify: `docs/index.html` (canvas plot region on narrow screens)

- [ ] **Step 1: Make the plot use full width on mobile**

The desktop `X(m)` maps the plot to the right 56–95% of the viewport. On narrow screens the text stacks above, so the plot should use the horizontal center instead. Update `X` in the canvas script:

```js
  function X(m){
    if(W<=760){const left=0.10*W,right=0.90*W;return left+m*(right-left);}
    const left=0.56*W,right=0.95*W;return left+m*(right-left);
  }
```

Also lower the plot band on mobile so it sits beside/below text. Update `drawHero`'s `baseY`:

```js
    const baseY = (W<=760 ? H*0.84 : H*0.66), amp=H*0.30;
```

- [ ] **Step 2: Verify mobile layout**

Run: `cd docs && python3 -m http.server 8000`; open devtools device toolbar (~390px wide), reload.
Expected: scenes stack (text block, then space); the plot spans the lower-center of the viewport using most of the width; readout legible; scrubbing still works by scrolling. Check a desktop width too (plot back on the right). Stop server.

- [ ] **Step 3: Commit**

```bash
git add docs/index.html
git commit -m "Responsive plot placement for narrow screens"
```

---

## Task 9: Final verification & polish

**Files:**
- Modify: `docs/index.html` (only if issues found)

- [ ] **Step 1: Full verification pass**

Run: `node specs/verify_hero_math.mjs` → Expected: `PASS`.
Run: `cd docs && python3 -m http.server 8000`; open `http://localhost:8000` and confirm ALL of:
- Scrub drives readout 0.97 → 0.05 smoothly across the 5 scenes.
- All five scenes legible; equations (scenes 2–4) render via KaTeX.
- Footer: GitHub link, simulation-report link both resolve; `simulation-report.html` opens the report; install one-liners and citation present.
- Reduced-motion shows static balanced end-state.
- No-JS shows full stacked text.
- Mobile width: columns stack, plot contained and legible.
- Browser console: no errors.

- [ ] **Step 2: Verify R build still ignores docs/ and specs/**

Run: `cd .. ; grep -E '\^docs\$|\^specs\$' .Rbuildignore`
Expected: both `^docs$` and `^specs$` present (so `R CMD check` ignores the page and plans). No code change needed if present.

- [ ] **Step 3: Commit any fixes**

```bash
git add -A
git commit -m "Final polish for WSGA landing page"
```

---

## Post-implementation (NOT part of automated execution — confirm with user)

These are outward-facing and must be done with the user, not autonomously:

1. **Version bump:** per `CLAUDE.md`, every PR bumps the unified version — invoke the `wsga-version-bump` skill (patch) before opening the PR. The landing page touches no package code, but the repo convention requires it.
2. **Open PR** `feature/landing-page` → `main`.
3. **After merge:** GitHub Pages already serves `main` `/docs`; the new `index.html` goes live at `alvarocarril.com/wsga/` automatically, with the report at `alvarocarril.com/wsga/simulation-report.html`. Verify both URLs after the Pages build completes.
