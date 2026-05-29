# WSGA landing page — design spec

**Date:** 2026-05-29
**Status:** Approved (brainstorming), pending implementation plan
**Author:** Alvaro Carril
**Scope:** A single static landing page for the `wsga` project, published via GitHub Pages at `alvarocarril.com/wsga/`, that explains the *intuition* of Weighted Subgroup Analysis to a broad academic audience.

---

## 1. Goal & audience

- **Audience:** broad academic (seminar / referee / curious-economist traffic). Primary goal is that the reader **grasps the idea** of WSGA. Software links are secondary.
- **Success:** a reader who has never heard of WSGA leaves understanding *why naive subgroup comparisons conflate the subgroup characteristic with correlated moderators, and how reweighting fixes it* — and can find the paper/packages if they want them.
- **Non-goals:** not a docs site, not a tutorial, not a pkgdown reference. No prose essay. The intuition is told in compact excerpts, not paragraphs.

## 2. Narrative spine

Abstract notation, mirroring the paper: subgroups `G ∈ {0,1}`, moderators `M`, treatment effect difference across subgroups. **5-scene arc**, each beat an *equation + plain-language* pair (technical and intuitive at once):

| # | Scene | Equation (accent) | Intuitive line (the star) | Animation `t` |
|---|-------|-------------------|---------------------------|---------------|
| 1 | **Hook** | — | Title, authors, *"Two subgroups, pulled into balance."* | `t=0` (apart) |
| 2 | **The problem** | `Δ = τ(G=1) − τ(G=0)` | *"Effects differ by subgroup G. But G is tangled with other moderators M. Is the effect different by G — or by what G is correlated with?"* | `t=0`, std-diff ≈ 0.97 highlighted |
| 3 | **The idea** | `wᵢ = f̂_pool(Mᵢ) / f̂_{G(i)}(Mᵢ)` | *"Up-weight the units that **look like** the other group."* (— and down-weight the ones that look like their own.) | `t: 0 → ~0.6`, dots resize |
| 4 | **Balance** | `std.diff(M) → 0` | *"Now M is distributed identically across subgroups. The remaining effect difference is attributable to G, holding M constant."* | `t=1`, curves coincide, std-diff ≈ 0.05 |
| 5 | **The honest caveat + links** | — | *"Causal reading still needs two assumptions: unobservables balanced too, and identification holding within each subgroup."* + CTAs | `t=1`, settled |

Exact equation/copy wording is finalizable during implementation; the arc, notation, and the equation+intuition structure are fixed. The caveat scene is **not optional** — it reflects the paper's conditional-unconfoundedness honesty.

## 3. Visual design — "Generative ink"

- **Palette:** paper `#f1ede3`, ink `#1d2433`, single red accent `#b03a2e`, muted text `#73767e`.
- **Type:**
  - Mono (IBM Plex Mono via CDN, fallback `ui-monospace, 'Courier New', monospace`) for kickers, labels, supporting clauses, plot annotations.
  - Serif (Georgia / Times stack — no CDN needed) for the large intuitive scene-heads.
  - Equations rendered with **KaTeX** (CDN) for crispness; styled small and quiet, left-ruled in red as a supporting accent (not the focal point).
- **Emphasis convention:** key intuitive phrases (e.g. *look like*) set in **red italic** within the serif scene-head.
- The aesthetic reference is creative-coding / gallery generative art (Molnár / Sol LeWitt / pen-plotter), not a product marketing page.

## 4. Layout — scene template

**Split, two-column, sticky.** Per scene:

- **Left column (~55–60%):** mono kicker (`03 · the idea`) → quiet red-ruled KaTeX equation → **large serif intuitive line (the star)** with red-italic emphasis → one muted mono supporting clause.
- **Right column (~40–45%):** a **contained, framed plot panel** (the live canvas), with a small `std. diff` readout. The plot is *not* full-bleed — it is one element balanced against the text.

The plot is **continuous across all 5 scenes** (single sticky canvas); only `t`, the readout, and dot sizes advance. Text columns scroll/cross-fade scene by scene over the sticky canvas.

## 5. The hero visualization (validated prototype exists)

Two subgroups' distributions over a moderator `M`, rendered as fine density curves (ink = G=0, red = G=1) with weight-sized dots along a baseline and dashed weighted-mean ticks.

- **Synthetic DGP:** two Gaussians with genuine common support — `μ₀=0.42, μ₁=0.58, σ=0.17`, N=160 per group. (Common support matters: well-separated groups can't be balanced — itself a real WSGA caveat.)
- **Weights (true IPW, mode 2 — balance both toward pooled):** `wᵢ = f_pool(Mᵢ) / f_{G(i)}(Mᵢ)`, using the **analytic** Gaussian densities (legitimate since the DGP is known), normalized to mean 1 within each group so effective N is preserved.
- **Live readout:** standardized difference in `M`, computed from weighted moments. Validated endpoints: **unweighted ≈ 0.97 → IPW-weighted ≈ 0.05** (matches the paper's 0.68→0.08 spirit). Max weight modest (~6–10), no clipping.
- **Known pitfall (resolved):** KDE-based densities bias the ratio via self-inclusion and floor the weighted std-diff at ~0.25. Must use analytic densities for the weights. Display curves may still use KDE for smoothness; the *weights and readout* must use analytic densities.

Prototype: `.superpowers/brainstorm/.../content/centerpiece-anim-v3.html` (gitignored; reference only).

## 6. Interaction — scroll-scrubbed

- The reader's **scroll position is the animation clock**: scrolling the 5-scene track drives `t ∈ [0,1]`, so they literally scrub the subgroups into balance.
- Implementation: a tall scroll track (≈5 × 100vh) with a sticky canvas; scroll progress → `t`. Eased so motion breathes.
- Optional idle drift is out of scope for v1 (can add later).

## 7. Accessibility & robustness

- **Scene text is real DOM** — readable without JS and by screen readers. The canvas is a progressive enhancement with an `aria-label`/`<figcaption>` summarizing the before/after.
- **`prefers-reduced-motion`:** skip scrubbing; render the balanced end-state statically with all scene text stacked and visible.
- **No-JS fallback:** all five scenes' text render as a normal stacked document.
- **Mobile:** canvas scales; columns stack (text above contained plot); scrub still works via vertical scroll.

## 8. Tech & build

- **Single self-contained `docs/index.html`.** Vanilla JS + Canvas 2D. **No build step, no framework, no bundler.** Only external calls: IBM Plex Mono + KaTeX via CDN (both have graceful fallbacks).
- Rationale: matches static Pages hosting, zero toolchain, easy to maintain, and canvas is the natural medium for the generative art.

## 9. Repo changes

1. **`docs/index.html`** → new landing page (becomes the `/wsga/` page).
2. **Relocate the simulation report:** current `docs/index.html` (the Quarto report) → **`docs/simulation-report.html`**. The report is generated from `simulation/simulation_report.html`; the "copy into docs" step's target filename changes from `index.html` to `simulation-report.html`.
3. **Scene 5 links** to `simulation-report.html`.
4. `specs/` added to `.Rbuildignore` (`^specs$`) — this spec lives at repo root, not under the published `docs/`.

## 10. Footer / CTAs (scene 5)

- GitHub repo (`github.com/acarril/wsga`)
- R install one-liner (`devtools::install_github("acarril/wsga")`) + Stata one-liner (`net install wsga` from the raw GitHub URL)
- Simulation report (`simulation-report.html`)
- Working-paper **citation as text** (Carril, Cazor, Gerardino, Litschig, Pomeranz — "Weighted Subgroup Analysis", working paper). **No PDF link** (none hosted; trivial to add later if a PDF is placed under `docs/`).

## 11. Verification

Serve `docs/` locally and confirm:
- Scrubbing drives the readout 0.97 → 0.05 smoothly across the 5 scenes.
- All five scenes legible; equations render (KaTeX) with fallback if CDN blocked.
- Footer links resolve; `simulation-report.html` loads.
- `prefers-reduced-motion` renders the static balanced end-state.
- Mobile layout: columns stack, plot contained, text readable.
- No-JS: all scene text present as stacked document.

## 12. Out of scope (v1)

- Idle/breathing drift animation.
- The DiD variant of the visualization (RD-flavored intuition only; the method generalizes, but the hero shows one balancing story).
- Interactive parameter controls (sliders for μ/σ/N).
- Hosting the paper PDF.
