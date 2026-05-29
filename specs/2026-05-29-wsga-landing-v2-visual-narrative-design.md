# WSGA landing page v2 — richer visual narrative (design spec)

**Date:** 2026-05-29
**Status:** Approved (brainstorming), pending implementation plan
**Builds on:** `specs/2026-05-29-wsga-landing-page-design.md` (v1, shipped in PR #48 / branch `feature/landing-page`)
**Scope:** Evolve the single-canvas hero from one visual idea (M-distributions reweighting) into a **dots-forward, multi-beat evolving stage** that gives every narrative scene its own visual support — and fix the scrub timing so balance is reached at scene 04.

---

## 1. Motivation (two user comments)

1. **Timing (easy):** the `std-diff` minimum currently lands at the very bottom (end of scene 05) because `progress()` maps scroll linearly across the whole 5-scene track. It should reach its minimum exactly as **scene 04 (Balance)** arrives; scene 05 holds.
2. **Visual storytelling (the substance):** today only the M-distribution reweighting is visual, which serves only scenes 03–04. Scenes 01, 02, 05 ride on text. In particular the **problem** (scene 02) — that subgroup effect differences are confounded by M — is never shown. The fix: a single *evolving stage* where dots are the primary actors and a live **effect-gap (Δ) readout** carries the causal story.

The chosen direction merges two explored options: **C** (dots aligning + reweighting = the *mechanism*) as the visual spine, and **B** (an effect gap Δ = the *problem* and *payoff*) as a compact live readout. They are complementary, not exclusive.

## 2. What stays the same

Generative-ink aesthetic, palette/type, split scene layout (text left / visual right), scroll-scrubbed interaction, `prefers-reduced-motion` + no-JS + mobile fallbacks, single self-contained `docs/index.html`, KaTeX + IBM Plex Mono via CDN, seeded determinism (mulberry32, **seed 42**, N=160 per group), analytic-density IPW weights (mode 2 toward pooled). The existing `WSGA` math object is extended, not replaced.

## 3. The evolving stage — three coordinated elements

The visual area holds three elements, all driven by scroll, all always showing true state (no dead/"recomputing" states):

1. **Dots** — individual units, primary actors. G0 = ink `#1d2433`, G1 = red `#b03a2e`. Dot radius encodes IPW weight.
2. **Density curves** — the *balance verdict*. Introduced at scene 02 (raw, separated); converge during reweighting. (Reuses the existing KDE display path; KDE is display-only — never used for weights.)
3. **Effect panel (Δ)** — a compact live mini-chart: a short vertical "effect τ" axis with two markers τ(G0) ink / τ(G1) red and the gap Δ bracketed between them, plus the existing `std-diff in M` number. Both markers and numbers animate continuously.

## 4. Scroll timeline — two phases (replaces single linear `t`)

`progress()` is replaced by a phased mapping over the 5-scene track. Define two scalars derived from scroll position:

- **`align ∈ [0,1]`** — Phase A, spans scene 01 → scene 02. Drives: dot-clouds gathering and aligning onto the M axis (at equal size), and the raw curves fading in. No reweighting.
- **`t ∈ [0,1]`** — Phase B, spans scene 02 → scene 04. Drives reweighting: dot radii, curve convergence, `std-diff` 0.97→0.04, and Δ 4pp→9pp. Eased.

`t` reaches 1 exactly as scene 04's top reaches the viewport; `t` holds at 1 through scene 05 (satisfies motivation #1). Concretely, anchor `align` to the scene-01→02 scroll span and `t` to the scene-02→04 scroll span, using each scene section's geometry (the existing per-scene `data-t` markers can encode anchor targets, or anchors are computed from section offsets).

Reduced-motion: render the final state — `align=1, t=1` (aligned, balanced, Δ=9pp) — statically, no scroll binding (as v1 did for `t=1`).

## 5. Beat-by-beat

| Scene | Dots | Curves | Effect panel (Δ) | M balance |
|---|---|---|---|---|
| **01 setup** | two faint clouds (2D, off-axis) | — | τ(G0), τ(G1) shown; **Δ = 4pp (naïve)** | — |
| **01→02** | gather + align onto M axis, **equal size** | **fade in** (raw, separated) | Δ = 4pp, flagged **"confounded?"** | `std-diff 0.97` appears |
| **02→04** | **resize** (reweight, exaggerated) | **slide into coincidence** | τ markers **diverge**, Δ **counts 4→9pp** live | `std-diff` 0.97 → 0.04 live |
| **04 balance** | balanced, settled | coincident | **Δ′ = 9pp**, labeled *"holding M constant"* | `std-diff 0.04` |
| **05 caveat** | held + faint **hidden-dimension cue** | held | Δ′ = 9pp, footnote *"…if unobservables balance too"* | held |

**Scene 05 hidden-dimension cue (in scope):** a subtle greyed axis or scatter of pale points the weighting does not touch, visually saying "an unobserved dimension remains" — giving scene 05 its own minimal visual. Keep it quiet (low contrast) so it doesn't compete with the settled balance.

**Dot reweighting must be exaggerated** relative to v1 (user note) so the size change reads clearly during Phase B.

**Visual juxtaposition (the "aha"):** as the M-distributions *converge*, the two τ markers *diverge* — "once the groups are comparable, the real (larger) effect difference emerges."

## 6. Synthetic effects (verified, seed 42, N=160)

Each unit has a per-unit treatment effect `τᵢ = α + β·Gᵢ + γ·Mᵢ`. Subgroup effects are **weighted means of `τᵢ` over the same units and IPW weights that drive the dots**, so Δ is genuinely recomputed from the one weight vector (not asserted). The naïve gap `Δ(t) = β + γ·(M̄₁(t) − M̄₀(t))`: balancing M removes the confound term, so weighted `Δ′ → β`. `γ < 0` means high-M units have smaller effects, so the raw imbalance *masks* the true gap (the Fujiwara "effect grows" story).

**Coefficients (verified to hit targets exactly on seed-42 data):**
- `α = 0.1895`, `β = 0.0924`, `γ = −0.3340`

**Realized values:**
- t=0 (raw): τ(G0) = 5.0pp, τ(G1) = 9.0pp, **Δ = 4.0pp**
- t=1 (weighted): τ(G0) = 2.3pp, τ(G1) = 11.3pp, **Δ = 9.0pp**

Displayed as percentage points (effects, e.g. "+9 pp"). `α` only sets marker height (cancels in Δ); β, γ set the gap story. If the DGP/seed ever changes, re-solve β, γ from the realized raw/weighted M-gaps.

## 7. Math object additions (`WSGA`)

Add to the existing object (alongside `G0,G1,w,kde,wmean,wsd,stdDiff`):
- effect coefficients `ALPHA, BETA, GAMMA`;
- `tau(p, g)` → `ALPHA + BETA*g + GAMMA*p.m` (per-unit effect);
- `tauMean(arr, g, t)` → weighted mean of `tau` over a group at reweighting level `t`;
- `effectGap(t)` → `tauMean(G1,1,t) − tauMean(G0,0,t)`.

`stdDiff` and `effectGap` both read the same `w(p,t)`.

## 8. Scene copy

Scene 01 gains the Δ-effect framing (it currently is "the setup"). The existing KaTeX equations (scenes 02–04) stay; scene 01 may gain a small `Δ = τ(G=1) − τ(G=0)` accent (the scene-02 equation already states this — keep one, avoid duplication). Exact wording finalizable in implementation; the beat structure and the equation/intuition pairing are fixed.

## 9. Implementation implications

- `progress()` → two-phase `{align, t}` mapping from section geometry.
- `drawHero` → gains: (a) the Phase-A dot gather/align from 2D clouds to the M axis, (b) the effect panel rendering (τ axis, two markers, Δ bracket + value), (c) exaggerated dot-radius mapping, (d) curves fading in at scene 02, (e) the scene-05 hidden-dimension cue.
- `WSGA` math gains the effect functions (§7).
- Node guard (`specs/verify_hero_math.mjs`) extends to assert **both** endpoints: `stdDiff` 0.97→~0.04 **and** `effectGap` ~0.040→~0.090 (tight, deterministic under seed 42).
- Single file, no build step, all v1 robustness paths preserved.

## 10. Verification

`node specs/verify_hero_math.mjs` → PASS asserting both metrics' endpoints. Headless-Chrome captures at: scene 01 (clouds + Δ=4pp), scene 02 (aligned + raw curves + confounded flag), mid Phase-B (resizing + diverging markers), scene 04 (balanced, std-diff 0.04, Δ=9pp), reduced-motion (final state static). Confirm scene-05 cue renders; mobile layout holds; no-JS shows all scene text; `t`=1 reached at scene 04 (not 05).

## 11. Out of scope (v2)

- Full B "morph" (a separate large effect-chart subsystem) — we chose the lightweight live Δ panel instead.
- DiD variant of the visualization.
- Interactive controls (sliders for coefficients/seed).
- Changing the v1 footer, links, report relocation, or page chrome.
