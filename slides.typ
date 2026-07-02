#import "touying/lib.typ": *
#import "@preview/pinit:0.1.4": *
#import "@preview/xarrow:0.3.0": xarrow
#import "@preview/cetz:0.4.1"
#import "@preview/mannot:0.3.0": *
#import "psi-slides-0.6.1.typ": *
#import "@preview/algorithmic:1.0.3"
#import algorithmic: algorithm

// color-scheme can be navy-red, blue-green, or pink-yellow
#show: psi-theme.with(aspect-ratio: "16-9",
                      color-scheme: "pink-yellow",
                             config-info(
                                title: [Koopmans functionals in practice],
                                subtitle: [minimisation, screening coefficients, automation, and more...],
                                author: [Edward Linscott],
                                // TODO: update with the real talk date
                                date: datetime(year: 2026, month: 5, day: 27),
                                location: [ICTP--MARVEL College],
                                references: [references.bib],
                             ),
                             config-common(show-strong-with-alert: false))

// Keep body bold pink (replaces touying's globally-applied alert rule);
// bib/footnote overrides below restore inherited colours there.
#show strong: set text(fill: rgb("#dc005a"))

#set footnote.entry(clearance: 0em)
// Citations use a footnote CSL style; drop the footnote number so a lone
// reference shows as plain grey text with no superscript / "1." marker.
#set footnote(numbering: (..) => "")
#show bibliography: it => {
  show strong: set text(fill: black)
  set text(0.6em)
  it
}
#show footnote.entry: it => {
  show strong: set text(fill: rgb("#888888"))
  it
}

// --- helpers ---------------------------------------------------------------
// accent colour matches the pink-yellow theme primary
#let accent = rgb("#dc005a")
#let twitter-blue = rgb("#1da1f2")
// theme colours used by slides ported from the 2025 PsiQuantum deck
#let primary = rgb("#dc005a")
#let secondary = rgb("#f0f500")
#let connector(label) = align(center, stack(spacing: 3pt, text(size: 0.7em, label), $arrow.r.long$))
// DFPT schematic: a supercell calculation = a primitive-cell calculation + a q-point sum
#let dfpt-supercell-diagram = cetz.canvas(length: 1.4cm, {
  import cetz.draw: *

  let cell = 1.6
  let dot = 0.22
  // two-atom basis at fractional (0, 0) and (0.5, 0.5) of each cell
  let basis = ((0, 0), (0.5, 0.5))

  // supercell (2x2), reference cell highlighted
  rect((0, 0), (cell, cell), fill: secondary, stroke: primary + 1.2pt)
  rect((cell, 0), (2 * cell, cell), stroke: primary + 1.2pt)
  rect((0, cell), (cell, 2 * cell), stroke: primary + 1.2pt)
  rect((cell, cell), (2 * cell, 2 * cell), stroke: primary + 1.2pt)
  for c in ((0, 0), (1, 0), (0, 1), (1, 1)) {
    for b in basis {
      circle(
        ((c.at(0) + b.at(0)) * cell, (c.at(1) + b.at(1)) * cell),
        radius: dot, fill: primary, stroke: none,
      )
    }
  }
  content((cell, -0.6), text(size: 0.8em)[supercell])

  content((2 * cell + 1.0, cell), text(weight: "bold", size: 1.8em)[$=$])

  // primitive cell
  let px = 2 * cell + 2.0
  rect((px, cell / 2), (px + cell, cell / 2 + cell), fill: secondary, stroke: primary + 1.2pt)
  for b in basis {
    circle((px + b.at(0) * cell, cell / 2 + b.at(1) * cell), radius: dot, fill: primary, stroke: none)
  }
  content((px + cell / 2, -0.6), text(size: 0.8em)[primitive cell])

  content((px + cell + 1.0, cell), text(weight: "bold", size: 1.8em)[$+$])

  // reciprocal space: q-point mesh
  let qx = px + cell + 2.0
  let qy = cell / 2
  let s = 0.9 * cell
  let ax = 1.5 * cell // axis length = 1.5 x the lattice edge length
  line((qx, qy), (qx + ax, qy), mark: (end: ">", fill: black), stroke: 1.2pt)
  line((qx, qy), (qx, qy + ax), mark: (end: ">", fill: black), stroke: 1.2pt)
  for c in ((0, 0), (1, 0), (0, 1), (1, 1)) {
    circle((qx + c.at(0) * s, qy + c.at(1) * s), radius: dot, fill: luma(140), stroke: none)
  }
  content((qx + s / 2, -0.6), text(size: 0.8em)[$bold(q)$ points])
})
// code listing helper
#let listing(path, lang: none, size: 0.5em) = block(
  fill: luma(245), inset: 7pt, radius: 3pt, width: 100%,
  text(size: size, raw(read(path), lang: lang, block: true)),
)
// show a slice of lines of a file as a code listing (start/end are 0-indexed, end exclusive)
#let listing-lines(path, start, end: none, lang: none, size: 0.5em) = {
  let lines = read(path).split("\n")
  let end = if end == none { lines.len() } else { end }
  block(
    fill: luma(245), inset: 7pt, radius: 3pt, width: 100%,
    text(size: size, raw(lines.slice(start, end).join("\n"), lang: lang, block: true)),
  )
}

#title-slide()

#slide[
  #align(center + horizon, image("figures/fig_en_curve_sl_annotated_zoom.svg", height: 90%))
]

== Core theory

$
  E_"Koopmans" [rho, {f_i}, {alpha_i}]
  = E_"DFT" [rho]
  + sum_i alpha_i
  (- underbrace(integral_0^(f_i) epsilon_i (f) dif f, "removes curvature")
   + underbrace(f_i eta_i, "restores linearity"))
$

#v(1em)

#pause
Differences to semi-local functionals:
#pause
- screening #pause
- different flavours #pause
- orbital-density dependence

@Dabo2010@Borghi2014@Colonna2019

= Screening

== Screening

In *Hartree-Fock* @Li2017
  $
    E_"ee"^"HF" = 1/2 sum_(i j) f_i f_j integral dif bold(r) dif bold(r)'
    (|psi_i (bold(r))|^2 |psi_j (bold(r)')|^2)/(bold(r) - bold(r)')
    - (psi_i^* (bold(r)) psi_j^* (bold(r)') psi_i (bold(r)') psi_j (bold(r)))/(bold(r) - bold(r)')
  #pause
    arrow.r (d^2 E_"ee"^"HF") / (d f_i^2) =^? 0
  $

#pause

#grid(columns: (auto, 1fr), column-gutter: 1em, align: top,
  [
    #image("figures/TjallingKoopmans1967.jpg", height: 40%)
  ],
  [
    #quote(block: true, attribution: [Koopmans' theorem, 1934@Koopmans1934])[
      In Hartree-Fock, the first ionisation energy of a system is equal to the negative of the highest occupied orbital's energy $-epsilon_"HOMO"^"HF"$.
    ]
    $ "IP" = E(N - 1) - E(N) = - epsilon_"HOMO"^"HF" $
    #pause
    ... exact only if the remaining orbitals do not relax.
  ],
)

#pagebreak()
*Koopmans functionals*
$
  E_"KI" [{rho_i}] = & E_"DFT" [rho]
  + sum_i (
    - integral_0^(f_i) epsilon_i (f) dif f
    + f_i integral_0^1 epsilon_i (f) dif f
  )
  #pause
  \ = & E_"DFT" [rho]
  + sum_i (
    - (E_"DFT" [rho] - mark(E_"DFT" [rho^(f_i arrow.r 0)], tag: #<eNm1>, color: #accent))
    + f_i (mark(E_"DFT" [rho^(f_i arrow.r 1)], tag: #<eNp1>, color: #accent) - mark(E_"DFT" [rho^(f_i arrow.r 0)], tag: #<eNm1b>, color: #accent))
  )
$

  #annot(<eNp1>, pos: bottom, dy: 2em)[total energy differences]
  #annot(<eNm1>, pos: bottom, dy: 2em)[cannot be evaluated directly]

#slide[
#align(center + horizon,
  image("figures/fig_pwl_DFT.svg", height: 100%)
)
]
#slide[
#align(center + horizon,
  image("figures/fig_pwl_uKI.svg", height: 100%)
)
]
#slide[
#align(center + horizon,
  image("figures/fig_pwl_alphaKI.svg", height: 100%)
)
]
== Screening: how to calculate the ideal $alpha$?

#align(center)[
  #only("1", image("figures/alpha_calc/fig_alpha_calc_step_0.svg", width: 40%))
  #only("2", image("figures/alpha_calc/fig_alpha_calc_step_1.svg", width: 40%))
  #only("3,4", image("figures/alpha_calc/fig_alpha_calc_step_2.svg", width: 40%))
  #only("5", image("figures/alpha_calc/fig_alpha_calc_step_3.svg", width: 40%))
  #only("6-", image("figures/alpha_calc/fig_alpha_calc_step_4.svg", width: 40%))
]

#only("4,5,6")[
  $
    lambda_(i i) (alpha) eq.triple
    lr(chevron.l phi_i | hat(h)^"DFT" + alpha hat(v)^"Koopmans" | phi_i chevron.r)
    = lr((dif E^"Koopmans") / (dif f_i) |)_(f_i = s)
  $
]
#only("7-")[Given this, how to work out the ideal $alpha$?]
#only("8-")[
  $ alpha_(n+1) = alpha_n (E_i (N - 1) - E(N) - lambda_(i i) (0)) / (lambda_(i i) (alpha_n) - lambda_(i i) (0)) $
]

== Screening

// Reveal each highlight (pink) in sync with its annotation: the mark colour is
// a plain value driven by the current subslide, so the equation re-renders each
// step (avoids touying reveal-wrappers inside the math fraction).
#slide(repeat: 5, self => {
  let ss = self.subslide
  let col(s) = if ss >= s { accent } else { black }
  [
    #v(8em)
    $
      alpha^(n+1) = mark(alpha^n, tag: #<alpha_guess>, color: #(col(1)))
      (mark(E_i (N - 1), tag: #<enm1>, color: #(col(2))) - mark(E(N), tag: #<en>, color: #(col(3))) - lambda_(i i) (0))
      / (mark(lambda_(i i) (alpha^n), tag: #<lama>, color: #(col(4))) - mark(lambda_(i i) (0), tag: #<lam0>, color: #(col(5))))
    $

    #annot(<alpha_guess>, pos: bottom, dy: 4.5em)[trial screening \ parameter]
    #if ss >= 2 { annot(<enm1>, pos: top, dy: -3.5em)[total energy with electron \ removed from orbital $i$] }
    #if ss >= 3 { annot(<en>, pos: top, dy: -1.5em)[total energy \ of neutral system] }
    #if ss >= 4 { annot(<lama>, pos: bottom, dy: 1em)[expectation value \ of $hat(H)^"Koopmans"$] }
    #if ss >= 5 { annot(<lam0>, pos: bottom, dy: 3.5em)[expectation value \ of $hat(H)^"DFT"$] }
    #v(5em)
  ]
})

== Screening via DFPT

How can we avoid explicit charged defect calculations in a supercell?

  #pause
  Reformulate in terms of DFPT@Colonna2019...
  $ alpha_i = 1 + (lr(chevron.l v_"pert"^i | Delta^i n chevron.r)) / (lr(chevron.l n_i | v_"pert"^i chevron.r)) $

  #pause
  ... in reciprocal space@Colonna2022
  $ alpha_(bold(0) i) = 1 + (sum_bold(q) lr(chevron.l v_("pert", bold(q))^(bold(0) i) | Delta_bold(q)^(bold(0) i) n chevron.r)) / (sum_bold(q) lr(chevron.l n_bold(q)^(bold(0) i) | v_("pert", bold(q))^(bold(0) i) chevron.r)) $

#pagebreak()
  #align(center + horizon, dfpt-supercell-diagram)


= Orbital-density-dependence
== Orbital-density-dependence

$
  & - integral_0^(f_i) epsilon_i (f) dif f
  + f_i integral_0^1 epsilon_i (f) dif f \
  =& E_"Hxc" [rho] + E_"Hxc" [rho - rho_i]
  + f_i (- E_"Hxc" [rho - rho_i] + E_"Hxc" [rho - rho_i + n_i])
$

#uncover("2-")[
  For filled orbitals with KI:
  $
    v_i^"KI" \/ alpha_i = - E_"H" [n_i]
    + E_"xc" [rho]
    - E_"xc" [rho - n_i]
    - integral dif bold(r)' thin v_"xc" (bold(r)', [rho]) thin n_i (bold(r)')
  $
]

== Orbital-density-dependence

Minimise $E[rho, rho_i]$ with CG minimisation:
#pause
/ outer loop: $phi_i^((n+1)) = phi_i^((n)) + Delta_i$
#pause
/ inner loop: $phi_i^((n+1)) = sum_j U_(i j) phi_j^((n))$
#pause

For more details see Borghi _et al._ (2015) @Borghi2015
#pause

- gives rise to a set of *minimising orbitals* (localised/variational) #pause
- diagonalising at the minimum gives rise to *diagonalising orbitals* (delocalised/canonical)

#pagebreak()

#grid(columns: (1fr, 1fr), column-gutter: 0em, align: center,
  rotate(90deg, image("figures/fig_nguyen_variational_orbital.png", height: 80%)),
  rotate(90deg, image("figures/fig_nguyen_canonical_orbital.png", height: 80%)),
  [variational],
  [canonical]
)
@Nguyen2018

== Importance of localisation

#align(center, image("figures/nguyen_bulk_limit.png", width: 55%))
@Nguyen2018

#v(-2em)

#only("2,3")[In the bulk limit for one cell $Delta E = E(N - delta N) - E(N)$]

#only("3")[Across all the cells $Delta E = 1/(delta N) (E(N - delta N) - E(N)) = - (dif E)/(dif N) = - epsilon_"HO"$]

#only("4-")[
  $
    lim_(n_i (r) arrow 0) v_i^"KI" \/ alpha_i =
    lim_(n_i (r) arrow 0)
    - E_"H" [n_i] + E_"xc" [rho] - E_"xc" [rho - n_i]
    - integral dif bold(r)' thin v_"xc" (bold(r)', [rho]) thin n_i (bold(r)') = 0
  $
]

#only("5-")[
#align(center, 
[*Solution: apply Koopmans correction to localised orbitals*])
]

== Orbital-density-dependence

Other features of orbital-density-dependence
#pause
- ODD functional means that we know $hat(H) lr(|phi_i chevron.r)$ for variational orbitals ${lr(|phi_i chevron.r)}$ but we don't know $hat(H)$ in general #pause
- Practically we can often use MLWFs #pause
- a natural generalisation in the direction of spectral functional theory@Ferretti2014

= Flavours
== Flavours of Koopmans functionals

#only("1")[$
  E_"Koopmans" [rho, {f_i}, {alpha_i}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i eta_i)
$]
#only("2")[$
  E_(#text(fill: accent)[KI]) [rho, {f_i}, {alpha_i}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i integral_0^1 epsilon_i (f) dif f)
$]
#only("3-")[$
  E_(#text(fill: accent)[KIPZ]) [rho, {f_i}, {alpha_i}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i { integral_0^1 epsilon_i (f) dif f - E_"Hxc" [n_i] })
$]

One degree of freedom: what should be the gradient of this linear term?

#only("2-")[
- the base functional $arrow$ "KI" (Koopmans integral). Enforces IP theorem. Does not affect energy/density!
]
#only("3-")[
- with a PZ correction $arrow$ "KIPZ"
]

#pause You might also see...

#pause
- "pKIPZ" = KIPZ Hamiltonian evaluated on the KI solution #pause
- "K" = an earlier iteration based off half-filling rather than integer endpoints (no longer used)


= To summarise...
== To summarise...

When performing a Koopmans calculation, you must decide...
#pause
- which flavour (KI, KIPZ)? #pause
- how are we treating the orbital-density-dependence?
  - how are we initialising our variational orbitals? (N.B. depends on the flavour!)
  - are we going to explicitly minimise the ODD? #pause
- how are we calculating the screening parameters? (finite differences, DFPT)

= Results
== Molecules

#set text(size: 0.9em)
Ionisation potentials of 100 molecules cf. CCSD(T)
#align(center, image("figures/colonna_2019_gw100_ip.jpeg", height: 29%))

Ultraviolet photoemission spectra
#align(center, image("figures/fig_nguyen_prl_spectra.png", height: 32%))

@Colonna2018@Nguyen2015

== Koopmans functionals: results for solids

#align(horizon, 
grid(columns: (35%, 55%), column-gutter: 2em, align: horizon + left,
  image("figures/fig_nguyen_prx_bandgaps.png", width: 100%),
  [
    #set text(size: 0.85em)
    Mean absolute error (eV) across prototypical semiconductors and insulators

    #v(1em)
    // style the KI and KIPZ columns (3 and 4) in accent bold
    #show table.cell.where(x: 3): set text(fill: accent, weight: "bold")
    #show table.cell.where(x: 4): set text(fill: accent, weight: "bold")
    #table(
      columns: 6,
      gutter: 0.5em,
      align: (left, center, center, center, center, center),
      stroke: none,
      table.header([], [PBE], [G#sub[0]W#sub[0]], [KI], [KIPZ], [QSG#sym.tilde[W]]),
      table.hline(),
      [$E_"gap"$], [2.54], [0.56], [0.27], [0.22], [0.18],
      [IP], [1.09], [0.39], [0.19], [0.21], [0.49],
    )
  ],
)
)
@Nguyen2018

// == Solids
// 
// #align(center)[
//   #set text(size: 0.8em)
//   // style the KI and KIPZ columns (4 and 5) in accent bold
//   #show table.cell.where(x: 4): set text(fill: accent, weight: "bold")
//   #show table.cell.where(x: 5): set text(fill: accent, weight: "bold")
//   #table(
//     columns: 7,
//     align: (left, ) + (center, ) * 6,
//     stroke: none,
//     inset: 4pt,
//     table.hline(), table.hline(y: 0),
//     table.header(
//       [], [PBE], [G#sub[0]W#sub[0]], [scG#sym.tilde("W")],
//       [KI\@[PBE,MLWFs]], [KIPZ\@PBE], [exp]
//     ),
//     table.hline(),
//     [$E_g$], [0.49], [1.06], [1.14], [1.16], [1.15], [1.17],
//     [$Gamma_(1v) arrow Gamma_(25'v)$], [11.97], [12.04], [], [11.97], [12.09], [12.5 ± 0.6],
//     [$X_(1v) arrow Gamma_(25'v)$], [7.82], [], [], [7.82], [], [7.75],
//     [$X_(4v) arrow Gamma_(25'v)$], [2.85], [2.99], [], [2.85], [2.86], [2.90],
//     [$L_(2'v) arrow Gamma_(25'v)$], [9.63], [9.79], [], [9.63], [9.74], [9.3 ± 0.4],
//     [$L_(1v) arrow Gamma_(25'v)$], [6.98], [7.18], [], [6.98], [7.04], [6.8 ± 0.2],
//     [$L_(3'v) arrow Gamma_(25'v)$], [1.19], [1.27], [], [1.19], [], [1.2 ± 0.2],
//     [$Gamma_(25'v) arrow Gamma_(15c)$], [2.48], [3.29], [], [3.17], [3.20], [3.35 ± 0.01],
//     [$Gamma_(25'v) arrow Gamma_(2'c)$], [3.28], [4.02], [], [3.95], [3.95], [4.15 ± 0.05],
//     [$Gamma_(25'v) arrow X_(1c)$], [0.62], [1.38], [], [1.28], [1.31], [1.13],
//     [$Gamma_(25'v) arrow L_(1c)$], [1.45], [2.21], [], [2.12], [2.13], [2.04 ± 0.06],
//     [$Gamma_(25'v) arrow L_(3c)$], [3.24], [4.18], [], [3.91], [3.94], [3.9 ± 0.1],
//     table.hline(stroke: 0.5pt),
//     [MSE], [0.35], [0.02], [], [0.01], [0.03], [],
//     [MAE], [0.44], [0.21], [], [0.14], [0.17], [],
//     table.hline(), table.hline(),
//   )
// ]
// @Shishkin2007@Hybertsen1986@Shishkin2007a@Madelung2004


// == Photocatalytic water-splitting
// 
// #align(horizon + center,
// 
//   image("figures/water_splitting.png", width: 70%)
// )
// #align(horizon + center,
//   image("figures/anatase_water_slab.png", width: 60%),
// )
// @Stojkovic2026
// 
// #pagebreak()
// $
//   epsilon_i^"aligned" &= lr((epsilon_i - lr(chevron.l V_H chevron.r)_("TiO"_2)))_("bulk TiO"_2"")
//   + lr((lr(chevron.l V_H chevron.r)_("TiO"_2) - lr(chevron.l V_H chevron.r)_"water"))_("TiO"_2\/"water")
//   + lr((lr(chevron.l V_H chevron.r)_"water" - V_"vac"))_("water"\/"vac")
// $
// #pause
// 
// #align(center, image("figures/band_offset.svg", width: 50%))
//   
// #pause
// #align(center, [*slab calculations only semi-local!*])
// #v(-2em)
// @Stojkovic2026
// 
// == Band alignment for water-splitting
// 
// #align(center)[
//   #only("1", image("figures/Water_BA_KI_1_potentials.svg", width: 55%))
//   #only("2", image("figures/Water_BA_KI_2_vacuum.svg", width: 55%))
//   #only("3", image("figures/Water_BA_KI_3_solvated.svg", width: 55%))
//   #only("4-", image("figures/Water_BA_KI_4_experiment.svg", width: 55%))
// ]
// @Stojkovic2026@Stojkovic_inprep

== Toy systems

For Hooke's atom (two electrons in a harmonic confining potential with Coulombic repulsion)

#v(1fr)
#grid(columns: (1fr, 1fr), column-gutter: 1em, align: center + horizon,
  image("figures/schubert_vxc.jpeg", width: 80%),
  uncover("2-", image("figures/schubert_vxc_integrated.jpeg", width: 80%)),
)
#v(1fr)
@Schubert2023

== Koopmans functionals: caveats

#pause
- restricted to systems with a non-zero band gap #pause
- empty state localisation in the bulk limit #pause
- can potentially break the crystal point group symmetry@Su2020

// == Resonance with other efforts
// 
// - Wannier transition-state method of Anisimov and Kozhevnikov @Anisimov2005
// - Optimally tuned hybrid functionals of Kronik, Pasquarello, and others (refer back to Leeor's talk on Wednesday) @Kronik2012@Wing2021
// - Ensemble DFT of Kronik and co-workers @Kraisler2013
// - Koopmans-Wannier of Wang and co-workers @Ma2016
// - Dielectric-dependent hybrid functionals of Galli and co-workers @Skone2016a
// - LOSC functionals of Yang and co-workers @Li2018


// === ported from the 2025 PsiQuantum deck (cost/scaling → AiiDA) =========
= Computational cost and scaling
== Computational cost and scaling
#align(center + horizon,
image("figures/timings/benchmark.svg", width: 80%)
)

#pagebreak()

The vast majority of the computational cost: determining screening parameters

$
  alpha_i = (chevron.l n_i|epsilon^(-1) f_"Hxc"|n_i chevron.r) / (chevron.l n_i|f_"Hxc"|n_i chevron.r)
$

#pause

- one screening parameter per orbital #pause
- must be computed #emph[ab initio] via *embarrassingly parallel* steps... #pause
  - $Delta$SCF@Nguyen2018@DeGennaro2022a: $cal(O)(N_"SC"^3) tilde cal(O)(N_bold(k)^3 N^3)$ #pause
  - DFPT@Colonna2018@Colonna2022: $cal(O)(N_bold(k)^2 N^3)$

== Machine-learned electronic screening

#slide[
  #align(
    center,
    grid(
      columns: 5,
      align: horizon,
      gutter: 1em,
      image("figures/orbital.emp.00191_cropped.png", height: 30%),
      $stretch(->)^("power spectrum decomposition")$,
      $vec(delim: "[", x_0, x_1, x_2, dots.v)$,
      $stretch(->)^("ridge regression")$,
      $alpha_i$,
    ),
  )

  $
    c^i_(n l m, k) & = integral dif bold(r) g_(n l) (r) Y_(l m)(theta,phi) n^i (
      bold(r) - bold(R)^i
    )
  $


  $
    p^i_(n_1 n_2 l,k_1 k_2) = pi sqrt(8 / (2l+1)) sum_m c_(n_1 l m,k_1)^(i *) c_(n_2 l m,k_2)^i
  $

  @Schubert2024
]

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 2,
      align: horizon + center,
      gutter: 1em,
      image("figures/water.png", height: 70%),
      image("figures/CsSnI3_disordered.png", height: 70%),

      "water", "CsSnI" + sub("3"),
    ),
  )
  @Schubert2024
]

The use-case

   #grid(columns: 8, column-gutter: 0.3em, row-gutter: 0.3em,
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        grid.cell(align: center + horizon, [...]),
        grid.cell(inset: 0.4em, align: center, fill: primary, colspan: 3, text(fill: white, "train", size: 1em, weight: "bold")),
        grid.cell(inset: 0.4em, align: center, fill: secondary, colspan: 5, text("predict", size: 1em, weight: "bold")),
  )

  #pause
  N.B. not a general model


#slide[
  #grid(
    columns: (1fr, 1fr),
    align: horizon + center,
    gutter: 1em,
    image(
      "figures/water_cls_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),
    image(
      "figures/CsSnI3_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),

    "water", "CsSnI" + sub("3"),
  )
  @Schubert2024
]

#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  @Schubert2024
]


// == Taking advantage of symmetries
// To compute screening parameters via DFPT...
// #algorithm(inset: 0.3em, indent: 1em, {
//   import algorithmic: *
//   Function("CalculateAlpha", ($n$,), {
//     For($bold(q) in "BZ"$,
//     {
//         For($bold(k) in "BZ"$, {Comment[Linear system $A x = b$ to obtain $Delta psi_(bold(k)+bold(q),v)(bold(r))$]})
//           Assign[$Delta rho^(0n)_(q)$][$sum_(bold(k)v)psi^*_(bold(k)v) (bold(r))Delta psi_(bold(k)+bold(q),v)(bold(r)) + c.c.$]
//           Assign[$Pi^((r))_(0 n, bold(q))$][$chevron.l Delta rho^(0 n)_(bold(q))|f_"Hxc"|rho^(0 n)_(bold(q)) chevron.r$]
//           Assign[$Pi^((u))_(0 n, bold(q))$][$chevron.l rho^(0 n)_bold(q)|f_"Hxc"|rho^(0 n)_bold(q) chevron.r$]
//     })
//     Return[$1 + sum_bold(q) Pi^((r))_(0 n, bold(q)) \/ sum_bold(q) Pi^((u))_(0 n, bold(q))$]
//   })
// })
// 
// #pagebreak()
// 
// #align(center,
//   image("figures/bz-to-ibz-outer.svg", height: 80%)
// )
// $bold(q) in "BZ" $ $arrow.r$ $bold(q) in "IBZ"(n)$ (the symmetry of the perturbation; lower than that of the primitive cell)
// #pagebreak()
// #align(center,
//   image("figures/bz-to-ibz-inner.svg", height: 80%)
// )
// $bold(k) in "BZ"$ $arrow.r$ $bold(k) in "IBZ"(bold(q))$ (can only use symmetries that leave $bold(q)$ invariant)
// 
// #align(horizon + center, image("figures/bz-to-ibz-speedup.svg", height: 100%))

= Running a Koopmans functional calculation
== The workflows

The general workflow:
- define/initialise a set of variational orbitals
- calculate the screening parameters ${alpha_i}$
- construct and diagonalise the Hamiltonian

== The workflows

#uncover("1-")[
  (a) finite difference calculations using a supercell
  #align(center, image("figures/supercell_workflow.svg", width: 100%))
]

#uncover("2-")[
  (b) DFPT in a primitive cell
  #align(center, image("figures/primitive_workflow.svg", width: 65%))
]

// == How do I run these calculations?
// 
// Complicated workflows mean that...
// #pause
// - lots of different codes that need to handshake #pause
// - lots of scope for human error #pause
// - reproducibility becomes difficult #pause
// - expert knowledge required
// 
// #pause Our solution...

#focus-slide()[
  #align(center, image("media/logos/koopmans_white_on_transparent.svg", width: 80%))
]

== The #raw("koopmans") package

#grid(columns: (55%, 45%), column-gutter: 1em, align: horizon + left,
  [
    - `v1.0` released in 2023@Linscott2023 (now `v1.2`)
    - implementations of Koopmans functionals
    - automated workflows
      - start-to-finish Koopmans calculations
      - Wannierisation
      - dielectric tensor
      - ...
    - built on top of ASE@Larsen2017
    - under the hood, calls #smallcaps[Quantum ESPRESSO]
    - does not require expert knowledge
  ],
  align(center)[
    #v(1.5em)
    #link("https://koopmans-functionals.org")[`koopmans-functionals.org`]
    #v(-0.5em)
    #image("figures/website_cropped.png", width: 95%)
  ],
)

== koopmans: the input file

#grid(columns: (1fr, 1fr), column-gutter: 1em,
  listing-lines("scripts/si.json", 0, end: 20, lang: "json", size: 0.85em),
  listing-lines("scripts/si.json", 20, lang: "json", size: 0.85em),
)

// == koopmans is scriptable
// 
// #listing("scripts/si.py", lang: "python", size: 0.75em)
// 
// #pause
// but don't get too used to it... #pause

#focus-slide()[
  🚧 koopmans `v2` is coming... 🚧
]

== 
@Huber2020
#v(-2em)
#align(center,
  [
  #grid(columns: 3, align: horizon, column-gutter: 0.5em,
    image("media/logos/koopmans_grey_on_transparent.svg", height: 3em),
    image("figures/handshake.png", height: 2em, alt: "handshake"),
    image("media/logos/aiida.svg", height: 3em)
  )
  ]
)

- remote compute
- parallel step execution
- provenance-tracking
- error recovery
- and more!

#pause
#align(center,
  image("figures/aiida-speed-up.svg", width: 70%)
)


// == Generic structure
// 
// #set text(size: 0.9em)
// #uncover("1-")[/ #raw("Workflow"):]
// #pad(left: 1.5em)[
//   #uncover("2-")[/ #raw("atoms"): an #raw("ASE") #raw("Atoms") object]
//   #uncover("3-")[/ #raw("calculations"): a list of #raw("ASE") calculators]
//   #uncover("4-")[/ #raw("kpoints"): a custom class containing $k$-point information]
//   #uncover("5-")[/ #raw("pseudopotentials"): a dictionary of pseudopotentials]
// ]
// 
// #uncover("6-")[We will see examples in the hands-on!]

== Take home messages

#align(center, grid(columns: 3, column-gutter: 1.5em, align: horizon,
  image("figures/colonna_2019_gw100_ip.jpeg", height: 23%),
  image("figures/fig_nguyen_prx_bandgaps.png", height: 23%),
  image("figures/supercell_workflow.svg", height: 23%),
))

- Koopmans functionals are more complicated than a simple semi-local DFT calculation, because of...
  - orbital-density-dependence
  - screening parameters
- Koopmans functionals are implemented in #smallcaps[Quantum ESPRESSO]
- the complexity of the workflows are handled by the #raw("koopmans") package
- keep your eyes peeled for `v2`!

== Take home messages

#align(center + horizon, image("figures/jctc.png", width: 100%))

#focus-slide()[#align(center, [`koopmans-functionals.org`])]

#focus-slide()[Thank you!]

== Acknowledgements

#align(center + horizon)[
  #grid(columns: 8, column-gutter: 0.4em, row-gutter: 0.4em, align: center,
    image("media/mugshots/nicola_colonna.png", height: 30%),
    image("media/mugshots/miki_bonacci.jpg", height: 30%),
    image("media/mugshots/aleksandr_poliukhin.jpg", height: 30%),
    image("media/mugshots/marija_stojkovic.jpg", height: 30%),
    image("media/mugshots/giovanni_cistaro.jpeg", height: 30%),
    image("media/mugshots/junfeng_qiao.jpeg", height: 30%),
    image("media/mugshots/yannick_schubert.jpg", height: 30%),
    image("media/mugshots/nicola_marzari.jpeg", height: 30%),
    [Nicola Colonna], [Miki Bonacci], [Aleksandr Poliukhin], [Marija Stojkovic], [Giovanni Cistaro], [Junfeng Qiao], [Yannick Schubert], [Nicola Marzari],
  )

  #v(1em)
  #grid(columns: 2, column-gutter: 2em, align: horizon + center,
    image("media/logos/SNF_logo_standard_web_color_pos_e.svg", height: 12%),
    image("media/logos/marvel_color_on_transparent.png", height: 12%),
  )

  #v(1em)
  Want to find out more? Go to #link("https://koopmans-functionals.org")[`koopmans-functionals.org`]

  #v(0.5em)
  slides available at #box(image("logos/github-favicon.png", height: 0.9em)) github/elinscott-talks
]

== References
#bibliography("references.bib")

// Spare slides
#focus-slide[Spare slides]

== Off-diagonal occupancies

#block(fill: luma(235), inset: 12pt, radius: 4pt, width: 100%)[
  *Recap from earlier*

  Key idea: construct a functional such that the _variational_ orbital energies
  $ epsilon_i^"Koopmans" = lr(chevron.l phi_i | H | phi_i chevron.r) = partial E_"Koopmans" \/ partial f_i $
  are...
  - independent of the corresponding occupancies $f_i$
  - equal to the corresponding total energy difference $E_i (N - 1) - E(N)$
]

#v(1em)
zero band gap $arrow$ occupancy matrix for variational orbitals is off-diagonal

== Automated Wannierisation
#slide()[
  Koopmans functionals rely heavily on Wannier functions...
  - to initialise the minmising orbitals, _or_
  - in place of the minimising orbitals entirely

#pause

#grid(
  columns: (2fr, 2fr, 3fr),
  align: center + horizon,
  gutter: 1em,
  image("figures/proj_disentanglement_fig1a.png", height: 45%),
  image("figures/new_projs.png", height: 45%),
  image("figures/target_manifolds_fig1b.png", height: 45%),

  text("projectability-based disentanglement") + cite(<Qiao2023>),
  text("use PAOs found in pseudopotentials"),
  text("parallel transport to separate manifolds") + cite(<Qiao2023a>),
)
]
