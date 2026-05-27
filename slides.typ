#import "touying/lib.typ": *
#import "@preview/pinit:0.1.4": *
#import "@preview/xarrow:0.3.0": xarrow
#import "@preview/cetz:0.3.3"
#import "@preview/mannot:0.3.0": *
#import "psi-slides-0.6.1.typ": *

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
#let connector(label) = align(center, stack(spacing: 3pt, text(size: 0.7em, label), $arrow.r.long$))
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

// ===========================================================================
== Core theory

$
  E_"Koopmans" [rho, #text(fill: accent)[$\{f_i\}$], #text(fill: accent)[$\{alpha_i\}$]]
  = E_"DFT" [rho]
  + sum_i #text(fill: accent)[$alpha_i$]
  (- underbrace(integral_0^(f_i) epsilon_i (f) dif f, "removes curvature")
   + underbrace(f_i #text(fill: accent)[$eta_i$], "restores linearity"))
$

#v(1em)

#pause
Differences to semi-local functionals:
- #pause different flavours
- #pause orbital-density dependence
- #pause screening

@Dabo2010@Borghi2014@Colonna2019

// ===========================================================================
== Flavours of Koopmans functionals

#only("1")[$
  E_"Koopmans" [rho, \{f_i\}, \{alpha_i\}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i eta_i)
$]
#only("2")[$
  E_(#text(fill: accent)[KI]) [rho, \{f_i\}, \{alpha_i\}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i integral_0^1 epsilon_i (f) dif f)
$]
#only("3-")[$
  E_(#text(fill: accent)[KIPZ]) [rho, \{f_i\}, \{alpha_i\}] = E_"DFT" [rho]
  + sum_i alpha_i (- integral_0^(f_i) epsilon_i (f) dif f + f_i \{ integral_0^1 epsilon_i (f) dif f - E_"Hxc" [n_i] \})
$]

One degree of freedom: what should be the gradient of this linear term?

- #pause the base functional $arrow$ "KI" (Koopmans integral). Enforces IP theorem. Does not affect energy/density!
- #pause with a PZ correction $arrow$ "KIPZ"

#pause You might also see...

- #pause "pKIPZ" = KIPZ Hamiltonian evaluated on the KI solution
- #pause "K" = an earlier iteration based off half-filling rather than integer endpoints (no longer used)

// ===========================================================================
== Orbital-density-dependence

$
  - integral_0^(f_i) epsilon_i (f) dif f
  + f_i integral_0^1 epsilon_i (f) dif f
  = E_"Hxc" [rho] + E_"Hxc" [rho - rho_i]
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

// ===========================================================================
== Orbital-density-dependence

#uncover("2-")[Initialise with MLWFs, then (optionally) solve with CG minimisation:]
#uncover("3-")[/ outer loop: $phi_i^((n+1)) = phi_i^((n)) + Delta_i$]
#uncover("4-")[/ inner loop: $phi_i^((n+1)) = sum_j U_(i j) phi_j^((n))$]

#uncover("4-")[For more details see @Borghi2015]

#uncover("5-")[Gives rise to a set of minimising orbitals (localised/variational)]

#uncover("6-")[Diagonalising at the minimum gives rise to diagonalising orbitals (delocalised/canonical)]

#v(1fr)
#grid(columns: (1fr, 1fr), column-gutter: 2em, align: center,
  uncover("5-", figure(rotate(90deg, reflow: true, image("figures/fig_nguyen_variational_orbital.png", height: 50%)), caption: [variational])),
  uncover("6-", figure(rotate(90deg, reflow: true, image("figures/fig_nguyen_canonical_orbital.png", height: 50%)), caption: [canonical])),
)
@Nguyen2018

// ===========================================================================
== Importance of localisation

#align(center, image("figures/nguyen_bulk_limit.png", width: 55%))
@Nguyen2018

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

// ===========================================================================
== Orbital-density-dependence

Other features of orbital-density-dependence
- #pause ODD functional means that we know $hat(H) lr(|phi_i chevron.r)$ for variational orbitals $\{lr(|phi_i chevron.r)\}$ but we don't know $hat(H)$ in general
- #pause Practically we can often use MLWFs
- #pause a natural generalisation in the direction of spectral functional theory (as discussed already by Andrea)@Ferretti2014

// ===========================================================================
== Screening

#set text(size: 0.9em)
- #pause In Hartree-Fock (the original "Koopmans' theorem"):
  $
    E_"ee"^"HF" = 1/2 sum_(i j) f_i f_j integral dif bold(r) dif bold(r)'
    (|psi_i (bold(r))|^2 |psi_j (bold(r)')|^2)/(bold(r) - bold(r)')
    - (psi_i^* (bold(r)) psi_j^* (bold(r)') psi_i (bold(r)') psi_j (bold(r)))/(bold(r) - bold(r)')
  $
  @Li2017
- #pause Account for screening post-hoc:
  $ (dif E)/(dif f_i) approx alpha_i (partial E)/(partial f_i) $
- #pause How to choose an appropriate value for $alpha_i$? Return to the original idea of Koopmans functionals:
  $ epsilon_i^"Koopmans" = E_i (N - 1) - E(N) $

// ===========================================================================
== Screening

#align(center)[
  #only("1", image("figures/alpha_calc/fig_alpha_calc_step_0.pdf", width: 45%))
  #only("2", image("figures/alpha_calc/fig_alpha_calc_step_1.pdf", width: 45%))
  #only("3,4", image("figures/alpha_calc/fig_alpha_calc_step_2.pdf", width: 45%))
  #only("5", image("figures/alpha_calc/fig_alpha_calc_step_3.pdf", width: 45%))
  #only("6-", image("figures/alpha_calc/fig_alpha_calc_step_4.pdf", width: 45%))
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

// ===========================================================================
== Screening

#v(5em)
$
  alpha^(n+1) = mark(alpha^n, color: #accent)
  (mark(E_i (N - 1), tag: #<enm1>, color: #accent) - mark(E(N), tag: #<en>, color: #accent) - lambda_(i i) (0))
  / (mark(lambda_(i i) (alpha^n), tag: #<lama>, color: #accent) - mark(lambda_(i i) (0), tag: #<lam0>, color: #accent))

  #annot(<enm1>, pos: top + left, dy: -3.5em)[total energy with electron \ removed from orbital $i$]
  #annot(<en>, pos: top, dy: -1.5em)[total energy \ of neutral system]
  #annot(<lama>, pos: bottom + left, dy: 1em)[expectation value \ of $hat(H)^"Koopmans"$]
  #annot(<lam0>, pos: bottom + right, dy: 3.5em)[expectation value \ of $hat(H)^"DFT"$]
$
#v(5em)

// ===========================================================================
== Screening via DFPT

How can we avoid explicit charged defect calculations in a supercell?

#uncover("2-")[
  Reformulate in terms of DFPT@Colonna2019...
  $ alpha_i = 1 + (lr(chevron.l v_"pert"^i | Delta^i n chevron.r)) / (lr(chevron.l n_i | v_"pert"^i chevron.r)) $
]

#uncover("3-")[
  ... in reciprocal space@Colonna2022
  $ alpha_(bold(0) i) = 1 + (sum_bold(q) lr(chevron.l v_("pert", bold(q))^(bold(0) i) | Delta_bold(q)^(bold(0) i) n chevron.r)) / (sum_bold(q) lr(chevron.l n_bold(q)^(bold(0) i) | v_("pert", bold(q))^(bold(0) i) chevron.r)) $

  (See Nicola Colonna's talk)
]

#uncover("4-")[N.B. even for the supercell, we can still reconstruct a band structure@DeGennaro2022]

// ===========================================================================
== To summarise...

When performing a Koopmans calculation, you must decide...
- #pause which flavour (KI, KIPZ)?
- #pause how are we treating the orbital-density-dependence?
  - how are we initialising our variational orbitals? (N.B. depends on the flavour!)
  - are we going to explicitly minimise the ODD?
- #pause how are we calculating the screening parameters? (finite differences, DFPT, ML...)

// ===========================================================================
== Koopmans functionals: results for molecules

#set text(size: 0.9em)
Ionisation potentials of 100 molecules cf. CCSD(T)
#align(center, image("figures/colonna_2019_gw100_ip.jpeg", height: 22%))

Ultraviolet photoemission spectra
#align(center, image("figures/fig_nguyen_prl_spectra.png", height: 30%))

@Colonna2018@Nguyen2015

// ===========================================================================
== Koopmans functionals: results for solids

#grid(columns: (35%, 60%), column-gutter: 1em, align: horizon + left,
  image("figures/fig_nguyen_prx_bandgaps.png", width: 100%),
  [
    #set text(size: 0.85em)
    Mean absolute error (eV) across prototypical semiconductors and insulators

    #v(1em)
    #table(
      columns: 6,
      align: (left, center, center, center, center, center),
      stroke: none,
      table.hline(),
      table.header([], [PBE], [G#sub[0]W#sub[0]], text(fill: accent, weight: "bold")[KI], text(fill: accent, weight: "bold")[KIPZ], [QSG$tilde("W")$]),
      table.hline(),
      [$E_"gap"$], [2.54], [0.56], text(fill: accent, weight: "bold")[0.27], text(fill: accent, weight: "bold")[0.22], [0.18],
      table.hline(stroke: 0.5pt),
      [IP], [1.09], [0.39], text(fill: accent, weight: "bold")[0.19], text(fill: accent, weight: "bold")[0.21], [0.49],
      table.hline(),
    )
  ],
)
@Nguyen2018

// ===========================================================================
== Koopmans functionals: results for solids

#align(center)[
  #set text(size: 0.5em)
  #let kihdr = text.with(fill: accent, weight: "bold")
  #table(
    columns: 7,
    align: (left, ) + (center, ) * 6,
    stroke: none,
    inset: 4pt,
    table.hline(), table.hline(y: 0),
    table.header(
      [], [PBE], [G#sub[0]W#sub[0]#super[a]], [scG$tilde("W")$#super[b]],
      kihdr[KI\@[PBE,MLWFs]], kihdr[KIPZ\@PBE], [exp#super[c]],
    ),
    table.hline(),
    [$E_g$], [0.49], [1.06], [1.14], kihdr[1.16], kihdr[1.15], [1.17],
    [$Gamma_(1v) arrow Gamma_(25'v)$], [11.97], [12.04], [], kihdr[11.97], kihdr[12.09], [12.5 ± 0.6],
    [$X_(1v) arrow Gamma_(25'v)$], [7.82], [], [], kihdr[7.82], kihdr[], [7.75],
    [$X_(4v) arrow Gamma_(25'v)$], [2.85], [2.99], [], kihdr[2.85], kihdr[2.86], [2.90],
    [$L_(2'v) arrow Gamma_(25'v)$], [9.63], [9.79], [], kihdr[9.63], kihdr[9.74], [9.3 ± 0.4],
    [$L_(1v) arrow Gamma_(25'v)$], [6.98], [7.18], [], kihdr[6.98], kihdr[7.04], [6.8 ± 0.2],
    [$L_(3'v) arrow Gamma_(25'v)$], [1.19], [1.27], [], kihdr[1.19], kihdr[], [1.2 ± 0.2],
    [$Gamma_(25'v) arrow Gamma_(15c)$], [2.48], [3.29], [], kihdr[3.17], kihdr[3.20], [3.35 ± 0.01],
    [$Gamma_(25'v) arrow Gamma_(2'c)$], [3.28], [4.02], [], kihdr[3.95], kihdr[3.95], [4.15 ± 0.05],
    [$Gamma_(25'v) arrow X_(1c)$], [0.62], [1.38], [], kihdr[1.28], kihdr[1.31], [1.13],
    [$Gamma_(25'v) arrow L_(1c)$], [1.45], [2.21], [], kihdr[2.12], kihdr[2.13], [2.04 ± 0.06],
    [$Gamma_(25'v) arrow L_(3c)$], [3.24], [4.18], [], kihdr[3.91], kihdr[3.94], [3.9 ± 0.1],
    table.hline(stroke: 0.5pt),
    [MSE], [0.35], [0.02], [], kihdr[0.01], kihdr[0.03], [],
    [MAE], [0.44], [0.21], [], kihdr[0.14], kihdr[0.17], [],
    table.hline(), table.hline(),
  )
  #set text(size: 0.7em)
  #align(left)[#super[a]@Shishkin2007 for $E_g$ and @Hybertsen1986 for the transitions; #super[b]@Shishkin2007a; #super[c]@Madelung2004]
]

// ===========================================================================
== Koopmans functionals: results for toy systems

For Hooke's atom (two electrons in a harmonic confining potential with Coulombic repulsion)

#v(1fr)
#grid(columns: (1fr, 1fr), column-gutter: 1em, align: center + horizon,
  image("figures/schubert_vxc.jpeg", width: 70%),
  uncover("2-", image("figures/schubert_vxc_integrated.jpeg", width: 70%)),
)
#v(1fr)
@Schubert2023

// ===========================================================================
== Koopmans functionals: caveats

- #pause restricted to systems with a non-zero band gap
- #pause empty state localization in the bulk limit
- #pause can potentially break the crystal point group symmetry@Su2020

// ===========================================================================
== The workflows

The general workflow:
- define/initialize a set of variational orbitals
- calculate the screening parameters $\{alpha_i\}$
- construct and diagonalize the Hamiltonian

// ===========================================================================
== The workflows

#uncover("1-")[
  (a) finite difference calculations using a supercell
  #align(center, image("figures/supercell_workflow.pdf", width: 100%))
]

#uncover("2-")[
  (b) DFPT in a primitive cell
  #align(center, image("figures/primitive_workflow.pdf", width: 65%))
]

// ===========================================================================
== How do I run these calculations?

Complicated workflows mean that...
- #pause lots of different codes that need to handshake
- #pause lots of scope for human error
- #pause reproducibility becomes difficult
- #pause expert knowledge required

#pause Our solution...

// ===========================================================================
== The #raw("koopmans") package

#align(center, image("figures/koopmans_grey_on_transparent.png", height: 18%))

#grid(columns: (55%, 45%), column-gutter: 1em, align: horizon + left,
  [
    #set text(size: 0.9em)
    - v1.0 released earlier this year@Linscott2023
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
    #link("https://koopmans-functionals.org")[`koopmans-functionals.org`]
    #v(0.5em)
    #image("figures/website_cropped.png", width: 95%)
  ],
)

// ===========================================================================
== koopmans: the input file

#grid(columns: (1fr, 1fr), column-gutter: 1em,
  listing-lines("scripts/si.json", 0, end: 20, lang: "json"),
  listing-lines("scripts/si.json", 20, lang: "json"),
)

// ===========================================================================
== koopmans is scriptable

#listing("scripts/si.py", lang: "python", size: 0.62em)

// ===========================================================================
== Generic structure

#set text(size: 0.9em)
#uncover("1-")[/ #raw("Workflow"):]
#pad(left: 1.5em)[
  #uncover("2-")[/ #raw("atoms"): an #raw("ASE") #raw("Atoms") object]
  #uncover("3-")[/ #raw("calculations"): a list of #raw("ASE") calculators]
  #uncover("4-")[/ #raw("kpoints"): a custom class containing $k$-point information]
  #uncover("5-")[/ #raw("pseudopotentials"): a dictionary of pseudopotentials]
]

#uncover("6-")[We will see examples in the hands-on!]

// ===========================================================================
== Take home messages

#align(center, grid(columns: 3, column-gutter: 1.5em, align: horizon,
  image("figures/colonna_2019_gw100_ip.jpeg", height: 18%),
  image("figures/fig_nguyen_prx_bandgaps.png", height: 18%),
  image("figures/supercell_workflow.pdf", height: 18%),
))

- Koopmans functionals are more complicated than a simple semi-local DFT calculation, because of...
  - orbital-density-dependence
  - screening parameters
- Koopmans functionals are implemented in #smallcaps[Quantum ESPRESSO]
- the complexity of the workflows are handled by the #raw("koopmans") package

// ===========================================================================
== Take home messages

#align(center + horizon, image("figures/jctc.png", width: 100%))

// ===========================================================================
== Acknowledgements

#align(center)[
  #set text(size: 0.8em)
  #grid(columns: 4, column-gutter: 1em, row-gutter: 0.5em, align: center,
    image("figures/nicola_marzari.jpg", height: 28%),
    image("figures/nicola_colonna2.png", height: 28%),
    image("figures/riccardo_degennaro.jpg", height: 28%),
    image("figures/yannick_schubert.jpg", height: 28%),
    [Nicola Marzari], [Nicola Colonna], [Riccardo De Gennaro], [Yannick Schubert],
  )

  #v(1em)
  #grid(columns: 2, column-gutter: 2em, align: horizon,
    image("logos/SNF_logo_standard_web_color_pos_e.png", height: 12%),
    image("figures/marvel_trimmed.png", height: 14%),
  )

  #v(1em)
  Want to find out more? Go to #link("https://koopmans-functionals.org")[`koopmans-functionals.org`]

  #v(0.5em)
  Follow #box(image("figures/Twitter_Bird.png", height: 0.8em)) #text(fill: twitter-blue)[\@ed_linscott] for updates | Slides available at #box(image("logos/github-favicon.png", height: 0.9em)) github/elinscott
]

// ===========================================================================
== References
#bibliography("references.bib")

// ===========================================================================
// Spare slides
// ===========================================================================
#focus-slide[Spare slides]

// ===========================================================================
== Koopmans functionals: off-diagonal occupancies

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

// ===========================================================================
== Learning the screening parameters

#align(center)[
  #grid(columns: 5, column-gutter: 1.2em, align: horizon + center,
    [#image("figures/orbital.emp.00191_cropped.png", width: 3cm) \ $rho_i (bold(r))$],
    connector[power spectrum \ decomposition],
    $ vec(x_0, x_1, x_2, dots.v) $,
    connector[ML model],
    $ alpha_i $,
  )
]
@Schubert2022

$
  c_(n l m, k = "orbital")^i & = integral dif bold(r) thin g_(n l) (r) Y_(l m) (theta, phi) rho^i (bold(r) - bold(R)^i) \
  p_(n_1 n_2 l, k_1 k_2)^i & = pi sqrt(8/(2 l + 1)) sum_m c_(n_1 l m, k_1)^(i *) c_(n_2 l m, k_2)^i
$

// ===========================================================================
// TODO: this spare slide's running header wrongly shows the *next* slide's
// title ("Resonance with other efforts"). It's a touying heading-location
// quirk tied to the "layout did not converge" warning, which destabilises
// page numbers; the slide content itself is correct. Investigate separately.
== Learning the screening parameters: results

#align(center)[
  #grid(columns: (1fr, 1fr), column-gutter: 1em, align: center + horizon,
    image("figures/CsSnI3_calc_vs_pred_Edward.png", height: 55%),
    image("figures/convergence_analysis_Edward.png", height: 55%),
  )

  #v(0.5em)
  loss of accuracy of the band gap of $tilde 0.02$ eV

  (cf. when calculating screening parameters _ab initio_)

  speedup of $70 times$
]
@Schubert2022

// ===========================================================================
== Resonance with other efforts

- Wannier transition-state method of Anisimov and Kozhevnikov @Anisimov2005
- Optimally tuned hybrid functionals of Kronik, Pasquarello, and others (refer back to Leeor's talk on Wednesday) @Kronik2012@Wing2021
- Ensemble DFT of Kronik and co-workers @Kraisler2013
- Koopmans-Wannier of Wang and co-workers @Ma2016
- Dielectric-dependent hybrid functionals of Galli and co-workers @Skone2016a
- LOSC functionals of Yang and co-workers @Li2018
