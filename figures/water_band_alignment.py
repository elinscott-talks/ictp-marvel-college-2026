#!/usr/bin/env python3
"""Band-edge alignment diagram for TiO2 surfaces (KI results only).

A polished re-make of ``Water_BA_full_with_error`` in the layout sketched for
the ICTP--MARVEL 2026 talk:

  * two environment groups -- "vacuum" and "solvated";
  * each group contains the three surfaces (anatase, rutile 1, rutile 2);
  * the valence-band edge (VBM) is a bar rising from the bottom of the axes;
  * the conduction-band edge (CBM) is a bar hanging from the top (vacuum level);
  * the empty space between the two bars is the band gap;
  * error bars (T-caps) sit on each band edge;
  * dashed lines mark the important (water-redox) potentials;
  * solid red lines mark the experimental band edges.

Energies are quoted relative to the vacuum level (eV), so VBM = -IP and
CBM = -EA and both numbers are negative.

To use: edit ``DATA`` / ``POTENTIALS`` / the colours below, then run

    python scripts/water_band_alignment.py

which writes ``figures/Water_BA_KI.pdf`` and ``figures/Water_BA_KI.svg``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
from matplotlib import rcParams


# --------------------------------------------------------------------------- #
#  Data structures -- edit these                                              #
# --------------------------------------------------------------------------- #
@dataclass
class Edge:
    """A single band edge relative to the vacuum level (eV) and its error."""

    value: float
    err: float = 0.0


@dataclass
class Surface:
    """The valence- and conduction-band edges of one surface."""

    vbm: Edge                       # valence-band maximum   (= -IP)
    cbm: Edge | None = None         # conduction-band minimum (= -EA); if None,
                                    # derived to preserve the vacuum band gap
    experiment: list[float] = field(default_factory=list)  # exp. VBM(s) -> red

    @property
    def gap(self) -> float:
        return self.cbm.value - self.vbm.value


# Order of the surfaces within each group, and the colour of each one.
SURFACES = ["anatase 101", "rutile 110 p", "rutile 110 c"]
# matplotlib's default qualitative palette (tab10): 0=blue, 1=orange, 2=green
_TAB10 = plt.get_cmap("tab10").colors
COLOURS = {
    "anatase 101":  _TAB10[1],   # orange
    "rutile 110 p": _TAB10[0],   # blue
    "rutile 110 c": _TAB10[2],   # green
}

# DATA[environment][surface] = Surface(...)
#   Vacuum has no error bars (err = 0). The solvated CBM is left as None and is
#   reconstructed below by assuming the band gap is unchanged on solvation.
DATA: dict[str, dict[str, Surface]] = {
    "vacuum": {
        "anatase 101":  Surface(vbm=Edge(-8.59), cbm=Edge(-4.08)),
        "rutile 110 p": Surface(vbm=Edge(-8.38), cbm=Edge(-4.34)),
        "rutile 110 c": Surface(vbm=Edge(-8.38), cbm=Edge(-4.34)),
    },
    "solvated": {
        "anatase 101":  Surface(vbm=Edge(-7.5, 0.10), experiment=[-7.50, -7.64]),
        "rutile 110 p": Surface(vbm=Edge(-8.0, 0.10), experiment=[-7.85]),
        "rutile 110 c": Surface(vbm=Edge(-7.9, 0.30), experiment=[-7.85]),
    },
}

# Reconstruct any missing CBM by holding the band gap fixed between
# environments: CBM = VBM + (vacuum gap), with the VBM's uncertainty.
for _surf in SURFACES:
    _gap = DATA["vacuum"][_surf].gap
    for _env in DATA:
        _s = DATA[_env][_surf]
        if _s.cbm is None:
            _s.cbm = Edge(_s.vbm.value + _gap, _s.vbm.err)

# Important reference potentials -> dashed lines: (value vs vacuum, label).
# Plain (sans-serif) Unicode text -- no mathtext/LaTeX.
POTENTIALS = [
    (-4.44, "H⁺/H₂"),     # H⁺/H₂
    (-5.67, "O₂/H₂O"),    # O₂/H₂O
]


# --------------------------------------------------------------------------- #
#  Plot styling -- tweak to taste                                             #
# --------------------------------------------------------------------------- #
ENVIRONMENTS = list(DATA.keys())
Y_LIMITS = (-12.5, 0.0)           # energy axis range (eV)
BAR_WIDTH = 0.82                  # width of each surface bar (group units)
GROUP_GAP = 1.4                   # blank space between the two groups
CB_ALPHA = 0.55                   # transparency of the conduction-band bars
EDGE_COLOUR = "#33373B"
EXP_COLOUR = "#D62728"
DASH_COLOUR = "#404040"           # dark grey
OUT_STEM = "figures/Water_BA_KI"  # written as .pdf and .svg


def _x_position(env_index: int, surf_index: int) -> float:
    """x-coordinate of a given surface bar."""
    return env_index * (len(SURFACES) + GROUP_GAP) + surf_index


def _error_bar(ax, x, edge, bar_side, colour, *, cap=0.11, lw=1.6):
    """Draw a band-edge error bar in two tones: white where it overlaps the
    filled bar, and the bar's own colour where it sticks out past the edge.

    ``bar_side`` is "below" for a valence edge (bar fills below it) or "above"
    for a conduction edge (bar fills above it).
    """
    e, err = edge.value, edge.err
    if bar_side == "below":            # valence bar fills below the edge
        inside_end, outside_end = e - err, e + err
    else:                              # conduction bar fills above the edge
        inside_end, outside_end = e + err, e - err
    for end, c in ((inside_end, "white"), (outside_end, colour)):
        ax.plot([x, x], [e, end], color=c, lw=lw, zorder=4,
                solid_capstyle="butt")
        ax.plot([x - cap, x + cap], [end, end], color=c, lw=lw, zorder=4,
                solid_capstyle="butt")


def make_figure(*, show_vacuum=True, show_solvated=True,
                show_experiment=True) -> plt.Figure:
    rcParams.update({
        "font.size": 12,
        "font.family": "sans-serif",
        "svg.fonttype": "none",   # keep text as text in the SVG
    })

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ymin, ymax = Y_LIMITS

    x_left = _x_position(0, 0) - BAR_WIDTH
    x_right = _x_position(len(ENVIRONMENTS) - 1, len(SURFACES) - 1) + BAR_WIDTH

    # --- dashed reference potentials -------------------------------------- #
    for value, label in POTENTIALS:
        ax.axhline(value, ls=(0, (6, 4)), lw=1.2, color=DASH_COLOUR, zorder=1)
        ax.text(x_left, value - 0.1, label, color=DASH_COLOUR,
                va="top", ha="left", fontsize=10)

    # --- the bars (drawn per environment, optionally hidden) -------------- #
    show = {"vacuum": show_vacuum, "solvated": show_solvated}
    group_centres: dict[str, float] = {}
    for ei, env in enumerate(ENVIRONMENTS):
        xs = [_x_position(ei, si) for si in range(len(SURFACES))]
        group_centres[env] = sum(xs) / len(xs)
        if not show.get(env, True):
            continue
        for si, surf in enumerate(SURFACES):
            x = xs[si]
            data = DATA[env][surf]
            colour = COLOURS[surf]
            half = BAR_WIDTH / 2.0

            # valence-band bar: bottom of axes up to the VBM
            ax.bar(x, data.vbm.value - ymin, bottom=ymin, width=BAR_WIDTH,
                   color=colour, linewidth=0, zorder=2)
            # conduction-band bar: vacuum level down to the CBM
            ax.bar(x, ymax - data.cbm.value, bottom=data.cbm.value,
                   width=BAR_WIDTH, color=colour, linewidth=0,
                   alpha=CB_ALPHA, zorder=2)

            # error bars: white where they overlap the bar, bar-colour beyond
            for edge, side in ((data.vbm, "below"), (data.cbm, "above")):
                if edge.err > 0:
                    _error_bar(ax, x, edge, side, colour)

            # experimental band edge(s) -> solid red line(s)
            if show_experiment:
                for exp in data.experiment:
                    ax.hlines(exp, x - half - 0.04, x + half + 0.04,
                              color=EXP_COLOUR, lw=2.4, zorder=5)

            # surface name written vertically inside the valence bar
            ax.text(x, ymin + 0.3, surf, rotation=90, va="bottom", ha="center",
                    fontsize=11, color="white", fontweight="bold", zorder=3)

    # --- axes cosmetics --------------------------------------------------- #
    ax.set_ylim(ymin, ymax)
    # x_left / x_right are symmetric about the bars, so equal margins give the
    # same amount of space to the left of vacuum and the right of solvated
    ax.set_xlim(x_left - 0.1, x_right + 0.1)
    ax.set_ylabel("Energy relative to the vacuum (eV)", fontstyle="italic")
    ax.set_xticks(list(group_centres.values()))
    ax.set_xticklabels(list(group_centres.keys()), fontsize=15)
    ax.tick_params(axis="x", length=0, pad=12)
    ax.tick_params(axis="y", length=4)

    # minimal legend for the experimental marker (only once it is shown)
    if show_experiment:
        ax.plot([], [], color=EXP_COLOUR, lw=2.4, label="experiment")
        ax.legend(loc="upper right", frameon=False, fontsize=10)

    fig.tight_layout()
    return fig


def main() -> None:
    # cumulative reveal steps for the talk; every version shares an identical
    # canvas (no per-figure tight bbox) so they overlay cleanly on a slide.
    steps = [
        ("1_potentials", dict(show_vacuum=False, show_solvated=False,
                              show_experiment=False)),
        ("2_vacuum",     dict(show_vacuum=True,  show_solvated=False,
                              show_experiment=False)),
        ("3_solvated",   dict(show_vacuum=True,  show_solvated=True,
                              show_experiment=False)),
        ("4_experiment", dict(show_vacuum=True,  show_solvated=True,
                              show_experiment=True)),
    ]
    for suffix, flags in steps:
        fig = make_figure(**flags)
        for ext in ("pdf", "svg"):
            path = f"{OUT_STEM}_{suffix}.{ext}"
            fig.savefig(path)
            print(f"wrote {path}")
        plt.close(fig)

    # also keep a canonical full figure under the plain stem
    fig = make_figure()
    for ext in ("pdf", "svg"):
        fig.savefig(f"{OUT_STEM}.{ext}")
        print(f"wrote {OUT_STEM}.{ext}")
    plt.close(fig)

    print("\nband gaps (CBM - VBM):")
    for env in ENVIRONMENTS:
        for surf in SURFACES:
            print(f"  {env:8s}  {surf:9s}  {DATA[env][surf].gap:5.2f} eV")


if __name__ == "__main__":
    main()
