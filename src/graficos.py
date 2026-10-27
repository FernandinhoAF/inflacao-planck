"""Estilo visual e funções de desenho reutilizadas no notebook."""

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap

# Paleta categórica validada para daltonismo (ordem fixa: azul, laranja, verde-água)
CORES_CONJUNTOS = {
    "Planck": "#2a78d6",
    "Planck + BK15": "#eb6834",
    "Planck + BK15 + BAO": "#1baf7a",
}
TINTA = "#2b2b29"        # texto e marcas neutras
TINTA_SUAVE = "#6b6a63"  # rótulos secundários
GRADE = "#e6e5e0"

ROTULO_NS = r"índice espectral $n_s$"
ROTULO_R = r"razão tensor-escalar $r$"

# Ordem fixa da paleta categórica (usada quando há várias séries, ex.: valores de n)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

# Escala divergente para γ: violeta (γ < 0) → cinza neutro (γ = 0) → laranja (γ > 0)
CMAP_GAMMA = LinearSegmentedColormap.from_list("gamma", ["#4b2fb0", "#8a8880", "#c75a1d"])
CINZA_DADOS = "#a8a69e"  # contornos neutros quando a cor está reservada para γ


def aplicar_estilo():
    plt.rcParams.update({
        "axes.prop_cycle": cycler(color=PALETA),
        "figure.dpi": 110,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.edgecolor": TINTA_SUAVE,
        "axes.labelcolor": TINTA,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlecolor": TINTA,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRADE,
        "grid.linewidth": 0.8,
        "xtick.color": TINTA_SUAVE,
        "ytick.color": TINTA_SUAVE,
        "legend.frameon": False,
        "lines.linewidth": 2,
    })


def contornos(ax, densidade, cor, rotulo=None, ns_lim=(0.94, 1.0), r_lim=(0.0, 0.25), preencher=True):
    """Desenha as regiões de 68% (mais escura) e 95% (mais clara)."""
    NS, R, D = densidade.grade(ns_lim, r_lim)
    n68, n95 = densidade.nivel(0.68), densidade.nivel(0.95)
    if preencher:
        ax.contourf(NS, R, D, levels=[n95, n68], colors=[cor], alpha=0.25)
        ax.contourf(NS, R, D, levels=[n68, np.inf], colors=[cor], alpha=0.55)
    ax.contour(NS, R, D, levels=[n95, n68], colors=[cor], linewidths=1.2)
    if rotulo:
        ax.fill_between([], [], color=cor, alpha=0.5, label=rotulo)  # entrada da legenda
    ax.set_xlim(*ns_lim)
    ax.set_ylim(*r_lim)
    ax.set_xlabel(ROTULO_NS)
    ax.set_ylabel(ROTULO_R)
