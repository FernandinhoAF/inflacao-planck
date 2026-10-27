"""
Ferramentas estatísticas para amostras com peso (como as cadeias do Planck).

Três ideias principais:

1. Estatísticas ponderadas: cada amostra tem um peso, então média,
   desvio e quantis precisam levar o peso em conta.

2. Estimativa de densidade por kernel (KDE): transforma a nuvem de
   pontos em uma "superfície de probabilidade" suave. Como r ≥ 0
   por definição, usamos o método da reflexão para não perder
   probabilidade na borda r = 0.

3. Regiões de credibilidade (68% e 95%): a menor região que contém
   68% (ou 95%) da probabilidade. A partir delas, dizemos se a
   previsão de um modelo é compatível com os dados.
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


# ---------------------------------------------------------------------------
# 1. Estatísticas ponderadas
# ---------------------------------------------------------------------------

def media_ponderada(x, pesos) -> float:
    return float(np.average(x, weights=pesos))


def desvio_ponderado(x, pesos) -> float:
    m = np.average(x, weights=pesos)
    return float(np.sqrt(np.average((x - m) ** 2, weights=pesos)))


def quantil_ponderado(x, pesos, q: float) -> float:
    """Valor abaixo do qual está uma fração q do peso total."""
    ordem = np.argsort(x)
    x, pesos = np.asarray(x)[ordem], np.asarray(pesos)[ordem]
    acumulado = np.cumsum(pesos) - 0.5 * pesos
    return float(np.interp(q, acumulado / pesos.sum(), x))


def tamanho_efetivo(pesos) -> float:
    """Quantas amostras independentes, com peso igual, a cadeia 'vale'."""
    pesos = np.asarray(pesos, dtype=float)
    return float(pesos.sum() ** 2 / (pesos**2).sum())


def resumo(df: pd.DataFrame, colunas: list[str]) -> pd.DataFrame:
    """Tabela com média, desvio e intervalo de 95% de cada parâmetro."""
    linhas = []
    for c in colunas:
        linhas.append({
            "parâmetro": c,
            "média": media_ponderada(df[c], df["peso"]),
            "desvio": desvio_ponderado(df[c], df["peso"]),
            "2,5%": quantil_ponderado(df[c], df["peso"], 0.025),
            "97,5%": quantil_ponderado(df[c], df["peso"], 0.975),
        })
    return pd.DataFrame(linhas).set_index("parâmetro")


# ---------------------------------------------------------------------------
# 2. Densidade de probabilidade no plano (nₛ, r)
# ---------------------------------------------------------------------------

class DensidadeNsR:
    """KDE ponderada de (nₛ, r) com reflexão na borda r = 0.

    Por que reflexão? A KDE "espalha" cada ponto com uma gaussiana. Perto
    de r = 0 metade dessa gaussiana cairia em r < 0, que é impossível,
    e a densidade ficaria artificialmente baixa justamente onde moram os
    modelos com r pequeno. Espelhando os pontos em r → −r e somando as
    duas metades, a probabilidade "devolvida" volta para r ≥ 0.
    """

    def __init__(self, ns, r, pesos, max_amostras: int = 20_000, semente: int = 42,
                 fator_banda: float | None = None):
        ns, r, pesos = map(lambda a: np.asarray(a, dtype=float), (ns, r, pesos))
        # Reamostragem proporcional ao peso: deixa o cálculo rápido sem
        # distorcer a distribuição (amostras pesadas aparecem mais vezes).
        rng = np.random.default_rng(semente)
        n_amostras = min(max_amostras, len(ns))
        idx = rng.choice(len(ns), size=n_amostras, replace=True, p=pesos / pesos.sum())
        self.ns, self.r = ns[idx], r[idx]
        self.kde = gaussian_kde(np.vstack([self.ns, self.r]), bw_method=fator_banda)
        self._niveis = None
        self._grades = {}

    def __call__(self, ns, r):
        """Densidade de probabilidade em pontos (ns, r) com r ≥ 0."""
        ns, r = np.broadcast_arrays(np.asarray(ns, float), np.asarray(r, float))
        pontos = np.vstack([ns.ravel(), r.ravel()])
        espelho = np.vstack([ns.ravel(), -r.ravel()])
        dens = self.kde(pontos) + self.kde(espelho)
        dens[r.ravel() < 0] = 0.0
        return dens.reshape(ns.shape)

    def grade(self, ns_lim=(0.94, 1.0), r_lim=(0.0, 0.25), pontos=120):
        """Avalia a densidade numa grade regular (para desenhar contornos).

        O resultado fica guardado: a mesma grade não é recalculada duas vezes.
        """
        chave = (tuple(ns_lim), tuple(r_lim), pontos)
        if chave not in self._grades:
            NS, R = np.meshgrid(np.linspace(*ns_lim, pontos), np.linspace(*r_lim, pontos))
            self._grades[chave] = (NS, R, self(NS, R))
        return self._grades[chave]

    # ---------------------------------------------------------------------
    # 3. Regiões de credibilidade
    # ---------------------------------------------------------------------
    def _calibrar(self):
        """Calibra os níveis de densidade usando as próprias amostras.

        Truque padrão: para cada amostra calculamos a densidade nela. O
        nível que deixa 68% das amostras "acima" delimita a região de 68%.
        Isso dispensa integrar numericamente numa grade.
        """
        if self._niveis is None:
            k = min(len(self.ns), 8_000)  # 8 mil amostras já estimam bem os níveis
            self._niveis = np.sort(self(self.ns[:k], np.abs(self.r[:k])))[::-1]
        return self._niveis

    def nivel(self, massa: float) -> float:
        """Densidade que delimita a região que contém `massa` da probabilidade."""
        d = self._calibrar()
        return float(d[min(int(massa * len(d)), len(d) - 1)])

    def credibilidade(self, ns, r):
        """Menor região de credibilidade que contém o ponto (em %).

        Ex.: 40 → o ponto está dentro da região de 68%;
             90 → fora da de 68%, mas dentro da de 95%;
             99 → fora da região de 95% (modelo em tensão com os dados).
        """
        d = self._calibrar()
        dens_ponto = np.atleast_1d(self(ns, r))
        # fração das amostras com densidade MAIOR que a do ponto
        frac = 1.0 - np.searchsorted(d[::-1], dens_ponto, side="right") / len(d)
        return 100 * frac if frac.size > 1 else float(100 * frac[0])


def classificar(credibilidade_pct) -> str:
    if credibilidade_pct <= 68:
        return "compatível (dentro de 68%)"
    if credibilidade_pct <= 95:
        return "compatível (dentro de 95%)"
    return "em tensão (fora de 95%)"
