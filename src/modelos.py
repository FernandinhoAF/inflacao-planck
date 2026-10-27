"""
Previsões de modelos de inflação cósmica com um campo escalar.

O potencial estudado é o binomial do artigo (Antunes, Benetti, Bittencourt
e Franco, Physics of the Dark Universe 48, 101867, 2025):

    V(φ) / α = φ^n · (1 + γ · φ^(m − n))

com φ medido em unidades da massa de Planck. Com γ = 0 ele vira o
potencial monomial V ∝ φ^n.

Para cada potencial calculamos, na aproximação de rolamento lento
("slow-roll"), as duas grandezas que o satélite Planck mede:

    nₛ = 1 − 6ε + 2η      (índice espectral escalar)
    r  = 16ε              (razão tensor-escalar)

em que ε e η são os parâmetros de rolamento lento avaliados no instante
em que a escala de referência (pivô) cruzou o horizonte, N* e-folds
antes do fim da inflação.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


@dataclass(frozen=True)
class PotencialPolinomial:
    """V(φ)/α = φ^n (1 + γ φ^(m−n)).

    n     : expoente do termo dominante
    gamma : intensidade do segundo termo (γ = 0 é o caso monomial)
    m     : expoente do segundo termo (por padrão m = n + 1)
    """

    n: float
    gamma: float = 0.0
    m: float | None = None

    @property
    def expoente_m(self) -> float:
        return self.n + 1 if self.m is None else self.m

    # --- o potencial e suas derivadas -----------------------------------
    def V(self, phi):
        n, m, g = self.n, self.expoente_m, self.gamma
        return phi**n + g * phi**m

    def dV(self, phi):
        n, m, g = self.n, self.expoente_m, self.gamma
        return n * phi ** (n - 1) + g * m * phi ** (m - 1)

    def d2V(self, phi):
        n, m, g = self.n, self.expoente_m, self.gamma
        return n * (n - 1) * phi ** (n - 2) + g * m * (m - 1) * phi ** (m - 2)

    # --- parâmetros de rolamento lento -----------------------------------
    def epsilon(self, phi):
        """ε_V = ½ (V'/V)²  — mede a inclinação do potencial."""
        return 0.5 * (self.dV(phi) / self.V(phi)) ** 2

    def eta(self, phi):
        """η_V = V''/V  — mede a curvatura do potencial."""
        return self.d2V(phi) / self.V(phi)

    def phi_maximo(self, teto: float = 200.0) -> float:
        """Maior φ em que o potencial ainda é válido (V > 0 e V' > 0).

        Com γ < 0 o segundo termo acaba "virando" o potencial para baixo;
        a inflação só faz sentido antes desse ponto de máximo.
        """
        n, m, g = self.n, self.expoente_m, self.gamma
        if g >= 0:
            return teto
        # V' = 0  ⇔  φ^(m−n) = −n / (γ m)
        phi_topo = (-n / (g * m)) ** (1.0 / (m - n))
        return min(teto, phi_topo * (1 - 1e-9))


def fim_da_inflacao(pot: PotencialPolinomial) -> float:
    """φ_e: valor do campo quando ε_V = 1 (fim da inflação).

    Varremos φ de cima para baixo e pegamos a primeira vez que ε
    atinge 1; depois refinamos com o método de Brent.
    """
    phi_max = pot.phi_maximo()
    grade = np.geomspace(1e-3, phi_max, 4000)[::-1]
    eps = pot.epsilon(grade) - 1.0
    cruzou = np.nonzero(np.diff(np.sign(eps)) != 0)[0]
    if len(cruzou) == 0:
        return np.nan
    i = cruzou[0]
    return brentq(lambda p: pot.epsilon(p) - 1.0, grade[i + 1], grade[i])


def numero_de_efolds(pot: PotencialPolinomial, phi_e: float, phi: float) -> float:
    """N(φ) = ∫_{φ_e}^{φ} V / V' dφ  (equação 2.10 do artigo, em M_pl = 1)."""
    valor, _ = quad(lambda p: pot.V(p) / pot.dV(p), phi_e, phi, limit=200)
    return valor


def phi_no_pivo(pot: PotencialPolinomial, N_estrela: float) -> float:
    """φ*: valor do campo N* e-folds antes do fim da inflação."""
    phi_e = fim_da_inflacao(pot)
    if not np.isfinite(phi_e):
        return np.nan
    phi_max = pot.phi_maximo()
    if numero_de_efolds(pot, phi_e, phi_max) < N_estrela:
        return np.nan  # o potencial não comporta N* e-folds de inflação
    return brentq(
        lambda p: numero_de_efolds(pot, phi_e, p) - N_estrela,
        phi_e,
        phi_max,
        xtol=1e-10,
    )


def previsao(n: float, N_estrela: float = 55.0, gamma: float = 0.0, m: float | None = None) -> dict:
    """Calcula (nₛ, r) previstos para um potencial, resolvendo tudo numericamente."""
    pot = PotencialPolinomial(n=n, gamma=gamma, m=m)
    phi = phi_no_pivo(pot, N_estrela)
    if not np.isfinite(phi):
        return {"n": n, "m": pot.expoente_m, "gamma": gamma, "N": N_estrela,
                "phi": np.nan, "ns": np.nan, "r": np.nan}
    eps, eta = pot.epsilon(phi), pot.eta(phi)
    return {
        "n": n,
        "m": pot.expoente_m,
        "gamma": gamma,
        "N": N_estrela,
        "phi": phi,
        "ns": 1 - 6 * eps + 2 * eta,
        "r": 16 * eps,
    }


# ---------------------------------------------------------------------------
# Fórmulas analíticas do artigo — usadas para conferir o cálculo numérico
# ---------------------------------------------------------------------------

def monomial_analitico(n: float, N_estrela: float) -> tuple[float, float]:
    """Equação (3.4) do artigo, com φ*² = (n² + 4nN*)/2."""
    nsm1 = -2 * (n + 2) / (n + 4 * N_estrela)
    r = 16 * n / (n + 4 * N_estrela)
    return 1 + nsm1, r


def binomial_primeira_ordem(n: float, N_estrela: float, gamma: float, caso: str) -> tuple[float, float]:
    """Expansões em primeira ordem em γ do artigo.

    caso = "m=n+1": equações (4.15) e (4.16)
    caso = "m=n+2": equações (4.28) e (4.29)
    """
    A = n + 4 * N_estrela
    ns0, r0 = monomial_analitico(n, N_estrela)
    if caso == "m=n+1":
        dns = 4 * np.sqrt(2) * ((2 + n) * n**1.5 + (1 - n) * A**1.5) / (3 * np.sqrt(n) * A**2)
        dr = -32 * np.sqrt(2) * np.sqrt(n) * (n**1.5 - A**1.5) / (3 * A**2)
    elif caso == "m=n+2":
        dns = 6 - 3 * n - n * (2 + n) * (2 * np.sqrt(2) + n) / A**2
        dr = 8 * n * (3 + n * (2 * np.sqrt(2) + n) / A**2)
    else:
        raise ValueError("caso deve ser 'm=n+1' ou 'm=n+2'")
    return ns0 + dns * gamma, r0 + dr * gamma
