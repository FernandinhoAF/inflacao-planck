"""Testes automáticos: rodam com `pytest` na raiz do projeto."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cadeias_planck import carregar_cadeias  # noqa: E402
from estatistica import DensidadeNsR, media_ponderada, quantil_ponderado  # noqa: E402
from modelos import binomial_primeira_ordem, monomial_analitico, previsao  # noqa: E402


# --- modelos ----------------------------------------------------------------

@pytest.mark.parametrize("n", [0.5, 2 / 3, 1, 2, 3])
@pytest.mark.parametrize("N", [50, 55, 60])
def test_monomial_numerico_igual_ao_analitico(n, N):
    p = previsao(n, N)
    ns, r = monomial_analitico(n, N)
    assert p["ns"] == pytest.approx(ns, abs=1e-8)
    assert p["r"] == pytest.approx(r, abs=1e-8)


def test_reproduz_valores_do_artigo():
    # n = 1/2, N* = 50 — ponto usado na figura 3 do artigo
    p = previsao(0.5, 50)
    assert p["ns"] == pytest.approx(0.975062, abs=1e-6)
    assert p["r"] == pytest.approx(0.0399002, abs=1e-6)


@pytest.mark.parametrize("caso,m_menos_n,gamma", [("m=n+1", 1, 1e-3), ("m=n+2", 2, 1e-5)])
@pytest.mark.parametrize("n", [0.5, 1, 2])
def test_binomial_concorda_com_primeira_ordem_para_gamma_pequeno(caso, m_menos_n, gamma, n):
    for g in (gamma, -gamma):
        p = previsao(n, 55, gamma=g, m=n + m_menos_n)
        ns, r = binomial_primeira_ordem(n, 55, g, caso)
        assert p["ns"] == pytest.approx(ns, abs=2e-5)
        assert p["r"] == pytest.approx(r, rel=2e-3)


# --- estatística ------------------------------------------------------------

def test_estatisticas_ponderadas():
    x = np.array([1.0, 2.0, 3.0])
    w = np.array([1.0, 0.0, 1.0])
    assert media_ponderada(x, w) == pytest.approx(2.0)
    assert quantil_ponderado(np.arange(1000.0), np.ones(1000), 0.5) == pytest.approx(499.5, abs=1)


def test_credibilidade_de_uma_gaussiana():
    # Gaussiana 2D longe da borda r = 0: o ponto a 1σ do centro fica na
    # fronteira da região de 39,3% (= 1 − e^{−1/2}).
    rng = np.random.default_rng(0)
    ns = rng.normal(0.965, 0.004, 40_000)
    r = rng.normal(0.5, 0.02, 40_000)
    dens = DensidadeNsR(ns, r, np.ones_like(ns))
    assert dens.credibilidade(0.965, 0.5) < 5
    assert dens.credibilidade(0.965 + 0.004, 0.5) == pytest.approx(39.3, abs=4)
    assert dens.credibilidade(0.965 + 0.004 * 3, 0.5) > 95


# --- leitura das cadeias ------------------------------------------------------

def test_leitura_por_nome_e_burn_in(tmp_path):
    (tmp_path / "teste.paramnames").write_text("a  A\nns  n_s\nr*  r\n")
    dados = np.column_stack([np.ones(10), np.zeros(10), np.arange(10), np.arange(10) + 0.5, np.arange(10) * 2])
    np.savetxt(tmp_path / "teste_1.txt", dados)
    df = carregar_cadeias(tmp_path, "teste", ["r", "ns"], fracao_burn_in=0.3)
    assert len(df) == 7  # 30% descartado
    assert df["ns"].iloc[0] == 3.5 and df["r"].iloc[0] == 6.0

    (tmp_path / "teste.properties.ini").write_text("burn_removed = T\n")
    assert len(carregar_cadeias(tmp_path, "teste", ["ns"])) == 10
