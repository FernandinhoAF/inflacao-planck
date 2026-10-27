"""
Leitura das cadeias de Monte Carlo (MCMC) públicas do satélite Planck.

O que é uma "cadeia": a colaboração Planck explorou o espaço de parâmetros
cosmológicos com um algoritmo de Monte Carlo. Cada linha do arquivo é uma
combinação de parâmetros compatível com os dados; quanto mais vezes uma
região aparece, mais provável ela é. O formato (CosmoMC) é:

    coluna 0      -> peso (quantas vezes a linha "conta")
    coluna 1      -> −log(verossimilhança)
    colunas 2...  -> parâmetros, na ordem do arquivo <raiz>.paramnames

Em vez de usar o número da coluna "na mão" (o que quebra fácil),
procuramos cada parâmetro pelo nome no arquivo .paramnames.
"""

from pathlib import Path

import numpy as np
import pandas as pd


def ler_nomes_dos_parametros(arquivo_paramnames: Path) -> list[str]:
    """Lê o .paramnames e devolve os nomes na ordem das colunas.

    Parâmetros derivados vêm marcados com '*' no fim do nome; removemos a marca.
    """
    nomes = []
    for linha in Path(arquivo_paramnames).read_text().splitlines():
        if linha.strip():
            nomes.append(linha.split()[0].rstrip("*"))
    return nomes


def burn_in_ja_removido(pasta: Path, raiz: str) -> bool:
    """O Planck indica no <raiz>.properties.ini se o 'burn-in' já foi removido.

    Burn-in é o trecho inicial da cadeia, antes de o algoritmo "esquentar";
    ele deve ser descartado. Se o arquivo não disser nada, assumimos que não.
    """
    props = Path(pasta) / f"{raiz}.properties.ini"
    if not props.exists():
        return False
    for linha in props.read_text().splitlines():
        chave, _, valor = linha.partition("=")
        if chave.strip().lower() == "burn_removed":
            return valor.strip().upper() in {"T", "TRUE", "1"}
    return False


def carregar_cadeias(
    pasta: str | Path,
    raiz: str,
    parametros: list[str],
    fracao_burn_in: float = 0.3,
) -> pd.DataFrame:
    """Junta todas as cadeias <raiz>_1.txt, <raiz>_2.txt, ... em uma tabela.

    Devolve as colunas: 'cadeia', 'peso' e os parâmetros pedidos.
    """
    pasta = Path(pasta)
    nomes = ler_nomes_dos_parametros(pasta / f"{raiz}.paramnames")
    faltando = [p for p in parametros if p not in nomes]
    if faltando:
        raise KeyError(f"Parâmetros ausentes em {raiz}: {faltando}")
    colunas = [0] + [2 + nomes.index(p) for p in parametros]

    descartar = 0.0 if burn_in_ja_removido(pasta, raiz) else fracao_burn_in

    arquivos = sorted(pasta.glob(f"{raiz}_[0-9]*.txt"), key=lambda a: int(a.stem.rsplit("_", 1)[1]))
    if not arquivos:
        raise FileNotFoundError(f"Nenhuma cadeia {raiz}_N.txt em {pasta}")

    partes = []
    for arq in arquivos:
        dados = np.loadtxt(arq, usecols=colunas, ndmin=2)
        inicio = int(len(dados) * descartar)
        df = pd.DataFrame(dados[inicio:], columns=["peso", *parametros])
        df.insert(0, "cadeia", int(arq.stem.rsplit("_", 1)[1]))
        partes.append(df)
    return pd.concat(partes, ignore_index=True)
