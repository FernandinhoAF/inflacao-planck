"""
Extrai das cadeias brutas do Planck 2018 só o que a análise usa.

O pacote completo do Planck tem 11 GB. Este script lê as cadeias do
modelo ΛCDM + r e guarda apenas o peso, nₛ e r de cada amostra num
arquivo compacto (dados/planck_ns_r.csv.gz). Assim o repositório fica
leve e qualquer pessoa consegue rodar o notebook.

Uso:
    python scripts/preparar_dados.py --origem caminho/para/dados_planck

onde `dados_planck` contém as pastas copiadas de
COM_CosmoParams_fullGrid_R3.01/base_r/ (veja dados/README.md).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

RAIZ_PROJETO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ_PROJETO / "src"))

from cadeias_planck import carregar_cadeias, ler_nomes_dos_parametros  # noqa: E402

# nome amigável -> (subpasta, raiz dos arquivos)
CONJUNTOS = {
    "Planck": (
        "plikHM_TTTEEE_lowl_lowE_lensing",
        "base_r_plikHM_TTTEEE_lowl_lowE_lensing",
    ),
    "Planck + BK15": (
        "plikHM_TTTEEE_lowl_lowE_BK15_lensing",
        "base_r_plikHM_TTTEEE_lowl_lowE_BK15_lensing",
    ),
    "Planck + BK15 + BAO": (
        "plikHM_TTTEEE_lowl_lowE_BK15_lensing",
        "base_r_plikHM_TTTEEE_lowl_lowE_BK15_lensing_post_BAO",
    ),
}

# r: razão tensor-escalar na escala k = 0,05 Mpc⁻¹ (a mesma do N* = 55 do artigo)
PARAMETROS_DESEJADOS = ["ns", "r"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--origem", required=True, type=Path, help="pasta com as subpastas do base_r")
    parser.add_argument("--destino", type=Path, default=RAIZ_PROJETO / "dados" / "planck_ns_r.csv.gz")
    args = parser.parse_args()

    tabelas = []
    for nome, (subpasta, raiz) in CONJUNTOS.items():
        pasta = args.origem / subpasta
        if not (pasta / f"{raiz}.paramnames").exists():
            print(f"[aviso] {nome}: não encontrei {raiz} em {pasta} — pulando")
            continue
        disponiveis = ler_nomes_dos_parametros(pasta / f"{raiz}.paramnames")
        parametros = [p for p in PARAMETROS_DESEJADOS if p in disponiveis]
        df = carregar_cadeias(pasta, raiz, parametros)
        df.insert(0, "conjunto", nome)
        tabelas.append(df)
        print(f"{nome:22s} {len(df):>8,d} amostras  ({df['cadeia'].nunique()} cadeias)  colunas: {parametros}")

    if not tabelas:
        sys.exit("Nenhum conjunto encontrado. Confira o caminho passado em --origem.")

    tudo = pd.concat(tabelas, ignore_index=True)
    args.destino.parent.mkdir(parents=True, exist_ok=True)
    tudo.to_csv(args.destino, index=False, float_format="%.6g")
    tamanho_mb = args.destino.stat().st_size / 1e6
    print(f"\nSalvo em {args.destino} ({tamanho_mb:.1f} MB)")


if __name__ == "__main__":
    main()
