# Dados

`planck_ns_r.csv.gz` é um extrato das cadeias de Monte Carlo públicas da missão Planck (2018),
modelo ΛCDM + r. Cada linha é uma amostra da cadeia:

| coluna | significado |
|---|---|
| `conjunto` | combinação de dados: `Planck`, `Planck + BK15`, `Planck + BK15 + BAO` |
| `cadeia` | número da cadeia de origem (o Planck roda 4 cadeias independentes) |
| `peso` | quantas vezes a amostra "conta" |
| `ns` | índice espectral escalar |
| `r` | razão tensor-escalar em k = 0,05 Mpc⁻¹ |

## Como recriar este arquivo a partir dos dados originais

1. Baixe o pacote completo (11 GB) do Planck Legacy Archive ou do espelho da NASA/IRSA:
   <https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_CosmoParams_fullGrid_R3.01.zip>
2. Dentro do zip, em `COM_CosmoParams_fullGrid_R3.01/base_r/`, copie as pastas
   `plikHM_TTTEEE_lowl_lowE_lensing` e `plikHM_TTTEEE_lowl_lowE_BK15_lensing`
   para uma pasta `dados_planck/` na raiz do projeto.
3. Rode:

   ```bash
   python scripts/preparar_dados.py --origem dados_planck
   ```

**Crédito:** baseado em observações obtidas com o Planck (https://www.esa.int/Planck), uma missão científica
da ESA com instrumentos e contribuições financiados diretamente pelos Estados-membros da ESA, pela NASA e pelo Canadá.
