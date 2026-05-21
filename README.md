[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18571398.svg)](https://doi.org/10.5281/zenodo.18571398)


# One Sided Attribute Chart to Monitor Weibull Mean:

Este repositório disponibiliza os scripts em **R** utilizados em minha pesquisa do mestrado em Engenharia de Produção para desenhar limites de um **np-chart** (contagem por atributos) a partir de um alvo de falso alarme (ARL0) e mapear limites para **pontos de alerta em tempo delineados por Weibull**.  
Também inclui a avaliação de desempenho sob cenários fora de controle (ARL1)

---

## Requisitos

- R (recomendado: 4.2+)
- Pacotes:
  - `dplyr`
  - `tidyr`
  - `openxlsx`

Instalação (se necessário):
```r
install.packages(c("dplyr", "tidyr", "openxlsx"))
