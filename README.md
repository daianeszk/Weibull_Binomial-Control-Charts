# Weibull Binomial Control Charts: Versões unilateral e bilateral

Este repositório disponibiliza os scripts em **R** utilizados em minha pesquisa do mestrado em Engenharia de Produção para desenhar limites de um **np-chart** (contagem por atributos) a partir de um alvo de falso alarme (ARL0) e mapear limites para **pontos de alerta em tempo delineados por Weibull**.  
Também inclui a avaliação de desempenho sob cenários fora de controle (ARL1)

Scripts incluídos:

- `UnilateralNP_Weibull_MeanDecrease.R` — desenho **unilateral** (sinaliza quando `N ≥ Y`)
- `BilateralNP_Weibull_MeanDecrease.R` — desenho **bilateral** (LCL e UCL, com `alpha_total/2` em cada cauda)

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
