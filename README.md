# Southern Atlantic cloud-forest demography

R code and analysis-ready data for the manuscript **“Hierarchically structured tree demographic variation among southern Atlantic cloud-forest areas.”**

The workflow analyzes three tree inventories from 80 permanent plots in four cloud-forest areas.

## Run the analysis

Install the required packages once:

```r
install.packages(c(
  "dplyr", "readr", "stringr", "FactoMineR", "lme4", "lmerTest",
  "performance", "broom.mixed", "scales", "ggplot2", "patchwork", "ragg"
))
```

From the repository root, run:

```r
source("R/run_all.R", encoding = "UTF-8")
```

The workflow writes tables and figures locally to `results/`.

Contact: Pedro Higuchi — `higuchip@gmail.com`
