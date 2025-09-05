# Snakemake workflow for [**Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis**](https://www.biorxiv.org/content/10.1101/2025.04.15.648932v1.full)


### setup

- workflow was implemented and last executed successfully with<br>
  **R v4.5.1 with Bioc 3.20, and Python v3.12.11 with Snakemake v9.5.1**
- R version and library have to be specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `.Rprofile` is used for handling and printing command line arguments
- `logs/` capture `.Rout` files from `R CMD BATCH` executions
- `data/` contains any synthetic and real data
- intermediate results are generated in `outs/` 
- visualizations are generated in `plts/`

### workflow

- `<x>` denotes a wildcard, namely: `bin`, `cor`rection, `loc`al or global,  
  `sim`ulation, `sta`tistic

- `00-get_sim.R` 
  - **out:** simulated count matrix, `data/sim/00-raw/{sim}.rds"`,<br>
   simulated ground truth DE genes: `"data/sim/00-truth/{sim}.rds"`

- `01-bulkDEA.R` 
  - **in:** `00-raw/{sim}.rds`
  - **out:** bulk-level differential expression analysis `outs/bulkDEA/{sim}.rds`
  - `edgeR`-based DEA

- `02-pbDEA.R`
  - **in:** `00-raw/{sim}.rds`,
  - **out:** pseudobulk differential state analysis `outs/pbDEA/{sim},{size}.rds`
  - perform pseudobulk-level differential state analysis
  - source different number of random combination and sample size of "code/02-pbDEA-{size}.R"


- `03-bbhw.R`
  - **in:**: `00-raw/{sim}.rds`, `outs/bulkDEA/{sim}.rds`, `outs/pbDEA/{sim},{size}.rds`
  - **out:**:  `outs/bbhw/{sim},{size},{bin},{cor},{loc}.rds`
  - perform bulk-based hypothesis weighing with different `bin` methods and `cor`rection methods with global or `loc`al coorrection


- `04-sta.R`
  - **in:**: `bbhw/{sim},{size},{bin},{cor},{loc}.rds`, `data/sim/00-truth/{sim}.rds`
  - **out:**: `sta/{sta},{sim},{size},{bin},{cor},{loc}.rds`
  - compute precision/recall/F1 scores
  - source method from on of `04-sta-<sta>.R`, but in this case we computed all scores in one function `04-sta-metrics.R`


