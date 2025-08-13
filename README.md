# Snakemake workflow for **Bulk-based hypothesis weighing to increase power in single-cell differential expression analysis**

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
  - **out:** simulated count matrix, `sim_dat+"00-raw/{sim}.rds"`,<br>
   simulated ground truth DE genes: `"sim_dat+"00-truth/{sim}.rds"`
  - two different simulation tools: `muscat`, `splatter`
