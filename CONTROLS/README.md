# MURANUCI mock-community controls

`Controls.Rmd` compares the expected composition of Community, DNA and Gut mock-community controls with species-level Kraken2/Bracken results. It generates Supplementary Figure S1 and supporting CSV and Excel tables. It reads existing reports and does not run Kraken2 or Bracken.

## Data availability and required inputs

The expected compositions and the control-label mapping are defined in the Rmd. Original report files and generated outputs are not distributed in this repository. Cloning the repository alone is not sufficient to execute the analysis.

Paths below are relative to the repository root (`results`):

| Input or helper | Purpose |
|---|---|
| `TAX-TAXPROFILER/bracken/db2/*.report_bracken_species.txt` | Kraken2/Bracken reports containing the control samples |
| `TAX-TAXPROFILER/R/import_kraken.R` | Shared report importer |

The importer expects tab-separated reports without a header, with six columns: percentage, clade read count, directly assigned read count, taxonomic rank, NCBI taxonomy identifier and indented taxon name. Filenames must match the naming convention and control mapping used in `Controls.Rmd`. This analysis does not require the private patient-identifier mapping.

## Dependencies

Use the project environment described in the [repository README](../README.md#software-requirements) and restore its recorded packages with `renv::restore()`.

The script uses `tidyverse`, `dplyr`, `scales`, `DT` and `WriteXLS`. The shared importer requires `doParallel` and `foreach` and starts 20 worker processes. Excel export through `WriteXLS` requires Perl. Rendering requires `rmarkdown`, `knitr` and Pandoc.

## Run the analysis

Open `MURANUCI_results.Rproj` or start R in the repository root with the `renv` environment active. With the input reports available, run:

```r
rmarkdown::render("CONTROLS/Controls.Rmd")
```

Rendering creates `Controls.html` beside the Rmd. For interactive chunk execution, use `CONTROLS/` as the working directory so that the relative input and helper paths resolve correctly. Existing outputs with the same names are overwritten.

## Analysis scope

The script extracts species-level records, harmonises the taxon-name synonym handled in the code and joins observed results to the expected sample-by-species combinations. Missing observations are represented as zero reads and zero observed abundance. It calculates detection, differences in percentage points and observed-to-expected recovery ratios.

The comparison tables and heatmap include taxa with an expected abundance of at least 0.1%. Lower-abundance taxa remain in the `expected_composition` object but are excluded from these comparisons and exports. The heatmap uses a linear scale and displays expected and observed abundances for each control replicate.

## Generated outputs

All outputs are written to `CONTROLS/`:

| File | Contents |
|---|---|
| `controls_expected_observed_heatmap.png` | Supplementary Figure S1: expected versus observed abundance |
| `Supplementary_mock_controls.csv` | Detailed comparison, including detection status and recovery ratios |
| `controls_table_long.csv` | One row per control replicate and expected species |
| `controls_table_long.xlsx` | Excel version of the long-format table |
| `controls_table_wide.csv` | Expected and observed percentages in separate columns for each control replicate |
| `controls_table_wide.xlsx` | Excel version of the wide-format table |

The exported table field `WGS_percent` stores the observed percentage extracted from the Kraken2/Bracken report; it is not a separate isolate-WGS measurement.
