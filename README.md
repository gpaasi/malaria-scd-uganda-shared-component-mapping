# Spatial co-clustering of malaria and sickle cell disease burden across health facility catchments in Uganda using Bayesian shared-component disease mapping

DOI Release License: MIT  
R >= 4.3  
renv  
Issues  
Pull requests  

## Summary

This repository implements the analytic workflow for the paper, **“Spatial co-clustering of malaria and sickle cell disease burden across health facility catchments in Uganda: a Bayesian shared-component disease mapping analysis.”** The pipeline preprocesses modeled higher-level health facility catchments and routine DHIS2 admissions data, quantifies global and local spatial autocorrelation, fits Bayesian shared-component disease-mapping models, compares baseline, extended, and covariate-adjusted specifications, and derives hotspot, stability, and posterior-uncertainty summaries.

The workflow is designed to examine whether facility-recorded malaria and sickle cell disease (SCD) admissions share residual spatial structure across modeled service geographies in Uganda. It moves from descriptive mapping to global and local spatial dependence, then to joint Bayesian modeling of shared and disease-specific spatial components, and finally to hotspot interpretation and uncertainty classification.

Key outputs are cleaned analytic catchment datasets, global and local autocorrelation statistics, Bayesian model summaries, posterior risk surfaces, hotspot typologies, uncertainty maps, and manuscript-ready tables and figures.

## Folder layout

.
├─ data/  
│  ├─ raw/                          # place source inputs here  
│  ├─ interim/                      # intermediate analytic files  
│  └─ processed/                    # final outputs land here  
├─ R/  
│  ├─ 01_q2_step1_preprocess_catchments_and_outcomes.R  
│  ├─ 02_q2_step2_spatial_weights_and_global_moran.R  
│  ├─ 03_q2_step3_bivariate_local_lisa_and_multivariate_geary.R  
│  ├─ 04_q2_step4_bayesian_shared_component_model.R  
│  ├─ 05_q2_step5_extended_shared_component_model_comparison.R  
│  ├─ 06_q2_step6_hotspot_typology_stability_and_prioritisation.R  
│  ├─ 07_q2_step7_covariate_adjusted_shared_component_model.R  
│  └─ 08_q2_step8_uncertainty_and_residual_hotspot_interpretation.R  
├─ figures/                         # manuscript and supplementary figures  
├─ tables/                          # manuscript and supplementary tables  
├─ manuscript/  
│  ├─ main/                         # main manuscript files  
│  └─ supplementary/                # supplementary files  
├─ docs/  
│  ├─ README_workflow_notes.md  
│  ├─ METHODS_notes.md  
│  └─ RESULTS_notes.md  
├─ outputs/                         # optional exported model objects and summaries  
├─ renv/                            # reproducible R environment, if used  
├─ .gitignore  
├─ CITATION.cff  
├─ LICENSE  
├─ README.md  
└─ renv.lock  

## Required inputs

Place under `data/raw/` unless you edit paths in the scripts.

### Essential

- `catchments_traveltime_voronoi1_60min_clean4.*`  
  Cleaned higher-level health facility catchment polygons used as the spatial unit of analysis

- `q2_catchment_level_dataset.*`  
  Catchment-level analysis dataset containing cumulative malaria and SCD admissions, expected counts, and covariates

- `facility_locations.*`  
  Health facility point locations used for spatial linkage and map outputs

- `uganda_regions.*`  
  Regional boundary layer used for contextual covariates and mapping

- `uganda_districts.*`  
  District boundary layer for reference mapping and descriptive summaries

### Optional

- `adj_shared_risk_surface.*`  
  Saved posterior summaries of the adjusted shared component

- `adj_shared_exceedance_surface.*`  
  Catchment-level exceedance probabilities for adjusted hotspot classification

- `lisa_clusters.*`  
  Local spatial autocorrelation outputs used in residual hotspot interpretation

- `stability_scores.*`  
  Catchment-level hotspot stability metrics used in uncertainty-aware prioritisation

## Software

- R 4.3 or newer  
- GDAL, GEOS, PROJ  
- Recommended workflow management with `renv`

### Main R packages

- `sf`
- `terra`
- `spdep`
- `dplyr`
- `tidyr`
- `readr`
- `ggplot2`
- `tmap`
- `INLA`
- `stringr`
- `purrr`

Install packages as needed, or restore the project environment with:

```r
renv::restore()
```

## Run the pipeline

### Step 1. Preprocess catchments and outcomes

```r
source("R/01_q2_step1_preprocess_catchments_and_outcomes.R")
```

Prepares catchment-level malaria and SCD admission data, expected counts, and core covariates for analysis.

### Step 2. Spatial weights and global Moran statistics

```r
source("R/02_q2_step2_spatial_weights_and_global_moran.R")
```

Constructs spatial-neighbour structures and computes global autocorrelation statistics for malaria, SCD, and joint patterns.

### Step 3. Bivariate local LISA and multivariate Geary analysis

```r
source("R/03_q2_step3_bivariate_local_lisa_and_multivariate_geary.R")
```

Calculates bivariate Moran’s I, local bivariate LISA, and multivariate Geary’s C, and produces local cluster outputs.

### Step 4. Basic Bayesian shared-component model

```r
source("R/04_q2_step4_bayesian_shared_component_model.R")
```

Fits the baseline shared-component disease-mapping model for malaria and SCD admissions.

### Step 5. Extended model and model comparison

```r
source("R/05_q2_step5_extended_shared_component_model_comparison.R")
```

Fits the extended shared-component model with disease-specific spatial terms and compares model performance.

### Step 6. Hotspot typology, stability, and prioritisation

```r
source("R/06_q2_step6_hotspot_typology_stability_and_prioritisation.R")
```

Builds hotspot typologies, stability metrics, and prioritisation summaries from spatial and Bayesian outputs.

### Step 7. Covariate-adjusted shared-component model

```r
source("R/07_q2_step7_covariate_adjusted_shared_component_model.R")
```

Fits the fully adjusted shared-component model including linked service population, facility class, ownership, and region.

### Step 8. Posterior uncertainty and residual hotspot interpretation

```r
source("R/08_q2_step8_uncertainty_and_residual_hotspot_interpretation.R")
```

Summarises posterior uncertainty, constructs uncertainty classes, and classifies high-confidence, moderate-confidence, and uncertain residual hotspots.

## Key outputs

Typical outputs include:

- cleaned catchment-level analysis datasets
- neighbour lists and spatial weights objects
- global Moran and Geary statistics
- local LISA cluster classifications
- baseline, extended, and adjusted Bayesian model summaries
- posterior shared and disease-specific spatial effects
- exceedance-probability maps
- hotspot typology and stability tables
- uncertainty-class and residual hotspot confidence maps
- manuscript-ready figures and tables

Example output locations:

- `data/processed/`
- `figures/`
- `tables/`
- `outputs/`

## Manuscript linkage

This repository supports the paper:

**Spatial co-clustering of malaria and sickle cell disease burden across health facility catchments in Uganda: a Bayesian shared-component disease mapping analysis**

Main manuscript files should be placed in `manuscript/main/` and supplementary files in `manuscript/supplementary/`.

## Reproducibility

The analysis is organised as an ordered stepwise workflow. Scripts are intended to be run sequentially from Step 1 through Step 8. Intermediate objects saved by earlier steps are used by later steps for model fitting, hotspot interpretation, and uncertainty analysis.

For reproducibility:

- keep all raw inputs in `data/raw/`
- write derived files to `data/interim/` or `data/processed/`
- avoid editing outputs manually
- use `renv` to lock the R package environment
- record all path changes if adapting the workflow to a new machine

## Troubleshooting

### File or path errors
Check that all required source files are in `data/raw/` and that script paths match your local project structure.

### CRS mismatch
Confirm that all spatial layers use a consistent coordinate reference system before neighbour construction, overlay, or mapping.

### Missing neighbours
Some sparse or isolated polygons may require careful handling when building contiguity-based neighbour lists.

### INLA model issues
If the Bayesian models fail to run, check:
- data completeness
- expected-count construction
- factor reference levels
- graph structure for spatial random effects
- local INLA installation

### Large precision estimates
Very large disease-specific spatial precision estimates can occur when those structured terms are heavily shrunk toward zero. These should be interpreted cautiously.

## Citing

Please cite this repository using `CITATION.cff` and cite the main methodological and substantive sources used in the analysis, including:

- Bayesian shared-component disease mapping methods  
- spatial autocorrelation methods including Moran’s I and LISA  
- multivariate Geary statistics  
- INLA methodology for Bayesian spatial modeling  
- Uganda malaria and sickle cell epidemiology sources used in the manuscript

## License

See `LICENSE`.

Recommended approach:
- code under the MIT License
- manuscript text, tables, and figures under a Creative Commons license if desired
- restricted health data should **not** be deposited publicly unless data-sharing permissions allow it

## Contact

**George Paasi**  
Makerere University College of Health Sciences, School of Medicine, Kampala, Uganda  
Mbale Clinical Research Institute, Mbale, Uganda  
Email: georgepaasi8@gmail.com  
ORCID: 0000-0001-6360-0589  
Postal address: P.O. Box 1966, Mbale, Uganda

Co-authors: Grace Ndeezi, Ruth Namazzi, Ezekiel Mupere, Ian Guyton Munabi, Sarah Kiguli, Richard Idro, and Peter Olupot-Olupot.
