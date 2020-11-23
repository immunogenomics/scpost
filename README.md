
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scPOST

<!-- badges: start -->

<!-- badges: end -->

Simulation of single-cell datasets for power analyses that estimate
power to detect cell state frequency shifts between conditions (e.g. an
expansion of a cell state in disease vs. healthy).

## Installation

Install the current version of scPOST from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/scpost")
```

### Installation notes:

  - You may need to install the latest version of devtools (because of
    the recent GitHub change from “master” to “main” terminology, which
    can cause previous versions of `install_github` to fail).

# Usage/Demos

## Quick start

Check out this [short
tutorial](https://github.com/immunogenomics/scpost/blob/main/vignettes/raFib_tutorial.ipynb)
that outlines a typical workflow.

## Basic simulation workflow

scPOST uses a workflow that comprises 3 general steps:

Step 1: Estimates variance parameters from real data

Step 2: Simulation of a new dataset based on estimated parameters in
Step 1

Step 3: Association testing with MASC to detect cell state frequency
shifts between conditions

![Workflow](docs/PowerFig1.pdf)

The following code performs these steps with a pre-loaded dataset
provided in the scpost package

### Step 1: Parameter estimation

``` r
raFib_freqEstimates <- estimateFreqVar(meta = ra_FibObj$meta, clusCol = 'clusOnlyFib', sampleCol = 'sample', logCov = TRUE)
```

``` r
raFib_pcEstimates <- estimatePCVar(pca = ra_FibObj$embeddings, npcs = 20, meta = ra_FibObj$meta, clusCol = 'clusOnlyFib',
                                   sampleCol = 'sample', batchCol = 'batch')
```

### Step 2 and 3 together:

In order to run many simulations, we recommend creating a parameter
table that lists the parameters you plan to use for each simulation. We
provide an example here.

``` r
# Set up parameter table for multiple simulations
set.seed(23)

# Set the number of samples, number of cells per sample, and create batch structure
ncases <- 17
nctrls <- 3
nbatches <- 4
batchStructure <- distribSamples(ncases = ncases, nctrls = nctrls, nbatches = nbatches)
ncells <- rep(250, times = ncases + nctrls)
names(ncells) <- batchStructure$sample_names

# Retrieve estimates
meanFreqs <- raFib_freqEstimates$meanFreq
cfcov <- raFib_freqEstimates$cfcov
centroids <- raFib_pcEstimates$centroids
pc_cov_list <- raFib_pcEstimates$pc_cov_list
batch_vars <- raFib_pcEstimates$batch_vars
sample_vars <- raFib_pcEstimates$sample_vars

# Set scaling parameters (realistic)
b_scale <- 1
s_scale <- 1
cf_scale <- 1

# Set cluster names to induce fc, magnitude of fc, and resolution to re-cluster
clus <- "clus0"
fc <- 5
resolution <- 0.6
mc.cores <- 1

save_path <- file.path(getwd(), "scpostSims/unbalanced/")
reps <- 1:5

# Create parameter table for running multiple simulations
params <- expand.grid(
    rep = reps, ncases = ncases, nctrls = nctrls, nbatches = nbatches, b_scale = b_scale, s_scale = s_scale,
    cf_scale = cf_scale, clus = clus, fc = fc, res_use = resolution, save_path = save_path
)
params$cond_induce = "cases"
params$seed <- sample(.Machine$integer.max, size = nrow(params))
```

Once you have the parameter table, you can perform many simulations
(note, a save directory is required so that you can save the results of
your simulations):

``` r
suppressWarnings({
    lapply(seq(nrow(params)), function(x){
            simDataset.withMASC(
                save_path = params[x, 'save_path'],
                rep = params[x, 'rep'],
                seed = params[x, 'seed'],
                ncases = params[x, 'ncases'],
                nctrls = params[x, 'nctrls'],
                nbatches = params[x, 'nbatches'],
                batchStructure = batchStructure,
                ncells = ncells,
                centroids = centroids,
                pc_cov_list = pc_cov_list,
                batch_vars = batch_vars,
                b_scale = params[x, 'b_scale'],
                sample_vars = sample_vars,
                s_scale = params[x, 's_scale'],
                cfcov = cfcov,
                cf_scale = params[x, 'cf_scale'],
                meanFreqs = meanFreqs,
                clus = params[x, 'clus'],
                fc = params[x, 'fc'],
                cond_induce = params[x, 'cond_induce'],
                res_use = params[x, 'res_use'], 
                mc.cores = 1,
                clusterData = TRUE,
                returnPCs = FALSE
            )
    })
})
```
