
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

## Tutorials

  - Get started with this [short
    tutorial](https://github.com/immunogenomics/scpost/blob/main/vignettes/GettingStarted_Tutorial.ipynb)
    that runs through a typical scPOST workflow.
  - If you just want the simulated dataset for analyses other than
    differential abundance testing, you check out this
    [tutorial](https://github.com/immunogenomics/scpost/blob/main/vignettes/RetrievingSimulations_Tutorial.ipynb)
    that shows how to retrieve your simulated datasets.

## Basic simulation workflow

scPOST uses a workflow that comprises 3 general steps

![Workflow](https://github.com/immunogenomics/scpost/blob/main/docs/images/PowerFig1.png)

The following code performs these steps with a pre-loaded dataset
(ra\_FibObj) provided in the scpost package

### Step 1: Parameter estimation

``` r
raFib_freqEstimates <- estimateFreqVar(
    meta = ra_FibObj$meta, 
    clusCol = 'clusOnlyFib', 
    sampleCol = 'sample', 
    logCov = TRUE
)
```

``` r
raFib_pcEstimates <- estimatePCVar(
    pca = ra_FibObj$embeddings, 
    npcs = 20, 
    meta = ra_FibObj$meta, 
    clusCol = 'clusOnlyFib',
    sampleCol = 'sample', 
    batchCol = 'batch'
)
```

### Step 2 and 3 together: Dataset simulation and association testing with MASC

In order to run many simulations, we recommend creating a parameter
table that lists the parameters you plan to use for each simulation. We
provide an example here.

``` r
# Set up parameter table for multiple simulations
set.seed(23)

# Set the number of samples, number of cells per sample, and create batch structure
ncases <- 10
nctrls <- 10
nbatches <- 4
batchStructure <- distribSamples(ncases = ncases, nctrls = nctrls, nbatches = nbatches)
ncells <- rep(250, times = ncases + nctrls)
names(ncells) <- batchStructure$sample_names

params <- createParamTable(
    nreps = 5,
    clus = "clus0",
    fc = 5,
    ncases = ncases,
    nctrls = nctrls,
    nbatches = nbatches,
    b_scale = 1,
    s_scale = 1,
    cf_scale = 1,
    res_use = 0.6,
    cond_induce = "cases",
    save_path = file.path(getwd(), "scpostSims/gettingStarted/")
)
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
