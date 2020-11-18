#' Estimate the mean and variation in each cell state's frequency across samples
#'
#' Given a metadata table containing each cell's sample identity and assigned cell
#' state, this function will return a list containing: a vector containing the mean frequency
#' of each cell state (in linear space) and a variance-covariance matrix containing the covariances 
#' between each cell state (log space by default). Both of these elements are used as direct 
#' inputs to our dataset generation functions.
#'
#' @param meta A metadata table that should include named columns containing: Cell IDs, 
#' cell state assignments, and sample identities for each cell.
#' @param clusCol The name of the metadata column containing the cell state cluster assignments for each cell.
#' @param sampleCol The name of the metadata column containing the sample identities of each cell.
#' @param logCov Boolean designating whether you want the variance-covariance matrix in log space (default)
#' or in linear space. The log space matrix is used as input to our simulation functions; only use
#' linear space for visualization.
#'
#' @return Returns a list containing: a vector containing the mean frequency
#' of each cell state (in linear space) and a variance-covariance matrix containing the covariances 
#' between each cell state (log space by default). The mean frequency vector will be transformed into 
#' log space during data simulation.
#'
#' @importFrom dplyr filter group_by tally select
#' @importFrom tidyr complete
#' @importFrom gtools mixedsort
#' @importFrom rlang .data
#' @export
estimateFreqVar <- function(meta, clusCol, sampleCol, logCov = TRUE){
    nclus <- meta[, clusCol] %>% unique %>% length
    sample_meta <- meta[, c(clusCol, sampleCol)]
    colnames(sample_meta) <- c('clus', 'sample')
    sample_list <- sample_meta[, 'sample'] %>% unique
    
    # Calculate frequencies of each cell state in each sample (gets converted to proportions)
    tot_freqTable <- do.call('rbind', sapply(sample_list, function(x){
        sample_meta %>% dplyr::filter(sample == x) %>% dplyr::group_by(.data$clus) %>% 
            dplyr::tally(name = "freq") %>% tidyr::complete(.data$clus, fill = list(freq = 0)) %>% 
            dplyr::select(-.data$clus)
    }))
    tot_freqTable <- tot_freqTable + 1
    rownames(tot_freqTable) <- sample_list
    colnames(tot_freqTable) <- paste0('clus', gtools::mixedsort(meta[, clusCol] %>% unique))
    tot_propTable <- sweep(tot_freqTable, 1, rowSums(tot_freqTable), '/')
    
    # Calculate the variance-covariance matrix for cell states across samples
    if(logCov == TRUE){
        tot_cfcov <- stats::cov(log(tot_propTable))
    } else if(logCov == FALSE){
        tot_cfcov <- stats::cov(tot_propTable)
    }
    
    # Calculate the mean frequency of each cell state across samples
    meanFreq <- tot_propTable %>% colMeans
    names(meanFreq) <- colnames(tot_propTable)
    
    return(list(meanFreq = meanFreq, cfcov = tot_cfcov))
}

#' Estimate the sources of variation in principal component (PC) space
#'
#' Given principal component (PC) embeddings and a metadata table containing cell state assignments, 
#' sample identities, and batch identities, this function uses principal variance component analysis (PVCA) 
#' to deconvolute the sources of variation in each PC. For PVCA, this function fits a linear mixed
#' effects model with sample and batch fit as random effects. We recommend that datasets have multiple
#' samples and batches in order to accurately estimate the variance contributions of sample and batch.
#'
#' @inheritParams estimateFreqVar
#' @param pca The PC embeddings from a PCA (cell x PC). Every PC should contain a value for each cell.
#' @param npcs The number of PCs whose variance will be deconvoluted. As this function's output is used 
#' as direct input for our simulation functions, this will determine how many PCs your generated datasets
#' will contain (e.g. npcs = 20 here will mean that the generated datasets will contain 20 PCs).
#' @param meta A metadata table that should include named columns containing: Cell IDs, 
#' cell state assignments, sample identities, and batch identities for each cell.
#' @param batchCol The name of the metadata column containing batch identities of each cell.
#' @param parallel Boolean determining whether the linear models will be fit in parallel (highly recommend)
#' if possible. Computation time is determined by number of cells in the data you are fitting.
#' @param mc.cores Number of cores used if models are being fit in parallel.
#' @param save The path the models will be saved to if a directory is provided. While you can save the 
#' output of this function (variance estimates), you can save the fitted models if you wish tocheck their
#' fitted values individually.
#'
#' @return Returns a list containing the following elements: a matrix containing the mean PC values for each
#' cell state, a list where each element is the residual variance-covariance matrix for each cell state,
#' a matrix containing the batch-associated variance for each cell state in each PC, a matrix containing the
#' sample-associated variance for each cell state in each PC. 
#'
#' @importFrom gtools mixedsort
#' @importFrom lme4 lmer lmerControl .makeCC
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom parallel mclapply
#' @export
estimatePCVar <- function(pca, npcs, meta, clusCol, sampleCol, batchCol, parallel = FALSE, mc.cores = 1, save = NULL){
    colnames(pca) <- paste0("PC", 1:ncol(pca))
    # Create strings containig model formulas
    models <- paste(paste0("PC", 1:npcs), "~ 1 + (1|sample) + (1|batch)")
    
    pca_model_meta <- cbind.data.frame(clus = meta[, clusCol],
                                       batch = meta[, batchCol],
                                       sample = meta[, sampleCol],
                                       pca[, 1:npcs])
    colnames(pca_model_meta) <- c("clus", "batch", "sample", paste0("PC", 1:npcs))
    clus_apply <- gtools::mixedsort(pca_model_meta[, "clus"] %>% unique)
    
    # Fit models
    if(parallel == TRUE){
        suppressWarnings({
            fitted_mods <- mclapply(clus_apply, function(x){
                mclapply(models, function(y){
                    lme4::lmer(formula = y, data = pca_model_meta %>% dplyr::filter(.data$clus == x), REML = TRUE,
                               control = lme4::lmerControl(optimizer = "nloptwrap", check.nlev.gtr.1 = "ignore",
                                                           check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))
                }, mc.preschedule = TRUE, mc.cores = mc.cores)
            })
        })
    } else{
        suppressWarnings({
            fitted_mods <- lapply(clus_apply, function(x){
                lapply(models, function(y){
                    lme4::lmer(formula = y, data = pca_model_meta %>% dplyr::filter(.data$clus == x), REML = TRUE,
                               control = lme4::lmerControl(optimizer = "nloptwrap", check.nlev.gtr.1 = "ignore",
                                                           check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))
                })
            })
        })
    }
    
    if(!is.null(save)){
        saveRDS(fitted_mods, save)
    }
    
    # Retrieve estimates from models
    mod_vars <- getModVars(mods = fitted_mods, npcs = npcs, clusNames = gtools::mixedsort(meta[, clusCol] %>% unique))
    return(mod_vars)
}

#' Retrieve variables from models fitted in the "estimatePCVar" function
#'
#' Given fitted models and the names for each cell state, retrieve for each cell state: the centroids 
#' (fixed effect), residual variance-covariance matrices (residuals), batch-associated variances (random effect variance),
#' and sample associated variances (random effect variance).
#'
#' @inheritParams estimatePCVar
#' @param mods A list of lme4 odels that this function will extract the variables from 
#' @param clusNames A list of unique names of the cell state cluster assignments 
#'
#' @return Returns a list containing the following elements: a matrix containing the mean PC values for each
#' cell state, a list where each element is the residual variance-covariance matrix for each cell state,
#' a matrix containing the batch-associated variance for each cell state in each PC, a matrix containing the
#' sample-associated variance for each cell state in each PC. 
#'
#' @importFrom lme4 VarCorr
getModVars <- function(mods, npcs, clusNames){
    # Get centroids
    centroids <- do.call(rbind, lapply(mods, function(x){
        a <- lapply(x, function(y){
            return(y@beta[1])
        })
    })) %>% data.frame
    centroids <- apply(centroids, 2, as.numeric)
    rownames(centroids) <- paste0("clus", clusNames)
    colnames(centroids) <- paste0("PC", 1:npcs, "_mean") 
    
    # Get residual variance
    pc_cov_list <- lapply(mods, function(x){
        resids <- do.call(rbind, lapply(x, function(y){
            return(stats::residuals(y))
        }))
        resids <- resids %>% t
        colnames(resids) <- paste0("PC", 1:npcs)
        return(stats::cov(resids))
    })
    
    # Get batch taus
    batch_vars <- do.call(rbind, lapply(mods, function(x){
        a <- lapply(x, function(y){
            b <- as.data.frame(lme4::VarCorr(y))
            return(b[2, "vcov"])
        })
    }))
    colnames(batch_vars) <- paste0("PC", 1:npcs, "_batchVar")
    batch_vars <- batch_vars %>% data.frame
    batch_vars <- apply(batch_vars, 2, as.numeric)
    
    # Get donor taus
    sample_vars <- do.call(rbind, lapply(mods, function(x){
        a <- lapply(x, function(y){
            b <- as.data.frame(lme4::VarCorr(y))
            return(b[1, "vcov"])
        })
    }))
    colnames(sample_vars) <- paste0("PC", 1:npcs, "_donorVar")
    sample_vars <- sample_vars %>% data.frame
    sample_vars <- apply(sample_vars, 2, as.numeric)
    
    return(list(centroids = centroids, pc_cov_list = pc_cov_list, batch_vars = batch_vars, sample_vars = sample_vars))
}