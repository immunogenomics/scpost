#' Simulate a dataset and save the results (metadata table and p-values from MASC analysis)
#'
#' This function simulates a dataset and then performs MASC analysis on the generated data. The function saves the results
#' as a list containing a metadata table for all simulated cells, as well as the MASC analysis results for all of the
#' simulated cell states. The simulated PC locations for all cells can optionally be saved.
#' 
#' @param save_path The name of the directory the results will be saved to.
#' @param rep A numeric value representing the replicate number of the simulated dataset.
#' @param seed A numeric value representing the seed that will be set before simulating the dataset.
#' @param ncases The number of cases.
#' @param nctrls The number of controls.
#' @param nbatches The number of batches that samples will be distributed into.
#' @param batchStructure The structure of the study design in which cases and controls are split into batches. These
#' structures are output by the "distributeSample" functions (which can then be modified if specific structure is desired).
#' If this parameter is kept as NULL, this function will automatically create a batchStructure with the "distributeSamples"
#' function. 
#' @param ncells A vector containing the number of cells that will be simulated for each sample. The vector must be 
#' the same length as the number of total samples (ncases + nctrls). 
#' @param centroids The mean PC values for each cell state. These are obtained as output from the "estimatePCVar" function.
#' @param pc_cov_list A list containing the residual variance-covariance matrices for each cell state. These are obtained
#' as output from the "estimatePCVar" function.
#' @param batch_vars A matrix containing the batch_associated variance for each cell state in each PC. These are obtained
#' as output from the "estimatePCVar" function.
#' @param b_scale The magnitude of batch-associated gene expression variation the simulated dataset will exhibit. Setting b_scale
#' = 1 will result in realistic levels of batch-associated variation (as derived from parameter estimation of the input dataset).
#' Increasing b_scale results in higher variation, while decreasing b_scale results in lower variation.
#' @param sample_vars A matrix containing the sample_associated variance for each cell state in each PC. These are obtained
#' as output from the "estimatePCVar" function.
#' @param s_scale The magnitude of sample-associated gene expression variation the simulated dataset will exhibit. Setting s_scale
#' = 1 will result in realistic levels of sample-associated variation (as derived from parameter estimation of the input dataset).
#' Increasing s_scale results in higher variation, while decreasing s_scale results in lower variation.
#' @param cfcov The cell state frequency variance-covariance across samples. This is obtained as output from the "estimateFreqVar"
#' function
#' @param cf_scale The magnitude of cell state frequency variation that cell states will exhibit across samples in the simulated
#' dataset. Setting cf_scale = 1 will result in realistic levels of cell state frequency variation (as derived from parameter 
#' estimation of the input dataset). Increasing cf_scale results in higher variation, while decreasing cf_scale results in lower variation.
#' @param meanFreqs A vector containing the mean frequencies of cell states (linear space) across samples from the original input 
#' dataset. This vector is obtained as output from the "estimateFreqVar" function.
#' @param clus The name of the cluster in which a fold change will be induced.
#' @param fc The magnitude of the fold change that will be induced in the chosen cluster. If no fold change is desired, set
#' fc = 1.
#' @param cond_induce The condition you wish to induce a fold change in. Setting cond_induce = "cases" will induce a fold 
#' change into cases, while setting cond_induce = "ctrls" will induce a fold change into controls.
#' @param res_use The resolution that will be used for clustering (Louvain method) the simulated dataset. 
#' @param mc.cores The number of cores that will be used for simulating the dataset and clustering.
#' @param clusterData Boolean determining whether the simulated dataset will be clustered. 
#' @param returnPCs Boolean determining whether the function will also return the simulated PC locations for reach cell.
#' @param null_mod The right-hand side of the formula that will be used as the null model in MASC analysis
#' @param full_mod The right-hand side of the formula that will be used as the full model in MASC analysis
#' @param adj_method The p-value correction method that will be used via the "p.adjust" function
#' @param verbose Print out time dataset was simulated at
#'
#' @return This function returns NULL and instead saves the results to the directory designated in "save_path". If returnPCs = FALSE, 
#' the saved results will be a list containing the metadata table for the simulated cells. The metadata table contains a dataframe with 
#' the following columns: "cellstate" which refers to the assigned cell state during simulation, "sample" which refers to the sample the 
#' cell was assigned to, "batch" which refers to the sample the cell was assigned to, and "condition" which refers to the condition the 
#' cell was assigned to (case or control). If clusterData = TRUE, the metadata table will also contain a column "new_clus", which refers 
#' to the new cluster assignments for the simulated cells. If returnPCs = TRUE, the saved results will contain the metadata table for 
#' the simulated cells and a matrix containing the simulated PC locations for each cell.
#'
#' @export
simDataset.withMASC <- function(save_path, rep = 1, seed = 1, ncases, nctrls, nbatches, batchStructure = NULL, ncells, centroids, 
                                pc_cov_list, batch_vars, b_scale = 1, sample_vars, s_scale = 1, cfcov, cf_scale = 1, meanFreqs, clus, 
                                fc = 1, cond_induce = "cases", res_use = 1.2, mc.cores = 1, clusterData = TRUE, returnPCs = FALSE,
                                null_mod = "1 + (1|batch) + (1|sample)", full_mod = "condition + (1|batch) + (1|sample)",
                                adj_method = "bonferroni", verbose = TRUE){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    
    data <- simPCs(
        ncases = ncases,
        nctrls = nctrls,
        nbatches = nbatches,
        batchStructure = batchStructure,
        ncells = ncells,
        centroids = centroids,
        pc_cov_list = pc_cov_list,
        batch_vars = batch_vars,
        b_scale = b_scale,
        sample_vars = sample_vars,
        s_scale = s_scale,
        cfcov = cfcov,
        cf_scale = cf_scale,
        meanFreqs = meanFreqs,
        clus = clus,
        fc = fc, 
        cond_induce = cond_induce,
        res_use = res_use,
        mc.cores = mc.cores,
        clusterData = clusterData,
        returnPCs = returnPCs
    )
    
    res <- getPvals.MASC(
        meta = data[["meta"]],
        clusterCol = "new_clus",
        null_mod = null_mod,
        full_mod = full_mod,
        mc.cores = mc.cores,
        adj_method = adj_method
    )
    data[["res"]] <- res
    
    filename <- paste('res_', clus, '_', fc, 'fc_', ncases, 'cases_', nctrls, 'ctrls_', nbatches, 'batches_',
                      mean(ncells), 'ncells_', b_scale, 'bscale_', s_scale, 'sscale_', cf_scale, 'cfscale_', res_use, 'reso_',
                      'rep', rep, '_seed', seed, '.rds', sep = "")
    saveRDS(object = data, file = file.path(save_path, filename))
    
    if(verbose){
        return(message(paste("Simulated dataset at", Sys.time())))
    }
}    
 

#' Simulate a dataset and save the results
#'
#' This function simulates a dataset and then saves the results as a list containing a metadata table for all simulated cells. The simulated 
#' PC locations for all cells can optionally be saved.
#' 
#' @inheritParams simDataset.withMASC
#'
#' @return This function returns NULL and instead saves the results to the directory designated in "save_path". If returnPCs = FALSE, 
#' the saved results will be a list containing the metadata table for the simulated cells. The metadata table contains a dataframe with 
#' the following columns: "cellstate" which refers to the assigned cell state during simulation, "sample" which refers to the sample the 
#' cell was assigned to, and "condition" which refers to the condition the cell was assigned to (case or control). If clusterData = TRUE, 
#' the metadata table will also contain a column "new_clus", which refers to the new cluster assignments for the simulated cells. 
#' If returnPCs = TRUE, the saved results will contain the metadata table for the simulated cells and a matrix containing the simulated 
#' PC locations for each cell.
#'
#' @export
simDataset.base <- function(save_path, rep = 1, seed = 1, ncases, nctrls, nbatches, batchStructure = NULL, ncells, centroids, 
                            pc_cov_list, batch_vars, b_scale = 1, sample_vars, s_scale = 1, cfcov, cf_scale = 1, meanFreqs, clus, fc = 1, 
                            cond_induce = "cases", res_use = 1.2, mc.cores = 1, clusterData = TRUE, returnPCs = FALSE, verbose = TRUE){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    
    data <- simPCs(
        ncases = ncases,
        nctrls = nctrls,
        nbatches = nbatches,
        batchStructure = batchStructure,
        ncells = ncells,
        centroids = centroids,
        pc_cov_list = pc_cov_list,
        batch_vars = batch_vars,
        b_scale = b_scale,
        sample_vars = sample_vars,
        s_scale = s_scale,
        cfcov = cfcov,
        cf_scale = cf_scale,
        meanFreqs = meanFreqs,
        clus = clus,
        fc = fc, 
        cond_induce = cond_induce,
        res_use = res_use,
        mc.cores = mc.cores,
        clusterData = clusterData,
        returnPCs = returnPCs
    )
    
    filename <- paste('data_', clus, '_', fc, 'fc_', ncases, 'cases_', nctrls, 'ctrls_', nbatches, 'batches_',
                      mean(ncells), 'ncells_', b_scale, 'bscale_', s_scale, 'sscale_', cf_scale, 'cfscale_', res_use, 'reso_', 
                      'rep', rep, '_seed', seed, '.rds', sep = "")
    saveRDS(object = data, file = file.path(save_path, filename))
    
    if(verbose){
        message(paste("Simulated dataset at", Sys.time()))
    }
}