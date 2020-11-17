#' Simulate principal component locations for simulated cells
#'
#' This function generates the principal component (PC) location for simulated cell. This function allows
#' for control over cell state frequency variation, as well as batch-associated and sample-associated gene
#' expression variation via (cf_scale), (b_scale), and (s_scale) respectively.
#' 
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
#' @param returnPCs Boolean determining whether the function will also return the simulated PC locations for reach cell
#'
#' @return If returnPCs = FALSE, this function returns a list containing the metadata table for the simulated cells. The
#' metadata table contains a dataframe with the following columns: "cellstate" which refers to the assigned cell state during
#' simulation, "sample" which refers to the sample the cell was assigned to, and "condition" which refers to the condition
#' the cell was assigned to (case or control). If clusterData = TRUE, the metadata table will also contain a column "new_clus", which
#' refers to the new cluster assignments for the simulated cells. If returnPCs = TRUE, this function returns a list containing: 
#' the metadata table for the simulated cells and a matrix containing the simulated PC locations for each cell.
#'
#' @import Seurat
#' @importFrom utils getFromNamespace
#' @export
simPCs <- function(ncases, nctrls, nbatches, batchStructure = NULL, ncells, centroids, pc_cov_list, batch_vars,
                   b_scale = 1, sample_vars, s_scale = 1, cfcov, cf_scale = 1, meanFreqs, clus, fc = 1, cond_induce = "cases",
                   res_use = 1.2, mc.cores = 1, clusterData = TRUE, returnPCs = FALSE){
    ndimensions <- ncol(centroids)
    nclusters <- nrow(centroids)
    # The "true" cell-states that simulated cells will be assigned to
    cellstates <- rownames(centroids)
    
    # Create a default batchstructure if none is provided
    if(is.null(batchStructure)){
        batchStructure <- distributeSamples(ncases = ncases, nctrls = nctrls, nbatches = nbatches)
    }
    # Determine which batches have at least one sample
    numSamplesBatch <- lapply(batchStructure[["batches"]], function(x){
        case_length <- length(x[["cases"]])
        ctrl_length <- length(x[["ctrls"]])
        return(case_length + ctrl_length)
    })
    nonEmptyBatches <- names(which(numSamplesBatch > 0))
    
    # Generate frequency distributions for each sample
    cf_sigma <- cfcov * cf_scale
    all_freqs <- generateFreqs(batchStructure = batchStructure, log_prior = log(meanFreqs), clus = clus, fc = fc,
                               cond_induce = cond_induce, cf_sigma = cf_sigma)
    
    # Genereate new dataset
    batch_shifts <- generateBatchShifts(batch_list = batchStructure[["batches"]], batch_vars = batch_vars, b_scale = b_scale)
    sample_shifts <- generateSampleShifts(sample_list = batchStructure[["sample_names"]], sample_vars = sample_vars, s_scale = s_scale)

    sim_pcs <- mclapply(nonEmptyBatches, function(x){
        generateSamplesInBatches(batchStructure = batchStructure, batch = x, case_freqs = all_freqs[["case_freqs"]],
                                 ctrl_freqs = all_freqs[["ctrl_freqs"]], cfcov = cfcov, cf_scale = cf_scale, centroids = centroids,
                                 batch_shifts = batch_shifts, sample_shifts = sample_shifts, pc_cov_list = pc_cov_list,
                                 ncells = ncells)
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, length(nonEmptyBatches)))
    
    new_pcs <- do.call("rbind", lapply(sim_pcs, "[[", 1)) %>% as.matrix
    meta <- do.call("rbind", lapply(sim_pcs, "[[", 2)) %>% data.frame
    meta$cellstate <- factor(meta$cellstate)
    
    # Cluster generated cells
    if(clusterData == TRUE){
        snn_ref <- BuildSNN(data.use = new_pcs, k.param = 30, prune.SNN = 1/15, nn.eps = 0, method = "annoy",
                        n.trees = 50, search.k = -1, mc.cores = mc.cores)
        new_clusters <- utils::getFromNamespace("RunModularityClustering", "Seurat")(SNN = snn_ref, modularity = 1, resolution = res_use,
                                                                                     algorithm = 1,n.start = 20, n.iter = 20, 
                                                                                     random.seed = 100, print.output = FALSE, 
                                                                                     temp.file.location = NULL, edge.file.name = NULL) 
        meta$new_clus <- factor(new_clusters)
    }
    
    if(returnPCs == TRUE){
        obj <- list(meta = meta, new_pcs = new_pcs)
    } else{
        obj <- list(meta = meta)
    }
    
    return(obj)
}

#' Simulate the samples in a single designated batch
#'
#' This function generates the principal component (PC) location for all simulated cells in a single designated
#' batch.
#'
#' @inheritParams simPCs
#' @param batch The name of the batch that will be simulated
#' @param case_freqs The cell state frequency distributions for all case samples. These are received as output from the
#' "generateFreqs" function
#' @param ctrl_freqs The cell state frequency distributions for all control samples. These are received as output from the
#' "generateFreqs" function
#' @param batch_shifts The simulated linear shifts in each batch for each cell state. These are received as output from the
#' "generateBatchShifts" function. 
#' @param sample_shifts The simulated linear shifts in each sample for each cell state. These are received as output from the
#' "generateSampleShifts" function. 
#'
#' @return Returns a metadata table and the simulated PC locations for all simulated cells in the designated batch
#' @importFrom MASS mvrnorm
generateSamplesInBatches <- function(batchStructure, batch = "batch1", case_freqs, ctrl_freqs, cfcov, cf_scale, centroids,
                                     batch_shifts, sample_shifts, pc_cov_list, ncells){
    ndimensions <- ncol(centroids)
    nclusters <- nrow(centroids)
    cases <- batchStructure$batches[[batch]]$cases
    ctrls <- batchStructure$batches[[batch]]$ctrls
    pcs_cases <- NULL
    pcs_ctrls <- NULL
    cellstate_cases <- NULL
    cellstate_ctrls  <- NULL
    ncells_cases <- 0
    ncells_ctrls <- 0
    
    if(length(cases) != 0){
        sim_cases <- lapply(cases, function(sample){
            adj_centroids <- centroids + batch_shifts[[batch]] + sample_shifts[[sample]]
            cellstate_assignments <- rep(1:nclusters, times = round(case_freqs[[sample]] * ncells[sample], 0))
            
            sim_pcs <- sapply(cellstate_assignments, function(x){
                MASS::mvrnorm(1, mu = adj_centroids[x, ] %>% as.numeric, Sigma = pc_cov_list[[x]])
            }) %>% data.frame %>% t
            
            return(list(sim_pcs = sim_pcs, cellstate_assignments = cellstate_assignments))
        })
        pcs_cases <- do.call("rbind", lapply(sim_cases, "[[", 1))
        cellstate_cases <- do.call("c", lapply(sim_cases, "[[", 2))
        ncells_cases <- unlist(lapply(sim_cases, function(x){x[[2]] %>% length}))
    } 
    if(length(ctrls) != 0){
        sim_ctrls <- lapply(ctrls, function(sample){
            adj_centroids <- centroids + batch_shifts[[batch]] + sample_shifts[[sample]]
            cellstate_assignments <- rep(1:nclusters, times = round(ctrl_freqs[[sample]] * ncells[sample], 0))
            
            sim_pcs <- sapply(cellstate_assignments, function(x){
                MASS::mvrnorm(1, mu = adj_centroids[x, ] %>% as.numeric, Sigma = pc_cov_list[[x]])
            }) %>% data.frame %>% t
            
            return(list(sim_pcs = sim_pcs, cellstate_assignments = cellstate_assignments))
        })
        pcs_ctrls <- do.call("rbind", lapply(sim_ctrls, "[[", 1))
        cellstate_ctrls <- do.call("c", lapply(sim_ctrls, "[[", 2))
        ncells_ctrls <- unlist(lapply(sim_ctrls, function(x){x[[2]] %>% length}))
    } 
    
    pcs_batch <- rbind(pcs_cases, pcs_ctrls)
    cellstates <- c(cellstate_cases, cellstate_ctrls) - 1
    
    meta_batch <- data.frame(
        cellstate = cellstates,
        sample = c(rep(cases, times = ncells_cases), rep(ctrls, times = ncells_ctrls)),
        condition = rep(c("case", "ctrl"), times = c(length(cellstate_cases), length(cellstate_ctrls)))
    )
    meta_batch$batch <- rep(batch, nrow(meta_batch))
    return(list(pcs_batch = pcs_batch, meta_batch = meta_batch))
}

#' Simulate batch shifts for each cell state
#'
#' This function generates linear shifts in each batch for each cell state. These shifts contribute to the adjusted centroid
#' that is used for simulating PC locations
#'
#' @inheritParams simPCs
#' @param batch_list A named list of the batch names e.g. batch1, batch 2, etc.
#'
#' @return Returns a matrix containing the batch shift for each cell state in each PC
#'
#' @importFrom MASS mvrnorm
generateBatchShifts <- function(batch_list, batch_vars, b_scale){
    ndimensions <- ncol(batch_vars)
    nclusters <- nrow(batch_vars)
    batch_sigma <- batch_vars * b_scale
    
    batch_shifts <- lapply(seq(length(batch_list)), function(x){
        shift <- sapply(seq(nclusters), function(x){
            MASS::mvrnorm(1, mu = rep(0, ndimensions), Sigma = diag(batch_sigma[x, ]))
        }) %>% t
        return(shift)
    })
    names(batch_shifts) <- names(batch_list)
    return(batch_shifts)
}

#' Simulate sample shifts for each cell state
#'
#' This function generates linear shifts in each samplefor each cell state. These shifts contribute to the adjusted centroid
#' that is used for simulating PC locations
#'
#' @inheritParams simPCs
#' @param sample_list A named list of the sample names e.g. sample1, sample2, etc.
#'
#' @return Returns a matrix containing the sample shift for each cell state in each PC
#'
#' @importFrom MASS mvrnorm
generateSampleShifts <- function(sample_list, sample_vars, s_scale){
    ndimensions <- ncol(sample_vars)
    nclusters <- nrow(sample_vars)
    sample_sigma <- sample_vars * s_scale
    sample_shifts <- lapply(seq(length(sample_list)), function(x){
        shift <- sapply(seq(nclusters), function(x){
            MASS::mvrnorm(1, mu = rep(0, ndimensions), Sigma = diag(sample_sigma[x, ]))
        }) %>% t
        return(shift)
    })
    names(sample_shifts) <- sample_list
    return(sample_shifts)
}