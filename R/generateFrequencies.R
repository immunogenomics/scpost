#' Generate cell state frequency distributions for samples
#'
#' Given a batchStructure and baseline frequency distribution (in log space), this function will
#' generate a cell state frequency distribution for each sample. This function also allows users to induce
#' a designated fold change (fc) into either case samples or control samples, and control the magnitude
#' of covariance that cell states have with each other via cf_sigma (e.g. increase or decrease the cell state frequency
#' variation across samples). The magnitude of the fold change will be the ratio of case to control cells (e.g. inducing 
#' a fold change of 2 in cases will result in there being, on average, 2 times more case cells than control cells of that
#' cluster
#'
#' @param batchStructure The structure of the study design in which cases and controls are split into batches. These
#' structures are output by the "distributeSample" functions (which can then be modified if specific structure is desired).
#' @param log_prior A named vector containing the mean frequencies of the prototype dataset's cell states (log space). The
#' "estimateFreqVar" function returns this a mean frequency vector in linear space (can be transformed into log space via the 
#' "log" function).
#' @param clus The name of the cluster in which a fold change will be induced.
#' @param fc The magnitude of the fold change that will be induced in the chosen cluster. If no fold change is desired, set
#' fc = 1.
#' @param cond_induce The condition you wish to induce a fold change in. Setting cond_induce = "cases" will induce a fold 
#' change into cases, while setting cond_induce = "ctrls" will induce a fold change into controls.
#' @param cf_sigma A matrix containing the covariance between cell states. This matrix is received as output from 
#' the "estimateFreqVar" function
#'
#' @return Returns a list containing: a list of cell state frequencies for all case samples and a list
#' of cell state frequencies for all control samples
#'
#' @importFrom MASS mvrnorm
generateFreqs <- function(batchStructure, log_prior, clus, fc = 1, cond_induce = "cases", cf_sigma){
    if(cond_induce == "cases"){
        case_freqs <- lapply(batchStructure[["cases"]], function(x){
            getInducedFreqs(log_prior = log_prior, clus = clus, fc = fc, cf_sigma = cf_sigma) %>% prop.table
        })
        ctrl_freqs <- lapply(batchStructure[["ctrls"]], function(x){
            exp(MASS::mvrnorm(n = 1, mu = log_prior, Sigma = cf_sigma)) %>% prop.table
        })
    } else if(cond_induce == "ctrls"){
        case_freqs <- lapply(batchStructure[["cases"]], function(x){
            exp(MASS::mvrnorm(n = 1, mu = log_prior, Sigma = cf_sigma)) %>% prop.table
        })
        ctrl_freqs <- lapply(batchStructure[["ctrls"]], function(x){
            getInducedFreqs(log_prior = log_prior, clus = clus, fc = fc, cf_sigma = cf_sigma) %>% prop.table
        })
    }
    names(case_freqs) <- batchStructure[["cases"]]
    names(ctrl_freqs) <- batchStructure[["ctrls"]]
    return(list(case_freqs = case_freqs, ctrl_freqs = ctrl_freqs))
}

#' Generate a cell state frequency distribution with induced frequencies
#'
#' Given a vector of cell state frequencies (log space), this function will induce a specified fold change (fc)
#' into the specified cluster (clus)
#'
#' @inheritParams generateFreqs
#'
#' @return Returns a generated vector in which the designated cluster has the designated fold change induced
#'
#' @importFrom MASS mvrnorm
getInducedFreqs <- function(log_prior, clus, fc = 1, cf_sigma){
    baseline <- exp(MASS::mvrnorm(n = 1, mu = log_prior, Sigma = cf_sigma)) %>% prop.table
    return(induceFC.capped(baseline = baseline, clus = clus, fc = fc))
}

#' Induce a fold change in an element in a cluster
#'
#' Given a baseline distribution, induce a fold change in a cluster
#'
#' @inheritParams generateFreqs
#' @param baseline A named vector containing the frequencies of cell states
#'
#' @return Returns a vector in which the designated cluster has the designated fold change induced
#'
induceFC.capped <- function (baseline, clus, fc = 1){
    baseline <- baseline %>% prop.table
    clus_idx <- which(names(baseline) == clus)
    clusxfc <- baseline[clus_idx] * fc
    if (clusxfc > 1) {
        clusxfc <- 1
    }
    otherclusters <- (baseline[-clus_idx]/sum(baseline[-clus_idx])) * 
        (1 - clusxfc)
    return(insertElems(otherclusters, pos = clus_idx, elems = clusxfc))
}

