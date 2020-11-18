#' Distribute samples equally into batches in equal numbers of cases and controls
#'
#' Given the number of cases, number of controls, and number of batches in a planned study design
#' this function distributes samples equally into batches, in equal numbers of cases and controls
#' based off the "split" function.
#'
#' @param ncases The number of cases.
#' @param nctrls The number of controls.
#' @param nbatches The number of batches that samples will be distributed into.
#'
#' @return Returns a list containing the following elements: a list containing the batch
#' structure of the study (info on which samples are cases/controls and which batch they
#' were placed in), a vector containing the names of all samples, a vector containing the
#' names of case samples, and a vector containing the names of control samples. This output is used directly
#' as input into our simulation functions (batchStructure parameter).
#'
#' @export
distribSamples <- function(ncases, nctrls, nbatches){
    sample_names <- paste0("sample", 1:sum(ncases, nctrls))
    case_idx <- sample(length(sample_names), ncases)
    case_names <- sample_names[case_idx]
    ctrl_names <- sample_names[-case_idx]
    case_split <- suppressWarnings(split(case_names, factor(1:nbatches)))
    ctrl_split <- suppressWarnings(split(ctrl_names, factor(1:nbatches)))
    batches <- list()
    for(i in seq(nbatches)){
        comb <- list(case_split[[i]], ctrl_split[[i]])
        names(comb) <- c("cases", "ctrls")
        batches[[i]] <- comb
    }
    names(batches) <- paste0("batch", seq(nbatches))
    return(list(batches = batches, sample_names = sample_names, cases = case_names, ctrls = ctrl_names))
}

#' Distribute each sample into its own batch
#'
#' Given the number of cases, the number of controls, and number of batches, this function will
#' distribute a sample into its own batch. The number of bathes MUST be equal to the total
#' number of samples (nbatch = ncases + nctrls)
#'
#' @inheritParams distribSamples
#'
#' @return Returns a list containing the following elements: a list containing the batch
#' structure of the study (info on which samples are cases/controls and which batch they
#' were placed in), a vector containing the names of all samples, a vector containing the
#' names of case samples, and a vector containing the names of control samples. This output is used directly
#' as input into our simulation functions (batchStructure parameter).
#'
#' @export
distribSamplePerBatch <- function(ncases, nctrls, nbatches){
    if(nbatches != (ncases + nctrls)){
        stop("Number of batches is not equal to total number of samples")
    }
    sample_names <- paste0("sample", 1:sum(ncases,  nctrls))
    case_idx <- sample(length(sample_names), ncases)
    case_names <- sample_names[case_idx]
    ctrl_names <- sample_names[-case_idx]
    case_split <- split(case_names, factor(1:ncases))
    ctrl_split <- split(ctrl_names, factor(1:nctrls))
    batches <- list()
    for(i in 1:ncases){
        comb <- list(case_split[[i]], NULL)
        names(comb) <- c("cases", "ctrls")
        batches[[i]] <- comb
    }
    for(j in 1:nctrls){
        batch_ind <- ncases + j
        comb <- list(NULL,  ctrl_split[[j]])
        names(comb) <- c("cases", "ctrls")
        batches[[batch_ind]] <- comb
    }
    names(batches) <- paste0("batch", seq(nbatches))
    return(list(batches = batches, sample_names = sample_names, cases = case_names, ctrls = ctrl_names))
}

#' Split each sample into a specified number of subsamples before distributing into batches
#'
#' Given the number of cases, the number of controls, and number of batches, this function will
#' first split each sample (regardless of case-control status) into a specified number of 
#' equally sized subsamples (numSubsamples). These subsamples will then be distributed into 
#' different batches. Each batch will contain cells from multiple samples (equal to numSubsamples)
#'
#' @inheritParams distribSamples
#' @param numSubsamples The number of subsamples that each sample will be divided into
#'
#' @return Returns a list containing the following elements: a list containing the batch
#' structure of the study (info on which samples are cases/controls and which batch they
#' were placed in), a vector containing the names of all samples, a vector containing the
#' names of case samples, and a vector containing the names of control samples. This output is used directly
#' as input into our simulation functions (batchStructure parameter).
#'
#' @export
distribSampleSplit <- function(ncases, nctrls, nbatches, numSubsamples){
    sample_names <- paste0("sample", 1:sum(ncases, nctrls))
    case_idx <- sample(length(sample_names), ncases)
    case_names <- sample_names[case_idx]
    ctrl_names <- sample_names[-case_idx]
    case_samples <- rep(case_names, numSubsamples)
    ctrl_samples <- rep(ctrl_names, numSubsamples)
    case_split  <- split(sample(case_samples), factor(1:nbatches))
    ctrl_split <- split(sample(ctrl_samples), factor(1:nbatches))
    
    batches <- list()
    for(i in 1:nbatches){
        comb <- list(case_split[[i]], ctrl_split[[i]])
        names(comb) <- c("cases", "ctrls")
        batches[[i]] <- comb
    }
    names(batches) <- paste0("batch", seq(nbatches))
    return(list(batches = batches, sample_names = sample_names, cases = case_names, ctrls = ctrl_names))
}

