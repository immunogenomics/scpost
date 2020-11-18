#' Test clusters for an expansion (or depletion) between conditions
#'
#' Given the results of a simulated dataset (meta), this function will test all of the new clusters
#' for differential abundance (expansion or depletion) between conditions.
#' 
#' @param meta A metadata table received as output from simulating datasets (via simDataset.base or simDataset.withMASC). The
#' metadata table should contain the following columns: "cellstate" which refers to the assigned cell state during simulation, 
#' "sample" which refers to the sample the cell was assigned to, "batch" which refers to the batch the cell was assigned to, and 
#' "condition" which refers to the condition the cell was assigned to (case or control).
#' @param clusterCol The name of the metadata column containing the new cell state clusters. Our functions designate this column 
#' "new_clus" in the saved results.
#' @param null_mod The right-hand side of the formula that will be used as the null model in MASC analysis
#' @param full_mod The right-hand side of the formula that will be used as the full model in MASC analysis
#' @param mc.cores The number of cores that will be used for simulating the dataset and clustering. The computation time will
#' depend on the number of simulated cells.
#' @param adj_method The p-value correction method that will be used via the "p.adjust" function
#'
#' @return Returns a dataframe containing the calculated MASC p-values for each new cluster, as well as
#' corrected p-values (default Bonferroni corrected)
#'
#' @importFrom parallel mclapply
#' @export
getPvals.MASC <- function(meta, clusterCol, null_mod, full_mod, mc.cores = 1, adj_method = "bonferroni"){
    cluster_list <- meta[ ,clusterCol] %>% levels
    mascRes <- mclapply(cluster_list, function(x){
        mascCalc(
            meta = meta,
            clusterCol = clusterCol,
            clusterName = x,
            null_mod = null_mod,
            full_mod = full_mod
        )
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, length(cluster_list)))
    
    res <- data.frame(
        new_clus = cluster_list,
        masc_pval = unlist(mascRes) %>% as.numeric
    )
    res$masc_adj <- stats::p.adjust(res$masc_pval, method = adj_method)
    return(res)
}

#' Perform a test for differential abundance for a cell state cluster with MASC
#'
#' Given a metadata table, this function will test a desginated cell state cluster (clusterName) for differential abundance (expansion
#' or depletion between conditions)
#' with MASC
#' 
#' @inheritParams getPvals.MASC
#' @param clusterName The name of the cell state cluster that will be tested forr differential abundance.
#'
#' @return Returns a p-value (numeric) obtained from MASC analysis
#'
#' importFrom dplyr group_by tally pull
#' importFrom lme4 glmer glmerControl .makeCC
#' @export
mascCalc <- function(meta, clusterCol, clusterName, null_mod, full_mod){
    meta$cluster.onehot <- ifelse(meta[[clusterCol]] == clusterName, 1, 0)
    meta_frame <- meta %>% dplyr::group_by(.data$condition, .data$batch, .data$sample, .data$cluster.onehot) %>% 
        dplyr::tally(name = "weights")
    weights <- meta_frame %>% dplyr::pull(.data$weights)
    null.fmla <- paste("cluster.onehot ~", null_mod)
    full.fmla <- paste("cluster.onehot ~", full_mod)
    
    null.model <- lme4::glmer(formula = null.fmla, data = meta_frame, family = "binomial", nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa", 
                                                           check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)), 
                        weights = weights)
    full.model <- lme4::glmer(formula = full.fmla, data = meta_frame, family = "binomial", nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa", 
                                                           check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)),
                        weights = weights)
    lrt <- stats::anova(null.model, full.model)
    return(lrt[["Pr(>Chisq)"]][[2]])
}