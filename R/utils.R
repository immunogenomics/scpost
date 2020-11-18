#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#  Insert an element or vector of elements to designated positions in a vector
#'
#' Implementation from Tutur Qhuhuit that inserts an element/s into a vector
#'
#' @param vect The vector you are inserting elements into
#' @param pos The position you want to insert an element/s into
#' @param elems The element/s you want to insert into a vector
#'
#' @return A vector containing both the original and inserted elements
insertElems <- function(vect, pos, elems){
  l <- length(vect)
  k <- 0
  for (i in 1:length(pos)){
      if(pos[i] == 1){
          vect <- c(elems[k+1], vect)
      } else if(pos[i] == (l+1)){
          vect <- c(vect, elems[k+1])
      } else{
          vect <- c(vect[1:(pos[i] - 1 + k)], elems[k+1], vect[(pos[i] + k):(l + k)])
      }
      k <- k + 1
  }
  return(vect)
}

#  Determine if all elements of a vector are equal (from Hadley Wickham)
#'
#' Simple implementation from Hadley Wickham that checks if all elements
#' of a vector are equal
#'
#' @param x Vector containing numeric elements
#' @param tol The tolerated acceptable error in determining whether elements are equal
#'
#' @return A boolean representing if all elements of a vector are equal
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

#' Retrieve the parameters from a simulation
#'
#' Given a list of files that contain the results from simulating a dataset (via simDataset.base or simDataset.withMASC),
#' this function will retrieve the parameters used to simulate the dataset
#'
#' @param fileList A list of filenames that refer to the results from simulating a dataset
#'
#' @return Returns a dataframe containing the parameters used to simulate a dataset
#'
#' @importFrom gtools mixedsort
#' @export
getParsFromResFile <- function(fileList){
    rep <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[12]
    })
    clus <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[2]
    })
    ind_fc <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[3]
        return(as.numeric(unlist(strsplit(a, 'f'))[1]))
    }) 
    ncases <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[4]
        return(as.numeric(unlist(strsplit(a, 'c'))[1]))
    })
    nctrls <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[5]
        return(as.numeric(unlist(strsplit(a, 'c'))[1]))
    })
    nbatches <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[6]
        return(unlist(strsplit(a, 'b'))[1])
    }) 
    ncells <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[7]
        return(unlist(strsplit(a, 'n'))[1])
    })
    bscale <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[8]
        return(unlist(strsplit(a, 'b'))[1])
    }) 
    sscale <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[9]
        return(unlist(strsplit(a, 's'))[1])
    }) 
    cfscale <- sapply(fileList, function(x){
        a <- unlist(strsplit(x, '_'))[10]
        return(unlist(strsplit(a, 'c'))[1])
    }) 
    reso <- sapply(fileList, function(x){
        l <- unlist(strsplit(x, '_'))[11]
        return(unlist(strsplit(l, 'r'))[1])
    }) 
    nsamples <- as.numeric(ncases) + as.numeric(nctrls)
    
    tbl <- data.frame(
        rep = rep,
        clus = factor(clus, levels = gtools::mixedsort(clus %>% unique)),
        ind_fc = factor(ind_fc, levels = gtools::mixedsort(ind_fc %>% unique)),
        ncases = factor(ncases, levels = gtools::mixedsort(ncases %>% unique)),
        nctrls = factor(nctrls, levels = gtools::mixedsort(nctrls %>% unique)),
        nsamples = factor(nsamples, levels = gtools::mixedsort(nsamples %>% unique)),
        nbatches = factor(nbatches, levels = gtools::mixedsort(nbatches %>% unique)),
        ncells = factor(ncells, levels = gtools::mixedsort(ncells %>% unique)),
        bscale = factor(bscale, levels = gtools::mixedsort(bscale %>% unique)),
        sscale = factor(sscale, levels = gtools::mixedsort(sscale %>% unique)),
        cfscale = factor(cfscale, levels = gtools::mixedsort(cfscale %>% unique)),
        reso = factor(reso, levels = gtools::mixedsort(reso %>% unique)),
        row.names = NULL
    )
    return(tbl)
}

#' Retrieve the power results from simulations
#'
#' Given results, this function will calculate the number of simulations in which at least one cell state
#' cluster had a significant p-value frorm MASC analysis, meaning the cluster was expanded or depleted between
#' conditions.
#'
#' @param resFiles A list of filenames that refer to the results from simulating a dataset
#' @param resTables The MASC results tables obtained from MASC analysis on simulated datasets
#' @param threshold The p-value threshold that determines significance
#' @param z Determines the size of the confidence interval (CI) where z represents the % CI
#' @param stratByClus Boolean that determines whether the power results will be calculated as the aggregate of all
#' clusters or stratified by cluster.
#'
#' @return Returns a dataframe containing the power calculations from the simulated datasets
#'
#' @importFrom dplyr group_by tally summarise select mutate
#' @export
getPowerFromRes <- function(resFiles, resTables, threshold = 0.05, z = 1.96, stratByClus = FALSE){
    parTable <- getParsFromResFile(fileList = resFiles)
    pvals <- data.frame(
        min_masc_unadj = sapply(resTables, function(x){
            x[, "masc_pval"] %>% min(na.rm = TRUE)
        }),
        min_masc_adj = sapply(resTables, function(x){
            x[, "masc_adj"] %>% min(na.rm = TRUE)
        }),
        parTable
    )
    if(stratByClus == TRUE){
        powerVals <- pvals %>% dplyr::group_by(.data$clus, .data$ind_fc, .data$ncases, .data$nctrls, .data$nsamples, .data$ncells, 
                                               .data$bscale, .data$sscale, .data$cfscale) %>%
            dplyr::tally(name = "trials") 
        masc_power <- pvals %>% dplyr::group_by(.data$clus, .data$ind_fc, .data$ncases, .data$nctrls, .data$nsamples, .data$ncells, 
                                                .data$bscale, .data$sscale, .data$cfscale) %>%
            dplyr::summarise(masc_power = sum(.data$min_masc_adj < threshold), .groups = "drop") %>% dplyr::select(.data$masc_power)
        powerTable <- cbind.data.frame(powerVals, masc_power = masc_power) %>% dplyr::mutate(masc_power = masc_power / .data$trials)
    } else{
        powerVals <- pvals %>% dplyr::group_by(.data$ind_fc, .data$ncases, .data$nctrls, .data$nsamples, .data$ncells, 
                                               .data$bscale, .data$sscale, .data$cfscale) %>%
            dplyr::tally(name = "trials") 
        masc_power <- pvals %>% dplyr::group_by(.data$ind_fc, .data$ncases, .data$nctrls, .data$nsamples, .data$ncells, 
                                                .data$bscale, .data$sscale, .data$cfscale) %>%
            dplyr::summarise(masc_power = sum(.data$min_masc_adj < threshold), .groups = "drop") %>% dplyr::select(.data$masc_power)
        powerTable <- cbind.data.frame(powerVals, masc_power = masc_power) %>% dplyr::mutate(masc_power = masc_power / .data$trials)
    }
    powerTable$masc_power_ci <- apply(powerTable, 1, function(x){
        calcBinomCI(p = as.numeric(x['masc_power']), n = as.numeric(x['trials']), z = z)
    })
    return(powerTable)
}

#' Calculate a binomial confidence interval
#'
#' Given an observed proportion (p) and a number of trials (n), calculate the
#' half-width of a (z)% binomial confidence interval
#'
#' @param p The observed proportion
#' @param n The total number of samples/trials used to calculate the proportion
#' @param z The quantile of a standard normal distribution
#'
#' @return A singular numeric value representing the half-width of the CI
#' 
calcBinomCI <- function(p, n, z = 1.96) {
   return(z * sqrt((p * (1 - p))/n))
}