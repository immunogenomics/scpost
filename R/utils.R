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
#' @export
calcBinomCI <- function(p, n, z = 1.96) {
   return(z * sqrt((p * (1 - p))/n))
}

#' Build a shared nearest-neighbor graph
#'
#' Returns a shared nearest-neighbor graph by first creating a nearest-neighbor
#' graph with Annoy and then compute shared nearest-neighbor with Seurat's cpp implementation
#'
#' @param data.use Data used to create nearest-neighbor graph. 
#' @param k.param Number of nearest-neighbors to calculate.
#' @param prune.SNN Cutoff for removing edges signaling low shared overlap of neighbors.
#' @param nn.eps Error bound for using RANN to calculate nearest-neighbors.
#' @param method Method to use for creating nearest-neighbor graph. Annoy tends to be faster, but does not give
#' exact solution
#' @param n.trees Number of trees used in Annoy search query. Higher value gives higher precision
#' @param search.k Number of nodes an Annoy search query touches. -1 invokes n.trees * n search
#' @param mc.cores Number of cores to use in a parallel format
#'
#' @return Returns an snn object containing the shared nearest-neighborr graph
#' 
#' @importFrom RANN nn2
#' @importFrom utils getFromNamespace
#' @import Seurat
BuildSNN <- function(data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0, method = "annoy", n.trees = 50,
                     search.k = -1, mc.cores = 1){
    switch(EXPR = method, rann = {
        knn <- RANN::nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    }, annoy = {
        annoyObj <- AnnoyBuildObj(data.use = data.use, n.trees = n.trees)
        knn <- AnnoyGetNN(obj = annoyObj, query = data.use, k.param = k.param, 
                          search.k = -1, get.distance = TRUE, mc.cores = mc.cores)
    })
    nn.ranked <- knn$nn.idx
    snn_res <- utils::getFromNamespace("ComputeSNN", "Seurat")(nn_ranked = nn.ranked, prune = prune.SNN)
    return(snn_res)
}

#' Build an Annoy object 
#'
#' @inheritParams BuildSNN
#'
#' @return Returns an annoy object, as documented by RcppAnnoy
#'
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
AnnoyBuildObj <- function(data.use, n.trees = 50){
    l <- ncol(x = data.use)
    annoyObj <- new(Class = RcppAnnoy::AnnoyEuclidean, l)
    for(i in seq(nrow(x = data.use))){
        annoyObj$addItem(i - 1, data.use[i, ])
    }
    annoyObj$build(n.trees)
    return(annoyObj)
}

#' Find nearest-neighbors of elements in an Annoy Object 
#'
#' Given an Annoy object, this function quickly finds the nearest-neighbors of each element
#'
#' @inheritParams BuildSNN
#' @param obj Annoy obj used to search for nearest-neighbors
#' @param query The data used to create the Annoy object, containing elements you want to find nearest-
#' neighbors of
#' @param get.distance Include distances 
#'
#' @return Returns a list containing the nearest-neighbors of each element and the distances between
#' elements and their nearest-neighbors
#'
#' @importFrom parallel mclapply
AnnoyGetNN <- function(obj, query, k.param = 30, search.k = -1, get.distance = TRUE, mc.cores = 1){
    n <- nrow(x = query)
    nn.idx <- mclapply(X = 1:n, function(x){
        nn <- obj$getNNsByVectorList(query[x, ], k.param, search.k, get.distance)
        list(nn$item + 1, nn$distance)
    }, mc.cores = mc.cores)
    idx <- do.call(rbind, lapply(nn.idx, "[[", 1))
    if(get.distance){
        dist <- do.call(rbind, lapply(nn.idx, "[[", 2))
    } else{
        dist <- matrix(nrow = n, ncol = k.param)
    }
    return(list(nn.idx = idx, nn.dists = dist))
}