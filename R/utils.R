
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
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

calcBinomCI <- function(p, n, z = 1.96) {
   return(z * sqrt((p * (1 - p))/n))
}

AnnoyBuildObj <- function(data.use, metric = "euclidean", n.trees = 50){
    l <- ncol(x = data.use)
    annoyObj <- switch(EXPR = metric, 
                       euclidean = new(Class = RcppAnnoy::AnnoyEuclidean, l), 
                       cosine = new(Class = RcppAnnoy::AnnoyAngular, l), 
                       manhattan = new(Class = RcppAnnoy::AnnoyManhattan, l), 
                       hamming = new(Class = RcppAnnoy::AnnoyHamming, l), 
                       stop("Enter valid distance metric"))
    for(i in seq(nrow(x = data.use))){
        annoyObj$addItem(i - 1, data.use[i, ])
    }
    annoyObj$build(n.trees)
    return(annoyObj)
}

AnnoyGetNN <- function(obj, query, k.param = 30, search.k = -1, get.distance = TRUE, mc.cores = 2){
    n <- nrow(x = query)
    cosine_dist <- methods::is(obj, "Rcpp_AnnoyAngular")
    nn.idx <- mclapply(X = 1:n, function(x){
        nn <- obj$getNNsByVectorList(query[x, ], k.param, search.k, get.distance)
        if(cosine_dist){
            nn$distance <- 0.5 * (res$distance * res$distance)
        }
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

BuildSNN <- function(data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0, method = "annoy", n.trees = 50,
                     search.k = -1, metric = "euclidean", mc.cores = 1){
    switch(EXPR = method, rann = {
        knn <- RANN::nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    }, annoy = {
        annoyObj <- AnnoyBuildObj(data.use = data.use, metric = metric, n.trees = n.trees)
        knn <- AnnoyGetNN(obj = annoyObj, query = data.use, k.param = k.param, 
                          search.k = -1, get.distance = TRUE, mc.cores = mc.cores)
    })
    nn.ranked <- knn$nn.idx
    snn_res <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    return(snn_res)
}