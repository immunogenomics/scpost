#' List of metadata and PC embeddings for the RA dataset after running Harmony integration
#'
#' @format 
#'   meta: data.frame of 5,265 rows (cells) and 10 columns. "preHarmClus" refers to
#'   cluster assignments from clustering the data before Harmony, while "harmClus"
#'   refers to cluster assignments from clustering the data after Harmony
#'   embeddings: data.frame of 5,265 row (cells) and 20 columns (harmonized PCs). Obtained from Harmony integration
#' 
#' @source \url{https://www.immport.org/shared/study/SDY998}
"ra_HarmObj"

#' List of metadata and PC embeddings for the RA dataset before running Harmony integration
#'
#' @format 
#'   meta: data.frame of 5,265 rows (cells) and 10 columns. "preHarmClus" refers to
#'   cluster assignments from clustering the data before Harmony, while "harmClus"
#'   refers to cluster assignments from clustering the data after Harmony
#'   embeddings: data.frame of 5,265 row (cells) and 20 columns (PCs).
#' 
#' @source \url{https://www.immport.org/shared/study/SDY998}
"ra_PreHarmObj"

#' List of metadata and PC embeddings for only the fibroblasts from the RA dataset before running Harmony integration. 
#' Thie was used for Fig 4a
#'
#' @format 
#'   meta: data.frame of 5,265 rows (cells) and 10 columns. "clusAllCells" refers to
#'   cluster assignments from clustering all of the RA data, while "clusOnlyFib"
#'   refers to cluster assignments from clustering only the fibroblasts (this was used in Fig 4a)
#'   embeddings: data.frame of 5,265 row (cells) and 20 columns (PCs).
#' 
#' @source \url{https://www.immport.org/shared/study/SDY998}
"ra_FibObj"
