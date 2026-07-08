#' Aggregate tile-level embeddings into a patient-level matrix
#'
#' Collapses a list of per-patient tile x feature embedding matrices into a
#' single patients x features matrix, by taking the column mean or median
#' across tiles within each patient.
#'
#' @param embedding_list Named list, one element per patient, each a tiles x
#'   features numeric matrix.
#' @param method \code{"mean"} (default) or \code{"median"} pooling across tiles.
#'
#' @return A patients x features matrix, row-named by the names of
#'   \code{embedding_list}.
#' @export
aggregate_embeddings <- function(embedding_list, method = "mean") {
  stopifnot(method %in% c("mean", "median"))
  agg_fn <- if (method == "mean") colMeans else function(x) apply(x, 2, median)
  Z <- do.call(rbind, lapply(embedding_list, agg_fn))
  rownames(Z) <- names(embedding_list)
  Z
}

#' Center and scale a patient-level embedding matrix
#'
#' Thin wrapper around \code{\link[base]{scale}} (center and unit-variance
#' scale each column), used to preprocess embeddings before mutual
#' information estimation.
#'
#' @param Z A patients x features numeric matrix.
#' @return The scaled matrix.
#' @keywords internal
#' @noRd
preprocess_embeddings <- function(Z) {
  scale(Z, center = TRUE, scale = TRUE)
}
