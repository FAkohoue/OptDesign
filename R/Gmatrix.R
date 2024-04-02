# Define the function to generate the inverse of the random factor variance-covariance matrix G
#
#' \code{Gmatrix} Generates the inverse of the random factor variance-covariance matrix G by
#' Cholesky decomposition.
#'
#' @import Matrix
#'
#' @param ng number of levels for the random factor. Required if Amat is not provided.
#' @param VarG variance of the random factor. Default value is 1
#' @param G G matrix, usually the numerator relationship matrix. Its dimension should be ntrt x ntrt
#'
#' @return the inverse of the variance-covariance matrix G provided
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @author
#' Salvador Gezan & Lazarus Mramba
#'
#' @examples
#' # Example 1: G-inverse from a matrix of 30 independent effects
#' Ginv <- Gmatrix(ng=30, VarG=1)
#' Ginv[1:5,1:5]
#'
#' # Example 2: G-inverse obtained from a pedigree with 30 genotypes
#' data(ped30fs)
#' head(ped30fs)
#' Amat <- makeA(ped30fs)
#' Amat <- as.matrix(Amat[-c(1:5),-c(1:5)])   # Only offspring, eliminating the 5 parents.
#' Amat[1:5,1:5]
#' Ginv <- Gmatrix(VarG=1, G=Amat)
#' Ginv[1:5,1:5]

#Gmatrix <- function(ng=0, VarG=1, G=NULL){
#  if (is.null(G)){
#    Ginv <- round((1/VarG)*Matrix::Diagonal(ng),10)
#  } else {
#    G <- VarG*as.matrix(G)
#    Ginv <- round(chol2inv(chol(as.matrix(G))),10)
#  }
#  return(Ginv=Ginv)
#}

Gmatrix <- function(ng = 0, VarG = 1, G = NULL) {
  if (is.null(G)) {
    Ginv <- round((1 / VarG) * Matrix::Diagonal(ng), 10)
  } else {
    G <- VarG * as.matrix(G)
    Ginv <- round(chol2inv(chol(as.matrix(G))), 10)
  }
  return(Ginv = Ginv)
}
