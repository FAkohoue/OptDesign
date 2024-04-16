# Define the function to evaluate the efficieny of the optimization
#
#' Generate variance-covariance matrix of treatments, traces and/or log(determinants) and other
#' related matrices for and RCBD with both random block and treatment effects.
#'
#' \code{VarCov_bRtR} generates the variance-covariance matrix (information matrix) of treatment
#' effects and calculates traces and/or log(determinants) for a linear mixed model with
#' both random block and treatment effects.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' or incomplete block designs (IBD) where 'Treatment' is the column of effects of interest.
#' @param criteria indicates the optimization criteria to report. It can be 'A' for A-optimal or
#' 'D' for D-optimal criteria. Default is 'A'.
#' @param s2Bl variance of the blocks effect. Default is 0.1
#' @param Ginv a variance-covariance matrix from pedigree or molecular data previously generated
#' @param Rinv an inverse of the error variance-covariance matrix.
#' @param K an intermediate matrix calculated from the original layout that was previously obtained.
#'
#' @return either a trace value or log of determinant of the variance-covariance matrix
#' (information matrix) for treatments, together with the K matrix to use in future runs.
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016) Generating experimental designs for spatially and
#' genetically correlated data using mixed models, Submitted to Australian and New Zealand
#' Journal of Statistics.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#'
#'
#' @export
#' @seealso \code{\code{\link{p_rep_design} \link{Rinv}}, \code{\link{Ginv}}

VarCov_bRtR_alp <- function(matdf, criteria="A", s2Bl=0.1, Ginv=NULL, Rinv=NULL, K=NULL) {
  # Checking if Rinv and Ginv are not provided
  if(is.null(Rinv)){
    stop('Rinv was not provided')
  }
  if(is.null(Ginv)){
    stop('Ginv was not provided')
  }

  # Obtaining Z matrix
  Z.block <- Matrix::sparse.model.matrix(~as.factor(matdf[,"Block"]) - 1)
  Z.trt <- Matrix::sparse.model.matrix(~as.factor(matdf[,"Treatment"]) - 1)
  Z <- cbind(Z.block,Z.trt)
  Z <- as.matrix(Z)

  # Obtaining Rinv
  Rinv <- as.matrix(Rinv)

  # Definining the block diagonal G and obtaining Ginv
  bb <- length(unique(matdf[,"Block"]))  # Number of blocks
  ng <- length(unique(matdf[,"Treatment"]))  # Number of treatments
  Ginv.block <- Gmatrix(ng=bb, VarG=s2Bl)   # Independent random effects with variance s2Bl
  Ginv.trt <- Matrix::drop0(round(Ginv,7))
  Ginv <- bdiag(Ginv.block, Ginv.trt)

  if(is.null(K)){

    # Obtaining X matrix: it is only the mean
    X <- matrix(data=1, nrow=nrow(matdf), ncol=1)

    # Obtaining C22 (for Random Blocks and Treatments)
    C11 <- Matrix::crossprod(as.matrix(X), as.matrix(Rinv)) %*% as.matrix(X) # It can be optimized with vector of 1
    C11inv <- 1/C11
    k1 <- Rinv %*% as.matrix(X)
    k2 <- Matrix::tcrossprod(as.matrix(C11inv), as.matrix(X))
    k3 <- k2 %*% Rinv
    K <- k1 %*% k3
    K <- as(K, "sparseMatrix")
    temp0 <- Matrix::crossprod(Z, Rinv) %*% Z + Ginv - Matrix::crossprod(Z, K) %*% Z

    # Rounding K matrix
    K <- round(K,7)

  } else {
    temp0 <- t(Z) %*% Rinv %*% Z + Ginv - t(Z) %*% K %*% Z
  }

  C22 <- solve(temp0)
  C22.trt <- C22[-c(1:bb),-c(1:bb)]
  C22.trt <- as(C22.trt, "sparseMatrix")  # This is the M(lambda) matrix for all random effects

  # Calculating Optimum Criteria (over matrix C22 of ALL random effects)
  if(criteria == "A"){
    OptimScore <- sum(Matrix::diag(C22.trt))
    return(list(OptimScore=OptimScore,K=K))
  }
  if(criteria == "D"){
    OptimScore <- log(Matrix::det(C22.trt))
    return(list(OptimScore=OptimScore,K=K))
  }

}
