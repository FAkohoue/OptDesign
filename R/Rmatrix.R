# Define the function to generate the inverse of the error (spatial or temporal) variance-covariance matrix R
#'
#' \code{Rmatrix} Generates the inverse of the error (spatial or temporal) variance-covariance
#' matrix R based on an autorregressive of order 1 homogeneous error structure (AR1) with or
#' without nugget. Best guesses of variance component parameters (rhox, rhoy, VarE and nugget)
#' need to be provided according to what is expected in the particular study.
#'
#' @import Matrix
#'
#' @param matdf an experimental design (layout) where 'Treatment' is the column of effects of interest
#' @param rhox spatial correlation between experimental units along the rows. Default value is 0.
#' @param rhoy spatial correlation between experimental units along the columns. Default value is 0.
#' @param VarE variance of the residuals. Default value is 1.
#' @param nugget spatial nugget error. Default value is 0.
#' @param regular a logical statement, if FALSE, a different (slower) algorithm for the AR1
#' error structure is called. Default is FALSE
#'
#' @return the inverse of the error variance-covariance matrix for homogeneuos AR1 structure in
#' sparce form.
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and genetically
#' correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' # Example: Unimproved regular-grid RCB designs
#' matdf <- rcbd(nblock=6, ntrt=9, rb=3, cb=3)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#' Rinv1 <- Rmatrix(matdf=matdf, VarE=0.9, regular=TRUE)
#' Rinv[1:5,1:5]
#' Rinv2 <- Rmatrix(matdf=matdf, VarE=0.9, rhox=0.6, rhoy=0.6, regular=TRUE)
#' Rinv2[1:5,1:5]
#' Rinv3<- Rmatrix(matdf=matdf, VarE=0.9, rhox=0.6, rhoy=0.6, regular=FALSE)
#' Rinv3[1:5,1:5]
#'
#' @export

Rmatrix <- function(matdf, VarE=1, rhox=0, rhoy=0, nugget=0, regular=FALSE) {

  # Debugging message to show the input parameters
  cat("Input Parameters - VarE:", VarE, "rhox:", rhox, "rhoy:", rhoy, "nugget:", nugget, "regular:", regular, "\n")

  #s2e <- 1-VarG-nugget
  s2e <- VarE
  matdf <- matdf[order(matdf[,"Row"], matdf[,"Col"]),]

  if(rhoy == 0 & rhox == 0) {
    Rinv <- round((1/(s2e+nugget)) * Matrix::Diagonal(nrow(matdf)), 10)
  }

  # Irregular Experiment
  if(regular == FALSE) {
    N <- nrow(matdf)
    R <- Matrix::Diagonal(N)

    # Debugging message to show progress
    cat("Calculating spatial correlation matrix for irregular experiment...\n")

    for(i in 1:(N-1)) {
      x1 <- matdf[,"Col"][i]
      y1 <- matdf[,"Row"][i]
      for (j in (i+1):nrow(matdf)){
        x2 <- matdf[,"Col"][j]
        y2 <- matdf[,"Row"][j]
        R[i,j] <- (rhox^abs(x2 - x1)) * (rhoy^abs(y2 - y1))
      }
    }

    R <- R + nugget * Matrix::Diagonal(N)
    R <- as.matrix(round(s2e * R, 10))
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    R <- as(R, "sparseMatrix")
    Rinv <- round(chol2inv(chol(R)), 10)
    Rinv <- as(Rinv, "sparseMatrix")
  }

  # Regular Experiment
  if(regular == TRUE) {
    N <- nrow(matdf)
    Tr <- max(matdf[,"Row"])
    Tc <- max(matdf[,"Col"])
    sigx <- Matrix::Diagonal(Tc)
    sigx <- rhox^abs(row(sigx) - col(sigx))
    sigy <- Matrix::Diagonal(Tr)
    sigy <- rhoy^abs(row(sigy) - col(sigy))
    R <- s2e * kronecker(sigy, sigx) + nugget * Matrix::Diagonal(N)
    R <- round(R, 7)
    R <- as(R, "sparseMatrix")
    Rinv <- round(chol2inv(chol(R)), 10)
    Rinv <- as(Rinv, "sparseMatrix")
  }

  Matrix::drop0(Rinv)
}
