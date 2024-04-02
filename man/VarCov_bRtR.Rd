\name{VarCov_bRtR}
\alias{VarCov_bRtR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to evaluate the efficieny of design optimization.
%%  ~~function to do ... ~~
}
\description{generates the variance-covariance matrix (information matrix) of treatment effects and calculates traces and/or log(determinants) for a linear mixed model with both random block and treatment effects.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VarCov_bRtR <- function(matdf, criteria = "A", s2Bl = 0.1, Ginv = NULL, Rinv = NULL, K = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matdf}{An experimental design (layout) based on a randomized complete block designs (RCBD) or incomplete block designs (IBD) where 'Treatment' is the column of effects of interest.
%%     ~~Describe \code{x} here~~
}
  \item{criteria}{Indicates the optimization criteria to report. It can be 'A' for A-optimal or 'D' for D-optimal criteria. Default is 'A'.
%%     ~~Describe \code{x} here~~
}
  \item{s2Bl}{Variance of the blocks effect. Default is 0.1.
%%     ~~Describe \code{x} here~~
}
  \item{Ginv}{A variance-covariance matrix from pedigree or molecular data previously generated.
%%     ~~Describe \code{x} here~~
}
  \item{Rinv}{An inverse of the error variance-covariance matrix.
%%     ~~Describe \code{x} here~~
}
  \item{K}{An intermediate matrix calculated from the original layout that was previously obtained.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{either a trace value or log of determinant of the variance-covariance matrix (information matrix) for treatments, together with the K matrix to use in future runs.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Mramba, L.K. and Gezan, S.A. (2016) Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
%% ~put references to the literature/web site here ~
}
\author{Lazarus Mramba & Salvador Gezan
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{p_rep_design}}, \code{\link{Gmatrix}}, \code{\link{Rmatrix}}, \code{\link{SwapMethods}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

VarCov_bRtR <- function(matdf, criteria="A", s2Bl=0.1, Ginv=NULL, Rinv=NULL, K=NULL) {
  # Checking if Rinv and Ginv are not provided
  if(is.null(Rinv)){
    stop('Rinv was not provided')
  }
  if(is.null(Ginv)){
    stop('Ginv was not provided')
  }

  # Obtaining Z matrix
  Z.block <- Matrix::sparse.model.matrix(~as.factor(matdf[,"Rep"]) - 1)
  Z.trt <- Matrix::sparse.model.matrix(~as.factor(matdf[,"Treatment"]) - 1)
  Z <- cbind(Z.block,Z.trt)
  Z <- as.matrix(Z)

  # Obtaining Rinv
  Rinv <- as.matrix(Rinv)

  # Definining the block diagonal G and obtaining Ginv
  bb <- length(unique(matdf[,"Rep"]))  # Number of blocks
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
