\name{Rmatrix}
\alias{Rmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to generate the inverse of the error (spatial or temporal) variance-covariance matrix R.
%%  ~~function to do ... ~~
}
\description{Generates the inverse of the error (spatial or temporal) variance-covariance matrix R based on an autorregressive of order 1 homogeneous error structure (AR1) with or without nugget. Best guesses of variance component parameters (rhox, rhoy, VarE and nugget) need to be provided according to what is expected in the particular study.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rmatrix(matdf, VarE=1, rhox=0, rhoy=0, nugget=0, regular=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matdf}{An experimental design (layout) where 'Treatment' is the column of effects of interest.
%%     ~~Describe \code{x} here~~
}
 \item{VarE}{Variance of the residuals. Default value is 1.
%%     ~~Describe \code{x} here~~
}
 \item{rhox}{Spatial correlation between experimental units along the rows. Default value is 0.
%%     ~~Describe \code{x} here~~
}
 \item{rhoy}{Spatial correlation between experimental units along the columns. Default value is 0.
%%     ~~Describe \code{x} here~~
}
 \item{nugget}{Spatial nugget error. Default value is 0.
%%     ~~Describe \code{x} here~~
}
 \item{regular}{A logical statement, if FALSE, a different (slower) algorithm for the AR1 error structure is called. Default is FALSE
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The inverse of the error variance-covariance matrix for homogeneuos AR1 structure in sparce form.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
%% ~put references to the literature/web site here ~
}
\author{Lazarus Mramba & Salvador Gezan
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{p_rep_design}}, \code{ \link{Gmatrix}}, \code{\link{SwapMethods}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

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
