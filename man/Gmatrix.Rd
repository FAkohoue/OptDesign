\name{Gmatrix}
\alias{Gmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to generate the inverse of the random factor variance-covariance matrix G.
%%  ~~function to do ... ~~
}
\description{This function generates the inverse of the random factor variance-covariance matrix G by Cholesky decomposition.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Gmatrix(ng = 0, VarG = 1, G = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ng}{Number of levels for the random factor. Required if Amat is not provided.
%%     ~~Describe \code{x} here~~
}
  \item{VarG}{Variance of the random factor. Default value is 1.
%%     ~~Describe \code{x} here~~
}
  \item{G}{G matrix, usually the numerator relationship matrix. Its dimension should be ntrt x ntrt.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The inverse of the variance-covariance matrix G provided.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
%% ~put references to the literature/web site here ~
}
\author{ Salvador Gezan & Lazarus Mramba
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{p_rep_design}}, \code{\link{Rmatrix}}, \code{\link{SwapMethods}}, \code{\link{alpha_lattice}}, \code{\link{chk_alpha_lattice}}, \code{\link{VarCov_bRtR}}, \code{\link{VarCov_bRtR_alp}}, \code{\link{SwapMethods_alp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
Gmatrix <- function(ng = 0, VarG = 1, G = NULL) {
  if (is.null(G)) {
    Ginv <- round((1 / VarG) * Matrix::Diagonal(ng), 10)
  } else {
    G <- VarG * as.matrix(G)
    Ginv <- round(chol2inv(chol(as.matrix(G))), 10)
  }
  return(Ginv = Ginv)
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
