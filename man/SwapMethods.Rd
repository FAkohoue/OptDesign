\name{SwapMethods}
\alias{SwapMethods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to randomly swap pairs of treatments for optimization.
%%  ~~function to do ... ~~
}
\description{Randomly selects a block and a single or several pairs of treatments, swaps them and creates a new experimental layout.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SwapMethods(matdf, pairs = 1, swapmethod = "within")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matdf}{An experimental design (layout) based on a randomized complete block designs (RCBD) or incomplete block designs (IBD).
%%     ~~Describe \code{x} here~~
}
  \item{pairs}{Number of pairs of treatments to be swapped in the run.
%%     ~~Describe \code{x} here~~
}
  \item{swapmethod}{Selected method to be used. This can be "within" which the default (that should be used for randomized complete block designs) or "across" or "any" which can be used for incomplete blocks or unbalanced designs.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{an new experimental design with swapped pairs of treatments with the 4 columns: Row, Col, Rep, Treatment
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Mramba, L.K. and Gezan, S.A. (2015), Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and Environmental Statistics, 150, 1-32.
%% ~put references to the literature/web site here ~
}
\author{Lazarus Mramba & Salvador Gezan
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{p_rep_design}}, \code{\link{Gmatrix}}, \code{\link{Rmatrix}}, \code{\link{alpha_lattice}}, \code{\link{chk_alpha_lattice}}, \code{\link{VarCov_bRtR}}, \code{\link{VarCov_bRtR_alp}}, \code{\link{SwapMethods_alp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

SwapMethods <- function(matdf, pairs=1, swapmethod="within") {

  trt <- max(matdf[,"Treatment"])
  gsize <- pairs*2
  stopifnot(gsize %% 2 == 0)
  if (gsize > trt) {
    stop("Number of swaps is larger than the number of Treatments")
  }

  # Swapping any blocks any pairs
  if (swapmethod == "any") {
    list <- c('within', 'across')
    sel <- sample(list, 1)
    swapmethod <- sel
  }

  mat <- as.data.frame(matdf, row.names = NULL)

  # Swapping only within a block for any pair
  if (swapmethod == "within") {
    b1 <- sample(mat$Rep, 1, replace = TRUE)  # Selects a block at random

    # Get the indices of the treatments within the selected block
    b1_indices <- which(mat$Rep == b1)
    # Sample gsize indices without replacement from the treatments in the selected block
    g1_indices <- sample(b1_indices, gsize, replace = FALSE)

    # Perform the swaps
    temp <- mat$Treatment[g1_indices]
    temp[seq(1, gsize, 2)] <- mat$Treatment[g1_indices[seq(2, gsize, 2)]]
    temp[seq(2, gsize, 2)] <- mat$Treatment[g1_indices[seq(1, gsize, 2)]]

    # Replace the treatments in the selected block with the swapped ones
    mat$Treatment[g1_indices] <- temp
  }

  # Swapping across blocks for any pairs
  if (swapmethod == "across") {
    blks <- sample(unique(mat$Rep), 2, replace = FALSE)

    # Get the indices of the treatments within the selected blocks
    blk1_indices <- which(mat$Rep == blks[1])
    blk2_indices <- which(mat$Rep == blks[2])

    # Sample half gsize indices without replacement from each block
    g1_indices <- sample(blk1_indices, gsize/2, replace = FALSE)
    g2_indices <- sample(blk2_indices, gsize/2, replace = FALSE)

    # Perform the swaps
    temp1 <- mat$Treatment[g1_indices]
    temp2 <- mat$Treatment[g2_indices]
    mat$Treatment[g1_indices] <- temp2
    mat$Treatment[g2_indices] <- temp1
  }

  # Check if 'Row' and 'Col' columns exist before sorting
  if
  if ("Row" %in% names(mat) && "Col" %in% names(mat)){
    return(mat[order(mat[,"Row"], mat[,"Col"]), ])
  } else {
    return(mat)}
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
