\name{SwapMethods_alp}
\alias{SwapMethods_alp}
\title{Function to Randomly Swap Pairs of Treatments in an Alpha Lattice Design}
\description{
This function randomly selects a block and swaps one or more pairs of treatments to optimize an Alpha lattice experimental design with repeated checks.
The function supports within-block and across-block swaps to improve spatial randomization.
}
\usage{
SwapMethods_alp(matdf, pairs = 1, swapmethod = "within")
}
\arguments{
  \item{matdf}{A data frame representing an experimental design layout based on a randomized complete block design (RCBD) or an incomplete block design (IBD).
  The data frame should contain at least the columns: \code{Block}, \code{Treatment}, and optionally \code{Row} and \code{Col}.}

  \item{pairs}{An integer specifying the number of treatment pairs to swap in the run. Each pair consists of two treatments that will be exchanged in the layout.}

  \item{swapmethod}{A character string specifying the swap method. Options are:
    \describe{
      \item{\code{"within"}}{Swaps occur only within a single block (recommended for RCBD).}
      \item{\code{"across"}}{Swaps occur between two different blocks (recommended for unbalanced designs).}
      \item{\code{"any"}}{Randomly selects between \code{"within"} and \code{"across"} for each swap.}
    }
  The default is \code{"within"}.}
}
\details{
This function is designed to improve the spatial distribution of treatments within an experimental design. The goal is to minimize local clustering of treatments while maintaining the structure of the original design.

\strong{Swapping Mechanism:}
\itemize{
  \item For \code{"within"}, the function randomly selects one block and swaps \code{pairs} of treatments within that block.
  \item For \code{"across"}, the function selects two blocks and swaps treatments between them.
  \item For \code{"any"}, the function randomly chooses between \code{"within"} and \code{"across"} at each iteration.
}

\strong{Assumptions:}
\itemize{
  \item The input \code{matdf} should contain a column \code{Block} that identifies the block assignments.
  \item The column \code{Treatment} should be numeric or character identifiers for the treatments.
  \item If the columns \code{Row} and \code{Col} are present, the returned object will be sorted accordingly.
}
}
\value{
A data frame with the same structure as \code{matdf}, but with the specified swaps applied.
The returned data frame contains:
  \itemize{
    \item \code{Row}: (Optional) Row number in the field layout.
    \item \code{Col}: (Optional) Column number in the field layout.
    \item \code{Block}: Block assignment.
    \item \code{Treatment}: Treatment identifier, modified with swapped values.
  }
}
\references{
Mramba, L.K. and Gezan, S.A. (2015). "Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data Using Mixed Models," Journal of Agricultural, Biological and Environmental Statistics, 150, 1-32.
}
\author{
Lazarus Mramba \email{lazarus.mramba@example.com} and Salvador Gezan \email{salvador.gezan@example.com}
}
\seealso{
\code{\link{p_rep_design}}, \code{\link{alpha_lattice}}, \code{\link{SwapMethods}}
}
\examples{

# Simulated example dataset
set.seed(123)
example_df <- data.frame(
  Block = rep(1:4, each = 5),
  Treatment = sample(1:20, 20, replace = FALSE),
  Row = rep(1:5, times = 4),
  Col = rep(1:4, each = 5)
)

# Swap 2 pairs within a block
new_design <- SwapMethods_alp(example_df, pairs = 2, swapmethod = "within")

# Swap across blocks
new_design_across <- SwapMethods_alp(example_df, pairs = 2, swapmethod = "across")

# Swap using a random method
new_design_any <- SwapMethods_alp(example_df, pairs = 2, swapmethod = "any")

# View results
head(new_design)

SwapMethods_alp <- function(matdf, pairs = 1, swapmethod = "within") {

  trt <- max(matdf[,"Treatment"])
  gsize <- pairs * 2
  stopifnot(gsize %% 2 == 0)
  if (gsize > trt) {
    stop("Number of swaps is larger than the number of Treatments")
  }

  # Randomly selecting a swap method if swapmethod == "any"
  if (swapmethod == "any") {
    swapmethod <- sample(c("within", "across"), 1)  # Selects randomly from 'within' or 'across'
  }

  mat <- as.data.frame(matdf, row.names = NULL)

  # Swapping only within a block for any pair
  if (swapmethod == "within") {
    b1 <- sample(mat$Block, 1, replace = TRUE)  # Selects a block at random

    # Get the indices of the treatments within the selected block
    b1_indices <- which(mat$Block == b1)

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
    blks <- sample(unique(mat$Block), 2, replace = FALSE)

    # Get the indices of the treatments within the selected blocks
    blk1_indices <- which(mat$Block == blks[1])
    blk2_indices <- which(mat$Block == blks[2])

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
  if ("Row" %in% names(mat) && "Col" %in% names(mat))
  {
    return(mat[order(mat[,"Row"], mat[,"Col"]), ])
  } else {
    return(mat)
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
