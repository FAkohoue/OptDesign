# Define the function to randomly swap pairs of treatments for optimization
#'
#'
#' \code{SwapMethods} Randomly selects a block and a single or several pairs of treatments, swaps them
#' and creates a new experimental layout.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' or incomplete block designs (IBD).
#' @param pairs number of pairs of treatments to be swapped in the run.
#' @param swapmethod is the selected method to be used. This can be "within" which the default
#' (that should be used for randomized complete block designs) or "across" or "any" which can be used for
#' incomplete blocks or unbalanced designs.
#'
#' @return an new experimental design with swapped pairs of treatments with the 4 columns:
#' Row, Col, Rep, Treatment
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2015), Optimal Randomized Complete Block Designs for Spatially
#' and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and
#' Environmental Statistics, 150, 1-32.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' # Example: Regular-grid experiment with independent random effects
#' matdf <- rcbd(nblock=2, ntrt=9, rb=3, cb=3)  # Original design
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping within a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="within")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping within 2 pairs
#' newmatdf <- SwapMethods(matdf=matdf, pairs=2, swapmethod="within")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping across a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="across")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping any (random for within or across) for a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="any")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' @export
#' @seealso \code{\link{rcbd}}
#'
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
  if ("Row" %in% names(mat) && "Col" %in% names(mat)) {
    return(mat[order(mat[,"Row"], mat[,"Col"]), ])
  } else {
    return(mat)
  }
}
