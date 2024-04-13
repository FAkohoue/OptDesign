#' Create an Alpha lattice design for agricultural experiments
#'
#' @param trt_names: Vector of unique treatment names
#' @param nrep: Number of replicates
#' @param nblock: Number of incomplete blocks
#' @param ntrt: Number of treatments
#' @param rb: Number of rows per block
#' @param cb: Number of columns per block
#' @param trial: Trial identifier (optional, can be used to differentiate between multiple trials)
#'
#' @return A matrix of field plan generated successfully
#' @export
#'
#' @examples
#' # Test the alpha_lattice_design function
# <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
#               'AB', 'BC', 'CD', 'DE', 'EF', 'FG', 'GH', 'HI', 'IJ', 'JK', 'KL', 'LM',
#               'BA', 'CB', 'DC', 'ED', 'FE', 'GF', 'HG', 'IH', 'JI', 'KJ', 'LK', 'ML')
#design <- alpha_lattice(trt_names = trt_names, nreplicates = 3, nblock = 3, ntreatments = 36, rb = 4, cb = 3, trial = 'ZKP')

# View the generated design
#print(design)
#'
#'
alpha_lattice <- function(trt_names, nrep, nblock, ntrt, rb, cb, trial) {

  # Calculate the total number of plots
  plots_per_block <- rb * cb
  total_plots <- nrep * nblock * plots_per_block

  # Check if the total number of treatment allocations is equal to the total plots
  if (ntrt * nrep != total_plots) {
    stop("The total number of treatment allocations must be equal to the total number of plots.")
  }

  # Check if the number of treatments can be divided among the blocks within a replicate
  if (ntrt %% nblock != 0) {
    stop("The number of treatments must be divisible by the number of blocks.")
  }

  cat("Generating block, row, column, and plot numbers...\n")

  # Generating rep numbers
  Rep_num <- rep(1:nrep, each = nblock * plots_per_block)

  # Generating block numbers
  Block_num <- rep(rep(1:nblock, each = plots_per_block), times = nrep)

  # Generating row numbers
  row_num <- rep(rep(1:rb, each = cb), times = nrep * nblock)

  # Generating column numbers
  col_num <- rep(rep(1:cb, each = rb), times = nrep * nblock)

  # Create the basic structure of the design
  design <- data.frame(
    Plot = 1:total_plots,
    Row = row_num,
    Col = col_num,
    Block = Block_num,
    Rep = Rep_num,
    Treatment = rep(NA, total_plots)
  )

  # Calculate block size
  block_size <- ntrt / nblock

  # Assign treatments to plots within each replicate
  trt_ID <- as.character(trt_names)
  for (rep in 1:nrep) {
    rep_indices <- which(design$Rep == rep)
    treatments <- sample(trt_ID, replace = FALSE)
    for (blk in 1:nblock) {
      block_indices <- rep_indices[((blk - 1) * block_size + 1):(blk * block_size)]
      design$Treatment[block_indices] <- treatments[((blk - 1) * block_size + 1):(blk * block_size)]
    }
  }

  # Generating plot numbers and adding trial
  trial_num <- rep(trial, length.out = total_plots)

  # Generating field plan dataframe
  matdf <- data.frame(design, Trial = trial_num)

  cat("Field plan generated successfully.\n")

  return(matdf)
}
