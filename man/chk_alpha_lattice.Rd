\name{chk_alpha_lattice}
\alias{chk_alpha_lattice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Alpha lattice design with checks included in each block for agricultural experiments
%%  ~~function to do ... ~~
}
\description{Create an Alpha lattice design with checks included in each block for agricultural experiments.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
chk_alpha_lattice(trt_names = NULL, nrep, nblock, ntrt, rb, cb, chk, trial = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt_names}{Vector of unique treatment names
%%     ~~Describe \code{x} here~~
}
  \item{nrep}{Number of replicates
%%     ~~Describe \code{x} here~~
}
  \item{nblock}{Number of incomplete blocks
%%     ~~Describe \code{x} here~~
}
  \item{ntrt}{Number of treatments
%%     ~~Describe \code{x} here~~
}
  \item{rb}{Number of rows per block
%%     ~~Describe \code{x} here~~
}
  \item{cb}{Number of columns per block
%%     ~~Describe \code{x} here~~
}
  \item{chk}{Vector of unique check treatment names
%%     ~~Describe \code{x} here~~
}
  \item{trial}{Trial identifier (optional, can be used to differentiate between multiple trials)
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix of field plan generated successfully
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Félicien Akohoue
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Gmatrix}}, \code{\link{Rmatrix}}, \code{\link{SwapMethods}}, \code{\link{p_rep_design}}, \code{\link{alpha_lattice}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Define the control treatments
controls <- c('Control1', 'Control2', 'Control3', 'Control4')

# Call the function with control treatments
trt_names <- as.character(1:100)

chk_alpha_lattice(
  trt_names = trt_names,
  nrep = 3,
  nblock = 5,
  ntrt = 100,
  rb = 6,
  cb = 4,
  chk = controls,
  trial = 'ZKP'
)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
chk_alpha_lattice <- function(trt_names = NULL, nrep, nblock, ntrt, rb, cb, chk, trial = NULL) {

  # Ensure the controls are a character vector
  controls <- as.character(chk)

  # Number of controls
  ncontrols <- length(chk)

  # Calculate the total number of plots
  plots_per_block <- rb * cb
  total_plots <- nrep * nblock * plots_per_block

  # Check if the total number of treatment allocations is equal to the total plots
  if ((ntrt * nrep) + (ncontrols * nblock * nrep) != total_plots) {
    stop("The total number of treatment allocations must be equal to the total number of plots.")
  }

  # Check if the number of non-control treatments can be divided among the blocks within a replicate
  if (mod(ntrt, nblock) != 0) {
    stop("The number of non-control treatments must be divisible by the number of blocks.")
  }

  cat("Generating block, row, column, and plot numbers...\n")

  # Generating rep numbers
  Rep_num <- rep(1:nrep, each = nblock * plots_per_block)

  # Generating block numbers
  Block_num <- rep(rep(1:nblock, each = plots_per_block), times = nrep)

  # Generating row numbers
  row_num <- rep(rep(1:rb, each = cb), times = nrep * nblock)

  # Generating column numbers
  col_num <- rep(rep(1:cb, times = rb), times = nrep * nblock)

  # Create the basic structure of the design
  design <- data.frame(
    Plot = 1:total_plots,
    Row = row_num,
    Col = col_num,
    Block = Block_num,
    Rep = Rep_num,
    Treatment = rep(NA, total_plots)
  )

  # Suppress warning messages for the following block of code
  suppressWarnings({
    # Calculate block size for non-control treatments
    block_size <- (ntrt / nblock)

    # Assign treatments to plots within each replicate
    trt_ID <- as.character(trt_names)
    for (rep in 1:nrep) {
      rep_indices <- which(design$Rep == rep)
      for (blk in 1:nblock) {
        block_indices <- rep_indices[((blk - 1) * plots_per_block + 1):(blk * plots_per_block)]

        # Shuffle the non-control treatments for each block
        shuffled_trt_ID <- sample(trt_ID)

        # Calculate the number of non-control treatments to be included in the block
        n_non_control <- ntrt / nblock

        # Assign non-control treatments to the block
        design$Treatment[block_indices] <- shuffled_trt_ID[(((blk-1) * n_non_control) + 1):(blk * n_non_control)]

        # Assign control treatments to the block (exactly once in each block)
        control_indices <- sample(block_indices, ncontrols, replace = FALSE)
        design$Treatment[control_indices] <- chk
      }
    }
  })

  # Generating plot numbers and adding trial
  trial_num <- rep(trial, length.out = total_plots)

  # Generating field plan dataframe
  matdf <- data.frame(design, Trial = trial_num)

  cat("Field plan generated successfully.\n")

  return(matdf)
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
