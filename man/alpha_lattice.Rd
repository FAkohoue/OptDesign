\name{alpha_lattice}
\alias{alpha_lattice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create an Alpha lattice design for agricultural experiments
%%  ~~function to do ... ~~
}
\description{The function is designed to create an Alpha lattice design for agricultural experiments. This function takes several parameters such as treatment names, number of replicates, number of incomplete blocks, number of treatments, rows and columns per block, and an optional trial identifier.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alpha_lattice(trt_names = NULL, nrep, nblock, ntrt, rb, cb, trial = NULL)
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

\seealso{\code{\link{Gmatrix}}, \code{\link{Rmatrix}}, \code{\link{SwapMethods}}, \code{\link{p_rep_design}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Test the alpha_lattice_design function
trt_names <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
               'AB', 'BC', 'CD', 'DE', 'EF', 'FG', 'GH', 'HI', 'IJ', 'JK', 'KL', 'LM',
               'BA', 'CB', 'DC', 'ED', 'FE', 'GF', 'HG', 'IH', 'JI', 'KJ', 'LK', 'ML')
alpha_lattice(trt_names = trt_names, nrep = 3, nblock = 3, ntrt = 36, rb = 4, cb = 3, trial = 'ZKP')

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
alpha_lattice <- function(trt_names, nrep, nblock, ntrt, rb, cb, trial) {

  # Calculate the total number of plots
  plots_per_block <- rb * cb
  total_plots <- nrep * nblock * plots_per_block

  # Check if the total number of treatment allocations is equal to the total plots
  if (ntrt * nrep != total_plots) {
      stop("The total number of treatment allocations must be equal to the total number of plots.")
  }

  # Check if the number of treatments can be divided among the blocks within a replicate
  if (mod(ntrt, nblock) != 0) {
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
