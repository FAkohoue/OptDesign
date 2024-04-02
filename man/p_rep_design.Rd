\name{p_rep_design}
\alias{p_rep_design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a field plan for agricultural experiments with treatments, controls, and p-reps.
%%  ~~function to do ... ~~
}
\description{ Generates a p-rep row column design for plant breeding field trials. The function takes various parameters to create a matrix representing the layout of unreplicated treatments,checks, and replicated treatments across experimental blocks. The function randomizes treatment assignments, ensures proper distribution of treatments within blocks, and corrects for any over-representation of replicated treatments. Finally, it generates a field plan dataframe with information about block, row, column, plot number, treatment assignment, and trial identifier.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\code{p_rep_design}{ Generates a p-rep row column design for plant breeding field trials. The function takes various parameters to create a matrix representing the layout of unreplicated treatments,checks, and replicated treatments across experimental blocks. The function randomizes treatment assignments, ensures proper distribution of treatments within blocks, and corrects for any over-representation of replicated treatments. Finally, it generates a field plan dataframe with information about block, row, column, plot number, treatment assignment, and trial identifier.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p_rep_design <- function(trt_data = NULL,
                         rep_trt_data = NULL,
                         check_data = NULL,
                         nblock = 0,
                         test_trt = c(0, 1),
                         p_rep_trt = c(0,1),
                         check_trt = c(0, 1),
                         rb = 0, cb = 0,
                         trial = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt_data}{Vector of unique treatment names.
%%     ~~Describe \code{x} here~~
}
  \item{rep_trt_data}{Vector of unique p-rep treatment names.
%%     ~~Describe \code{x} here~~
}
  \item{check_data}{Vector of unique check names.
%%     ~~Describe \code{x} here~~
}
  \item{nblock}{Number of blocks in the experimental design.
%%     ~~Describe \code{x} here~~
}
  \item{test_trt}{Vector containing the number of unreplicated test treatments.
%%     ~~Describe \code{x} here~~
}
  \item{p_rep_trt}{Vector containing the number of p-rep treatments and their replications.
%%     ~~Describe \code{x} here~~
}
  \item{check_trt}{Vector containing the number of control treatments and their replications.
%%     ~~Describe \code{x} here~~
}
  \item{rb}{Number of rows per block.
%%     ~~Describe \code{x} here~~
}
  \item{cb}{Number of columns per block.
%%     ~~Describe \code{x} here~~
}
  \item{trial}{Trial identifier (optional, can be used to differentiate between multiple trials).
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
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\code{\code{ \link{Gmatrix}}, \code{\link{Rmatrix}, \code{\link{SwapMethods}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
p_rep_design <- function(trt_data = NULL,
                         rep_trt_data = NULL,
                         check_data = NULL,
                         nblock = 0,
                         test_trt = c(0, 1),
                         p_rep_trt = c(0,1),
                         check_trt = c(0, 1),
                         rb = 0, cb = 0,
                         trial = NULL) {

  # Begin the process of generating the field plan and alert the user
  cat("Starting field plan generation...\n")

  # Extract parameters for main test, p-rep test, and control treatments from the provided vectors
  test <- test_trt[1]
  nrep_test <- test_trt[2]

  rep_test <- p_rep_trt[1]
  prep_test <- p_rep_trt[2]

  control <- check_trt[1]
  nrep_control <- check_trt[2]

  # Check if the design dimensions are consistent with the number of blocks, rows, and columns provided
  # The total number of plots within each block should match the product of rows and columns per block
  N1 <- rb * cb # Total number of plots in a single block
  N2 <- control + (test + rep_test * prep_test) / nblock # Number of plots per block based on treatment and control replications

  if (N1 != N2) {
    stop("Design has problems with its parameters. Verify the dimensions of the blocks versus the number of plots derived from treatments and controls.")
  }

  cat("Verifying input dimensions...\n")

  # Randomize the order of treatment, p-rep treatment, and control names to avoid systematic bias
  cat("Randomizing unreplicated treatment names...\n")
  trt_names <- as.character(trt_data)
  trt_names <- sample(trt_names, size = test, replace = FALSE)

  cat("Randomizing p-rep treatment names...\n")
  rep_test_names <- sample(as.character(rep_trt_data), replace = FALSE)

  cat("Randomizing ckeck names...\n")
  control_names <- sample(as.character(check_data), replace = FALSE)

  # Create matrices to organize the plot data for treatments, p-reps, and controls within each block
  cat("Creating plot matrices...\n")
  trt_plot <- matrix(rep(trt_names, nrep_test), ncol = nblock, byrow = TRUE)
  rep_test_plot <- matrix(rep(rep_test_names, prep_test), ncol = nblock, byrow = TRUE)
  control_plot <- matrix(rep(control_names, nrep_control), nrow = nblock, byrow = TRUE)

  # Within each block, shuffle the treatments, p-reps, and controls to ensure a random distribution
  for (i in 1:nblock) {

    # Ensure rep_test is not repeated within block
    rep_test_plot[, i] <- sample(rep_test_names, prep_test, replace = FALSE)

    # Ensure ckeck is not repeated within block
    control_plot[i, ] <- sample(control_names, control, replace = FALSE)
  }

  # Combine the treatment assignment matrices and randomly shuffle the treatments across all blocks to create the overall experimental layout
  new_mat <- rbind(t(control_plot), rep_test_plot, trt_plot)
  trts <- apply(new_mat, 2, sample)
  treatment <- matrix(trts, ncol = 1, byrow = FALSE)

  # Ensure that no p-rep treatment is repeated more than the allowed number of times across the entire experimental design
  cat("Correcting for any over-representation of p-rep treatments...\n")
  rep_test_counts <- table(treatment)[rep_test_names]

  # Loop to correct any over-representation by replacing excess occurrences with other p-rep treatments
  while (any(rep_test_counts > prep_test, na.rm = TRUE)) {
    idx_repeat <- which(rep_test_counts > prep_test)
    for (i in idx_repeat) {
      repeat_indices <- which(treatment == rep_test_names[i])
      to_replace <- sample(repeat_indices, rep_test_counts[i] - prep_test, replace = FALSE)

      # Sample replacement treatments only from rep_test_names
      replacement_treatments <- sample(rep_test_names, length(to_replace), replace = TRUE)

      treatment[to_replace] <- replacement_treatments
    }
    rep_test_counts <- table(treatment)[rep_test_names]
  }

  # Generating block, row, column, and plot numbers
  cat("Generating block, row, column, and plot numbers...\n")
  total_plots = length(treatment)
  plots_per_block = rb * cb
  number_of_blocks = ceiling(total_plots / plots_per_block)
  rep_num <- rep(1:number_of_blocks, each = plots_per_block)

  # Truncating rep numbers to total_plots length
  rep_num <- rep_num[1:total_plots]

  # Generating row numbers
  r <- rep(1:rb, each = cb)
  row_num <- r

  while(length(row_num) < total_plots) {
    r <- r + rb
    row_num <- c(row_num, r)
  }

  # Truncating row numbers to total_plots length
  row_num <- row_num[1:total_plots]

  # Generating column numbers
  col_num <- rep(rep(1:cb, each = rb), times = ceiling(total_plots / (cb * rb)))
  # Truncating column numbers to total_plots length
  col_num <- col_num[1:total_plots]

  # Generating plot numbers and adding trial
  plot_num <- 1:total_plots
  trial_num <- rep(trial, length.out = nrow(treatment))

  # Combining row and column numbers
  coord <- cbind(Row = rep(row_num, each = 1), Col = rep(col_num, each = 1))

  # Generating field plan dataframe
  matdf <- data.frame(coord, Rep = rep_num, Treatment = as.vector(treatment), Plot = plot_num, Trial = trial_num)

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
