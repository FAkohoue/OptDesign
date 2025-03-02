#' Define a function to create a field plan for agricultural experiments with treatments, controls, and p-reps
#'
#' @param trt_data: Vector of unique treatment names
#' @param rep_trt_data: Vector of unique p-rep treatment names
#' @param check_data: Vector of unique check names
#' @param nblock: Number of blocks in the experimental design
#' @param test_trt: Vector containing the number of unreplicated test treatments
#' @param p_rep_trt: Vector containing the number of p-rep treatments and their replications
#' @param check_trt: Vector containing the number of control treatments and their replications
#' @param rb: Number of rows per block
#' @param cb: Number of columns per block
#' @param trial: Trial identifier (optional, can be used to differentiate between multiple trials)
#'
#' @return A matrix of field plan generated successfully
#' @export
#'
#' @examples
#' trt_names <- as.character(1:100)
#' controls <- c('Control1', 'Control2', 'Control3', 'Control4')
#' reptrt <- as.character(101:110)
#' N1 <- 4 * 4
#' N2 <- 4*1 + (100*1 + 10*2)/10
#' p_rep_design(trt_data = trt_names,
#'           rep_trt_data = reptrt,
#'           check_data = controls,
#'           nblock = 10,
#'          test_trt = c(100,1),
#'          p_rep_trt = c(10,2),
#'          check_trt = c(4,10),
#'          rb = 4, cb = 4,
#'          trial = 'TOL')
#'
p_rep_design <- function(trt_data = NULL,
                         rep_trt_data = NULL,
                         check_data = NULL,
                         nblock = 0,
                         test_trt = c(0, 1),
                         p_rep_trt = c(0,1),
                         check_trt = c(0, 1),
                         rb = 0, cb = 0,
                         trial = NULL) {

  cat("Starting field plan generation...\n")

  # Validate input lengths
  if (length(trt_data) != test_trt[1]) {
    stop("Length of trt_data does not match the number of test treatments specified in test_trt[1].")
  }
  if (length(rep_trt_data) != p_rep_trt[1]) {
    stop("Length of rep_trt_data does not match the number of p-rep treatments specified in p_rep_trt[1].")
  }
  if (length(check_data) != check_trt[1]) {
    stop("Length of check_data does not match the number of controls specified in check_trt[1].")
  }

  # Extract parameters
  test <- test_trt[1]
  nrep_test <- test_trt[2]

  rep_test <- p_rep_trt[1]
  prep_test <- p_rep_trt[2]

  control <- check_trt[1]
  nrep_control <- check_trt[2]

  # Check design dimensions
  N1 <- rb * cb
  N2 <- control + (test * nrep_test + rep_test * prep_test) / nblock

  if (N1 != N2) {
    stop("Design parameters inconsistent. Verify block dimensions versus treatment allocations.")
  }

  cat("Verifying input dimensions...\n")

  # Randomize treatment names
  cat("Randomizing treatment names...\n")
  trt_names <- sample(trt_data)

  # Generate p-rep treatments ensuring exact replication
  cat("Assigning p-rep treatments...\n")
  rep_test_names <- sample(rep_trt_data)
  p_rep_vector <- rep(rep_test_names, each = prep_test)
  p_rep_shuffled <- sample(p_rep_vector)
  rep_test_plot <- matrix(p_rep_shuffled, nrow = prep_test, ncol = nblock)

  # Randomize controls and create control matrix
  cat("Assigning control treatments...\n")
  control_names <- sample(check_data)
  control_plot <- matrix(rep(control_names, nrep_control),
                         nrow = nblock, byrow = TRUE)
  for (i in 1:nblock) {
    control_plot[i, ] <- sample(control_plot[i, ])
  }

  # Create test treatment matrix
  cat("Assigning test treatments...\n")
  trt_plot <- matrix(rep(trt_names, nrep_test),
                     ncol = nblock, byrow = TRUE)

  # Combine components and shuffle
  cat("Combining components and randomizing...\n")
  new_mat <- rbind(t(control_plot), rep_test_plot, trt_plot)
  trts <- apply(new_mat, 2, sample)
  treatment <- as.vector(trts)

  # Generate field layout
  cat("Generating field layout...\n")
  total_plots <- length(treatment)
  plots_per_block <- rb * cb
  number_of_blocks <- ceiling(total_plots / plots_per_block)

  # Generate block IDs
  rep_num <- rep(1:number_of_blocks, each = plots_per_block)[1:total_plots]

  # Generate row numbers
  row_num <- rep(1:rb, each = cb, times = number_of_blocks)[1:total_plots]

  # Generate column numbers
  col_num <- rep(1:cb, times = rb * number_of_blocks)[1:total_plots]

  # Create final data frame
  matdf <- data.frame(
    Row = row_num,
    Col = col_num,
    Rep = rep_num,
    Treatment = treatment,
    Plot = 1:total_plots,
    Trial = if (!is.null(trial)) trial else "Trial"
  )

  cat("Field plan generated successfully.\n")
  return(matdf)
}
