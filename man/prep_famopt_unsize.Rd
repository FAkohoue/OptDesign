\name{prep_famopt_unsize}
\alias{prep_famopt_unsize}
\title{
Create a Partially Replicated (p-rep) Experimental Design Layout with unequal column size to adapt to irregular field dimensions
}
\description{
This function generates a partially replicated (p-rep) experimental design with unequal column sizes for field trials, ensuring that check treatments are balanced across blocks,
p-rep (partially replicated) and unreplicated treatments are distributed evenly, and similar families of treatments are not placed in adjacent positions.
It also supports an optional serpentine layout for plot numbering and can auto-correct layout dimensions.
}
\usage{
prep_famopt_unsize(check_treatments, check_families,
                  p_rep_treatments, p_rep_reps, p_rep_families,
                  unreplicated_treatments, unreplicated_families,
                  n_blocks, n_rows, n_cols = NULL,
                  cols_per_row,
                  serpentine = TRUE, seed = NULL,
                  attempts = 1000, warn_and_correct = TRUE, fix_rows = TRUE)
}
\arguments{
  \item{check_treatments}{A character vector of check treatment identifiers. Each check treatment appears exactly once in each block, ensuring full replication across blocks.}

  \item{check_families}{A character vector of the same length as \code{check_treatments}, specifying the family for each check treatment.}

  \item{p_rep_treatments}{A character vector of partially replicated (p-rep) treatment identifiers. These treatments have replication counts less than or equal to the number of blocks.}

  \item{p_rep_reps}{An integer vector indicating the number of replicates for each treatment in \code{p_rep_treatments}.
  Must be the same length as \code{p_rep_treatments}, and each value should be \eqn{\leq} \code{n_blocks} (each p-rep treatment can appear at most once per block).}

  \item{p_rep_families}{A character vector of the same length as \code{p_rep_treatments}, specifying the family for each p-rep treatment.}

  \item{unreplicated_treatments}{A character vector of unreplicated treatment identifiers. Each unreplicated treatment appears exactly once in the entire experiment.}

  \item{unreplicated_families}{A character vector of the same length as \code{unreplicated_treatments}, specifying the family for each unreplicated treatment.}

  \item{n_blocks}{An integer specifying the number of blocks (replicates). Each check treatment will appear in every block exactly once. \code{n_blocks} should be at least 1.}

  \item{n_rows}{An integer specifying the number of rows in the field layout.}

  \item{n_cols}{An integer specifying the number of columns in the field layout.}

  \item{cols_per_row}{A vector indicating the size of each columns according to the field reality. \code{cols_per_row} is required only when \code{n_cols} is set to \code{NULL}.}

  \item{serpentine}{A logical value indicating whether to use a serpentine plot layout.
  If \code{TRUE}, plot numbering (and filling of treatments) follows a snake-like pattern: left-to-right in the first row, right-to-left in the second row, and so on.
  If \code{FALSE}, plot numbering proceeds left-to-right for every row.}

  \item{seed}{An optional integer for setting the random seed to ensure reproducibility of the design.}

  \item{attempts}{An integer specifying the maximum number of iterations for optimizing family distribution across blocks.}

  \item{warn_and_correct}{A logical value indicating whether to automatically adjust \code{n_rows} and \code{n_cols}
  if the provided layout dimensions do not match the required number of plots.
  If \code{TRUE}, the function will issue a warning and adjust the dimensions to fit the exact number of plots needed.
  If \code{FALSE}, the function will stop execution and require the user to manually correct the input.}

  \item{fix_rows}{Logical indicating whether to automatically fix the number of rows and vary the number of columns when \code{warn_and_correct} is \code{TRUE}
  If \code{TRUE}, the function will adjust the adjust the number of columns to fit the exact number of plots needed while keeping the number of rows at the indicated value.
  If \code{FALSE}, both rows and columns can vary when applying \code{warn_and_correct}.}
}
\details{
\strong{Design Constraints and Features:}
\itemize{
  \item \strong{Balanced Check Treatments:} Each check treatment is placed exactly once in each block. The checks are randomized across blocks such that every block contains all check treatments exactly once, ensuring balance of checks across blocks.

  \item \strong{Partially Replicated (p-rep) Treatments:} P-rep treatments have replication counts less than or equal to the number of blocks. The design ensures no p-rep treatment appears more than once within the same block. Replicates of each p-rep treatment are spread across different blocks as evenly as possible (e.g., a treatment with 2 replicates in a 4-block design will be present in 2 of the 4 blocks, chosen at random).

  \item \strong{Unreplicated Treatments:} Each treatment listed in \code{unreplicated_treatments} is included exactly once in the entire experiment (one replicate total). These unreplicated entries are assigned to blocks such that they are spread across blocks rather than clustered. The function attempts to distribute unreplicated treatments across distinct blocks, so no block is overloaded with unreplicated entries while others have none (unless the number of unreplicated treatments is less than the number of blocks, in which case each unreplicated treatment still goes into a separate block).

  \item \strong{Avoid Adjacent Similar Families:} Within each block, treatments from the same family are not placed in adjacent plot positions. This means that when arranging the sequence of treatments within a block, the function ensures that no two consecutive plots (in the block's planting order) contain treatments from the same family. If a block's composition makes this impossible (e.g., too many treatments from one family in a single block), the function will issue a warning and attempt to minimize family adjacency.

  \item \strong{Spatial Optimization:} The randomization is done in a way that reduces clustering of treatments from the same family across the field. Each block is filled with a random order of treatments (respecting the above constraints), and blocks can be arranged contiguously in the field. The user-specified \code{serpentine} option can influence the neighbor relationships of plots.

  \item \strong{Serpentine Layout:} If \code{serpentine = TRUE}, the assignment of treatments to the field follows a serpentine path. This means that odd-numbered rows are filled left-to-right, while even-numbered rows are filled right-to-left. If \code{serpentine = FALSE}, each row is filled left-to-right in the same orientation.

  \item \strong{Reproducibility via Seed:} By specifying the \code{seed} parameter, users can obtain the same randomization and layout across runs.

  \item \strong{Automatic Layout Correction:} If the provided \code{n_rows} and \code{n_cols} do not exactly accommodate all the required plots, the function can adjust the layout. When \code{warn_and_correct = TRUE}, the function will automatically modify \code{n_rows} and/or \code{n_cols} to ensure that \code{n_rows * n_cols} equals the total number of plots required.

  \item \strong{Unequal column size:} If \code{cols_per_row} is defined, and \code{n_cols} is \code{NULL}, the function automatically allows for unequal column sizes.
}

\strong{Input Validation:}
The function checks that inputs are consistent:
\itemize{
  \item Verifies that the length of family vectors matches the corresponding treatment vectors.
  \item Ensures replication counts are valid (e.g., no p-rep replication exceeds \code{n_blocks}).
  \item Confirms that the total number of plots required matches the field dimensions (adjusting or warning as needed).
}

\strong{Output Format:}
The function returns a list with two components:
\enumerate{
  \item \strong{layout_matrix}: A matrix (of size \code{n_rows} x \code{n_cols}) representing the field layout. Each entry is the Treatment identifier assigned to that plot position.
  \item \strong{field_book}: A data frame (tibble) with columns \code{Plot}, \code{Block}, \code{Row}, \code{Column}, \code{Treatment}, and \code{Family}. Each row corresponds to a plot in the field, sorted by \code{Plot} order.
}
}

\value{
A list with two elements:
\describe{
   \item{layout_matrix}{A matrix of size \code{n_rows} x \code{n_cols} containing treatment identifiers arranged according to the design layout.}
   \item{field_book}{A data frame with columns \code{Plot}, \code{Block}, \code{Row}, \code{Column}, \code{Treatment}, \code{Family} describing each plot in the field layout.}
}
}

\references{Félicien Akohoue
%% ~put references to the literature/web site here ~
}
\author{Félicien Akohoue
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{## Create a p-rep design with 3 checks, 3 p-rep, and 6 unreplicated treatment in 10 blocks. Columns have unequal sizes.
check_treatments <- c("CheckA", "CheckB", "CheckC")
check_families <- c("FamCheck1", "FamCheck2", "FamCheck3")   # Families for checks
p_rep_treatments <- c("T1", "T2", "T3")                    # One partially replicated treatment
p_rep_reps <- c(8, 8, 8)                              # T1 will appear in 2 out of 3 blocks
p_rep_families <- c("Fam1", "Fam2", "Fam2")   # Family for T1
unreplicated_treatments <- c("U1", "U2", "U3", "U4", "U5", "U6")              # One unreplicated treatment
unreplicated_families <- c("Fam1", "Fam2", "Fam3", "Fam4", "Fam5", "Fam6")              # Family for U1
n_blocks <- 10
n_rows <- 5

cols_per_row <- c(10, 15, 12, 10, 13)
design <- prep_famopt_unsize(check_treatments, check_families,
                              p_rep_treatments, p_rep_reps, p_rep_families,
                              unreplicated_treatments, unreplicated_families,
                              n_blocks = n_blocks, n_rows = n_rows, n_cols = NULL,
                              cols_per_row = cols_per_row,
                              serpentine = TRUE, seed = NULL,
                              attempts = 1000, warn_and_correct = TRUE, fix_rows = TRUE)
#' # --- Create layout matrix ---
ncols <- max(design$field_book$Column)
layout_matrix <- matrix(NA, nrow = 5, ncol = ncols)
for (i in seq_len(nrow(combined_fieldbook))) {
  r <- design$field_book$Row[i]
  c <- design$field_book$Column[i]
  layout_matrix[r, c] <- design$field_book$Treatment[i]
}
print(layout_matrix)

#' # --- View the field book (first few entries) ---
head(design$field_book)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
prep_famopt_unsize <- function(check_treatments, check_families,
                               p_rep_treatments, p_rep_reps, p_rep_families,
                               unreplicated_treatments, unreplicated_families,
                               n_blocks, n_rows = NULL, n_cols = NULL,
                               cols_per_row = NULL,
                               serpentine = FALSE, seed = NULL,
                               attempts = 1000, warn_and_correct = TRUE, fix_rows = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  # Validate input lengths
  if (length(check_treatments) != length(check_families)) {
    stop("Length of check_families must match length of check_treatments.")
  }
  if (length(p_rep_treatments) != length(p_rep_reps) ||
      length(p_rep_treatments) != length(p_rep_families)) {
    stop("Lengths of p_rep_treatments, p_rep_reps, and p_rep_families must all match.")
  }
  if (length(unreplicated_treatments) != length(unreplicated_families)) {
    stop("Length of unreplicated_families must match length of unreplicated_treatments.")
  }
  if (any(p_rep_reps > n_blocks)) {
    stop("Each p-rep treatment's replication count must not exceed the number of blocks.")
  }

  # Create family lookup table
  family_lookup <- setNames(
    c(check_families, p_rep_families, unreplicated_families),
    c(check_treatments, p_rep_treatments, unreplicated_treatments)
  )

  # Calculate total treatments
  total_checks <- n_blocks * length(check_treatments)
  total_prep <- sum(p_rep_reps)
  total_unrep <- length(unreplicated_treatments)
  total_required <- total_checks + total_prep + total_unrep

  # Determine field layout
  if (!is.null(cols_per_row)) {
    n_rows <- length(cols_per_row)
    field_size <- sum(cols_per_row)
  } else {
    if (is.null(n_rows) || is.null(n_cols)) {
      stop("If cols_per_row is not provided, both n_rows and n_cols must be specified.")
    }
    field_size <- n_rows * n_cols
    cols_per_row <- rep(n_cols, n_rows)
  }

  # Adjust layout if needed
  if (field_size != total_required) {
    if (warn_and_correct && is.null(cols_per_row)) {
      warning(paste0("Field size (", field_size, ") does not match required (", total_required, "). Adjusting layout..."))
      if (fix_rows) {
        n_cols <- ceiling(total_required / n_rows)
      } else {
        n_rows <- ceiling(total_required / n_cols)
      }
      cols_per_row <- rep(n_cols, n_rows)
      field_size <- sum(cols_per_row)
    } else {
      stop(paste0("Field size (", field_size, " plots) does not match required (", total_required, ")."))
    }
  }

  # Assign p-rep treatments to blocks
  p_rep_assignments <- vector("list", length(p_rep_treatments))
  names(p_rep_assignments) <- p_rep_treatments
  available_blocks <- lapply(seq_along(p_rep_treatments), function(i) sample(seq_len(n_blocks)))

  for (i in seq_along(p_rep_treatments)) {
    treatment <- p_rep_treatments[i]
    reps <- p_rep_reps[i]
    if (reps > length(available_blocks[[i]])) {
      stop(paste0("Not enough blocks for p-rep treatment '", treatment, "'"))
    }
    assigned <- sample(available_blocks[[i]], reps)
    available_blocks[[i]] <- setdiff(available_blocks[[i]], assigned)
    p_rep_assignments[[treatment]] <- assigned
  }

  # Distribute unreplicated treatments
  unrep_shuffled <- sample(unreplicated_treatments)
  block_unrep_list <- split(unrep_shuffled, rep(seq_len(n_blocks), length.out = length(unrep_shuffled)))

  # Block assembly function with check randomization
  assemble_block <- function(check_trts, p_rep_trts, unrep_trts, block_num) {
    valid_block <- FALSE
    attempt <- 1
    n_checks <- length(check_trts)

    while (!valid_block && attempt <= attempts) {
      # Combine and shuffle all treatments
      all_trts <- sample(c(check_trts, p_rep_trts, unrep_trts))

      # Check family separation
      families <- family_lookup[all_trts]
      family_ok <- !any(families[-1] == families[-length(families)])

      # Check distribution
      check_pos <- which(!is.na(match(all_trts, check_trts)))
      spacing_ok <- if (n_checks > 1) {
        min_spacing <- floor(length(all_trts) / (n_checks + 1))
        all(diff(check_pos) >= min_spacing)
      } else TRUE

      valid_block <- family_ok && spacing_ok
      attempt <- attempt + 1
    }

    if (!valid_block) {
      warning("Could not find ideal arrangement for checks in block ", block_num)
    }

    data.frame(
      Treatment = all_trts,
      Family = families,
      Block = block_num,
      stringsAsFactors = FALSE
    )
  }

  # Build all blocks
  blocks <- vector("list", n_blocks)
  for (b in seq_len(n_blocks)) {
  p_rep_in_block <- names(p_rep_assignments)[sapply(p_rep_assignments, function(x) any(x == b))]
  unrep_in_block <- if (b <= length(block_unrep_list)) sample(block_unrep_list[[b]]) else character(0)

  blocks[[b]] <- assemble_block(
    check_trts = check_treatments,
    p_rep_trts = p_rep_in_block,
    unrep_trts = unrep_in_block,
    block_num = b
  )
}

  # Combine all blocks
  final_data <- do.call(rbind, blocks)
  final_data$Plot <- seq_len(nrow(final_data))

  # Assign row/column positions
  get_row_col <- function(plot_num, cols_per_row) {
    cum_counts <- cumsum(cols_per_row)
    row <- which(plot_num <= cum_counts)[1]
    col <- plot_num - ifelse(row == 1, 0, cum_counts[row - 1])
    c(row, col)
  }

  pos <- t(sapply(final_data$Plot, get_row_col, cols_per_row))
  final_data$Row <- pos[, 1]
  final_data$Column <- pos[, 2]

  # Apply serpentine pattern
  if (serpentine) {
    for (r in seq_along(cols_per_row)) {
      idx <- which(final_data$Row == r)
      if (mod(r, 2) == 0) {
        final_data$Column[idx] <- rev(final_data$Column[idx])
      }
    }
  }

  # Create layout matrix
  layout_matrix <- lapply(seq_len(n_rows), function(r) {
    row_data <- final_data[final_data$Row == r, ]
    row_data <- row_data[order(row_data$Column), ]
    row_data$Treatment
  })

  return(list(
    layout_matrix = layout_matrix,
    field_book = final_data
  ))
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
