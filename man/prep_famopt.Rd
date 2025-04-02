\name{prep_famopt}
\alias{prep_famopt}
\title{
Create a Partially Replicated (p-rep) Experimental Design Layout
}
\description{
This function generates a partially replicated (p-rep) experimental design for field trials, ensuring that check treatments are balanced across blocks,
p-rep (partially replicated) and unreplicated treatments are distributed evenly, and similar families of treatments are not placed in adjacent positions.
It also supports an optional serpentine layout for plot numbering and can auto-correct layout dimensions.
}
\usage{
prep_famopt(check_treatments, check_families, p_rep_treatments, p_rep_reps,
            p_rep_families, unreplicated_treatments, unreplicated_families,
            n_blocks, n_rows, n_cols, serpentine = FALSE, seed = NULL,
            attempts = 100, warn_and_correct = TRUE)
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

  \item{serpentine}{A logical value indicating whether to use a serpentine plot layout.
  If \code{TRUE}, plot numbering (and filling of treatments) follows a snake-like pattern: left-to-right in the first row, right-to-left in the second row, and so on.
  If \code{FALSE}, plot numbering proceeds left-to-right for every row.}

  \item{seed}{An optional integer for setting the random seed to ensure reproducibility of the design.}

  \item{attempts}{An integer specifying the maximum number of iterations for optimizing family distribution across blocks.}

  \item{warn_and_correct}{A logical value indicating whether to automatically adjust \code{n_rows} and \code{n_cols}
  if the provided layout dimensions do not match the required number of plots.
  If \code{TRUE}, the function will issue a warning and adjust the dimensions to fit the exact number of plots needed.
  If \code{FALSE}, the function will stop execution and require the user to manually correct the input.}
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
\examples{## Create a p-rep design with 2 checks, 1 p-rep, and 1 unreplicated treatment in 3 blocks.
 check_treatments <- c("CheckA", "CheckB")
 check_families <- c("FamCheck1", "FamCheck2")   # Families for checks
 p_rep_treatments <- c("T1")                     # One partially replicated treatment
 p_rep_reps <- c(2)                              # T1 will appear in 2 out of 3 blocks
 p_rep_families <- c("Fam1")                     # Family for T1
 unreplicated_treatments <- c("U1")              # One unreplicated treatment
 unreplicated_families <- c("Fam2")              # Family for U1
 n_blocks <- 3
 n_rows <- 3
 n_cols <- 3
 set.seed(42)
 design <- prep_famopt(check_treatments, check_families,
                              p_rep_treatments, p_rep_reps, p_rep_families,
                              unreplicated_treatments, unreplicated_families,
                              n_blocks, n_rows, n_cols, serpentine = TRUE, seed = NULL, attempts = 1000)
print(design)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
prep_famopt <- function(check_treatments, check_families,
                        p_rep_treatments, p_rep_reps, p_rep_families,
                        unreplicated_treatments, unreplicated_families,
                        n_blocks, n_rows, n_cols,
                        serpentine = FALSE, seed = NULL, attempts = 100, warn_and_correct = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (length(check_treatments) != length(check_families)) {
    stop("Length of check_families must match length of check_treatments.")
  }
  if (length(p_rep_treatments) != length(p_rep_reps) || length(p_rep_treatments) != length(p_rep_families)) {
    stop("Lengths of p_rep_treatments, p_rep_reps, and p_rep_families must all match.")
  }
  if (length(unreplicated_treatments) != length(unreplicated_families)) {
    stop("Length of unreplicated_families must match length of unreplicated_treatments.")
  }
  if (any(p_rep_reps > n_blocks)) {
    stop("Each p-rep treatment's replication count must not exceed the number of blocks.")
  }
  if (n_blocks < 1) {
    stop("n_blocks must be at least 1.")
  }
  
  total_checks <- n_blocks * length(check_treatments)
  total_prep <- sum(p_rep_reps)
  total_unrep <- length(unreplicated_treatments)
  total_required <- total_checks + total_prep + total_unrep
  
  field_size <- n_rows * n_cols
  if (field_size != total_required) {
    if (warn_and_correct) {
      warning(paste0("Field size (", n_rows, " rows x ", n_cols, " cols = ",
                     field_size, " plots) does not match required (", total_required, "). Adjusting layout dimensions..."))
      n_cols <- ceiling(total_required / n_rows)
      field_size <- n_rows * n_cols
    } else {
      stop(paste0("Provided field dimensions (", field_size, " plots) do not match the required (",
                  total_required, "). Adjust `n_rows` or `n_cols`.")) 
    }
  }
  
  # Create a lookup mapping for treatments to families
  family_lookup <- setNames(
    c(check_families, p_rep_families, unreplicated_families),
    c(check_treatments, p_rep_treatments, unreplicated_treatments)
  )
  
  # Initialize empty block list
  blocks <- vector("list", n_blocks)
  
  # Assign p-rep treatments to unique blocks
  p_rep_assignments <- vector("list", length(p_rep_treatments))
  names(p_rep_assignments) <- p_rep_treatments
  
  available_blocks <- vector("list", length(p_rep_treatments))
  for (i in seq_along(p_rep_treatments)) {
    available_blocks[[i]] <- sample(seq_len(n_blocks))  # Random block order
  }
  
  for (i in seq_along(p_rep_treatments)) {
    treatment <- p_rep_treatments[i]
    num_reps <- p_rep_reps[i]
    
    if (num_reps > length(available_blocks[[i]])) {
      stop(paste0("Error: Not enough unique blocks available for p-rep treatment '", treatment, "'.")) 
    }
    
    assigned_blocks <- sample(available_blocks[[i]], num_reps)  # Select unique blocks
    available_blocks[[i]] <- setdiff(available_blocks[[i]], assigned_blocks)  # Remove used blocks
    
    p_rep_assignments[[treatment]] <- assigned_blocks
  }
  
  # Shuffle and distribute unreplicated treatments across blocks
  unrep_treatments_shuffled <- sample(unreplicated_treatments)
  block_unrep_list <- split(unrep_treatments_shuffled, rep(seq_len(n_blocks), length.out = length(unrep_treatments_shuffled)))
  
  # Populate blocks
  for (b in seq_len(n_blocks)) {
    block_treatments <- c(check_treatments)  # Always include all checks in each block
    
    # Add p-rep treatments assigned to this block
    for (treatment in names(p_rep_assignments)) {
      if (is.element(b, p_rep_assignments[[treatment]])) {
        block_treatments <- c(block_treatments, treatment)
      }
    }
    
    # Add unreplicated treatments assigned to this block
    if (b <= length(block_unrep_list)) {
      block_treatments <- c(block_treatments, block_unrep_list[[b]])
    }
    
    # Assign correct family using the lookup vector
    block_families <- sapply(block_treatments, function(trt) family_lookup[trt])
    
    # Ensure no two adjacent treatments belong to the same family
    valid_order <- FALSE
    max_attempts <- attempts
    attempt <- 1
    
    while (!valid_order && attempt <= max_attempts) {
      block_order <- sample(seq_along(block_treatments))
      block_shuffled <- block_treatments[block_order]
      block_fam_shuffled <- block_families[block_order]
      
      if (!any(block_fam_shuffled[-1] == block_fam_shuffled[-length(block_fam_shuffled)])) {
        valid_order <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!valid_order) {
      warning(paste0("Could not completely avoid adjacent families in block ", b,
                     " after ", max_attempts, " attempts.")) 
    }
    
    blocks[[b]] <- data.frame(Treatment = block_shuffled, Family = block_fam_shuffled, Block = b)
  }
  
  final_data <- do.call(rbind, blocks)
  final_data$Plot <- seq_len(nrow(final_data))
  final_data$Row <- (final_data$Plot - 1) %/% n_cols + 1
  final_data$Column <- (final_data$Plot - 1) %% n_cols + 1
  
  # Apply serpentine layout if required
  if (serpentine) {
    for (r in seq_len(n_rows)) {
      if (r %% 2 == 0) {
        final_data$Column[final_data$Row == r] <- rev(final_data$Column[final_data$Row == r])
      }
    }
  }
  
  # Create layout matrix in correct serpentine order
  layout_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(nrow(final_data))) {
    r <- final_data$Row[i]
    c <- final_data$Column[i]
    
    # Apply serpentine logic directly to column position
    actual_col <- if (serpentine && (r %% 2 == 0)) {
      n_cols - c + 1  # reverse column if even row
    } else {
      c
    }
    
    layout_matrix[r, actual_col] <- final_data$Treatment[i]
  }
  
  return(list(layout_matrix = layout_matrix, field_book = final_data))
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
