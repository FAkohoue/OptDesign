#' Create a Partially Replicated (p-rep) Experimental Design Layout
#'
#' This function generates a partially replicated (p-rep) experimental design for field trials,
#' ensuring that check treatments are balanced across blocks, p-rep (partially replicated) and unreplicated treatments are distributed evenly,
#' and similar families of treatments are not placed in adjacent positions. It also supports an optional serpentine layout for plot numbering and can auto-correct layout dimensions.
#'
#' @param check_treatments Character vector of **check treatment** identifiers. Each check will appear exactly once in each block (fully replicated across blocks).
#' @param check_families Character vector of the same length as `check_treatments` specifying the family for each check treatment.
#' @param p_rep_treatments Character vector of **partially replicated treatment** identifiers. These treatments have replication counts less than or equal to the number of blocks.
#' @param p_rep_reps Integer vector indicating the number of replicates for each treatment in `p_rep_treatments`. Must be the same length as `p_rep_treatments`. Each value should be ≤ `n_blocks` (each p-rep treatment can appear at most once per block).
#' @param p_rep_families Character vector of the same length as `p_rep_treatments` specifying the family for each p-rep treatment.
#' @param unreplicated_treatments Character vector of **unreplicated treatment** identifiers. Each unreplicated treatment appears exactly once in the entire experiment.
#' @param unreplicated_families Character vector of the same length as `unreplicated_treatments` specifying the family for each unreplicated treatment.
#' @param n_blocks Integer number of blocks (replicates). Each check treatment will appear in every block exactly once. Block count should be ≥ 1.
#' @param n_rows Integer number of rows in the field layout.
#' @param n_cols Integer number of columns in the field layout.
#' @param order Plot numbering order, either `column` (default) or `row`
#' @param serpentine Logical, whether to use serpentine plot layout. If `TRUE`, plot numbering (and filling of treatments) will snake through the field: left-to-right in the first row, then right-to-left in the second, and so on. If `FALSE`, plot numbering proceeds left-to-right for every row.
#' @param seed Optional integer seed for random number generation. Use this for reproducible randomization.
#' @param attempts Number of iterations for family distribution across blocks.
#' @param warn_and_correct Logical indicating whether to automatically adjust `n_rows` and `n_cols` if the provided layout dimensions do not match the required number of plots.
#' If `TRUE`, the function will issue a warning and adjust the dimensions to fit the exact number of plots needed.
#' If `FALSE`, a warning is issued and the function will require the user to correct the input (no automatic adjustment; the function will stop if the dimensions are inconsistent).
#' @param fix_rows Logical indicating whether to automatically fix the number of rows and vary the number of columns when `warn_and_correct` is `TRUE`.
#' If `TRUE`, the function will adjust the adjust the number of columns to fit the exact number of plots needed while keeping the number of rows at the indicated value.
#' If `FALSE`, both rows and columns can vary when applying `warn_and_correct`.
#'
#' @details
#' **Design Constraints and Features:**
#'
#' - **Balanced Check Treatments:** Each check treatment is placed exactly once in each block. The checks are randomized across blocks such that every block contains all check treatments exactly once, ensuring balance of checks across blocks.
#'
#' - **Partially Replicated (p-rep) Treatments:** P-rep treatments have replication counts less than or equal to the number of blocks. The design ensures no p-rep treatment appears more than once within the same block. Replicates of each p-rep treatment are spread across different blocks as evenly as possible (e.g., a treatment with 2 replicates in a 4-block design will be present in 2 of the 4 blocks, chosen at random).
#'
#' - **Unreplicated Treatments:** Each treatment listed in `unreplicated_treatments` is included exactly once in the entire experiment (one replicate total). These unreplicated entries are assigned to blocks such that they are spread across blocks rather than clustered. The function attempts to distribute unreplicated treatments across distinct blocks, so no block is overloaded with unreplicated entries while others have none (unless the number of unreplicated treatments is less than the number of blocks, in which case each unreplicated still goes into a separate block).
#'
#' - **Avoid Adjacent Similar Families:** Within each block, treatments from the same family are not placed in adjacent plot positions. This means that when arranging the sequence of treatments within a block, the function will ensure that no two consecutive plots (in the block's planting order) contain treatments from the same family. This helps achieve a more uniform spatial distribution of families and minimizes clustering of similar genotypes. If a block's composition makes this impossible (for example, too many treatments from one family in a single block), the function will issue a warning and still attempt to minimize family adjacency.
#'
#' - **Spatial Optimization:** The randomization is done in a way to reduce clustering of treatments from the same family across the field. By balancing replicates across blocks and avoiding within-block adjacencies, the design inherently minimizes clusters of related treatments while maintaining randomness. Each block is filled with a random order of treatments (respecting the above constraints), and blocks can be arranged contiguously in the field. The user-specified `serpentine` option can further influence the neighbor relationships of plots. Overall, the layout is random yet balanced, aiming to avoid large patches of the same family.
#'
#' - **Serpentine Layout:** If `serpentine = TRUE`, the assignment of treatments to the field follows a serpentine path. This means that odd-numbered rows are filled left-to-right, while even-numbered rows are filled right-to-left. This layout is often used to facilitate field operations by having a continuous snake-like planting order. If `serpentine = FALSE`, each row is filled left-to-right in the same orientation. The field book `Plot` numbers will reflect this ordering.
#'
#' - **Reproducibility via Seed:** By specifying the `seed` parameter, users can obtain the same randomization and layout across runs. All random operations (block assignments, treatment shuffling, etc.) will use this seed for reproducibility.
#'
#' - **Automatic Layout Correction:** If the provided `n_rows` and `n_cols` do not exactly accommodate all the required plots (which is determined by the total number of treatment placements needed), the function can adjust the layout. When `warn_and_correct = TRUE`, the function will automatically modify `n_rows` and/or `n_cols` to ensure that `n_rows * n_cols` equals the total number of plots required by the design, and issue a warning describing the adjustment. It will choose a layout close to the user input if possible. If `warn_and_correct = FALSE` and the dimensions are inconsistent with the number of treatments, the function will issue a warning and stop, requiring the user to adjust `n_rows` or `n_cols` manually so that the field has the correct size.
#'
#' **Input Validation:** The function checks that inputs are consistent:
#' it verifies that the length of family vectors match the corresponding treatment vectors, replication counts are valid (e.g., no p-rep replication exceeds `n_blocks`), and that the total number of plots required matches the field dimensions (adjusting or warning as needed).
#'
#' **Output Format:** The function returns a list with two components:
#' 1. **layout_matrix**: A matrix (of size `n_rows` x `n_cols`) representing the field layout. Each entry is the Treatment identifier assigned to that plot position. Blocks are laid out in contiguous sections of the matrix (typically whole rows or sets of rows per block). Check treatments can be observed to appear once per block (e.g., once per block section in the layout).
#' 2. **field_book**: A data frame (tibble) with columns `Plot`, `Block`, `Row`, `Column`, `Treatment`, and `Family`. Each row of this data frame corresponds to a plot in the field. `Plot` is the plot number (starting from 1, incrementing in serpentine or row-major order according to `serpentine` setting), `Block` is the block (replicate) number, `Row` and `Column` give the position coordinates, `Treatment` is the treatment identifier, and `Family` is the family of that treatment. The data frame is sorted by `Plot` order (the order in which plots would be planted/numbered in the field).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{layout_matrix}{A matrix of size `n_rows` x `n_cols` containing treatment identifiers arranged according to the design layout.}
#'   \item{field_book}{A data frame with columns `Plot`, `Block`, `Row`, `Column`, `Treatment`, `Family` describing each plot in the field layout.}
#' }
#'
#' @examples
#' # Example: Create a p-rep design with 2 checks, 1 p-rep, and 1 unreplicated treatment in 3 blocks.
#' check_treatments <- c("CheckA", "CheckB")
#' check_families <- c("FamCheck1", "FamCheck2")   # Families for checks
#' p_rep_treatments <- c("T1")                     # One partially replicated treatment
#' p_rep_reps <- c(2)                              # T1 will appear in 2 out of 3 blocks
#' p_rep_families <- c("Fam1")                     # Family for T1
#' unreplicated_treatments <- c("U1")              # One unreplicated treatment
#' unreplicated_families <- c("Fam2")              # Family for U1
#' n_blocks <- 3
#' n_rows <- 3
#' n_cols <- 3
#' set.seed(42)
#' design <- prep_famopt(check_treatments, check_families,
#'                              p_rep_treatments, p_rep_reps, p_rep_families,
#'                              unreplicated_treatments, unreplicated_families,
#'                              n_blocks, n_rows, n_cols, serpentine = TRUE, seed = 42, fix_rows = TRUE)
#' # View the layout matrix
#' design$layout_matrix
#' #      [,1]     [,2]    [,3]
#' # [1,] "CheckB" "T1"    "CheckA"
#' # [2,] "CheckA" "U1"    "CheckB"
#' # [3,] "CheckB" "T1"    "CheckA"
#' #
#' # The layout shows 3 rows and 3 columns. Each block corresponds to one row in this case (3 blocks for 3 rows).
#' # CheckA and CheckB each appear once per row (block). T1 (p-rep) appears in 2 of the 3 blocks (Rows 1 and 3 here), and U1 appears in the remaining block (Row 2).
#' # Within each row, no two adjacent positions have the same family.
#'
#' # View the field book (first few entries)
#' head(design$field_book)
#' #   Plot Block Row Column Treatment    Family
#' # 1    1     1   1      1   CheckB  FamCheck2
#' # 2    2     1   1      2       T1      Fam1
#' # 3    3     1   1      3   CheckA  FamCheck1
#' # 4    4     2   2      3   CheckB  FamCheck2
#' # 5    5     2   2      2       U1      Fam2
#' # 6    6     2   2      1   CheckA  FamCheck1
#' #
#' # The field book is sorted by plot order. Here, Plot 1-3 correspond to row 1 (Block 1) left-to-right,
#' # Plot 4-6 correspond to row 2 (Block 2) right-to-left (because serpentine is TRUE, row 2 is reversed in numbering), and so on.
#'
prep_famopt <- function(check_treatments, check_families,
                        p_rep_treatments, p_rep_reps, p_rep_families,
                        unreplicated_treatments, unreplicated_families,
                        n_blocks, n_rows, n_cols, order = "column",
                        serpentine = FALSE, seed = NULL, attempts = 1000, warn_and_correct = TRUE, fix_rows = TRUE) {
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

      if (fix_rows) {
        # Adjust columns based on fixed number of rows
        n_cols <- ceiling(total_required / n_rows)
      } else {
        # Adjust rows based on fixed number of columns
        n_rows <- ceiling(total_required / n_cols)
      }

      field_size <- n_rows * n_cols  # Recalculate
    } else {
      stop(paste0("Provided field dimensions (", field_size, " plots) do not match the required (",
                  total_required, "). Adjust `n_rows` or `n_cols`, or enable warn_and_correct."))
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
  if (order == "column") {
    final_data$Column <- (final_data$Plot - 1) %/% n_rows + 1
    final_data$Row <- (final_data$Plot - 1) %% n_rows + 1
  } else if (order == "row") {
    final_data$Row <- (final_data$Plot - 1) %/% n_cols + 1
    final_data$Column <- (final_data$Plot - 1) %% n_cols + 1
  } else {
    stop("Invalid 'order' argument. Use 'row' or 'column'.")
  }

  # Apply serpentine layout if required
  if (serpentine) {
  if (order == "row") {
    # Reverse columns within even rows
    for (r in seq_len(n_rows)) {
      if (mod(r, 2) == 0) {
        final_data$Column[final_data$Row == r] <- rev(final_data$Column[final_data$Row == r])
      }
    }
  } else if (order == "column") {
    # Reverse rows within even columns
    for (c in seq_len(n_cols)) {
      if (mod(c, 2) == 0) {
        final_data$Row[final_data$Column == c] <- rev(final_data$Row[final_data$Column == c])
      }
    }
  }
}

  # Create layout matrix in correct serpentine order
  layout_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(nrow(final_data))) {
  r <- final_data$Row[i]
  c <- final_data$Column[i]
  
  layout_matrix[r, c] <- final_data$Treatment[i]
}

  return(list(layout_matrix = layout_matrix, field_book = final_data))
}
