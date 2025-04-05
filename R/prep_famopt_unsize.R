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
#' @param cols_per_row A vector indicating the size of each columns according to the field reality. cols_per_row is required only when n_cols is set to `NULL`.
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
#' - **Unequal column size** If `cols_per_row` is defined, and `n_cols` is `NULL`, the function automatically allows for unequal column sizes.
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
# Example: Create a p-rep design with 3 checks, 3 p-rep, and 6 unreplicated treatment in 10 blocks. Columns have unequal sizes.
# check_treatments <- c("CheckA", "CheckB", "CheckC")
# check_families <- c("FamCheck1", "FamCheck2", "FamCheck3")   # Families for checks
# p_rep_treatments <- c("T1", "T2", "T3")                    # One partially replicated treatment
# p_rep_reps <- c(8, 8, 8)                              # T1 will appear in 2 out of 3 blocks
# p_rep_families <- c("Fam1", "Fam2", "Fam2")   # Family for T1
# unreplicated_treatments <- c("U1", "U2", "U3", "U4", "U5", "U6")              # One unreplicated treatment
# unreplicated_families <- c("Fam1", "Fam2", "Fam3", "Fam4", "Fam5", "Fam6")              # Family for U1
# n_blocks <- 10
# n_rows <- 5
#
# cols_per_row <- c(10, 15, 12, 10, 13)
#' design <- prep_famopt_unsize(check_treatments, check_families,
#'                              p_rep_treatments, p_rep_reps, p_rep_families,
#'                              unreplicated_treatments, unreplicated_families,
#'                              n_blocks = n_blocks, n_rows = n_rows, n_cols = NULL,
#'                              cols_per_row = cols_per_row,
#'                              serpentine = TRUE, seed = NULL,
#'                              attempts = 1000, warn_and_correct = TRUE, fix_rows = TRUE)
#'
#' # --- Create layout matrix ---
# ncols <- max(design$field_book$Column)
# layout_matrix <- matrix(NA, nrow = 5, ncol = ncols)
# for (i in seq_len(nrow(combined_fieldbook))) {
#   r <- design$field_book$Row[i]
#   c <- design$field_book$Column[i]
#   layout_matrix[r, c] <- design$field_book$Treatment[i]
# }

#' design$layout_matrix
#' #      [,1]     [,2]     [,3] [,4]     [,5]     [,6]     [,7]     [,8]     [,9]     [,10]    [,11]    [,12]    [,13]    [,14] [,15]
#' # [1,] "CheckC" "CheckA" "T1" "T2"     "CheckB" "U1"     "T3"     "U6"     "CheckC" "CheckB" NA       NA       NA       NA    NA
#' # [2,] "T2"     "CheckB" "U5" "T3"     "T1"     "T3"     "CheckB" "U2"     "CheckA" "T2"     "CheckC" "T3"     "T1"     "T2"  "CheckA"
#' # [3,] "CheckA" "CheckC" "U4" "CheckA" "T1"     "CheckB" "CheckC" "CheckC" "CheckA" "T1"     "T2"     "CheckB" NA       NA    NA
#' # [4,] "T2"     "CheckC" "T3" "CheckB" "CheckC" "T2"     "T1"     "CheckA" "U3"     "T3"     NA       NA       NA       NA    NA
#' # [5,] "T1"     "CheckB" "T3" "CheckA" "T2"     "CheckB" "CheckA" "CheckC" "CheckC" "T3"     "T1"     "CheckA" "CheckB" NA    NA
#'
#'
#' # View the field book (first few entries)
# > head(design$field_book)
#         Treatment   Family Block Plot Row Column
# CheckC    CheckC FamCheck3     1    1   1      1
# CheckA    CheckA FamCheck1     1    2   1      2
# T1            T1      Fam1     1    3   1      3
# T2            T2      Fam2     1    4   1      4
# CheckB    CheckB FamCheck2     1    5   1      5
# U1            U1      Fam1     1    6   1      6
#' #
#' # The field book is sorted by plot order.
#'
#'
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
    p_rep_in_block <- names(p_rep_assignments)[sapply(p_rep_assignments, function(x) b %in% x)]
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
