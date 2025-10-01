#' Calculate Pairwise Jost's D Differentiation
#'
#' Calculates pairwise Jost's D differentiation indices between samples, 
#' optionally grouped by sample groups and per MOTU (Molecular Operational 
#' Taxonomic Unit).
#'
#' @param df Data frame containing MOTU abundance data
#' @param sample_names Vector of column names representing samples
#' @param sample_groups Optional vector of group assignments for samples. 
#'   Must be the same length as sample_names and in the same order.
#'   If NULL, each sample is treated as its own group
#' @param motu_col Name of the column containing MOTU identifiers (default: "MOTU")
#' @param rarefy Logical, whether to rarefy data within MOTU (default: TRUE)
#' @param rarefy_to Function to determine rarefaction depth (default: min)
#'
#' @return A list containing Jost's D matrices for each MOTU, and optionally 
#'   grouped results if sample_groups is provided
#'
#' @details Jost's D is a measure of differentiation that is independent of 
#'   within-population diversity. This function calculates it for each MOTU 
#'   separately and can group samples for analysis. The function rarefies data
#'   within each MOTU to standardize sampling effort before calculating 
#'   Jost's D this rarefaction uses the funcion vegan::rrarefy and 
#'   the specified rarefaction depth defined by the user with the parameter
#'   rarefy_to (default is the minimum sample sum across samples with non-zero
#'   counts within each MOTU).
#'
#' @examples
#' # Example with mock data
#' df <- data.frame(
#'   MOTU = rep(c("MOTU1", "MOTU2"), each = 5),
#'   sample1 = rpois(10, 5), 
#'   sample2 = rpois(10, 3),
#'   sample3 = rpois(10, 4)
#' )
#' result <- pairwise_djost(df, c("sample1", "sample2", "sample3"))
#'
#' @importFrom vegan vegdist
#' @export
pairwise_djost <-
    function(df,
             sample_names,
             sample_groups = NULL,
             motu_col = "MOTU",
             rarefy = TRUE,
             rarefy_to = min) {
        suppressMessages(library(vegan))

        if (is.null(sample_groups)) {
            grouped <- FALSE
            sample_groups <- sample_names
        } else {
            grouped <- TRUE
        }

        df <- df[rowSums(df[, sample_names]) > 0, ]
        # check if the MOTU column is present
        if (!(motu_col %in% colnames(df))) {
            stop("The MOTU column is not present in the dataframe")
        }
        djost_x_motu <- list()
        if (grouped) {
            djost_x_motu_grouped <- list()
        }
        for (motu in unique(df[, motu_col])) {
            motu_df <- df[df[, motu_col] == motu, sample_names, drop = F]

            if (dim(motu_df)[1] == 1) {
                valid_motu <- FALSE
            } else {
                valid_motu <- TRUE
                if (rarefy) {
                    motu_df <-
                        suppressWarnings(
                            t(rrarefy(
                                t(motu_df),
                                floor(rarefy_to(colSums(motu_df)[colSums(motu_df) > 0]))
                            ))
                        )
                }
                motu_df_rel <- motu_df
                motu_df_rel[, colSums(motu_df_rel) > 0] <-
                    t(t(motu_df_rel[, colSums(motu_df_rel) > 0]) /
                        colSums(motu_df_rel)[colSums(motu_df_rel) > 0])
            }


            dist_values <- c()
            for (i in 1:(length(sample_names) - 1)) {
                for (j in (i + 1):length(sample_names)) {
                    sample_a <- sample_names[i]
                    sample_b <- sample_names[j]
                    # Djost = (Ht - Hs)/(1 - Hs) * n/(n-1)
                    # n = populations = 2
                    # Ht total haplotipe diversity
                    # Hs mean within population diversity
                    if (valid_motu) {
                        if (sum(colSums(motu_df_rel[, c(sample_a, sample_b)]) %in% 0) > 0) {
                            djost <- NA
                        } else {
                            ht <- rowMeans(motu_df_rel[, c(sample_a, sample_b)])
                            ht <- 1 - sum(ht^2)
                            ha <- 1 - sum((motu_df_rel[, sample_a])^2)
                            hb <- 1 - sum((motu_df_rel[, sample_b])^2)
                            hs <- mean(c(ha, hb))
                            djost <- (ht - hs) / (1 - hs)
                        }
                    } else {
                        djost <- NA
                    }

                    dist_values <- c(dist_values, djost)
                }
            }
            dist_djost <- structure(
                dist_values,
                Size = length(sample_names), # Number of elements in the original dataset
                Labels = sample_names, # Optional: Labels for the elements
                Diag = FALSE, # Whether the diagonal is included (default is FALSE)
                Upper = FALSE, # Whether upper triangular values are included (default is FALSE)
                class = "dist" # Assign the "dist" class
            )
            djost_x_motu[[motu]] <- dist_djost

            if (grouped) {
                dist_values_grouped <- c()
                # Check if the sample groups are valid
                if (length(sample_groups) != length(sample_names)) {
                    stop("The length of sample_groups must match the number of samples")
                }
                for (i in 1:(length(unique(sample_groups)) - 1)) {
                    pop_a <- unique(sample_groups)[i]
                    pop_samples_a <- sample_names[sample_groups == pop_a]
                    pop_a_abund <- rowSums(motu_df[, pop_samples_a])
                    for (j in (i + 1):length(unique(sample_groups))) {
                        pop_b <- unique(sample_groups)[j]
                        pop_samples_b <- sample_names[sample_groups == pop_b]
                        pop_b_abund <- rowSums(motu_df[, pop_samples_b])
                        if (valid_motu) {
                            if (sum(pop_a_abund) == 0 ||
                                sum(pop_b_abund) == 0) {
                                djost <- NA
                            } else {
                                pop_a_abund <- pop_a_abund / sum(pop_a_abund)
                                pop_b_abund <- pop_b_abund / sum(pop_b_abund)
                                ht <- rowMeans(data.frame(a = pop_a_abund, b = pop_b_abund))
                                ht <- 1 - sum(ht^2)
                                ha <- 1 - sum(pop_a_abund^2)
                                hb <- 1 - sum(pop_b_abund^2)
                                hs <- mean(c(ha, hb))
                                djost <- (ht - hs) / (1 - hs)
                            }
                        } else {
                            djost <- NA
                        }
                        dist_values_grouped <- c(dist_values_grouped, djost)
                    }
                }
                dist_djost <- structure(
                    dist_values_grouped,
                    Size = length(unique(sample_groups)), # Number of elements in the original dataset
                    Labels = unique(sample_groups), # Optional: Labels for the elements
                    Diag = FALSE, # Whether the diagonal is included (default is FALSE)
                    Upper = FALSE, # Whether upper triangular values are included (default is FALSE)
                    class = "dist" # Assign the "dist" class
                )
                djost_x_motu_grouped[[motu]] <- dist_djost
            }
        }
        djost_x_motu <- lapply(djost_x_motu, function(d) {
            if (!is.null(names(d))) names(d) <- NULL
            d
        })
        result_array <- base::simplify2array(djost_x_motu)
        mean_distance_matrix <- apply(result_array, 1, mean, na.rm = TRUE)
        mean_distance_matrix <- structure(
            mean_distance_matrix,
            Size = length(sample_names), # Number of elements in the original dataset
            Labels = sample_names, # Optional: Labels for the elements
            Diag = FALSE, # Whether the diagonal is included (default is FALSE)
            Upper = FALSE, # Whether upper triangular values are included (default is FALSE)
            class = "dist" # Assign the "dist" class
        )
        if (grouped) {
            result_array_grouped <- simplify2array(djost_x_motu_grouped)
            mean_distance_matrix_grouped <- apply(result_array_grouped, 1, mean, na.rm = TRUE)
            mean_distance_matrix_grouped <- structure(
                mean_distance_matrix_grouped,
                Size = length(unique(sample_groups)), # Number of elements in the original dataset
                Labels = unique(sample_groups), # Optional: Labels for the elements
                Diag = FALSE, # Whether the diagonal is included (default is FALSE)
                Upper = FALSE, # Whether upper triangular values are included (default is FALSE)
                class = "dist" # Assign the "dist" class
            )

            output <-
                list(
                    "mean_matrix" = mean_distance_matrix,
                    "djost_x_motu" = djost_x_motu,
                    "mean_matrix_grouped" = mean_distance_matrix_grouped,
                    "djost_x_motu_grouped" = djost_x_motu_grouped
                )
        } else {
            output <-
                list(
                    "mean_matrix" = mean_distance_matrix,
                    "djost_x_motu" = djost_x_motu
                )
        }
        return(output)
    }
