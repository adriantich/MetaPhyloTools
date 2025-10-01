
#' Rarefy Data Within MOTUs
#'
#' Performs rarefaction of sample data within each MOTU (Molecular Operational 
#' Taxonomic Unit) to standardize sampling effort.
#'
#' @param df Data frame containing MOTU abundance data
#' @param sample_names Vector of column names representing samples
#' @param motu_col Name of the column containing MOTU identifiers (default: "MOTU")
#' @param rarefy Logical, whether to perform rarefaction (default: TRUE)
#' @param rarefy_to Function or numeric value to determine rarefaction depth 
#'   (default: min function)
#' @param rel_abund Logical, whether to convert to relative abundance after 
#'   rarefaction (default: TRUE)
#'
#' @return A data frame with rarefied abundance data
#'
#' @details This function rarefies data within each MOTU separately, which is 
#'   important for metabarcoding data where different MOTUs may have very 
#'   different abundance ranges.
#'
#' @examples
#' # Example with mock data
#' df <- data.frame(
#'   MOTU = rep(c("MOTU1", "MOTU2"), each = 3),
#'   sample1 = rpois(6, 10), 
#'   sample2 = rpois(6, 15)
#' )
#' result <- rarefy_within_motu(df, c("sample1", "sample2"))
#'
#' @importFrom vegan rrarefy
#' @export
rarefy_within_motu <- function(df,
                               sample_names,
                               motu_col = "MOTU",
                               rarefy = TRUE,
                               rarefy_to = min,
                               rel_abund = TRUE) {
    suppressMessages(library(vegan))
    if (is.numeric(rarefy_to)) {
        value <- rarefy_to
        rarefy_to <- function(x) {
            return(value)
        }
    } else if (!is.function(rarefy_to)) {
        stop("rarefy_to must be a function or a numeric value")
    }

    df <- df[rowSums(df[, sample_names]) > 0, ]
    # check if the MOTU column is present
    if (!(motu_col %in% colnames(df))) {
        stop("The MOTU column is not present in the dataframe")
    }

    for (x in df[, motu_col]) {
        motu_df <- df[df[, motu_col] == x, sample_names, drop = F]
        if (dim(motu_df)[1] > 1) {
            if (rarefy) {
                motu_df <-
                    suppressWarnings(
                        t(rrarefy(
                            t(motu_df),
                            floor(rarefy_to(colSums(motu_df)[colSums(motu_df) > 0]))
                        ))
                    )
            }
        }
        if (rel_abund) {
            if (dim(motu_df)[1] == 1) {
                motu_df_rel <- motu_df
                motu_df_rel[, motu_df_rel > 0] <- 1
            } else {
                motu_df_rel <- motu_df
                motu_df_rel[, colSums(motu_df_rel) > 0] <-
                    t(t(motu_df_rel[, colSums(motu_df_rel) > 0]) /
                        colSums(motu_df_rel)[colSums(motu_df_rel) > 0])
            }
        }

        df[df[, motu_col] == x, sample_names] <- motu_df
    }

    return(df)
}
