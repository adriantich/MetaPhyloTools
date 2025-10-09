#' Create haplotype network data for visualization
#'
#' This function processes MOTU (Molecular Operational Taxonomic Unit) data and metadata
#' to create a haplotype network and associated pie chart data for phylogenetic analysis.
#' It calculates mean abundances across samples grouped by a specified variable and
#' generates network data suitable for haplotype network visualization.
#'
#' @param motu_tab A data.frame containing MOTU abundance data. Rows represent MOTUs/haplotypes,
#'   columns represent samples. Must include columns specified by \code{id_col} and \code{seq_col}.
#' @param metadata A data.frame containing sample metadata with grouping information.
#'   Must include columns specified by \code{grouping_col} and \code{sample_col}.
#' @param grouping_col Character string. Name of the column in \code{metadata} that contains
#'   the grouping variable (e.g., "treatment", "location").
#' @param sample_col Character string. Name of the column in \code{metadata} that contains
#'   sample identifiers matching column names in \code{motu_tab}.
#' @param seq_col Character string. Name of the column in \code{motu_tab} that contains
#'   DNA sequences for haplotype analysis.
#' @param id_col Character string. Name of the column in \code{motu_tab} that contains
#'   unique identifiers for each MOTU/haplotype.
#' @param size_correction Numeric. Factor to correct haplotype frequencies for visualization.
#'   Default is 1 (no correction).
#'
#' @return A list containing two elements:
#'   \item{network}{A haplotype network object created by \code{pegas::haploNet()}, or NULL if
#'     only one haplotype is found. Contains network structure with frequency attributes.}
#'   \item{pieseas}{A matrix with haplotypes as rows and groups as columns, containing
#'     mean abundances for pie chart visualization.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Groups samples by the specified grouping variable
#'   \item Calculates mean abundances for each MOTU across samples within each group
#'   \item Filters out MOTUs with zero total abundance
#'   \item Converts DNA sequences to DNAbin format for phylogenetic analysis
#'   \item Creates haplotypes and constructs a haplotype network
#'   \item Attaches frequency information to the network for visualization
#' }
#'
#' @note
#' \itemize{
#'   \item Requires packages \code{ape} and \code{pegas}
#'   \item DNA sequences should be aligned and of equal length
#'   \item Sample identifiers in metadata must match column names in motu_tab
#'   \item If only one haplotype is found, network will be NULL with a warning
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result <- haploNet_data(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "treatment",
#'   sample_col = "sample_id",
#'   seq_col = "sequence",
#'   id_col = "motu_id",
#'   size_correction = 100
#' )
#' 
#' # Access results
#' network <- result$network
#' pie_data <- result$pieseas
#' }
#'
#' @export
haploNet_data <- function(motu_tab, metadata, grouping_col, sample_col, seq_col, id_col, size_correction = 1) {
    suppressMessages(library(ape))
    suppressMessages(library(pegas))

    pieseas <- NULL
    groups <- unique(metadata[[grouping_col]])
    # important, pieseas can't be a data.frame, it has to be a matrix
    for (case in groups) {
        cols <- metadata[[sample_col]][metadata[[grouping_col]] == case]
        if (length(cols) < 1) {
            next
        }
        if (sum(!(cols %in% colnames(motu_tab))) > 0) {
            warning(paste(
                "Sample columns",
                paste(cols[!cols %in% colnames(motu_tab)],
                    collapse = ", "
                ),
                "not found in motu_tab. Skipping the cases.",
                "\nCheck the sample_col and grouping_col",
                "arguments or the motu_tab samples."
            ))
            cols <- cols[cols %in% colnames(motu_tab)]
            if (length(cols) < 1) {
                warning(paste("No valid sample columns found for case", case, ". Skipping this case."))
                next
            }
        }
        df_case <- data.frame(
            col = rowMeans(motu_tab[, cols]),
            row.names = motu_tab[, id_col]
        )
        colnames(df_case) <- case
        if (is.null(pieseas)) {
            pieseas <- df_case
        } else {
            pieseas <- cbind(pieseas, df_case)
        }
    }
    if (is.null(pieseas)) {
        stop("No valid samples found in the provided metadata and motu_tab. Please check your inputs.")
    }
    pieseas <- as.matrix(pieseas)
    mida <- rowSums(pieseas)
    pieseas <- pieseas[mida > 0, , drop = FALSE]
    rownames(pieseas) <- seq_len(nrow(pieseas))
    if (nrow(pieseas) == 0) {
        stop("No valid haplotypes found. Please check your input data.")
    } else if (nrow(pieseas) == 1) {
        warning("Only one haplotype found. The haplotype network will be empty.")
        network <- NULL
    } else {
        # Convert sequences to dnaBIN and create haplotypes
        rr <- (cbind(as.character(motu_tab[mida > 0, id_col]), as.character(motu_tab[mida > 0, seq_col])))
        r <- t(sapply(strsplit(rr[, 2], ""), tolower))
        rownames(r) <- rr[, 1]
        y <- as.DNAbin(r)
        haplo <- haplotype(y)
        network <- haploNet(haplo, getProb = FALSE)
        mida <- mida[mida > 0]

        attr(network, "labels") <- mida # per usar les frequencies com a labels
        attr(network, "freq") <- mida / size_correction
    }


    return(list(network = network, pieseas = pieseas))
}
