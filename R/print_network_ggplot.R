#' Create haplotype network visualization using ggplot2
#'
#' This function creates a modern, customizable haplotype network plot using ggplot2.
#' It processes MOTU data and metadata to generate a network visualization with nodes
#' representing haplotypes and edges showing relationships between them.
#'
#' @param motu_tab A data.frame containing MOTU abundance data. Rows represent MOTUs/haplotypes,
#'   columns represent samples. Must include columns specified by \code{id_col} and \code{seq_col}.
#' @param output_file Character string. File path for the output PDF plot. Default is
#'   "haplotype_network.pdf". Currently not used in this function.
#' @param metadata A data.frame containing sample metadata with grouping information.
#'   Must include columns specified by \code{grouping_col} and \code{sample_col}.
#' @param grouping_col Character string. Name of the column in \code{metadata} that contains
#'   the grouping variable (e.g., "treatment", "location").
#' @param sample_col Character string. Name of the column in \code{metadata} that contains
#'   sample identifiers matching column names in \code{motu_tab}.
#' @param seq_col Character string. Name of the column in \code{motu_tab} that contains
#'   DNA sequences for haplotype analysis. Default is "sequence".
#' @param id_col Character string. Name of the column in \code{motu_tab} that contains
#'   unique identifiers for each MOTU/haplotype. Default is "id".
#' @param size_correction Numeric. Factor to correct haplotype frequencies for visualization.
#'   Default is 1 (no correction).
#' @param plot_title Character string. Main title to be displayed on the plot.
#'   Default is "Haplotype Network".
#' @param bg A vector of colors for the groups. If NULL, default colors will be used.
#' @param method Character string. Layout method for the network. Must be either "pegas"
#'   (uses pegas package layout) or "fruchtermanreingold" (uses Fruchterman-Reingold layout).
#'   Default is "pegas".
#' @param legend_pos Character string. Position of the legend. Default is "left_bottom".
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return A ggplot object containing the haplotype network visualization.
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Calls \code{haploNet_data()} to process input data and create network structure
#'   \item Uses \code{haplodata4ggplot()} to convert network data to ggplot-compatible format
#'   \item Creates the final plot using \code{haplo_ggplot()}
#' }
#'
#' The function supports two layout methods:
#' \itemize{
#'   \item \code{"pegas"}: Uses the default layout from the pegas package
#'   \item \code{"fruchtermanreingold"}: Uses the Fruchterman-Reingold force-directed layout
#' }
#'
#' @note
#' \itemize{
#'   \item Requires packages \code{ape}, \code{pegas}, and \code{ggplot2}
#'   \item This function depends on helper functions: \code{haploNet_data()}, 
#'     \code{haplodata4ggplot()}, and \code{haplo_ggplot()}
#'   \item The \code{output_file} parameter is included for consistency but not currently used
#'   \item DNA sequences should be aligned and of equal length
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot <- print_network_ggplot(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "treatment",
#'   sample_col = "sample_id",
#'   plot_title = "Treatment Comparison"
#' )
#' 
#' # With custom parameters
#' plot <- print_network_ggplot(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "location",
#'   sample_col = "sample_id",
#'   seq_col = "dna_sequence",
#'   id_col = "motu_id",
#'   size_correction = 100,
#'   method = "fruchtermanreingold",
#'   bg = c("#FF6B6B", "#4ECDC4", "#45B7D1"),
#'   legend_pos = "right"
#' )
#' 
#' # Display the plot
#' print(plot)
#' 
#' # Save the plot
#' ggsave("my_network.pdf", plot, width = 12, height = 10)
#' }
#'
#' @seealso
#' \code{\link{haploNet_data}}, \code{\link{haploNet_plot}}
#'
#' @export
print_network_ggplot <- function(
    motu_tab,
    output_file = "haplotype_network.pdf",
    metadata,
    grouping_col,
    sample_col,
    seq_col = "sequence",
    id_col = "id",
    size_correction = 1,
    plot_title = "Haplotype Network",
    bg = NULL,
    method = "pegas",
    legend_pos = "left_bottom",
    ...) {
    # motu <- 'HBLN_000000001'
    # motu_tab <- df[df[, motu_col] == motu, ]
    # output_file <- 'example_network.pdf'
    # metadata <- metadata
    # grouping_col <- 'codi.comunitat'
    # sample_col <- 'original_samples'
    # seq_col <- 'fin_nosample'
    # id_col <- 'id'
    # size_correction <- 1000
    # plot_title = "Haplotype Network"
    # bg = NULL
    # method= "pegas"
    # method= "fruchtermanreingold"

    if (!(method %in% c("pegas", "fruchtermanreingold"))) {
        stop("Method must be either 'pegas' or 'fruchtermanreingold'.")
    }

    haploNet_data_object <- haploNet_data(
        motu_tab = motu_tab,
        metadata = metadata,
        grouping_col = grouping_col,
        sample_col = sample_col,
        seq_col = seq_col,
        id_col = id_col,
        size_correction = size_correction
    )
    data4ggplot <- haplodata4ggplot(
        haploNet_data_object = haploNet_data_object,
        method = method
    )
    haplo_plot <- haplo_ggplot(
        data = data4ggplot,
        legend_pos = legend_pos
    )
    return(haplo_plot)
}