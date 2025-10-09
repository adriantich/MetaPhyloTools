#' Create and save a haplotype network plot to PDF
#'
#' This function is a convenient wrapper that combines haplotype network data generation
#' and plotting into a single step. It processes MOTU data and metadata to create a
#' haplotype network visualization and saves it directly to a PDF file.
#'
#' @param motu_tab A data.frame containing MOTU abundance data. Rows represent MOTUs/haplotypes,
#'   columns represent samples. Must include columns specified by \code{id_col} and \code{seq_col}.
#' @param output_file Character string. File path for the output PDF plot. Default is
#'   "haplotype_network.pdf".
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
#' @param bg A vector of colors for the pie chart segments. If NULL, default colors will be used.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return Invisible NULL. The function is called for its side effect of creating a PDF file.
#'
#' @details
#' This function is a streamlined workflow that:
#' \enumerate{
#'   \item Calls \code{haploNet_data()} to process the input data and create network structure
#'   \item Calls \code{haploNet_plot()} to generate and save the plot to PDF
#'   \item Displays a confirmation message with the output file path
#' }
#'
#' The resulting PDF contains:
#' \itemize{
#'   \item Haplotype network with nodes representing haplotypes
#'   \item Pie charts at each node showing group composition
#'   \item Node sizes proportional to haplotype frequencies
#'   \item Legend and title
#'   \item High-resolution output (15x15 inches) suitable for publication
#' }
#'
#' @note
#' \itemize{
#'   \item This is a convenience function that combines \code{haploNet_data()} and \code{haploNet_plot()}
#'   \item Requires packages \code{ape} and \code{pegas}
#'   \item For more control over the plotting process, use \code{haploNet_data()} and 
#'     \code{haploNet_plot()} separately
#'   \item DNA sequences should be aligned and of equal length
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' print_network(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "treatment",
#'   sample_col = "sample_id",
#'   output_file = "treatment_network.pdf"
#' )
#' 
#' # With custom parameters
#' print_network(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "location",
#'   sample_col = "sample_id",
#'   seq_col = "dna_sequence",
#'   id_col = "motu_id",
#'   size_correction = 1000,
#'   plot_title = "Geographic Distribution",
#'   bg = c("red", "blue", "green"),
#'   output_file = "geographic_network.pdf"
#' )
#' }
#'
#' @seealso
#' \code{\link{haploNet_data}}, \code{\link{haploNet_plot}}, \code{\link{print_network_ggplot}}
#'
#' @export

print_network <- function(
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


    haploNet_data_object <- haploNet_data(
        motu_tab = motu_tab,
        metadata = metadata,
        grouping_col = grouping_col,
        sample_col = sample_col,
        seq_col = seq_col,
        id_col = id_col,
        size_correction = size_correction
    )
    haploNet_plot(
        haploNet_data_object = haploNet_data_object,
        output_file = output_file,
        bg = bg,
        plot_title = plot_title
    )
    message(paste("Haplotype network saved to", output_file))
    return(invisible())
}