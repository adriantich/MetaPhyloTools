#' Plot haplotype network with pie charts
#'
#' This function creates a publication-ready haplotype network plot with pie charts
#' showing the relative abundance or frequency of haplotypes across different groups.
#' The plot includes a legend and title, and is saved as a PDF file.
#'
#' @param haploNet_data_object A list object returned by \code{haploNet_data()} containing
#'   network and pie chart data. Must include elements:
#'   \itemize{
#'     \item{\code{network}}{A haplotype network object from \code{pegas::haploNet()}}
#'     \item{\code{pieseas}}{A matrix with pie chart data for each haplotype}}
#' @param output_file Character string. File path for the output PDF plot (e.g., "network_plot.pdf").
#' @param bg A vector of colors for the pie chart segments. Length should match the number
#'   of groups in the pie chart data. Used for both pie segments and legend.
#' @param plot_title Character string. Main title to be displayed at the top of the plot.
#'
#' @return No return value. The function creates a PDF file at the specified location.
#'
#' @details
#' The function creates a multi-panel plot with:
#' \itemize{
#'   \item Main haplotype network with pie charts at each node
#'   \item Node sizes proportional to haplotype frequencies
#'   \item Title at the top of the plot
#'   \item Legend at the bottom showing group categories (currently hardcoded as "Inside"/"Outside")
#' }
#'
#' Plot specifications:
#' \itemize{
#'   \item PDF dimensions: 15 x 15 inches
#'   \item High-resolution output suitable for publication
#'   \item Customizable colors and title
#'   \item Error handling for plot generation issues
#' }
#'
#' @note
#' \itemize{
#'   \item Requires the \code{pegas} package for network plotting
#'   \item Legend labels are currently hardcoded as "Inside" and "Outside"
#'   \item The function assumes a two-group comparison (modify legend for more groups)
#'   \item If plotting fails, an error message is displayed with details
#' }
#'
#' @examples
#' \dontrun{
#' # Create network data first
#' network_data <- haploNet_data(
#'   motu_tab = my_motu_data,
#'   metadata = my_metadata,
#'   grouping_col = "location",
#'   sample_col = "sample_id",
#'   seq_col = "sequence",
#'   id_col = "motu_id"
#' )
#' 
#' # Create the plot
#' haploNet_plot(
#'   haploNet_data_object = network_data,
#'   output_file = "haplotype_network.pdf",
#'   bg = c("red", "blue"),
#'   plot_title = "Haplotype Network Analysis"
#' )
#' }
#'
#' @export
haploNet_plot <- function(haploNet_data_object, output_file, bg, plot_title) {
    network <- haploNet_data_object$network
    pieseas <- haploNet_data_object$pieseas
    pdf(output_file, width = 15, height = 15)
    par(oma = c(2, 2, 2, 2), mar = c(15, 15, 15, 15))
    par(fig = c(0, 1, 0.1, 0.9))
    tryCatch(
        {
            plot(network,
                size = attr(network, "freq"), threshold = 0, labels = F,
                pie = pieseas, bg = bg,
                fast = T, show.mutation = 0
            )
        },
        error = function(e) {
            cat("Custom error: The haplotype network plot could not be generated.\n")
            cat("Reason:", e$message, "\n")
            stop("Plot error.")
        }
    )
    par(new = T)
    plot.new()
    par(fig = c(0, 1, 0, 1))
    legend("top", bty = "n", text.font = 2, cex = 1.5, legend = plot_title)
    par(fig = c(0, 1, 0, 0.4), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = T, xpd = NA)
    legend("bottom",
        legend = c("Inside", "Outside"),
        fill = bg,
        cex = 1.5,
        text.font = 2,
        bty = "n",
        horiz = T,
        xpd = NA
    )
    dev.off()
}
