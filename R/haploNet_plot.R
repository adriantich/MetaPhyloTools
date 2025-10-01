
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
