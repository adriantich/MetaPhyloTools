
#' Create Haplotype Network Plots with ggplot2
#'
#' Creates haplotype network visualizations using ggplot2 with customizable 
#' legend positioning and scatter pie charts for nodes.
#'
#' @param data A list containing 'nodes' and 'connections' data frames for 
#'   the network. The 'nodes' data frame should contain x, y coordinates and 
#'   size information. The 'connections' data frame should specify network edges.
#' @param legend_pos Position for the legend. Options: "left_bottom", 
#'   "right_bottom", "top_left", "top_right", "left", "right", "top", 
#'   "bottom", "none"
#'
#' @return A ggplot object representing the haplotype network
#'
#' @details This function creates publication-ready haplotype network plots 
#'   with nodes represented as pie charts (using scatterpie) and customizable 
#'   legend positioning to avoid overlap with the network. The aim of this 
#'   function is to solve a problem often encountered when plotting haplotype
#'   networks with ggplot2, where the legend can overlap with the network 
#'   itself.
#'
#' @examples
#' # Example with mock network data
#' nodes <- data.frame(
#'   x = runif(5), 
#'   y = runif(5), 
#'   size = runif(5, 1, 3),
#'   pop1 = runif(5),
#'   pop2 = runif(5)
#' )
#' connections <- data.frame(from = c(1,2,3), to = c(2,3,4))
#' data <- list(nodes = nodes, connections = connections)
#' plot <- haplo_ggplot(data)
#'
#' @importFrom ggplot2 ggplot geom_segment theme_void coord_equal
#' @importFrom scatterpie geom_scatterpie
#' @export
haplo_ggplot <- function(data, legend_pos = "left_bottom") {
    suppressMessages(library(scatterpie))
    suppressMessages(library(ggplot2))
    nodes <- data$nodes
    connections <- data$connections
    if (!(legend_pos %in% c("left_bottom", "right_bottom", "top_left", "top_right", "left", "right", "top", "bottom", "none"))) {
        stop("legend_pos must be one of 'left_bottom', 'right_bottom', 'top_left', 'top_right', 'left', 'right', 'top', 'bottom', or 'none'.")
    } else {
        legend_yes <- TRUE
        if (legend_pos == "left_bottom") {
            x0 <- min(nodes$x)*0.999
            y0 <- min(nodes$y)*0.999
        } else if (legend_pos == "right_bottom") {
            x0 <- max(nodes$x)*1.001
            y0 <- min(nodes$y)*0.999
        } else if (legend_pos == "top_left") {
            x0 <- min(nodes$x)*0.999
            y0 <- max(nodes$y)*1.001
        } else if (legend_pos == "top_right") {
            x0 <- max(nodes$x)*1.001
            y0 <- max(nodes$y)*1.001
        } else if (legend_pos == "left") {
            x0 <- min(nodes$x)*0.999
            y0 <- (max(nodes$y) + min(nodes$y)) / 2
        } else if (legend_pos == "right") {
            x0 <- max(nodes$x)*1.001
            y0 <- (max(nodes$y) + min(nodes$y)) / 2
        } else if (legend_pos == "top") {
            x0 <- (max(nodes$x) + min(nodes$x)) / 2
            y0 <- max(nodes$y)*1.001
        } else if (legend_pos == "bottom") {
            x0 <- (max(nodes$x) + min(nodes$x)) / 2
            y0 <- min(nodes$y)*0.999
        } else if (legend_pos == "none") {
            x0 <- NULL
            y0 <- NULL
            legend_yes <- FALSE
        }
    }
    if (legend_yes) {
        circles_df <- data.frame(
            x = nodes$x,
            y = nodes$y,
            radius = nodes$size / 2
        )
        legend_coords <- resolve_overlap(
            circles_df,
            x = x0,
            y = y0,
            r = max(circles_df$radius)
        )
    }

    haplot <-  
        ggplot() +
        geom_point(data = nodes, aes(x = x, y = y), size = 0.1, color = "black") +
        geom_segment(data = connections,
            aes(x = x, y = y, xend = x.1, yend = y.1),
            color = "black", size = 0.5
        ) +
        geom_scatterpie(data = nodes, aes(x = x, y = y, r = size / 2), cols = colnames(nodes)[-c(1:4)]) +
        coord_equal()
    if (legend_yes) {
        haplot <- haplot +
            geom_scatterpie_legend(circles_df$radius ,
                x = legend_coords$x, y = legend_coords$y
            )
    }
    return(haplot)
}



resolve_overlap <- function(circles_df, x, y, r, dist = 0.1, fixed_dist = FALSE) {
    if (fixed_dist) {
        extra_dist <- dist
    } else {
        extra_dist <- r * (1 + dist)
    }
    dists <- sqrt((circles_df$x - x)^2 + (circles_df$y - y)^2)
    overlaps <- dists < (circles_df$radius + r * 1.1)
    i <- 0
    while (any(overlaps) && i != nrow(circles_df)) {
        i <- i + 1
        dist_from_center <- which.min(dists - circles_df$radius)
        scale <- (circles_df$radius[dist_from_center] + extra_dist * 1.1) / dists[dist_from_center]
        x <- circles_df$x[dist_from_center] + (x - circles_df$x[dist_from_center]) * scale
        y <- circles_df$y[dist_from_center] + (y - circles_df$y[dist_from_center]) * scale
        dists <- sqrt((circles_df$x - x)^2 + (circles_df$y - y)^2)
        overlaps <- dists < (circles_df$radius + r * 1.1)
    }

    return(data.frame(x = x, y = y, r = r))
}
