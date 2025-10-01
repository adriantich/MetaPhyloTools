
#' Calculate Percolation Threshold for Ecological Networks
#'
#' Calculates the percolation threshold for ecological networks by identifying 
#' the critical distance at which the network fragments into small components.
#' This function implements the percolation analysis described in EDENetworks.
#' 
#'
#' @param diss A distance matrix representing dissimilarities between nodes
#'
#' @return A list containing:
#' \describe{
#'   \item{graph}{The original igraph object}
#'   \item{S_list}{Vector of S values for each threshold}
#'   \item{critical_distance}{The critical distance at percolation threshold}
#' }
#'
#' @details Percolation threshold is defined as follows in the original reference:
#' from
#' 
#' Arnaud-Haond et al. EDENetworks: Ecological and Evolutionary Networks.
#' 
#' Percolation threshold: When links are removed from a connected network, it
#' eventually fragments into small components. The point where this happens 
#' is called the percolation threshold. More accurately, this is the point
#' where the so-called giant component (whose size is of the order of the 
#' network size) disappears and there is no long-range connectivity; even 
#' before the percolation threshold small disconnected fragments will appear,
#' yet a substantial fraction of nodes belongs to the giant component.
#' The precise location of this percolation point is made using the definition
#' classically proposed for finite systems (Stauffer and Aharony, 1994) by 
#' calculating the average proposed for finite systems (Stauffer and 
#' Aharony, 1994) by calculating the average size of the clusters excluding
#' the largest one:
#' 
#' S = sum{s<Smax}(s^2ns)
#' 
#' as a function of the last distance value removed, thr, and identifying the 
#' critical distance with the one at which <S>* has a maximum. N is the total
#' number of nodes not included in the largest cluster and ns is the number 
#' of clusters containing s nodes.
#'
#' @references 
#' Arnaud-Haond et al. EDENetworks: Ecological and Evolutionary Networks. 
#' \url{https://www.researchgate.net/profile/Sophie-Arnaud-Haond/publication/268433406_EDENetworks_Ecological_and_Evolutionary_Networks/links/54f705a30cf28d6dec9c7ad8/EDENetworks-Ecological-and-Evolutionary-Networks.pdf}
#'
#' @examples
#' # Create example distance matrix
#' dist_mat <- as.matrix(dist(iris[1:10, 1:4]))
#' result <- edenetwork_percolation(dist_mat)
#' print(result$critical_distance)
#'
#' @importFrom igraph graph_from_adjacency_matrix components delete_edges E
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

edenetwork_percolation <- function(diss) {
    # from
    # https://www.researchgate.net/profile/Sophie-Arnaud-Haond/publication/268433406_EDENetworks_Ecological_and_Evolutionary_Networks/links/54f705a30cf28d6dec9c7ad8/EDENetworks-Ecological-and-Evolutionary-Networks.pdf
    # Percolation threshold: When links are removed from a connected network, it
    # eventually fragments into small components. The point where this happens is called
    # the percolation threshold. More accurately, this is the point where the so-called giant
    # component (whose size is of the order of the network size) disappears and there is
    # no long-range connectivity; even before the percolation threshold small disconnected
    # fragments will appear, yet a substantial fraction of nodes belongs to the giant
    # component
    # The precise location of this percolation point is made using the definition classically
    # proposed for finite systems (Stauffer and Aharony, 1994) by calculating the average
    # size of the clusters excluding the largest one:
    # S = sum{s<Smax}(s^2ns)
    # as a function of the last distance value removed, thr, and identifying the critical
    # distance with the one at which <S>* has a maximum. N is the total number of nodes
    # not included in the largest cluster and ns is the number of clusters containing s
    # nodes.
    suppressMessages(library(igraph))
    # Create a graph from the distance matrix
    g <- graph_from_adjacency_matrix(as.matrix(diss), mode = "undirected", weighted = TRUE, diag = FALSE)

    # Initialize an empty list to store S for each threshold
    S_list <- list()

    # Sort the edge weights (distances) in descending order
    thresholds <- sort(E(g)$weight, decreasing = TRUE)

    # Create a progress bar
    pb <- txtProgressBar(min = 0, max = length(thresholds), style = 3)

    # For each threshold
    for (i in seq_along(thresholds)) {
        # Update the progress bar
        setTxtProgressBar(pb, i)

        thr <- thresholds[i]
        # Remove edges with weight greater than the threshold
        g_thr <- delete_edges(g, E(g)[weight > thr])

        # Identify the clusters
        clusters <- components(g_thr)

        # Exclude the largest cluster
        clusters_size <- clusters$csize
        clusters_size <- clusters_size[clusters_size < max(clusters_size)]
        clusters_size <- table(clusters_size)

        # Calculate S
        S <- sum(as.numeric(names(clusters_size))^2 * (clusters_size))

        # Store S in the list
        S_list[[as.character(thr)]] <- S
    }

    # Close the progress bar
    close(pb)

    # Identify the threshold at which S has a maximum
    critical_distance <- as.numeric(names(which.max(unlist(S_list))))

    summary_list <- list(graph = g, S_list = unlist(S_list), critical_distance = critical_distance)
    return(summary_list)
}
