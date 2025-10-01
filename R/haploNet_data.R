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
            if(length(cols) < 1) {
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
