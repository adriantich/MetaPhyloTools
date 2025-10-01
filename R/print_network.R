

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