hmap <- function(dataMatrix, orderCol = T, orderRow = T, dendroLineSize = 0.5, 
                 fontSize = 20, colorPalette = "Spectral", scaleName = "value", distMethod = "euclidean", 
                   clustMethod = "complete", revColors=F, main="", tileColor=NA) {
  
  data_m <- tibble::rownames_to_column(dataMatrix) %>% reshape2::melt()
  
  # Cluster rows
  if (orderRow) {
    dd.row <- as.dendrogram(hclust(dist(dataMatrix, method = distMethod), method = clustMethod))
    row.ord <- order.dendrogram(dd.row)
    ordered_row_names <- row.names(dataMatrix[row.ord, ])
    data_m$rowname <- factor(data_m$rowname, levels = ordered_row_names)
  }
  
  # Cluster columns
  if (orderCol) {
    dd.col <- as.dendrogram(hclust(dist(t(dataMatrix), method = distMethod), 
                                   method = clustMethod))
    col.ord <- order.dendrogram(dd.col)
    ordered_col_names <- colnames(dataMatrix[, col.ord])
    data_m$variable <- factor(data_m$variable, levels = ordered_col_names)
  }
  
  heat_plot <- ggplot2::ggplot(data_m, ggplot2::aes(x = variable, y = rowname, fill = value)) + 
    ggplot2::geom_tile(color=tileColor) + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(axis.line = ggplot2::element_line(size = 0),
                   text = ggplot2::element_text(size = fontSize),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
    ggplot2::scale_y_discrete(position = "right") + 
    ggplot2::xlab("") + 
    ggplot2::ylab("") + 
    ggplot2::ggtitle(main) +
    ggplot2::scale_fill_distiller(palette = colorPalette, name = scaleName, direction=ifelse(revColors,-1,1)) +
    guides(size = FALSE) 
  # theme(plot.margin = margin(0,0,0,4, "cm"))
  
  final_plot <- heat_plot
  
  if (orderRow) {
    dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")
    dendro_row <- cowplot::axis_canvas(heat_plot, axis = "y", coord_flip = TRUE) + 
      ggplot2::geom_segment(data = ggdendro::segment(dendro_data_row), ggplot2::aes(y = -y, 
                                                                                    x = x, xend = xend, yend = -yend), size = dendroLineSize) + 
      ggplot2::coord_flip()
    final_plot <- cowplot::insert_yaxis_grob(final_plot, dendro_row, grid::unit(0.2, 
                                                                                "null"), position = "left")
  }
  
  if (orderCol) {
    dendro_data_col <- ggdendro::dendro_data(dd.col, type = "rectangle")
    dendro_col <- cowplot::axis_canvas(heat_plot, axis = "x") + 
      ggplot2::geom_segment(data = ggdendro::segment(dendro_data_col), ggplot2::aes(x = x, y = y, xend = xend, yend = yend), size = dendroLineSize)
    final_plot <- cowplot::insert_xaxis_grob(final_plot, dendro_col, grid::unit(0.2, 
                                                                                "null"), position = "top")
  }
  
  cowplot::ggdraw(final_plot)
  
}
