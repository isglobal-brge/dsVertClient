ds.vertPCAplot <- function(data_names, variable_names, plot_xaxis = 1, plot_yaxis = 2, datasources = NULL) {
  principal_components <- ds.vertPCA(data_names, variable_names, datasources)$principal_components

  pca_plot <- plot(principal_components[, plot_xaxis], principal_components[, plot_yaxis], xlab = "PC1", ylab = "PC2", main = "PCA Plot")

  return(pca_plot)
}
