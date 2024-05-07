ds.vertCor <- function(data_names, variable_names, datasources = NULL) {
  # standard, if NULL look for connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  U_combined <- NULL

  study_names  <- names(datasources)
  for(i in 1:length(study_names)) {
    temp <-variable_names[[i]]
    assign("temp", temp, envir = .GlobalEnv)
    output <- datashield.aggregate(datasources[i], quote(vertCorDS(D, temp)))[[1]]
    U_combined <- cbind(U_combined, output)
  }


  # Perform SVD on the combined U matrix
  svd_final <- svd(U_combined)
  # Calculate and scale the correlation matrix
  correlation <- svd_final$v %*% diag(svd_final$d^2) %*% t(svd_final$v)
  correlation <- correlation / correlation[1, 1]
  return(correlation)
}
