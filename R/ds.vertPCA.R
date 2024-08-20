ds.vertPCA <- function(data_names, variable_names, datasources = NULL) {
  # standard, if NULL look for connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  U_combined <- NULL

  study_names  <- names(datasources)
  for(i in 1:length(study_names)) {
    temp <-var_names[[i]]
    assign("temp", temp, envir = .GlobalEnv)
    output <- datashield.aggregate(conns[i], quote(vertCorDS(D, temp)))[[1]]
    U_combined <- cbind(U_combined, output)
  }

  svd_result <- svd(U_combined)
  U <- svd_result$u
  D <- diag(svd_result$d)
  V <- svd_result$v
  principal_components <- U %*% D
  variance <- (svd_result$d^2) / sum(svd_result$d^2)
  return(list(principal_components = principal_components, variance = variance))
}
