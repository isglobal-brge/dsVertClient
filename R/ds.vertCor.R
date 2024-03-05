#' Client function for vertCor for finding correlation of vertically partitioned data in DataSHIELD
#'
#' @param variable_names names of variables to be used in correlation matrix
#' @param datasources conns for the defined study
#'
#' @return correlation matrix
#' @export
#'
#' @examples TODO
ds.vertCor <- function(variable_names, datasources = NULL) {
  # standard, if NULL look for connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  study_names  <- names(datasources)

  calltext <- call("vertCorDS", variable_names)

  output <- DSI::datashield.aggregate(datasources, calltext)

  U_combined <- NULL
  #print(output)

  for(i in 1:length(study_names)){
    result <- output
    #print(output)
    U_combined <- cbind(U_combined, result)
  }


  # Perform SVD on the combined U matrix
  svd_final <- svd(U_combined)

  # Calculate and scale the correlation matrix
  correlation <- svd_final$v %*% diag(svd_final$d^2) %*% t(svd_final$v)
  correlation <- correlation / correlation[1, 1]

  return(correlation)
}
