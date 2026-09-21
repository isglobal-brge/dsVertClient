root <- Sys.getenv("DSVERT_VALIDATION_ROOT", "/workspace/rc1")
out <- file.path(root, "dsVertClient/inst/validation/v1.3.0")
library(dsVertClient)
status <- ds.vertMethodStatus()
write.csv(status, file.path(out, "method_status_v130.csv"), row.names = FALSE)
print(names(status))
print(table(status$status))
print(table(status$release_contract))
server <- read.dcf(file.path(root, "dsVert/DESCRIPTION"))
endpoints <- trimws(unlist(strsplit(server[1, c("AggregateMethods", "AssignMethods")], ",", fixed = TRUE)))
summary <- list(rows = nrow(status), by_status = as.list(table(status$status)),
  by_release_contract = as.list(table(status$release_contract)), server_endpoint_count = length(endpoints),
  R = R.version.string, session = capture.output(sessionInfo()))
jsonlite::write_json(summary, file.path(out, "registry_counts_v130.json"), pretty = TRUE, auto_unbox = TRUE)
