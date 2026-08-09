# Retained implementations remain executable only inside the package's legacy
# regression suite. The installed binding is the constant FALSE and has no
# option/environment/request override.
.dsvert_test_replace_client_binding <- function(name, value) {
  namespace <- asNamespace("dsVertClient")
  was_locked <- bindingIsLocked(name, namespace)
  if (was_locked) unlockBinding(name, namespace)
  assign(name, value, envir = namespace)
  if (was_locked) lockBinding(name, namespace)
  invisible(value)
}

.dsvert_test_replace_client_binding(
  ".dsvert_quarantine_test_mode", function() TRUE)
