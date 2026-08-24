test_that("every original dsVertClient analysis export remains available", {
  original_exports <- c(
    "ds.getIdentityPks", "ds.isPsiAligned", "ds.psiAlign", "ds.vert.align",
    "ds.vert.chisq", "ds.vert.chisq_cross", "ds.vert.confint",
    "ds.vert.contrast", "ds.vert.cor", "ds.vert.cox", "ds.vert.coxph",
    "ds.vert.desc", "ds.vert.fisher", "ds.vert.gee", "ds.vert.glm",
    "ds.vert.glmm", "ds.vert.ipw", "ds.vert.is_aligned",
    "ds.vert.lasso", "ds.vert.lasso_1step", "ds.vert.lasso_cv",
    "ds.vert.lasso_iter", "ds.vert.lasso_proximal", "ds.vert.lmm",
    "ds.vert.lr", "ds.vert.mi", "ds.vert.multinom", "ds.vert.nb",
    "ds.vert.ordinal", "ds.vert.pca", "ds.vert.wald", "ds.vertChisq",
    "ds.vertChisqCross", "ds.vertConfint", "ds.vertContrast",
    "ds.vertCor", "ds.vertCox", "ds.vertCoxDiscreteNonDisclosive",
    "ds.vertCoxProfileNonDisclosive", "ds.vertDesc", "ds.vertFisher",
    "ds.vertGEE", "ds.vertGLM", "ds.vertGLMM", "ds.vertIPW",
    "ds.vertLASSO", "ds.vertLASSO1Step", "ds.vertLASSOCV",
    "ds.vertLASSOIter", "ds.vertLASSOProximal", "ds.vertLMM",
    "ds.vertLMM.k3", "ds.vertLR", "ds.vertMI", "ds.vertMultinom",
    "ds.vertMultinomJoint", "ds.vertMultinomJointNewton",
    "ds.vertNBFullRegTheta", "ds.vertOrdinal",
    "ds.vertOrdinalJointNewton", "ds.vertPCA", "ds.vertWald"
  )

  exports <- getNamespaceExports("dsVertClient")
  expect_setequal(intersect(original_exports, exports), original_exports)
  expect_true(all(vapply(original_exports, function(name) {
    exists(name, envir = asNamespace("dsVertClient"), mode = "function",
           inherits = FALSE)
  }, logical(1L))))

  maturity <- ds.vertMethodStatus()
  status <- maturity$status[match(original_exports, maturity$method)]
  expect_false(anyNA(status))
  expect_true(all(status == "promoted"))
})
