compile_gpu <- function(..., real_type = "float", gpu = TRUE) {
  ## Avoid weird errors by making sure these versions are ok here;
  ## later we can drop this and rely on the DESCRIPTION.
  loadNamespace("odin.dust", versionCheck = list(op = ">=", version = "0.2.28"))
  loadNamespace("lostturnip", versionCheck = list(op = ">=", version = "0.1.2"))

  path <- system.file("odin/m4_2.R", package = "monkeyspam", mustWork = TRUE)
  odin.dust::odin_dust_(path, real_type = real_type, gpu = gpu, ...)
}
