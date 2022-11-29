##' Compile model for GPU
##'
##' @title Compile model for GPU
##'
##' @param ... Additional arguments to [odin.dust::odin_dust()]
##'
##' @param real_type Name of real type to use, probably "float"
##'
##' @param gpu Configuration, passed through to [odin.dust::odin_dust()]
##'
##' @export
compile_gpu <- function(..., real_type = "float", gpu = TRUE) {
  ## Avoid weird errors by making sure these versions are ok here;
  ## later we can drop this and rely on the DESCRIPTION.
  loadNamespace("odin.dust", versionCheck = list(op = ">=", version = "0.2.28"))
  loadNamespace("lostturnip", versionCheck = list(op = ">=", version = "0.1.2"))

  path <- system.file("odin/model.R", package = "mpoxspam", mustWork = TRUE)
  odin.dust::odin_dust_(path, real_type = real_type, gpu = gpu, ...)
}
