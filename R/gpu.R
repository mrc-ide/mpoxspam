compile_gpu <- function(..., real_type = "float", gpu = TRUE) {
  path <- system.file("odin/m4_2.R", package = "monkeyspam", mustWork = TRUE)
  odin.dust::odin_dust_(path, real_type = real_type, gpu = gpu, ...)
}
