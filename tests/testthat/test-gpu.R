test_that("can run the gpu model on the cpu", {
  skip_if_not_installed("odin.dust")
  skip_if_not_installed("lostturnip")

  ## Looks like next step:
  ## - get dust, odin.dust to generate C++14 configuration
  ## - just remove the grepl check in odin.dust?
  ## - work out why we get some conflicting unpacking of newI, newIseed and time - this looks like the compare?

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = FALSE)
  expect_equal(gen$public_methods$name(), "m4_2")

  pars <- reference_pars()

  mod_cpu <- gen$new(pars, 0, 5, seed = 1L)
  mod_gpu <- gen$new(pars, 0, 5, seed = 1L, gpu_config = 0)

  end <- 140
  info <- mod_cpu$info()

  expect_equal(mod_gpu$info(), mod_cpu$info())

  res_cpu <- mod_cpu$run(end)
  res_gpu <- mod_gpu$run(end)

  expect_equal(res_gpu, res_cpu)
})
