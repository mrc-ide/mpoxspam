test_that("can run the gpu model on the cpu", {
  skip_if_not_installed("odin.dust")
  skip_if_not_installed("lostturnip")

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = FALSE)
  expect_equal(gen$public_methods$name(), "model")

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


test_that("Can run the gpu compare on the cpu", {
  skip_if_not_installed("odin.dust")
  skip_if_not_installed("lostturnip")

  gen <- compile_gpu(gpu = FALSE, gpu_generate = TRUE, verbose = FALSE)

  pars <- reference_pars()
  data <- filter_data(pars$dt)

  np <- 10
  mod_cpu <- gen$new(pars, 0, np, seed = 1L)
  mod_gpu <- gen$new(pars, 0, np, seed = 1L, gpu_config = 0)
  mod_cpu$set_data(dust::dust_data(data, "time_end"))
  mod_gpu$set_data(dust::dust_data(data, "time_end"))

  for (i in seq_len(nrow(data))) {
    y_gpu <- mod_gpu$run(data$day_end[[i]])
    y_cpu <- mod_cpu$run(data$day_end[[i]])
    expect_equal(y_gpu, y_cpu)
    expect_equal(mod_gpu$compare_data(),
                 mod_cpu$compare_data())
  }
})
