evaluate_pgf <- function(f, x) {
  vapply(f, function(fi) fi(x), numeric(1))
}
