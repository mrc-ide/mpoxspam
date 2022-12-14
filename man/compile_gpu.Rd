% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpu.R
\name{compile_gpu}
\alias{compile_gpu}
\title{Compile model for GPU}
\usage{
compile_gpu(..., real_type = "float", gpu = TRUE)
}
\arguments{
\item{...}{Additional arguments to \code{\link[odin.dust:odin_dust]{odin.dust::odin_dust()}}}

\item{real_type}{Name of real type to use, probably "float"}

\item{gpu}{Configuration, passed through to \code{\link[odin.dust:odin_dust]{odin.dust::odin_dust()}}}
}
\description{
Compile model for GPU
}
