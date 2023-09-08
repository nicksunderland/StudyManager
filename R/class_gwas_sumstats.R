#' Title
#'
#' @slot parameter numeric.
#'
#' @return object
#' @export
#' @importFrom methods callNextMethod
#'
GWASsumstats <- setClass(
  Class = "GWASsumstats",
  contains = "Study",
  slots = list(
    parameter = "numeric"
  ),
  prototype = list(
    parameter = numeric()
  )
)

setMethod(
  f = "initialize",
  signature = "GWASsumstats",
  definition = function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@data_files <- create_data_files(.Object)
    return(.Object)
  }
)

