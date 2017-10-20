#' @import Methods
#' @import h5
NULL

#' A class for loom
#'
#' @exportClass loom
#'
loom <- setClass(
  Class = 'loom',
  contains = 'H5File'
)