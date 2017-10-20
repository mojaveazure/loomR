#' @import h5
#' @importFrom methods setClass setMethod setGeneric callNextMethod
NULL

#' A class for loom
#'
#' @name loom-class
#' @rdname loom-class
#' @exportClass loom
#'
loom <- setClass(
  Class = 'loom',
  slots = c(
    version = 'ANY'
  ),
  contains = 'H5File'
)

#' @importFrom utils packageVersion
#'
setMethod(
  f = 'initialize',
  signature = 'loom',
  definition = function(.Object, name, mode = 'a') {
    .Object <- callNextMethod(
      .Object,
      name = name,
      mode = mode
    )
    .Object@version <- packageVersion(pkg = 'loom')
    return(.Object)
  }
)
