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

#' Validate a loom object
#'
#' @param object A loom object
#'
#' @return TRUE if a valid loom object
#'
validateLoom <- function(object) {
  # A loom file is a specific HDF5
  # We need a dataset in /matrix that's a two-dimensional dense matrix
  root.datasets <- list.datasets(.Object = object, path = '/', recursive = FALSE)
  if (length(x = root.datasets) != 1) {
    stop("There can only be one dataset at the root of the loom file")
  }
  if (root.datasets != '/matrix') {
    stop("The root dataset must be called '/matrix'")
  }
  dim.matrix <- object[root.datasets]@dim # Rows x Columns
  # There must be groups called '/col_attrs', '/row_attrs', and '/layers'
  required.groups <- c('/col_attrs', '/row_attrs', '/layers')
  root.groups <- list.groups(.Object = object, path = '/', recursive = FALSE)
  group.msg <- paste0(
    "There can only be three groups in the loom file: '",
    paste(required.groups, collapse = "', '"),
    "'"
  )
  if (length(x = root.groups) != 3) {
    stop(group.msg)
  }
  if (!all(required.groups %in% root.groups)) {
    stop(group.msg)
  }
  vapply(
    X = required.groups[1:2],
    FUN = function(group) {
      if (length(x = list.groups(.Object = object, path = group, recursive = FALSE)) > 0) {
        stop(paste("Group", group, "cannot have subgroups"))
      }
      if (length(x = list.attributes(.Object = object[group])) > 0) {
        stop(paste("Group", group, "cannot have subattributes"))
      }
      for (dataset in list.datasets(.Object = object, path = group)) {
        break
      }
    }
  )
}
