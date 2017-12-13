# Validate a loom object
#
# @param object A loom object
#
# @return None, errors out if object is an invalid loom connection
#
# @seealso \code{\link{loom-class}}
#
validateLoom <- function(object) {
  if (!inherits(x = object, what = 'loom')) {
    stop("No need to validate a non-loom object")
  }
  # A loom file is a specific HDF5
  # We need a dataset in /matrix that's a two-dimensional dense matrix
  root.datasets <- list.datasets(object = object, path = '/', recursive = FALSE)
  if (length(x = root.datasets) != 1) {
    stop("There can only be one dataset at the root of the loom file")
  }
  if (root.datasets != 'matrix') {
    stop("The root dataset must be called '/matrix'")
  }
  # There must be groups called '/col_attrs', '/row_attrs', and '/layers'
  required.groups <- c('row_attrs', 'col_attrs', 'layers')
  dim.matrix <- object[[root.datasets]]$dims # Columns x Rows
  names(x = dim.matrix) <- required.groups[c(2, 1)]
  root.groups <- list.groups(object = object, path = '/', recursive = FALSE)
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
  unlist(x = sapply(
    X = required.groups[1:2],
    FUN = function(group) {
      if (length(x = list.groups(object = object[[group]], recursive = FALSE)) > 0) {
        stop(paste("Group", group, "cannot have subgroups"))
      }
      if (length(x = list.attributes(object = object[[group]])) > 0) {
        stop(paste("Group", group, "cannot have subattributes"))
      }
      for (dataset in list.datasets(object = object[[group]])) {
        if (object[[paste(group, dataset, sep = '/')]]$dims != dim.matrix[group]) {
          stop(paste("All datasets in group", group, "must be of length", required.groups[group]))
        }
      }
    }
  ))
  for (dataset in list.datasets(object = object[['/layers']])) {
    if (any(object[[paste('layers', dataset, sep = '/')]]$dims != dim.matrix)) {
      stop(paste("All datasets in '/layers' must be", dim.matrix[1], 'by', dim.matrix[2]))
    }
  }
}
