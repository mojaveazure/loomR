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

# A function to determine if a vector is a vector and not a list
#
# @param x An object
#
# @return TRUE if 'x' is a vector or a factor, otherwise FALSE
#
is.actual_vector <- function(x) {
  return((is.vector(x = x) || is.factor(x = x)) && !is.list(x = x))
}

# Check additions to /matrix
#
# @param x A list of vectors to add to /matrix
# @param n The number of genes needed in each cell
#
# @return 'x' as a list of vectors
#
check.matrix_data <- function(x, n) {
  # Coerce x into a list, where each
  # entry in the list is a new cell added
  if (is.actual_vector(x = x)) {
    x <- list(x)
  } else if (is.matrix(x = x) || is.data.frame(x = x)) {
    x <- as.list(x = x)
  }
  if (!is.list(x = x)) {
    stop("Matrix data must be a list of vectors")
  }
  # Ensure that each entry in the list is a vector of length n
  vector.check <- vapply(
    X = x,
    FUN = is.actual_vector,
    FUN.VALUE = logical(length = 1L)
  )
  if (!all(vector.check)) {
    stop('Each new cell added must be represented as a vector')
  }
  # Ensure each new cell added has data for the number of genes present
  for (i in 1:length(x = x)) {
    cell.add <- x[[i]]
    if (length(x = cell.add) > n) {
      stop(paste(
        "Cannot add genes to a loom file, the maximum number of genes allowed is",
        n
      ))
    } else if (length(x = cell.add) < n) {
      cell.add[(length(x = cell.add) + 1):n] <- NA
    }
    x[[i]] <- cell.add
  }
  return(x)
}

# Get the number of cells being added to /matrix
#
# @param x A list of vectors to add to /matrix
#
# @return The number of cells in x
#
nCells.matrix_data <- function(x) {
  return(length(x = x))
}

# Add missing cells to data added to /matrix
#
# @param x A list of vectors to add to /matrix
# @param n The number of genes each cell needs
# @param m2 The number of cells being added to the loom file
#
# @return 'x' with the proper number of cells
#
addCells.matrix_data <- function(x, n, m2) {
  if (length(x = x) < m2) {
    x[(length(x = x) + 1):m2] <- list(rep.int(x = NA, times = n))
  }
  return(x)
}

# Check additions to /layers
#
# @param x A list of matrices to add layers in /layers
# @param n The number of genes needed for each layer
# @param layers.names Names found in /layers
#
# @return 'x' as a list of matricies with 'n' rows for each layer present in /layers
#
check.layers <- function(x, n, layers.names) {
  # Coerce x into a list of matricies
  if (is.null(x = x)) {
    x <- list()
  } else if (is.matrix(x = x) || is.data.frame(x = x)) {
    x <- list(as.matrix(x = x))
  }
  if (!is.list(x = x)) {
    stop("Layers data must be a list of matricies")
  }
  # Ensure we have enough layer additions for each layer
  # Manage named lists, taking only those with the same name as in /layers
  # or is named with an empty string
  if (!is.null(x = names(x = x))) {
    x.use <- which(x = names(x = x) %in% layers.names | names(x = x) == '')
    x <- x[x.use]
  }
  if (length(x = x) > length(x = layers.names)) {
    stop("Cannot add more layers than already present")
  } else if (length(x = x) < length(x = layers.names)) {
    x[(length(x = x) + 1):length(x = layers.names)] <- data.frame(rep.int(x = NA, times = n))
  }
  # Ensure that we have all genes needed
  for (i in 1:length(x = x)) {
    layer <- x[[i]]
    if (is.data.frame(x = layer)) {
      layer <- as.matrix(x = layer)
    } else if (is.actual_vector(x = layer)) {
      layer <- matrix(data = layer, ncol = 1)
    }
    if (!is.matrix(x = layer)) {
      stop("Layers data must be a list of matrices")
    }
    if (nrow(x = layer) > n) {
      stop(paste(
        "Cannot add genes to a loom file, the maximum number of genes allowed is",
        n
      ))
    } else if (nrow(x = layer) < n) {
      layer <- as.data.frame(x = layer)
      layer[(nrow(x = layer) + 1):n, ] <- NA
      layer <- as.matrix(x = layer)
    }
    x[[i]] <- layer
  }
  # Set names
  x.unnamed <- which(x = !(names(x = x) %in% layers.names))
  if (length(x = x.unnamed) == 0) {
    x.unnamed <- 1:length(x = x)
  }
  names.unused <- which(x = !(layers.names %in% names(x = x)))
  names(x = x)[x.unnamed] <- layers.names[names.unused]
  return(x)
}

# Get the number of cells being added to /layers
#
# @param x A list of matricies to add to /layers
#
# @return The number of cells within each matrix
#
nCells.layers <- function(x) {
  return(vapply(X = x, FUN = ncol, FUN.VALUE = integer(length = 1L)))
}

# Add missing cells to data added to /matrix
#
# @param x A list of matricies to add to /layers
# @param n The number of genes each cell needs
# @param m2 The number of cells being added to the loom file
#
# @return 'x' with the proper number of cells
#
addCells.layers <- function(x, n, m2) {
  layers.extend <- vapply(X = x, FUN = ncol, FUN.VALUE = integer(length = 1L))
  layers.extend <- which(x = layers.extend != m2)
  for (i in layers.extend) {
    layer <- x[[i]]
    layer.new <- matrix(nrow = n, ncol = m2)
    layer.new[, 1:ncol(x = layer)] <- layer
    x[[i]] <- layer.new
    gc(verbose = FALSE)
  }
  return(x)
}

# Check additions to /col_attrs
#
# @param x A list of vectors to add to /col_attrs
# @param attrs.names Names of attributes found in /col_attrs
#
# @return 'x' as a list of vectors for each attribute found in /col_attrs
#
check.col_attrs <- function(x, attrs.names) {
  # Coerce x into a list of vectors
  if (is.null(x = x)) {
    x <- list()
  } else if (is.actual_vector(x = x)) {
    x <- list(x)
  } else if (is.matrix(x = x) || is.data.frame(x = x)) {
    x <- as.list(x = x)
  }
  if (!is.list(x = x)) {
    stop("Attribute data must be a list of vectors")
  }
  # Ensure we have enough attribute additions for each col_attr
  # Manage named lists, taking only those with the same name as in /col_attrs
  # or is named with an empty string
  if (!is.null(x = names(x = x))) {
    x.use <- which(x = names(x = x) %in% attrs.names | names(x = x) == '')
    x <- x[x.use]
  }
  if (length(x = x) > length(x = attrs.names)) {
    stop("Cannot add more column attributes than already present")
  } else if (length(x = x) < length(x = attrs.names)) {
    x[(length(x = x) + 1):length(x = attrs.names)] <- NA
  }
  if (!all(vapply(X = x, FUN = is.actual_vector, FUN.VALUE = logical(length = 1L)))) {
    stop("Attribute data must be a list of vectors")
  }
  # Set names
  x.unnamed <- which(x = !(names(x = x) %in% attrs.names))
  if (length(x = x.unnamed) == 0) {
    x.unnamed <- 1:length(x = x)
  }
  names.unused <- which(x = !(attrs.names %in% names(x = x)))
  names(x = x)[x.unnamed] <- attrs.names[names.unused]
  return(x)
}

# Get the number of cells being added to /col_attrs
#
# @param x A list of vectors to add to /col_attrs
#
# @return The number of cells for each attribute
#
nCells.col_attrs <- function(x) {
  return(vapply(X = x, FUN = length, FUN.VALUE = integer(length = 1L)))
}

# Add missing cells to data added to /col_attrs
#
# @param x A list of vectors to add to /col_attrs
# @param m2 The number of cells being added to the loom file
#
# @return 'x' with the proper number of cells
#
addCells.col_attrs <- function(x, m2) {
  attrs.extend <- vapply(X = x, FUN = length, FUN.VALUE = integer(length = 1L))
  attrs.extend <- which(x = attrs.extend != m2)
  for (i in attrs.extend) {
    attr <- x[[i]]
    attr <- c(attr, rep.int(x = NA, times = m2 - length(x = attr)))
    x[[i]] <- attr
  }
  return(x)
}

# Create a progress bar
#
# @return A progress bar
#
#' @importFrom utils txtProgressBar
#
new.pb <- function() {return(txtProgressBar(style = 3, char = '='))}
