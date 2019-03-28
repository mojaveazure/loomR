# Generate chunk points
#
# @param data.size How big is the data being chunked
# @param chunk.size How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
chunkPoints <- function(data.size, chunk.size) {
  return(vapply(
    X = 1L:ceiling(data.size / chunk.size),
    FUN = function(i) {
      return(c(
        start = (chunk.size * (i - 1L)) + 1L,
        end = min(chunk.size * i, data.size)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}

# Convert sparse matrix pointers to column indicies
#
# A function to get the column (j) indices of a sparse matrix from the pointers (p)
#
# @param p A vector of sparse matrix pointers
#
# @return A vector of column (j) indices
# @author Josh O'Brien
#
# @examples
# dat <- c(0, 0, 1, 4, 0, 2, 0, 9, 0)
# smat <- Matrix::Matrix(data = dat, nrow = 3, sparse = TRUE)
# PointerToIndex(p = smat@p)
#
PointerToIndex <- function(p) {
  # From StackOverflow
  # https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix
  dp <- diff(x = p)
  j <- rep.int(x = seq_along(along.with = dp), times = dp)
  return(j)
}

# Try to convert a character vector to bytes
#
# @param x A character vector
#
# @return the value of said string in bytes
#
#' @importFrom stringr str_extract
#
charToBytes <- function(x) {
  x <- as.character(x = x)
  patterns <- c('^\\d+[A-Z]?b?$', '^\\d\\.?\\d*e[\\+-]\\d+$')
  x <- unlist(x = sapply(
    X = patterns,
    FUN = grep,
    USE.NAMES = FALSE,
    x = x,
    perl = TRUE,
    ignore.case = TRUE,
    value = TRUE
  ))
  if (length(x = x) < 1) {
    stop(paste(
      "No valid character strings passed\n",
      "values in 'x' must pass either the",
      paste(patterns, collapse = '\n or '),
      "regex pattern"
    ))
  }
  x <- tolower(x = x)
  bytes <- vapply(
    X = x,
    FUN = function(i) {
      exp <- FALSE
      mult.factor <- if (grepl(pattern = 'k', x = i)) {
        1e3
      } else if (grepl(pattern = 'm', x = i)) {
        1e6
      } else if (grepl(pattern = 'g', x = i)) {
        1e9
      } else {
        1
      }
      if (grepl(pattern = patterns[2], x = i)) {
        return(as.numeric(x = i))
      }
      dig <- as.numeric(x = str_extract(string = i, pattern = '\\d+'))
      return(dig * mult.factor)
    },
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  )
  return(bytes)
}

# Get HDF5 data types
#
# @param x An R object or string describing HDF5 datatype
#
# @return The corresponding HDF5 data type
#
# @ rdname getDtype
#
#' @import hdf5r
#
# @seealso \code\link{hdf5r::h5types}
#
getDtype <- function(x) {
  return(switch(
    EXPR = class(x = x),
    'numeric' = h5types$float,
    'integer' = h5types$int,
    'character' = H5T_STRING$new(size = max(nchar(x))),
    'logical' = H5T_LOGICAL$new(),
    stop(paste("Unknown data type:", class(x = x)))
  ))
}

# @describeIn getDtype A version of getDtype that works specifically for hdf5r types,
# useful for getDtype2(x = class(x = object[['dataset']]$get_type())[1])
#
getDtype2 <- function(x) {
  return(getDtype(x = switch(
    EXPR = x,
    'H5T_FLOAT' = numeric(),
    'H5T_INTEGER' = integer(),
    'H5T_STRING' = character(),
    'H5T_LOGICAL' = logical(),
    stop(paste("Unkown data type:", x))
  )))
}

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
  required.groups <- c('row_attrs', 'col_attrs', 'layers', 'row_graphs', 'col_graphs')
  dim.matrix <- object[['matrix']]$dims # Columns x Rows
  names(x = dim.matrix) <- required.groups[c(2, 1)]
  root.groups <- list.groups(object = object, path = '/', recursive = FALSE)
  group.msg <- paste0(
    paste("There can only be", length(x = required.groups), "groups in the loom file: '"),
    paste(required.groups, collapse = "', '"),
    "'"
  )
  reopen.msg <- paste(
    group.msg,
    "Reopen in 'r+' mode to automatically add missing groups",
    sep = '\n'
  )
  if (any(grepl(pattern = 'edges', x = root.groups))) {
    if (object$mode != 'r') {
      message("Moving edge groups to graph groups to conform to loom v2.0.1")
      edges <- grep(pattern = 'edges', x = root.groups, value = TRUE)
      for (group in edges) {
        graph <- gsub(pattern = 'edges', replacement = 'graphs', x = group)
        object$link_move_to(dst_loc = object, dst_name = graph, src_name = group)
      }
      root.groups <- list.groups(object = object, path = '/', recursive = FALSE)
    } else {
      object$close_all()
      stop(reopen.msg)
    }
  }
  if (length(x = root.groups) > length(x = required.groups)) {
    stop(group.msg)
  } else if (length(x = root.groups) < length(x = required.groups)) {
    if (all(root.groups %in% required.groups)) {
      if (object$mode != 'r') {
        missing.groups <- required.groups[!(required.groups %in% root.groups)]
        for (group in missing.groups) {
          object$create_group(name = group)
        }
        root.groups <- list.groups(object = object, path = '/', recursive = FALSE)
      } else {
        object$close_all()
        stop(reopen.msg)
      }
    } else {
      stop(group.msg)
    }
  }
  if (!all(required.groups %in% root.groups)) {
    stop(group.msg)
  }
  # Check row and column attributes
  for (group in required.groups[1:2]) {
    # No subgroups
    if (length(x = list.groups(object = object[[group]], recursive = FALSE)) > 0) {
      stop(paste("Group", group, "cannot have subgroups"))
    }
    # All datasets must have their first (last) dimmension equal to M(row) or N(column)
    for (dataset in list.datasets(object = object[[group]])) {
      dataset.dim <- object[[group]][[dataset]]$dims
      dataset.dim <- dataset.dim[length(x = dataset.dim)]
      if (dataset.dim != dim.matrix[group]) {
        print(dataset)
        print(object[[group]][[dataset]])
        print(dim.matrix)
        stop("All datasets in group ", group, " must be of length ", dim.matrix[group])
      }
    }
  }
  # Check row and column graphs
  graph.groups <- grep(pattern = 'graphs', x = required.groups, value = TRUE)
  graph.msg <- "There can only be three datasets in a graph: 'a', 'b', and 'w'"
  for (group in graph.groups) {
    group.datasets <- list.datasets(object = object[[group]], recursive = FALSE)
    if (length(x = group.datasets) > 0) {
      stop(paste("All datasets in", group, "must be in a graph group"))
    }
    graphs <- list.groups(object = object[[group]], full.names = TRUE, recursive = FALSE)
    for (graph in graphs) {
      graph.datasets <- list.datasets(object = object[[graph]], full.names = TRUE)
      if (length(x = graph.datasets) != 3) {
        stop(graph.msg)
      }
      if (!all(basename(path = graph.datasets) %in% c('a', 'b', 'w'))) {
        stop(graph.msg)
      }
      graph.lengths <- lapply(
        X = graph.datasets,
        FUN = function(dset) {
          return(object[[dset]]$dims)
        }
      )
      if (length(x = unique(x = graph.lengths)) != 1) {
        stop("All graph datasets must be the same length")
      }
    }
  }
  # Check layers
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
# @param n The number of features needed in each cell
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
  # Ensure each new cell added has data for the number of features present
  for (i in 1:length(x = x)) {
    cell.add <- x[[i]]
    if (length(x = cell.add) > n) {
      stop(paste(
        "Cannot add features to a loom file, the maximum number of features allowed is",
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
# @param n The number of features each cell needs
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
# @param n The number of features needed for each layer
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
  # Ensure that we have all features needed
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
        "Cannot add features to a loom file, the maximum number of features allowed is",
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
# @param n The number of features each cell needs
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
newPB <- function() {
  return(txtProgressBar(style = 3, char = '='))
}

# Break a vector into a list of consecutive values
#
# @param x An integer vector
# @param max.length Maximum length of consecutive values;
# set to NULL if no maximum is desired
#
# @return A list where each value is a run of consecutive values from x
#
#' @importFrom stats na.omit
#
break.consecutive <- function(x, max.length = NULL, min.length = NULL) {
  # Coerce x to integers, find unique values, sort it
  x <- sort(x = unique(x = na.omit(object = as.integer(x = x))))
  if (is.null(x = max.length)) {
    max.length <- length(x = x)
  }
  if (is.null(x = min.length)) {
    min.length <- 1L
  }
  # Functions to allocate a vector and increment a counter
  alloc <- function() {return(vector(mode = 'integer', length = max.length))}
  inc <- function(x) {return(x + 1L)}
  # Results list, counters, and holding vector, all preallocated for speed
  consecutive <- vector(mode = 'list', length = length(x = x))
  index <- 1L
  counter <- 0L
  run <- alloc()
  # For every value in x, if it's consecutive, add it to the holding vector
  # Otherwise, add the run to the list of consecutive values
  for (i in x) {
    if (counter == max.length || (counter != 0L && counter >= min.length && i != run[counter] + 1L)) {
      consecutive[[index]] <- run[1:counter]
      index <- inc(x = index)
      run <- alloc()
      counter <- 0L
    }
    counter <- inc(x = counter)
    run[counter] <- i
  }
  consecutive[[index]] <- run[1:counter]
  consecutive <- consecutive[1:index]
  return(consecutive)
}

# Calculate nUMI and nFeature
#
# @param object An object
# @param chunk.size Number of cells to chunk over
# @param is.expr Expression threshold for 'detected' feature. For most datasets, particularly UMI datasets, will be set to 0 (default).
# If not, when initializing, this should be set to a level based on pre-normalized counts (i.e. require at least 5 counts to be treated as expresesd).
# @param verbose Display a progress bar?
#
# @return Stores nUMI in \code{col_attrs/nUMI}; stores nFeature in \code{col_attrs/nFeature}; will overwrite any existing datasets at these locations
#
#' @importFrom Matrix rowSums
#
calc.umi <- function(
  object,
  chunk.size = 1000,
  is.expr = 0,
  verbose = TRUE
) {
  if (verbose) {
    message("Calculating number of UMIs per cell")
  }
  object$apply(
    name = 'col_attrs/nCount',
    FUN = rowSums,
    MARGIN = 2,
    chunk.size = chunk.size,
    dataset = 'matrix',
    overwrite = TRUE,
    verbose = verbose
  )
  if (verbose) {
    message("Calculating number of features expressed per cell")
    message("Using a threshold of ", is.expr, " for feature expression")
  }
  object$apply(
    name = 'col_attrs/nFeature',
    FUN = function(mat, is.expr) {
      return(rowSums(x = mat > is.expr))
    },
    MARGIN = 2,
    chunk.size = chunk.size,
    dataset = 'matrix',
    overwrite = TRUE,
    verbose = verbose,
    is.expr = is.expr
  )
  invisible(x = object)
}
