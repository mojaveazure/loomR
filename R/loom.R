#' @import hdf5r
#' @importFrom R6 R6Class
NULL

#' A class for loom files
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @aliases loom
#' @return Object of \code{\link{R6::R6Class}} to generate \code{loom} objects
#' @format An \code{\link{R6::R6Class}} object
#' @seealso \code{\link{hdf5r::H5File}}
#'
#' @field version Version of loomR object was created under
#' @field shape Shape of \code{/matrix} in columns (cells) by rows (genes)
#' @field chunksize Chunks set for this dataset in columns (cells) by rows (genes)
#' @field matrix The main data matrix, stored as columns (cells) by rows (genes)
#' @field layers Additional data matricies, the same shape as \code{/matrix}
#' @field col.attrs Extra information about cells
#' @field row.attrs Extra information about genes
#'
#' @section Methods:
#' \describe{
#'   \item{\code{add.layer(layer)}}{Add a data layer to this loom file, must be in column (cells) by row (genes) orientation}
#'   \item{\code{add.attribute(attribute, MARGIN)}}{
#'     Add extra information to this loom file where
#'     \code{attribute} is a named list where each element is a vector that is as long as one dimension of \code{/matrix} and
#'     \code{MARGIN} is either 1 for cells or 2 for genes
#'   }
#'   \item{\code{add.row.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 2)}}
#'   \item{\code{add.col.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 1)}}
#'   \item{\code{add.meta.data(meta.data)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 1)}}
#'   \item{\code{batch.scan(chunk.size, MARGIN, index.use, dataset.use, force.reset)}, \code{batch.next()}}{
#'     Scan a dataset in the loom file from \code{index.use[1]} to \code{index.use[2]}, iterating by \code{chunk.size}.
#'     \code{dataset.use} can be the name, not \code{group/name}, unless the name is present in multiple groups.
#'     Pass \code{MARGIN = 1} to iterate on cells or \code{MARGIN = 2} to iterate on genes for 'matrix' or any dataset in 'layers'.
#'     To force reevaluation of the iterator object, pass \code{force.reset = TRUE}.
#'     \code{MARGIN} does not need to be set for datasets in 'row_attrs' or 'col_attrs'.
#'     \code{chunk.size} defaults to \code{self$chunksize}, \code{MARGIN} defaults to 1,
#'     \code{index.use} defaults to \code{1:self$shape[MARGIN]}, \code{dataset.use} defaults to 'matrix'
#'   }
#'   \item{\code{apply(name, FUN, MARGIN, chunk.size, dataset.use, display.progress, ...)}}{
#'     Apply a function over a dataset within the loom file, stores the results in the loom file.
#'     \code{name} must be the full name of the dataset ('group/name').
#'     \code{apply} will always use the entire dataset specified in \code{dataset.use}
#'   }
#'   \item{\code{map(FUN, MARGIN, chunk.size, index.use, dataset.use, display.progress, expected, ...)}}{
#'     Map a function onto a dataset within the loom file, returns the result into R.
#'     The result will default to the shape of the dataset used; to change pass either 'vector' or 'matrix' to \code{expected}.
#'   }
#' }
#'
#' @importFrom iterators nextElem
#' @importFrom itertools hasNext ihasNext ichunk
#' @importFrom utils packageVersion txtProgressBar setTxtProgressBar
#'
#' @export
#'
loom <- R6Class(
  # The loom class
  # Based on the H5File class from hdf5r
  # Not clonable (no loom$clone method), doesn't make sense since all data is on disk, not in memory
  # Yes to portability, other packages can subclass the loom class
  # Class is locked, other fields and methods cannot be added
  classname = 'loom',
  inherit = hdf5r::H5File,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  # Public fields and methods
  # See above for documentation
  public = list(
    # Fields
    version = NULL,
    shape = NULL,
    chunksize = NULL,
    matrix = NULL,
    layers = NULL,
    col.attrs = NULL,
    row.attrs = NULL,
    # Methods
    initialize = function(
      filename = NULL,
      mode = c('a', 'r', 'r+', 'w', 'w-'),
      ...
    ) {
      # If the file exists, run validation steps
      do.validate <- file.exists(filename) && !(mode %in% c('w', 'w+'))
      super$initialize(filename = filename, mode = mode, ...)
      if (do.validate) {
        # Run the validation steps
        validateLoom(object = self)
        # Store /matrix and the shape of /matrix
        self$matrix <- self[['matrix']]
        self$shape <- self[['matrix']]$dims
        # Store the chunk size
        chunks <- h5attr(x = self, which = 'chunks')
        chunks <- gsub(pattern = '(', replacement = '', x = chunks, fixed = TRUE)
        chunks <- gsub(pattern = ')', replacement = '', x = chunks, fixed = TRUE)
        chunks <- unlist(x = strsplit(x = chunks, split = ','))
        self$chunksize <- as.integer(x = chunks)
        # Store version information
        self$version <- as.character(x = tryCatch(
          # Try getting a version
          # If it doesn't exist, can we write to the file?
          # If so, store the version as this version of loomR
          # Otherwise, store the version as NA_character_
          expr = h5attr(x = self, which = 'version'),
          error = function(e) {
            if (mode != 'r') {
              version <- packageVersion(pkg = 'loomR')
              h5attr(x = self, which = 'version') <- version
            } else {
              version <- NA_character_
            }
            return(version)
          }
        ))
        # Load layers
        private$load_layers()
        # Load attributes
        private$load_attributes(MARGIN = 1) # Cells (col_attrs)
        private$load_attributes(MARGIN = 2) # Genes (row_attrs)
      } else {
        # Assume new HDF5 file
        self$version <- as.character(x = packageVersion(pkg = 'loomR'))
      }
    },
    # Addding attributes and layers
    add.layer = function(layers) {
      if (self$mode == 'r') {
        stop("Cannot add a layer in read-only mode")
      }
      # Value checking
      if (!is.list(x = layers) || is.null(x = names(x = layers))) {
        stop("'layers' must be a named list")
      }
      if (is.null(x = self$shape)) {
        stop(private$err_msg)
      }
      # Add layers
      for (i in 1:length(x = layers)) {
        if (!is.matrix(x = layers[[i]])) {
          layers[[i]] <- as.matrix(x = layers[[i]])
        }
        do.transpose <- FALSE
        if (any(dim(x = layers[[i]]) != self$shape)) {
          if (all(rev(x = dim(x = layers[[i]])) == self$shape)) {
            do.transpose <- TRUE
          } else {
            stop(paste(
              "All layers must have",
              self$shape[1],
              "rows for cells and",
              self$shape[2],
              "columns for genes"
            ))
          }
        }
        if (do.transpose) {
          self[['layers', names(x = layers)[i]]] <- t(x = layers[[i]])
        } else {
          self[['layers', names(x = layers)[i]]] <- layers[[i]]
        }
      }
      self$flush()
      private$load_layers()
      invisible(x = self)
    },
    add.attribute = function(attribute, MARGIN) {
      if (self$mode == 'r') {
        stop("Cannot add attributes in read-only mode")
      }
      # Value checking
      if (!is.list(x = attribute) || is.null(x = names(x = attribute))) {
        stop("'attribute' must be a named list")
      }
      if (length(x = attribute) > 1) {
        for (i in attribute) {
          if (!is.vector(x = attribute)) {
            stop("All attributes must be one-dimensional vectors")
          }
        }
      }
      if (length(x = which(x = names(x = attribute) != '')) != length(x = attribute)) {
        stop("Not all attributes had names provided")
      }
      if (!MARGIN %in% c(1, 2)) {
        stop("'MARGIN' must be 1 or 2")
      }
      # Add the attributes as datasets for our MARGIN's group
      if (is.null(x = self$shape)) {
        stop(private$err_msg)
      }
      grp.name <- c('col_attrs', 'row_attrs')[MARGIN]
      grp <- self[[grp.name]]
      for (i in 1:length(x = attribute)) {
        if (length(attribute[[i]]) != self$shape[MARGIN])
          stop(paste(
            "All",
            switch(EXPR = MARGIN, '1' = 'cell', '2' = 'gene'),
            "attributes must be of length",
            self$shape[MARGIN]
          ))
        grp[[names(x = attribute)[i]]] <- attribute[[i]]
      }
      self$flush()
      gc(verbose = FALSE)
      # Load the attributes for this margin
      private$load_attributes(MARGIN = MARGIN)
      invisible(x = self)
    },
    add.row.attribute = function(attribute) {
      self$add.attribute(attribute = attribute, MARGIN = 2)
      invisible(x = self)
    },
    add.col.attribute = function(attribute) {
      self$add.attribute(attribute = attribute, MARGIN = 1)
      invisible(x = self)
    },
    add.meta.data = function(meta.data) {
      self$add.col.attribute(attribute = meta.data)
      invisible(x = self)
    },
    # Chunking functions
    batch.scan = function(
      chunk.size = NULL,
      MARGIN = 1,
      index.use = NULL,
      dataset.use = 'matrix',
      force.reset = FALSE
    ) {
      if (is.null(x = private$it) || !grepl(pattern = dataset.use, x = private$iter.dataset) || force.reset) {
        # Check the existence of the dataset
        private$iter.dataset <- grep(
          pattern = dataset.use,
          x = list.datasets(object = self),
          value = TRUE
        )
        if (length(x = private$iter.dataset) != 1) {
          stop(paste0("Cannot find dataset '", dataset.use, "' in the loom file"))
        }
        if (grepl(pattern = 'col_attrs', x = private$iter.dataset)) {
          MARGIN <- 1
        } else if (grepl(pattern = 'row_attrs', x = private$iter.dataset)) {
          MARGIN <- 2
        }
        # Check the margin
        if (!(MARGIN %in% c(1, 2))) {
          stop("MARGIN must be 1 (cells) or 2 (genes)")
        } else {
          private$iter.margin <- MARGIN
        }
        if (is.null(x = chunk.size)) {
          chunk.size <- self$chunksize[private$iter.margin]
        }
        # Set the indices to use
        index.use <- private$iter_range(index.use = index.use)
        # Setup our iterator
        private$it <- ihasNext(iterable = ichunk(
          iterable = index.use[1]:index.use[2],
          chunkSize = chunk.size
        ))
        private$iter.index <- c(index.use[1], ceiling(x = index.use[2] / chunk.size))
      }
      # Return the times we iterate, this isn't important, we only need the length of this vector
      return(private$iter.index[1]:private$iter.index[2])
    },
    batch.next = function(return.data = TRUE) {
      # Ensure that we have a proper iterator
      if (!'hasNext.ihasNext' %in% suppressWarnings(expr = methods(class = class(x = private$it)))) {
        stop("Please setup the iterator with self$batch.scan")
      }
      # Do the iterating
      if (hasNext(obj = private$it)) {
        # Get the indices for this chunk
        chunk.indices <- unlist(x = nextElem(obj = private$it))
        if (return.data) {
          # If we're working with a matrix dataset, ensure chunking on proper axis
          if (private$iter.dataset == 'matrix' || grepl(pattern = 'layers', x = private$iter.dataset)) {
            to.return <- switch(
              EXPR = private$iter.margin,
              '1' = self[[private$iter.dataset]][chunk.indices, ],
              '2' = self[[private$iter.dataset]][, chunk.indices]
            )
          } else {
            # Otherwise, iterating over an attribute (1 dimensional)
            to.return <- self[[private$iter.dataset]][chunk.indices]
          }
        } else {
          to.return <- chunk.indices
        }
        # Determine if we reset the iterator
        if (!hasNext(obj = private$it)) {
          private$reset_batch()
        }
        return(to.return)
      } else {
        # Just in case
        private$reset_batch()
        return(NULL)
      }
    },
    apply = function(
      name,
      FUN,
      MARGIN = 1,
      chunk.size = NULL,
      dataset.use = 'matrix',
      display.progress = TRUE,
      ...
    ) {
      if (self$mode == 'r') {
        stop("Cannot write to disk in read-only mode")
      }
      if (!inherits(x = FUN, what = 'function')) {
        stop("FUN must be a function")
      }
      # Check that we're storing our results properly
      results.basename <- basename(path = name) # Basename of our results
      results.dirname <- gsub(pattern = '/', replacement = '', x = dirname(path = name)) # Groupname of our results
      dirnames <- c('col_attrs', 'row_attrs', 'layers') # Allowed group names
      if (name %in% list.datasets(object = self)) {
        stop(paste("A dataset with the name", name, "already exists!"))
      }
      # Checks datset, index, and MARGIN
      # Sets full dataset path in private$iter.dataset
      # Sets proper MARGIN in private$iter.margin
      batch <- self$batch.scan(
        chunk.size = chunk.size,
        MARGIN = MARGIN,
        dataset.use = dataset.use,
        force.reset = TRUE
      )
      MARGIN <- private$iter.margin
      dataset.use <- private$iter.dataset
      # Ensure that our group name is allowed
      name.check <- which(x = dirnames == results.dirname)
      if (!any(name.check)) {
        private$reset_batch()
        stop(paste(
          "The dataset must go into one of:",
          paste(dirnames, collapse = ', ')
        ))
      }
      # Check that our group matches our iteration
      # ie. To store in col_attrs, MARGIN must be 1
      if (name.check %in% c(1, 2) && name.check != private$iter.margin) {
        private$reset_batch()
        stop(paste(
          "Iteration must be over",
          c('cells', 'genes')[name.check],
          paste0("(MARGIN = ", name.check, ")"),
          "to store results in",
          paste0("'", dirnames[name.check], "'")
        ))
      }
      # Check how we store our results
      dataset.matrix <- ('matrix' %in% private$iter.dataset || grepl(pattern = 'layers', x = private$iter.dataset))
      results.matrix <- name.check == 3
      # Get a connection to the group we're iterating over
      group <- self[[results.dirname]]
      if (display.progress) {
        pb <- txtProgressBar(char = '=', style = 3)
      }
      # Have to initialize the dataset differently than
      first <- TRUE
      for (i in 1:length(x = batch)) {
        # Get the indices we're iterating over
        chunk.indices <- self$batch.next(return.data = FALSE)
        # Get the data and apply FUN
        chunk.data <- if (dataset.matrix) {
          switch(
            EXPR = MARGIN,
            '1' = self[[dataset.use]][chunk.indices, ],
            '2' = self[[dataset.use]][, chunk.indices]
          )
        } else {
          self[[private$iter.datset]][chunk.indices]
        }
        chunk.data <- FUN(chunk.data, ...)
        # If this is the first iteration
        # Initialize the dataset within group, set first to FALSE
        if (first) {
          group[[results.basename]] <- chunk.data
          first <- FALSE
        } else {
          # If we're writign to a matrix
          # Figure out which way we're writing the data
          if (results.matrix) {
            switch(
              EXPR = MARGIN,
              '1' = group[[results.basename]][chunk.indices, ] <- chunk.data,
              '2' = group[[results.basename]][, chunk.indices] <- chunk.data
            )
          } else {
            # Just write to the vector
            group[[results.basename]][chunk.indices] <- chunk.data
          }
        }
        if (display.progress) {
          setTxtProgressBar(pb = pb, value = i / length(x = batch))
        }
      }
      # Clean up and allow chaining
      private$reset_batch()
      # Load layers and attributes
      private$load_layers()
      private$load_attributes(MARGIN = 1) # Cells (col_attrs)
      private$load_attributes(MARGIN = 2) # Genes (row_attrs)
      invisible(x = self)
    },
    map = function(
      FUN,
      MARGIN = 1,
      chunk.size = NULL,
      index.use = NULL,
      dataset.use = 'matrix',
      display.progress = TRUE,
      expected = NULL,
      ...
    ) {
      if (!inherits(x = FUN, what = 'function')) {
        stop("FUN must be a function")
      }
      # Checks datset, index, and MARGIN
      # Sets full dataset path in private$iter.dataset
      # Sets proper MARGIN in private$iter.margin
      batch <- self$batch.scan(
        chunk.size = chunk.size,
        MARGIN = MARGIN,
        index.use = index.use,
        dataset.use = dataset.use,
        force.reset = TRUE
      )
      MARGIN <- private$iter.margin
      dataset.use <- private$iter.dataset
      # Check how we store our results
      # And what the shape of our dataset is
      results.matrix <- TRUE
      dataset.matrix <- TRUE
      if (grepl(pattern = 'col_attrs', x = private$iter.dataset)) {
        results.matrix <- FALSE
        dataset.matrix <- FALSE
      } else if (grepl(pattern = 'row_attrs', x = private$iter.dataset)) {
        results.matrix <- FALSE
        dataset.matrix <- FALSE
      }
      if (!is.null(x = expected)) {
        results.matrix <- switch(
          EXPR = expected,
          'vector' = FALSE,
          'matrix' = TRUE,
          stop("'expected' must be one of 'matrix', 'vector', or NULL")
        )
      }
      # Determine the shape of our results
      index.use <- private$iter_range(index.use = index.use)
      # Create our results holding object
      if (results.matrix) {
        switch(
          EXPR = private$iter.margin,
          '1' = results <- matrix(
            nrow = length(x = index.use[1]:index.use[2]),
            ncol = self$shape[2]
          ),
          '2' = results <- matrix(
            nrow = self$shape[1],
            ncol = length(x = index.use[1]:index.use[2])
          )
        )
      } else {
        results <- vector(length = length(x = index.use[1]:index.use[2]))
      }
      if (display.progress) {
        pb <- txtProgressBar(char = '=', style = 3)
      }
      for (i in 1:length(x = batch)) {
        chunk.indices <- self$batch.next(return.data = FALSE)
        chunk.data <- if (dataset.matrix) {
          switch(
            EXPR = MARGIN,
            '1' = self[[dataset.use]][chunk.indices, ],
            '2' = self[[dataset.use]][, chunk.indices]
          )
        } else {
          self[[dataset.use]][chunk.indices]
        }
        chunk.data <- FUN(chunk.data, ...)
        if (results.matrix) {
          if (MARGIN == 1) {
            results[chunk.indices, ] <- chunk.data
          } else if (MARGIN == 2) {
            results[, chunk.indices] <- chunk.data
          }
        } else {
          results[chunk.indices] <- chunk.data
        }
        if (display.progress) {
          setTxtProgressBar(pb = pb, value = i / length(x = batch))
        }
      }
      private$reset_batch()
      return(results)
    },
    # Functions that modify `/matrix'
    add.cells = function(matrix.data, attributes.data = NULL, layers.data = NULL) {
      lengths <- vector(
        mode = 'integer',
        length = 1 + length(x = attributes.data) + length(x = layers.data)
      )
      lengths[1] <- length(x = matrix.data)
      attributes.end <- 1 + length(x = attributes.data)
      if (attributes.end != 1) {
        lengths[2:attributes.end] <- vapply(
          X = attributes.data,
          FUN = length,
          FUN.VALUE = integer(length = 1L)
        )
      }
      if (attributes.end != length(x = lengths)) {
        lengths[(attributes.end + 1):length(x = lengths)] <- vapply(
          X = layers.data,
          FUN = length,
          FUN.VALUE = integer(length = 1L)
        )
      }
      return(lengths)
    },
    add.loom = function() {}
  ),
  # Private fields and methods
  # @field err_msg A simple error message if this object hasn't been created with loomR::create or loomR::connect
  # @field it Iterator object for batch.scan and batch.next
  # @field iter.dataset Dataset for iterating on
  # @field iter.margin Margin to iterate over
  # @field iter.index # Index values for iteration
  # \describe{
  #   \item{\code{load_attributes(MARGIN)}}{Load attributes of a given MARGIN into \code{self$col.attrs} or \code{self$row.attrs}}
  #   \item{\code{load_layers()}}{Load layers into \code{self$layers}}
  #   \item{\code{reset_batch()}}{Reset the batch iterator fields}
  #   \item{\code{iter_range(index.use)}}{Get the range of indices for a batch iteration}
  # }
  private = list(
    # Fields
    err_msg = "This loom object has not been created with either loomR::create or loomR::connect, please use these functions to create or connect to a loom file",
    it = NULL,
    iter.dataset = NULL,
    iter.margin = NULL,
    iter.index = NULL,
    # Methods
    load_attributes = function(MARGIN) {
      attribute <- switch(
        EXPR = MARGIN,
        '1' = 'col_attrs',
        '2' = 'row_attrs',
        stop('Invalid attribute dimension')
      )
      group <- self[[attribute]]
      attributes <- unlist(x = lapply(
        X = names(x = group),
        FUN = function(x) {
          d <- list(group[[x]])
          names(x = d) <- x
          return(d)
        }
      ))
      if (MARGIN == 1) {
        self$col.attrs <- attributes
      } else if (MARGIN == 2) {
        self$row.attrs <- attributes
      }
    },
    load_layers = function() {
      self$layers <- unlist(x = lapply(
        X = names(x = self[['layers']]),
        FUN = function(n) {
          d <- c(self[['layers', n]])
          names(x = d) <- n
          return(d)
        }
      ))
    },
    reset_batch = function() {
      private$it <- NULL
      private$iter.dataset <- NULL
      private$iter.margin <- NULL
      private$iter.index <- NULL
    },
    iter_range = function(index.use) {
      if (is.null(private$iter.margin)) {
        stop("Batch processing hasn't been set up")
      }
      if (is.null(x = index.use)) {
        # If no index was provided, use entire range for this margin
        index.use <- c(1, self$shape[private$iter.margin])
      } else if (length(x = index.use) == 1) {
        # If one index was provided, start at one and go to index
        index.use <- c(1, index.use)
      } else {
        # Use index.use[1] and index.use[2]
        index.use <- c(index.use[1], index.use[2])
      }
      # Ensure the indices provided fit within the range of the dataset
      index.use[1] <- max(1, index.use[1])
      index.use[2] <- min(index.use[2], self$shape[private$iter.margin])
      # Ensure that index.use[1] is greater than index.use[2]
      if (index.use[1] > index.use[2]) {
        stop(paste0(
          "Starting index (",
          index.use[1],
          ") must be lower than the ending index (",
          index.use[2],
          ")"
        ))
      }
      return(index.use)
    }
  )
)

#' Create a loom object
#'
#' @param filename The name of the new loom file
#' @param data The data for \code{/matrix}, should be cells as rows and genes as columns
#' @param gene.attrs A named list of vectors with extra data for genes, each vector must be as long as the number of genes in \code{data}
#' @param cell.attrs A named list of vectors with extra data for cells, each vector must be as long as the number of cells in \code{data}
#' @param chunk.dims A one- or two-length integer vector of chunksizes for \code{/matrix}, defaults to 'auto' to automatically determine chunksize
#'
#' @return A connection to a loom file
#'
#' @importFrom utils packageVersion
#'
#' @seealso \code{\link{loom-class}}
#'
#' @export
#'
create <- function(
  filename,
  data,
  gene.attrs = NULL,
  cell.attrs = NULL,
  layers = NULL,
  chunk.dims = 'auto'
) {
  if (file.exists(filename)) {
    stop(paste('File', file, 'already exists!'))
  }
  if (!is.matrix(x = data)) {
    data <- as.matrix(x = data)
  }
  if (length(x = chunk.dims) > 2 || length(x = chunk.dims) < 1) {
    stop("'chunk.dims' must be a one- or two-length integer vector or 'auto'")
  } else if (length(x = chunk.dims == 1)) {
    if (!grepl(pattern = '^auto$', x = chunk.dims, perl = TRUE)) {
      chunk.dims <- rep.int(x = as.integer(x = chunk.dims), times = 2)
    }
  } else {
    chunk.dims <- as.integer(x = chunk.dims)
  }
  new.loom <- loom$new(filename = filename, mode = 'w-')
  # Create the matrix
  new.loom$create_dataset(
    name = 'matrix',
    robj = data,
    chunk_dims = chunk.dims
  )
  new.loom$matrix <- new.loom[['matrix']]
  new.loom$shape <- new.loom[['matrix']]$dims
  # Groups
  new.loom$create_group(name = 'layers')
  new.loom$create_group(name = 'row_attrs')
  new.loom$create_group(name = 'col_attrs')
  # Check for the existance of gene or cell names
  if (!is.null(x = colnames(x = data))) {
    new.loom$add.row.attribute(attribute = list('gene_names' = colnames(x = data)))
  }
  if (!is.null(x = rownames(x = data))) {
    new.loom$add.col.attribute(attribute = list('cell_names' = rownames(x = data)))
  }
  # Store some constants as HDF5 attributes
  h5attr(x = new.loom, which = 'version') <- new.loom$version
  h5attr(x = new.loom, which = 'chunks') <- paste0(
    '(',
    paste(new.loom[['matrix']]$chunk_dims, collapse = ', '),
    ')'
  )
  # Add layers
  if (!is.null(x = layers)) {
    new.loom$add.layer(layer = layers)
  }
  if (!is.null(x = gene.attrs)) {
    new.loom$add.row.attribute(attribute = gene.attrs)
  }
  if (!is.null(x = cell.attrs)) {
    new.loom$add.col.attribute(attribute = cell.attrs)
  }
  # Set last bit of information
  chunks <- new.loom[['matrix']]$chunk_dims
  chunks <- gsub(pattern = '(', replacement = '', x = chunks, fixed = TRUE)
  chunks <- gsub(pattern = ')', replacement = '', x = chunks, fixed = TRUE)
  chunks <- unlist(x = strsplit(x = chunks, split = ','))
  new.loom$chunksize <- as.integer(x = chunks)
  # Return the connection
  return(new.loom)
}

#' Validate a loom object
#'
#' @param object A loom object
#'
#' @return None, errors out if object is an invalid loom connection
#'
#' @seealso \code{\link{loom-class}}
#'
#' @export
#'
validateLoom <- function(object) {
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

#' Connect to a loom file
#'
#' @param filename The loom file to connect to
#' @param mode How do we connect to it? Pass 'r' for read-only or 'r+' for read/write
#'
#' @return A loom file connection
#'
#' @seealso \code{\link{loom-class}}
#'
#' @export
#'
connect <- function(filename, mode = "r") {
  if (!(mode %in% c('r', 'r+'))) {
    stop("'mode' must be one of 'r' or 'r+'")
  }
  new.loom <- loom$new(filename = filename, mode = mode)
  return(new.loom)
}

CreateLoomFromSeurat <- function(object, filename) {
  object.data=t(object@raw.data[rownames(object@data),object@cell.names])
  object.meta.data=object@meta.data
  row_attrs=list(); col_attrs=list()
  gene.names=colnames(object.data)
  object.meta.data$ident = object@ident
  object.meta.data$CellID = object@cell.names
  for(i in 1:ncol(object.meta.data)) {
    col_attrs[[colnames(object.meta.data)[i]]]=object.meta.data[,i]
  }
  row_attrs[["Gene"]]=gene.names
  create(filename,object.data,gene.attrs = row_attrs, cell.attrs = col_attrs)
}

#' Map a function or a series of functions over a loom file
#'
#' @param X A loom object
#' @param MARGIN The dimmension to map over, pass 1 for cells or 2 for genes
#' @param FUN A function to map to the loom file
#' @param chunk.size Chunk size to use, defaults to \code{loomfile$chunksize[MARGIN]}
#' @param index.use Indices of the dataset to use, defaults to \code{1:loomfile$shape[MARGIN]}
#' @param dataset.use Dataset to use, defauts to 'matrix'
#' @param display.progress Display a progress bar
#' @param expected Shape of expected results. Can pass either 'matrix' or 'vector'; defaults to shape of 'dataset.use'
#' @param ... Extra parameters for FUN
#'
#' @return The results of the map
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
map <- function(
  X,
  MARGIN = 1,
  FUN,
  chunk.size = NULL,
  index.use = NULL,
  dataset.use = 'matrix',
  display.progress = TRUE,
  expected = NULL,
  ...
) {
  if (!inherits(x = X, what = 'loom')) {
    stop("map only works on loom objects")
  }
  if (!inherits(x = FUN, what = 'function')) {
    stop("FUN must be a function")
  }
  # Check for existance of dataset
  if (!any(grepl(pattern = dataset.use, x = list.datasets(object = X)))) {
    stop(paste("Cannot find dataset", dataset.use, "in the loom file"))
  }
  # Figure out if we're returning a vector or matrix
  full.dataset <- grep(
    pattern = dataset.use,
    x = list.datasets(object = X),
    value = TRUE
  )
  results.matrix <- TRUE
  dataset.matrix <- TRUE
  if (grepl(pattern = 'col_attrs', x = full.dataset)) {
    MARGIN <- 1
    results.matrix <- FALSE
    dataset.matrix <- FALSE
  } else if (grepl(pattern = 'row_attrs', x = full.dataset)) {
    MARGIN <- 2
    results.matrix <- FALSE
    dataset.matrix <- FALSE
  }
  if (!is.null(x = expected)) {
    results.matrix <- switch(
      EXPR = expected,
      'vector' = FALSE,
      'matrix' = TRUE,
      stop("'expected' must be one of 'matrix', 'vector', or NULL")
    )
  }
  # Determine the shape of our results
  if (!(MARGIN %in% c(1, 2))) {
    stop("MARGIN must be either 1 (cells) or 2 (genes)")
  }
  if (is.null(x = index.use)) {
    index.use <- c(1, X$shape[MARGIN])
  } else if (length(x = index.use) == 1) {
    index.use <- c(1, index.use)
  }
  index.use[1] <- max(1, index.use[1])
  index.use[2] <- min(index.use[2], X$shape[MARGIN])
  batch <- X$batch.scan(
    chunk.size = chunk.size,
    MARGIN = MARGIN,
    index.use = index.use,
    dataset.use = dataset.use,
    force.reset = TRUE
  )
  # Create our results holding object
  if (results.matrix) {
    switch(
      EXPR = MARGIN,
      '1' = results <- matrix(
        nrow = length(x = index.use[1]:index.use[2]),
        ncol = X$shape[2]
      ),
      '2' = results <- matrix(
        nrow = X$shape[1],
        ncol = length(x = index.use[1]:index.use[2])
      )
    )
  } else {
    results <- vector(length = length(x = index.use[1]:index.use[2]))
  }
  if (display.progress) {
    pb <- txtProgressBar(char = '=', style = 3)
  }
  for (i in 1:length(x = batch)) {
    chunk.indices <- X$batch.next(return.data = FALSE)
    chunk.data <- if (dataset.matrix) {
      switch(
        EXPR = MARGIN,
        '1' = X[[dataset.use]][chunk.indices, ],
        '2' = X[[dataset.use]][, chunk.indices]
      )
    } else {
      X[[dataset.use]][chunk.indices]
    }
    chunk.data <- FUN(chunk.data, ...)
    if (results.matrix) {
      if (MARGIN == 1) {
        results[chunk.indices, ] <- chunk.data
      } else if (MARGIN == 2) {
        results[, chunk.indices] <- chunk.data
      }
    } else {
      results[chunk.indices] <- chunk.data
    }
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = i / length(x = batch))
    }
  }
  return(results)
}

#need to comment
#need to add progress bar
#but otherwise, pretty cool
#for paul to try :
# f <- connect("~/Downloads/10X43_1.loom")
# mean_var = map(f,f_list = c(mean,var),chunksize = 5000)
# nGene <- map(f, f_list = function(x) length(which(x>0)), MARGIN = 2)
map <- function(self, f_list = list(mean, var), MARGIN=1, chunksize=1000, selection) {
  n_func = length(f_list)
  if (n_func == 1) {
    f_list <- list(f_list)
  }
  if (MARGIN == 1) {
    results <- list()
    for (j in 1:n_func) {
      results[[j]] <- numeric(0)
    }
    rows_per_chunk <- chunksize
    ix <- 1
    while (ix <= self@shape[1]) {
      rows_per_chunk <- min(rows_per_chunk, self@shape[1] - ix + 1)
      chunk <- self["matrix"][ix:(ix + rows_per_chunk - 1), ]
      for (j in 1:n_func) {
        new_results <- apply(chunk, 1, FUN = f_list[[j]])
        results[[j]] <- c(results[[j]], new_results)
      }
      ix <- ix + chunksize
    }
  }
  if (MARGIN == 2) {
    results <- list()
    for (j in 1:n_func) {
      results[[j]] <- numeric(0)
    }
    cols_per_chunk <- chunksize
    ix <- 1
    while (ix <= self@shape[2]) {
      cols_per_chunk <- min(cols_per_chunk, self@shape[2] - ix + 1)
      chunk <- self["matrix"][, ix:(ix + cols_per_chunk - 1)]
      for (j in 1:n_func) {
        new_results <- apply(chunk, 2, FUN = f_list[[j]])
        results[[j]] <- c(results[[j]], new_results)
      }
      ix <- ix + chunksize
    }
  }
  if (n_func == 1) {
    results <- results[[1]]
  }
  return(results)
}
