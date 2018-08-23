#' @include internal.R
#' @import hdf5r
#' @import Matrix
#' @importFrom R6 R6Class
#' @importFrom WGCNA transposeBigData
NULL

#' A class for connections loom files
#'
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @aliases loom
#' @format An \code{\link{R6::R6Class}} object
#' @seealso \code{\link{loomR}}, \code{\link{hdf5r::H5File}}
#'
#' @field version Version of loomR object was created under
#' @field shape Shape of \code{/matrix} in genes (columns) by cells (rows)
#' @field matrix The main data matrix, stored as columns (cells) by rows (genes)
#' @field layers Additional data matricies, the same shape as \code{/matrix}
#' @field col.attrs Extra information about cells
#' @field row.attrs Extra information about genes
#'
#' @section Methods:
#' \describe{
#'   \item{\code{add.graph(a, b, w, name, MARGIN, overwrite)}, \code{add.graph.matrix(mat, name, MARGIN, overwrite)}}{
#'     Add a graph to the loom object; can add either in coorindate format (\code{add.graph}) or matrix format (\code{add.graph.matrix}).
#'     Stores graph in coordinate format as \code{[row, col]_graphs/name/a} (row indices),
#'     \code{[row, col]_graphs/name/b} (column indices), and \code{[row, col]_graphs/name/w} (values)
#'     \describe{
#'       \item{\code{a}}{Integer vector of row indices for graph, must be the same lengths as \code{b} and \code{w}}
#'       \item{\code{b}}{Integer vector of column indices for graph, must be the same lengths as \code{a} and \code{w}}
#'       \item{\code{w}}{Numeric vector of values for graph, must be the same lengths as \code{a} and \code{b}}
#'       \item{\code{mat}}{Graph provided as a matrix (sparse or dense) or data.frame}
#'       \item{\code{name}}{Name to store graph, will end up being \code{col_graphs/name} or \code{row_graphs/name}, depending on \code{MARGIN}}
#'       \item{\code{MARGIN}}{Store the graph in \code{row_graphs} (1) or \code{col_graphs} (2), defaults to 2}
#'       \item{\code{overwrite}}{Can overwrite existing graph?}
#'     }
#'   }
#'   \item{\code{add.layer(layer, chunk.size, overwrite)}}{
#'     Add a data layer to this loom file, must be the same dimensions as \code{/matrix}
#'     \describe{
#'       \item{\code{layer}}{A named list of matrices to be added as layers}
#'       \item{\code{chunk.size}}{Number of rows from each layer to stream at once, defaults to 1000}
#'       \item{\code{overwrite}}{If a layer already exists, overwrite with new data, defaults to FALSE}
#'     }
#'   }
#'   \item{\code{add.attribute(attribute, MARGIN, overwrite)}}{
#'     Add extra information to this loom file.
#'     \describe{
#'       \item{\code{attribute}}{A named list where the first dimmension of each element as long as one dimension of \code{/matrix}}
#'       \item{\code{MARGIN}}{Either 1 for genes or 2 for cells}
#'       \item{\code{overwrite}}{Can overwrite existing attributes?}
#'     }
#'   }
#'   \item{\code{add.row.attribute(attribute), add.col.attribute(attribute)}}{
#'     Add row or column attributes
#'   }
#'   \item{\code{get.attribute.df(MARGIN, attribute.names, row.names, col.names)}}{
#'     Get a group of row or column attributes as a data frame, will only return attributes that have one dimension
#'     \describe{
#'       \item{\code{MARGIN}}{Either '1' or '2' to get row- or column-attributes, respectively}
#'       \item{\code{attribute.names}}{A vector of attribute dataset basenames}
#'       \item{\code{row.names}}{Basename of the rownames dataset}
#'       \item{\code{col.names}}{Basename of the colnames dataset}
#'     }
#'   }
#'   \item{\code{get.graph(name, MARGIN)}}{
#'     Get a graph as a sparse matrix
#'     \describe{
#'       \item{\code{name}}{Name of the graph, can be either the basename or full name of the grpah}
#'       \item{\code{MARGIN}}{
#'         Load the graph from \code{row_graphs} (1) or \code{col_graphs} (2), defaults to 2.
#'         Ignored if full path to graph is passed to \code{name}
#'       }
#'     }
#'   }
#'   \item{\code{get.sparse(dataset.use, genes.use, cells.use, chunk.size, gene.names, cell.names, display.progress)}}{
#'     Get a sparse matrix representation of a dataset
#'     \describe{
#'       \item{\code{dataset.use}}{Dataset to be returned as a sparse matrix}
#'       \item{\code{genes.use}}{Vector of indices of genes to return}
#'       \item{\code{cells.use}}{Vector of indices of cells to return}
#'       \item{\code{chunk.size}}{Number of cells to iterate through at a time}
#'       \item{\code{gene.names}}{Dataset to gene names}
#'       \item{\code{cell.names}}{Dataset to cell names}
#'       \item{\code{display.progress}}{Display a progress bar}
#'     }
#'   }
#'   \item{\code{chunk.points(chunk.size = NULL, MARGIN = 2, dataset.use = 'matrix')}, \code{chunk.indices(chunk.size = NULL, MARGIN = 2, dataset.use = 'matrix')}}{
#'     Generate start/end points or indices of chunk size \code{chunk.size} for iteration
#'     \describe{
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by, defaults to \code{self[[dataset.use]]$chunk_dims}}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2), defaults to 2}
#'       \item{\code{dataset.use}}{Name of dataset to use, can be the name, not \code{group/name}, unless the name is present in multiple groups}
#'     }
#'   }
#'   \item{\code{scan(dataset.use = 'matrix', MARGIN = 2, chunk.size = NULL, index.use = NULL, display.progress = TRUE)}}{
#'     Scan through a dataset
#'     \describe{
#'       \item{\code{dataset.use}}{Name of dataset to use, can be the name, not \code{group/name}, unless the name is present in multiple groups}
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by, defaults to \code{self$chunksize}}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2), defaults to 2}
#'       \item{\code{index.use}}{Which specific values of \code{dataset.use} to use, defaults to \code{1:self$shape[MARGIN]} (all values)}
#'       \item{\code{display.progress}}{Display a progress bar}
#'     }
#'   }
#'   \item{\code{batch.scan(chunk.size, MARGIN, index.use, dataset.use, force.reset)}, \code{batch.next(return.data)}}{
#'     Scan a dataset in the loom file from \code{index.use[1]} to \code{index.use[2]}, iterating by \code{chunk.size}.
#'     \describe{
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by, defaults to \code{self$chunksize}}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2), defaults to 2}
#'       \item{\code{index.use}}{Which specific values of \code{dataset.use} to use, defaults to \code{1:self$shape[MARGIN]} (all values)}
#'       \item{\code{dataset.use}}{Name of dataset to use, can be the name, not \code{group/name}, unless the name is present in multiple groups}
#'       \item{\code{force.reset}}{Force a reset of the internal iterator}
#'       \item{\code{return.data}}{Return data for a given chunk, if \code{FALSE}, returns the indices across \code{MARGIN} for said chunk}
#'     }
#'   }
#'   \item{\code{apply(name, FUN, MARGIN, chunk.size, dataset.use, overwrite, display.progress, ...)}}{
#'     Apply a function over a dataset within the loom file, stores the results in the loom file. Will not make multidimensional attributes.
#'     \describe{
#'       \item{\code{name}}{Full name ('group/name') of dataset to store results to}
#'       \item{\code{FUN}}{Function to apply}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2), defaults to 2}
#'       \item{\code{index.use}}{Which specific values of \code{dataset.use} to use, defaults to \code{1:self$shape[MARGIN]} (all values)}
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by, defaults to \code{self$chunksize}}
#'       \item{\code{dataset.use}}{Name of dataset to use}
#'       \item{\code{overwrite}}{Overite \code{name} if already exists}
#'       \item{\code{display.progress}}{Display progress}
#'       \item{\code{max.size}}{Maximum HDF5 chunk size}
#'       \item{\code{chunk.dims}}{HDF5 chunk dimensions}
#'       \item{\code{display.progress}}{Display a progress bar}
#'       \item{\code{...}}{Extra parameters to pass to \code{FUN}}
#'     }
#'   }
#'   \item{\code{map(FUN, MARGIN, chunk.size, index.use, dataset.use, display.progress, expected, ...)}}{
#'     Map a function onto a dataset within the loom file, returns the result into R.
#'     \describe{
#'       \item{\code{FUN}}{}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2), defaults to 2}
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by, defaults to \code{self$chunksize}}
#'       \item{\code{index.use}}{Which specific values of \code{dataset.use} to use, defaults to \code{1:self$shape[MARGIN]} (all values)}
#'       \item{\code{dataset.use}}{Name of dataset to use}
#'       \item{\code{display.progress}}{Display a progress bar}
#'       \item{\code{...}}{Extra parameters to pass to \code{FUN}}
#'     }
#'   }
#'   \item{\code{add.cells(matrix.data, attributes.data = NULL, layers.data = NULL, do.transpose = TRUE,
#'   big.data = FALSE, display.progress = TRUE, ...)}}{
#'     Add cells to a loom file.
#'     \describe{
#'       \item{\code{matrix.data}}{A matrix of m2 cells (columns) with n rows (num genes, \code{self$shape[1]}),
#'       since this will be transposed (with \code{do.transpose = TRUE}), and to avoid transposing large datasets,
#'       a matrix of m2 rows and n columns can also be accepted (with \code{do.transpose = FALSE}).}
#'       \item{\code{attributes.data}}{A list where each entry is named for one of the datasets in \code{self[['col_attrs']]}; each entry is a vector of length m2.}
#'       \item{\code{layers.data}}{A list where each entry is named for one of the datasets in \code{self[['layers']]}; each entry is an n-by-m2 matrix where n is the number of genes in this loom file and m2 is the number of cells being added.}
#'       \item{\code{do.transpose}}{Logical, indicates whether the matrix and layer data should be transposed before added to the loom object.}
#'       \item{\code{big.data}}{Logical, indicates whether matrix and layer data should be transposed by chunks.}
#'       \item{\code{display.progress}}{Display a progress bar}
#'       \item{\code{...}}{Arguments to be passed to \code{WGCNA::transposeBigData}}
#'     }
#'   }
#'   \item{\code{add.loom(other, other.key, self.key, ...)}}{
#'     Add the contents of another loom file to this one.
#'     \describe{
#'       \item{\code{other}}{An object of class \code{loom} or a filename of another loom file}
#'       \item{\code{other.key}}{Row attribute in \code{other} to add by}
#'       \item{\code{self.key}}{Row attribute in this loom file to add by}
#'       \item{\code{...}}{Ignored for now}
#'     }
#'   }
#'   \item{\code{last.modified(x = 'self')}}{
#'     Get the time an object (group, dataset, file) was most recently modified
#'     \describe{
#'       \item{\code{x}}{A character vector of groups or datasets to get the timestamp of, defaults to 'self' or the file as a whole; pass \code{NULL} to get all timestamps present}
#'     }
#'   }
#'   \item{\code{get.changes.since(x)}}{
#'     Get changes since a given time
#'     \describe{
#'       \item{\code{x}}{An object of class \code{POSIXct}}
#'     }
#'   }
#' }
#'
#' @importFrom iterators nextElem
#' @importFrom pbapply pbapply pbsapply
#' @importFrom itertools hasNext ihasNext ichunk
#' @importFrom utils packageVersion setTxtProgressBar
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
    matrix = NULL,
    layers = NULL,
    col.attrs = NULL,
    row.attrs = NULL,
    # Methods
    initialize = function(
      filename = NULL,
      mode = c('a', 'r', 'r+', 'w', 'w-'),
      skip.validate = FALSE,
      ...
    ) {
      # If the file exists, run validation steps
      do.validate <- file.exists(filename) && !(mode %in% c('w', 'w+'))
      private$skipped.validation <- skip.validate
      super$initialize(filename = filename, mode = mode, ...)
      if (do.validate) {
        # Run the validation steps
        if (skip.validate) {
          warning("Skipping validation step, some fields are not populated")
        } else {
          validateLoom(object = self)
        }
        # Store /matrix and the shape of /matrix
        if (skip.validate) {
          if (getOption(x = 'verbose')) {
            warning("Not setting matrix field")
          }
        } else {
          self$matrix <- self[['matrix']]
        }
        self$update.shape()
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
              h5attr(x = self, which = 'version') <- as.character(x = version)
            } else {
              version <- NA_character_
            }
            return(version)
          }
        ))
        # Load layers
        private$load_layers()
        # Load attributes
        private$load_attributes(MARGIN = 1) # Genes (row_attrs)
        private$load_attributes(MARGIN = 2) # Cells (col_attrs
      } else {
        # Assume new HDF5 file
        self$version <- as.character(x = packageVersion(pkg = 'loomR'))
      }
    },
    finalizer = function() {
      self$close_all(close_self = TRUE)
    },
    load.fields = function() {
      private$load_layers()
      private$load_attributes(MARGIN = 1)
      private$load_attributes(MARGIN = 2)
    },
    update.shape = function() {
        self$shape <- rev(x = self[['matrix']]$dims)
    },
    # Addding attributes, layers, and graphs
    add.graph = function(a, b, w, name, MARGIN = 2, overwrite = FALSE) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Coerce datasets to proper types
      a <- as.integer(x = a)
      b <- as.integer(x = b)
      w <- as.double(x = w)
      # Check lengths of each vector
      graph.lengths <- vapply(
        X = list(a, b, w),
        FUN = length,
        FUN.VALUE = integer(length = 1L)
      )
      if (length(x = unique(x = graph.lengths)) > 1) {
        stop("'a', 'b', and 'w' must all be the same length")
      }
      # Check the name, automatically assign a group if provided in `name`
      if (dirname(path = name) != '.') {
        name.group <- dirname(path = name)
        MARGIN <- if (grepl(pattern = 'col', x = name.group)) {
          2
        } else if (grepl(pattern = 'row', x = name.group)) {
          1
        } else {
          cate("Unknown group:", name.group)
          cate("Using", MARGIN, "as default")
          MARGIN
        }
      }
      # Shorten `name` to the just the basename
      name <- basename(path = name)
      # Assign group
      group <- switch(
        EXPR = MARGIN,
        '1' = 'row',
        '2' = 'col',
        stop("'MARGIN' must be 1 or 2")
      )
      group <- paste0(group, '_graphs')
      # Check for existance of previous graph
      if (name %in% names(x = self[[group]])) {
        if (overwrite) {
          self[[group]]$link_delete(name = name)
        } else {
          stop(paste("A graph with the name", name, "exists already!"))
        }
      }
      # Create our datasets
      self[[group]]$create_group(name = name)
      self[[group]][[name]][['a']] <- a
      self[[group]][[name]][['b']] <- b
      self[[group]][[name]][['w']] <- w
      private$timestamp(x = file.path(group, name))
      invisible(x = self)
    },
    add.graph.matrix = function(mat, name, MARGIN = 2, overwrite = FALSE) {
      if (!inherits(x = mat, what = 'dgCMatrix')) {
        mat <- Matrix(
          data = mat,
          nrow = nrow(x = mat),
          ncol = ncol(x = mat),
          sparse = TRUE
        )
      }
      self$add.graph(
        a = mat@i,
        b = PointerToIndex(p = mat@p),
        w = mat@x,
        name = name,
        MARGIN = MARGIN,
        overwrite = overwrite
      )
      invisible(x = self)
    },
    add.layer = function(
      layers,
      chunk.size = 1000,
      overwrite = FALSE,
      display.progress = TRUE
    ) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      self$update.shape()
      # Value checking
      if (!is.list(x = layers) || is.null(x = names(x = layers))) {
        stop("'layers' must be a named list")
      }
      if (is.null(x = self$shape)) {
        stop(private$err_init)
      }
      # Add layers
      for (i in 1L:length(x = layers)) {
        if (!inherits(x = layers[[i]], what = c('Matrix', 'matrix'))) {
          layers[[i]] <- as.matrix(x = layers[[i]])
        }
        do.transpose <- FALSE
        shape.use <- rev(x = self$shape)
        if (any(dim(x = layers[[i]]) != shape.use)) {
          if (all(rev(x = dim(x = layers[[i]])) == shape.use)) {
            do.transpose <- TRUE
          } else {
            stop(paste(
              "All layers must have",
              shape.use[1],
              "rows for cells and",
              shape.use[2],
              "columns for genes"
            ))
          }
        }
        if (do.transpose) {
          layers[[i]] <- t(x = layers[[i]])
        }
        layer.name <- names(x = layers)[i]
        if (layer.name %in% list.datasets(object = self[['layers']])) {
          if (overwrite) {
            self[['layers']]$link_delete(name = layer.name)
          } else {
            stop(paste(
              "A layer with the name",
              layer.name,
              "already!"
            ))
          }
        }
        dtype <- guess_dtype(
          x = layers[[i]][1, 1],
          string_len = getOption(x = "loomR.string_len")
        )
        try(
          expr = {
            if (dtype$get_cset() != getOption(x = 'loomR.ascii')) {
              dtype$set_cset(cset = getOption(x = "loomR.ascii"))
            }
          },
          silent = TRUE
        )
        layer.shape <- dim(x = layers[[i]])
        layer.space <- H5S$new(
          type = 'simple',
          dims = layer.shape,
          maxdims = c(Inf, layer.shape[2])
        )
        chunk.dims <- guess_chunks(
          space_maxdims = layer.space$maxdims,
          dtype_size = dtype$get_size(),
          chunk_size = 4e9
        )
        gc(verbose = FALSE)
        self[['layers']]$create_dataset(
          name = layer.name,
          dtype = dtype,
          space = layer.space,
          chunk_dims = chunk.dims,
          gzip_level = 4
        )
        chunk.size <- chunk.dims[1]
        chunk.points <- chunkPoints(
          data.size = dim(x = layers[[i]])[1],
          chunk.size = chunk.size
        )
        if (display.progress) {
          catn(
            "Adding a layer to",
            names(x = layers)[i],
            "(layer",
            i,
            "of",
            paste0(length(x = layers), ')')
          )
          pb <- newPB()
        }
        for (col in 1:ncol(x = chunk.points)) {
          row.start <- chunk.points[1, col]
          row.end <- chunk.points[2, col]
          self[['layers']][[layer.name]][row.start:row.end, ] <- as.matrix(
            x = layers[[i]][row.start:row.end, ]
          )
          if (display.progress) {
            setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
          }
        }
        if (display.progress) {
          close(con = pb)
        }
      }
      self$flush()
      gc(verbose = FALSE)
      private$timestamp(x = file.path('layers', names(x = layers)))
      private$load_layers()
      invisible(x = self)
    },
    add.attribute = function(attribute, MARGIN, overwrite = FALSE) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      self$update.shape()
      # Value checking
      if (is.data.frame(x = attribute)) {
        attribute <- as.list(x = attribute)
      }
      is.actual.list <- is.list(x = attribute)
      if (!is.actual.list || is.null(x = names(x = attribute))) {
        stop("Attributes must be provided as a named list")
      }
      if (!MARGIN %in% c(1, 2)) {
        stop("'MARGIN' must be 1 or 2")
      }
      # length.use <- rev(x = self[['matrix']]$dims)[MARGIN]
      length.use <- self$shape[MARGIN]
      dim.msg <- paste(
        "At least one dimmension for each",
        switch(EXPR = MARGIN, '1' = 'gene', '2' = 'cell'),
        "attribute must be",
        length.use
      )
      for (i in 1:length(x = attribute)) {
        if (is.matrix(x = attribute[[i]]) || is.data.frame(x = attribute[[i]])) {
          margin.use <- which(x = dim(x = attribute[[i]]) == length.use)
          if (!length(x = margin.use)) {
            stop(dim.msg)
          }
          margin.use <- margin.use[1]
          attribute[[i]] <- switch(
            EXPR = margin.use,
            '1' = t(x = as.matrix(x = attribute[[i]])),
            '2' = as.matrix(x = attribute[[i]]),
            stop("All attributes must be one- or two-dimmensional")
          )
        } else {
          if (length(x = attribute[[i]]) != length.use) {
            stop(dim.msg)
          }
          attribute[[i]] <- as.vector(x = attribute[[i]])
        }
      }
      if (length(x = which(x = names(x = attribute) != '')) != length(x = attribute)) {
        stop("Not all attributes had names provided")
      }
      grp.name <- c('row_attrs', 'col_attrs')[MARGIN]
      grp <- self[[grp.name]]
      if (!overwrite) {
        names.fail <- which(x = names(x = attribute) %in% names(x = grp))
        if (length(x = names.fail) != 0) {
          stop(paste0(
            "The following names are already used for ",
            switch(EXPR = MARGIN, '1' = 'row', '2' = 'column'),
            " attributes: '",
            paste(names(x = attribute)[names.fail], collapse = ''),
            "'"
          ))
        }
      }
      # Add the attributes as datasets for our MARGIN's group
      for (i in 1:length(x = attribute)) {
        try(expr = grp$link_delete(name = names(x = attribute)[i]), silent = TRUE)
        grp[[names(x = attribute)[i]]] <- attribute[[i]]
      }
      self$flush()
      gc(verbose = FALSE)
      # Load the attributes for this margin
      private$load_attributes(MARGIN = MARGIN)
      private$timestamp(x = file.path(grp.name, names(x = attribute)))
      invisible(x = self)
    },
    add.row.attribute = function(attribute, overwrite = FALSE) {
      self$add.attribute(attribute = attribute, MARGIN = 1, overwrite = overwrite)
      invisible(x = self)
    },
    add.col.attribute = function(attribute, overwrite = FALSE) {
      self$add.attribute(attribute = attribute, MARGIN = 2, overwrite = overwrite)
      invisible(x = self)
    },
    add.meta.data = function(meta.data, overwrite = FALSE) {
      self$add.col.attribute(attribute = meta.data, overwrite = overwrite)
      invisible(x = self)
    },
    # Get data from the loom file
    get.attribute.df = function(
      MARGIN = 2,
      attribute.names = NULL,
      row.names = "gene_names",
      col.names = "cell_names"
    ) {
      # takes multiple row_attrs or col_attrs and combines them into a data.frame,
      # removing rows or columns that are entirely NAs.
      #
      # attribute.layer either "row" for row_attrs or "col" col_attrs
      # attribute.names name of rows or columns to combine into matrix
      # row.names either a character vector or name of an element in row_attrs
      # col.names either a character vector or name of an element in col_attrs
      # if (!attribute.layer %in% c("row", "col")) {
      #   stop("Invalid attribute.layer. Please select either 'row' or 'col'.")
      # }
      attribute.layer <- switch(
        EXPR = MARGIN,
        '1' = 'row',
        '2' = 'col',
        stop("'MARGIN' must be either '1' for row attributes or '2' for column attributes")
      )
      attribute.layer <- paste0(attribute.layer, "_attrs")
      # by default return all attributes
      if (is.null(attribute.names)) {
        attribute.names <- self[[attribute.layer]]$names
      }
      # check that attribute names are present
      if (!all(attribute.names %in% self[[attribute.layer]]$names)) {
        invalid.names <- attribute.names[which(!attribute.names %in% self[[attribute.layer]]$names)]
        stop(paste0("Invalid attribute.names: ", paste0(invalid.names, collapse = ", ")))
      }
      attr.paths <- paste0(attribute.layer, "/", attribute.names)
      # keep only the one-dimensional attributes
      dims <- sapply(attr.paths, function(x) length(self[[x]]$dims))
      data.lst <- lapply(attr.paths[dims == 1], function(x) data.frame(self[[x]][], stringsAsFactors = FALSE))
      combined.df <- Reduce(cbind, data.lst)
      colnames(combined.df) <- attribute.names[dims == 1]
      if (attribute.layer == "col_attrs") {
        rownames(combined.df) <- self[[paste0(attribute.layer, "/", col.names)]][]
      } else {
        rownames(combined.df) <- self[[paste0(attribute.layer, "/", row.names)]][]
      }
      # check if any row is all NAs
      rows.to.remove <- unname(obj = which(x = apply(
        X = combined.df,
        MARGIN = 1,
        FUN = function(x) {
          return(all(is.na(x = x)))
        }
      )))
      if (length(x = rows.to.remove) > 1) {
        combined.df <- combined.df[-rows.to.remove, ]
      }
      return(combined.df)
    },
    get.graph  = function(name, MARGIN = 2) {
      # Check the name, automatically assign a group if provided in `name`
      if (dirname(path = name) != '.') {
        name.group <- dirname(path = name)
        MARGIN <- if (grepl(pattern = 'col', x = name.group)) {
          2
        } else if (grepl(pattern = 'row', x = name.group)) {
          1
        } else {
          cate("Unknown group:", name.group)
          cate("Using", MARGIN, "as default")
          MARGIN
        }
      }
      # Shorten `name` to the just the basename
      name <- basename(path = name)
      # Assign group
      group <- switch(
        EXPR = MARGIN,
        '1' = 'row',
        '2' = 'col',
        stop("'MARGIN' must be 1 or 2")
      )
      group <- paste0(group, '_graphs')
      # Validate graph
      if (!name %in% list.groups(object = self[[group]], full.names = FALSE, recursive = FALSE)) {
        stop(paste("Cannot find a graph named", name, "in", group))
      }
      graph.datasets <- list.datasets(object = self[[group]][[name]], full.names = FALSE)
      if (length(x = graph.datasets) != 3) {
        stop(paste0(
          "Malformed graph: wrong number of components for graph (",
          length(x = graph.datasets),
          ")"
        ))
      }
      if (!all(graph.datasets %in% c('a', 'b', 'w'))) {
        stop("Malformed graph: unknown dataset names")
      }
      graph <- sparseMatrix(
        i = self[[group]][[name]][['a']][] + 1,
        j = self[[group]][[name]][['b']][],
        x = self[[group]][[name]][['w']][]
      )
      return(graph)
    },
    get.sparse = function(
      dataset.use,
      genes.use = NULL,
      cells.use = NULL,
      chunk.size = NULL,
      gene.names = NULL,
      cell.names = NULL,
      display.progress = TRUE
    ) {
      self$update.shape()
      genes.use <- if (is.null(x = genes.use)) {
        1L:self$shape[1]
      } else {
        sort(x = genes.use)
      }
      cells.use <- if (is.null(x = cells.use)) {
        1L:self$shape[2]
      } else {
        sort(x = cells.use)
      }
      datasets.avail <- c(
        '/matrix',
        list.datasets(object = self, path = '/layers', full.names = TRUE)
      )
      dataset.use <- if(self$exists(name= dataset.use)) {
        dataset.use
      } else {
        grep(pattern = dataset.use, x = datasets.avail, value = TRUE)
      }
      if (length(x = dataset.use) < 1) {
        stop("Cannot find requested dataset, 'get.sparse' only works for '/matrix' and datasets in '/layers'")
      } else if (length(x = dataset.use) > 1) {
        stop(private$err_ambiguous)
      }
      chunk.size <- ifelse(
        test = is.null(x = chunk.size),
        yes = self[[dataset.use]]$chunk_dims[1],
        no = chunk.size
      )
      chunk.points <- self$chunk.points(
        chunk.size = chunk.size,
        MARGIN = 2,
        dataset.use = dataset.use
      )
      myapply <- ifelse(test = display.progress, yes = pbapply, no = apply)
      i <- vector(mode = 'integer')
      p <- vector(mode = 'integer', length = 1L)
      dat <- myapply(
        X = chunk.points,
        MARGIN = 2,
        FUN = function(points) {
          chunk.indices <- points[1]:points[2]
          indices.use <- chunk.indices[chunk.indices %in% cells.use]
          indices.use <- indices.use - chunk.indices[1] + 1
          if (length(x = indices.use) < 1) {
            return(NULL)
          }
          chunk.data <- t(x = self[[dataset.use]][chunk.indices, ][indices.use, genes.use])
          i.chunk <- apply(
            X = chunk.data,
            MARGIN = 2,
            FUN = function(row) {
              return(which(x = row != 0))
            }
          )
          p.chunk <- vapply(X = i.chunk, FUN = length, FUN.VALUE = integer(length = 1L))
          names(x = p.chunk) <- names(x = i.chunk) <- NULL
          i <<- c(i, unlist(x = i.chunk))
          p <<- c(p, max(p) + cumsum(x = p.chunk))
          gc(verbose = FALSE)
          return(chunk.data[chunk.data != 0])
        }
      )
      dat <- as.vector(x = unlist(x = dat))
      gc(verbose = FALSE)
      data.return <- sparseMatrix(
        i = i,
        p = p,
        x = dat,
        dims = c(length(x = genes.use), length(x = cells.use))
      )
      if (!is.null(x = gene.names)) {
        if (!self$exists(name = gene.names)) {
          gene.names <- grep(
            pattern = gene.names,
            x = list.datasets(object = self, path = 'row_attrs', full.names = TRUE),
            value = TRUE
          )
        }
        if (length(x = gene.names) == 1) {
          rownames(x = data.return) <- self[[gene.names]][genes.use]
        } else {
          warning("Cannot find dataset ", gene.names)
        }
      }
      if (!is.null(x = cell.names)) {
        if (!self$exists(name = cell.names)) {
          cell.names <- grep(
            pattern = cell.names,
            x = list.datasets(object = self, path = 'col_attrs', full.names = TRUE),
            value = TRUE
          )
        }
        if (length(x = cell.names) == 1) {
          colnames(x = data.return) <- self[[cell.names]][cells.use]
        } else {
          warning("Cannot find dataset ", gene.names)
        }
      }
      return(data.return)
    },
    # Chunking functions
    chunk.points = function(chunk.size = NULL, MARGIN = 2, dataset.use = 'matrix') {
      if (!self$exists(name = dataset.use)) {
        dataset.use <- grep(
          pattern = dataset.use,
          x = list.datasets(
            object = self,
            path = '/',
            full.names = TRUE,
            recursive = TRUE
          ),
          value = TRUE
        )
      }
      if (length(x = dataset.use) < 1) {
        stop("unfound")
      } else if (length(x = dataset.use) > 1) {
        stop(private$err_ambiguous)
      }
      if (grepl(pattern = 'row_attrs', x = dataset.use)) {
        MARGIN <- 1
      } else if (grepl(pattern = 'col_attrs', x = dataset.use)) {
        MARGIN <- 2
      }
      if (!MARGIN %in% c(1, 2)) {
        stop("MARGIN")
      }
      if (is.null(x = chunk.size)) {
        chunk.size <- self[[dataset.use]]$chunk_dims[-MARGIN]
      }
      chunk.points <- chunkPoints(
        data.size = self[[dataset.use]]$dims[-MARGIN],
        chunk.size = chunk.size
      )
      return(chunk.points)
    },
    chunk.indices = function(chunk.size = NULL, MARGIN = 2, dataset.use = 'matrix') {
      return(apply(
        X = self$chunk.points(
          chunk.size = chunk.size,
          MARGIN = MARGIN,
          dataset.use = dataset.use
        ),
        MARGIN = 2,
        FUN = function(x) {
          return(x[1]:x[2])
        }
      ))
    },
    scan = function(
      check = FALSE,
      reset = FALSE
    ) {
      .NotYetImplemented()
      scan.check <- vapply(
        X = list(private$counter, private$chunk.indices, private$dataset.use),
        FUN = is.null,
        FUN.VALUE = logical(length = 1L)
      )
      if (any(scan.check)) {
        stop('init')
      } else if (check) {
        return(private$counter < length(x = private$chunk.indices))
      } else if (reset) {
        private$counter <- 0
        invisible(x = NULL)
      } else if (private$counter >= length(x = private$chunk.indices)) {
        return(NULL)
      } else {
        private$counter <- private$counter + 1
        data.return <- if (private$MARGIN == 1) {
          self[[private$dataset.use]][, private$chunk.indices[[private$counter]]]
        } else {
          self[[private$dataset.use]][private$chunk.indices[[private$counter]], ]
        }
        return(data.return)
      }
    },
    batch.scan = function(
      chunk.size = NULL,
      MARGIN = 2,
      index.use = NULL,
      dataset.use = 'matrix',
      force.reset = FALSE
    ) {
      .Deprecated(
        new = 'loom$map',
        package = 'loomR',
        msg = private$deprecated,
        old = 'loom$batch.scan'
      )
      if (is.null(x = private$it) || !grepl(pattern = dataset.use, x = private$iter.dataset) || force.reset) {
        # Check the existence of the dataset
        private$iter.dataset <- grep(
          pattern = dataset.use,
          x = list.datasets(object = self),
          value = TRUE
        )
        if (length(x = private$iter.dataset) != 1) {
          private$reset_batch()
          stop(paste0("Cannot find dataset '", dataset.use, "' in the loom file"))
        }
        if (grepl(pattern = 'row_attrs', x = private$iter.dataset)) {
          MARGIN <- 1
        }
        else if (grepl(pattern = 'col_attrs', x = private$iter.dataset)) {
          MARGIN <- 2
        }
        # Check the margin
        if (!(MARGIN %in% c(1, 2))) {
          private$reset_batch()
          stop("MARGIN must be 1 (genes) or 2 (cells)")
        } else {
          private$iter.margin <- MARGIN
        }
        if (is.null(x = chunk.size)) {
          chunk.size <- rev(x = self[[private$iter.dataset]]$chunk_dims)[private$iter.margin]
          # chunk.size <- rev(x = self$chunksize)[private$iter.margin]
        }
        private$iter.chunksize <- chunk.size
        # Set the indices to use
        index.use <- private$iter_range(index.use = index.use)
        # Setup our iterator
        private$it <- ihasNext(iterable = ichunk(
          iterable = index.use[1]:index.use[2],
          chunkSize = chunk.size
        ))
        private$iter.index <- index.use
      }
      # Return the times we iterate, this isn't important, we only need the length of this vector
      return(1:ceiling((private$iter.index[2] - private$iter.index[1]) / private$iter.chunksize))
    },
    batch.next = function(return.data = TRUE) {
      .Deprecated(
        new = 'loom$map',
        package = 'loomR',
        msg = private$deprecated,
        old = 'loom$batch.next'
      )
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
              '1' = self[[private$iter.dataset]][, chunk.indices], # Genes
              '2' = self[[private$iter.dataset]][chunk.indices, ] # Cells
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
      MARGIN = 2,
      index.use = NULL,
      chunk.size = NULL,
      dataset.use = 'matrix',
      overwrite = FALSE,
      display.progress = TRUE,
      max.size = '400mb',
      chunk.dims = NULL,
      dtype.use = NULL,
      ...
    ) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      if (!inherits(x = FUN, what = 'function')) {
        stop("FUN must be a function")
      }
      # Check that we're storing our results properly
      results.basename <- basename(path = name) # Basename of our results
      results.dirname <- gsub(pattern = '/', replacement = '', x = dirname(path = name)) # Groupname of our results
      dirnames <- c('row_attrs', 'col_attrs', 'layers') # Allowed group names
      if (name %in% list.datasets(object = self)) {
        if (overwrite) {
          self$link_delete(name = name)
        } else {
          stop(paste("A dataset with the name", name, "already exists!"))
        }
      }
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
      if (name.check %in% c(1, 2) && name.check != MARGIN) {
        private$reset_batch()
        stop(paste(
          "Iteration must be over",
          c('genes', 'cells')[name.check],
          paste0("(MARGIN = ", name.check, ")"),
          "to store results in",
          paste0("'", dirnames[name.check], "'")
        ))
      }
      # Check how we store our results
      dataset.matrix <- ('matrix' %in% dataset.use || grepl(pattern = 'layers', x = dataset.use))
      results.matrix <- name.check == 3
      # Ensure index.use is integers within the bounds of [1, self$shape[MARGIN]]
      if (is.null(x = index.use)) {
        index.use <- 1:self$shape[MARGIN]
      } else {
        # Filter index.use to values between 1 and self$shape[MARGIN]
        index.use <- sort(x = index.use)
        index.use <- as.integer(x = index.use)
        index.use <- index.use[index.use >= 1 & index.use <= self$shape[MARGIN]]
        index.use <- as.vector(x = na.omit(object = index.use))
        # If we still have values, figure out NAs, otherwise set index.use to NULL
        if (length(x = index.use) == 0) {
          warning("No values passed to 'index.use' fall within the data, using all values")
          index.use <- 1:self$shape[MARGIN]
        }
      }
      # Trial to get class of new dataset
      # Do a trial run to figure out the class of dataset
      na.use <- NA
      trial.use <- if (is.null(x = index.use)) {
        sample(x = 1:self$shape[MARGIN], size = 3, replace = FALSE)
      } else {
        sample(x = index.use, size = 3, replace = FALSE)
      }
      trial.use <- sort(x = trial.use)
      trial <- if (grepl(pattern = 'layers', x = dataset.use) || dataset.use == 'matrix') {
        switch(
          EXPR = MARGIN,
          '1' = self[[dataset.use]][, trial.use],
          '2' = self[[dataset.use]][trial.use, ]
        )
      } else {
        self[[dataset.use]][trial.use]
      }
      trial <- FUN(trial, ...)
      if (is.list(x = trial)) {
        trial <- unlist(x = trial)
      }
      trial <- as.vector(x = trial)
      class(x = na.use) <- class(x = trial)
      # Get a connection to the group we're iterating over
      group <- self[[results.dirname]]
      # Make results dataset
      if (is.null(x = dtype.use)) {
        dtype.use <- guess_dtype(
          x = na.use,
          string_len = getOption(x = "loomR.string_len")
        )
      }
      switch(
        EXPR = results.dirname,
        'layers' = {
          dims.use <- self[['matrix']]$dims
          max.dims <- c(Inf, dims.use[2])
          max.chunk.dims <- self[['matrix']]$dims
        },
        'row_attrs' = {
          max.dims <- dims.use <- self[['matrix']]$dims[2]
          max.chunk.dims <- self[['matrix']]$dims[2]
        },
        'col_attrs' = {
          dims.use <- self[['matrix']]$dims[1]
          max.dims <- Inf
          max.chunk.dims <- self[['matrix']]$dims[1]
        }
      )
      results.space <- H5S$new(
        type = 'simple',
        dims = dims.use,
        maxdims = max.dims
      )
      if (is.null(x = chunk.dims)) {
        mem.size <- charToBytes(x = max.size)
        if (mem.size > 4e9) {
          cate("HDF5 limits internal chunk sizes to 4 GB, setting to 4 GB")
          mem.size <- 4e9
        }
        chunk.dims <- guess_chunks(
          space_maxdims = results.space$maxdims,
          dtype_size = dtype.use$get_size(),
          chunk_size = mem.size
        )
        gc(verbose = FALSE)
      }
      chunk.dims <- pmin(chunk.dims, max.chunk.dims)
      group$create_dataset(
        name = results.basename,
        dtype = dtype.use,
        space = results.space,
        chunk_dims = as.integer(x = chunk.dims),
        gzip_level = 4
      )
      chunk.size = chunk.dims[1]
      chunk.points <- self$chunk.points(
        chunk.size = chunk.size,
        MARGIN = MARGIN,
        dataset.use = dataset.use
      )
      # Start the iteration
      if (display.progress) {
        catn("Writing results to", name)
        myapply <- pbapply
      } else {
        myapply <- apply
      }
      unique(x = myapply(
        X = chunk.points,
        MARGIN = 2,
        FUN = function(points, ...) {
          # Get the indices we're iterating over
          start <- points[1]
          end <- points[2]
          chunk.indices <- start:end
          indices.use <- chunk.indices[chunk.indices %in% index.use]
          indices.use <- indices.use - chunk.indices[1] + 1
          if (length(x = indices.use) >= 1) {
            # Get the data and apply FUN
            chunk.data <- if (dataset.matrix) {
              switch(
                EXPR = MARGIN,
                '1' = {
                  # Chunk genes
                  x <- self[[dataset.use]][, chunk.indices]
                  x[, indices.use]
                },
                '2' = {
                  # Chunk cells
                  x <- self[[dataset.use]][chunk.indices, ]
                  x[indices.use, ]
                }
              )
            } else {
              x <- self[[private$iter.datset]][chunk.indices]
              x[indices.use]
            }
            chunk.data <- FUN(chunk.data, ...)
            if (results.matrix) {
              # If we're writign to a matrix
              # Figure out which way we're writing the data
              switch(
                EXPR = MARGIN,
                '1' = {
                  chunk.full <- matrix(
                    nrow = nrow(x = chunk.data),
                    ncol = length(x = chunk.indices)
                  )
                  chunk.full[, indices.use] <- chunk.data
                  group[[results.basename]][, chunk.indices] <- chunk.full
                  # group[[results.use]][, chunk.indices] <- chunk.full
                },
                '2' = {
                  chunk.full <- matrix(
                    nrow = length(x = chunk.indices),
                    ncol = ncol(x = chunk.data)
                  )
                  chunk.full[indices.use, ] <- chunk.data
                  group[[results.basename]][chunk.indices, ] <- chunk.full
                  # group[[results.use]][chunk.indices, ] <- chunk.full
                }
              )
            } else {
              # Just write to the vector
              chunk.full <- vector(length = length(x = chunk.indices))
              chunk.full[indices.use] <- chunk.data
              group[[results.basename]][chunk.indices] <- chunk.full
              # group[[results.use]][chunk.indices] <- chunk.full
            }
          }
          gc(verbose = FALSE)
          return(NULL)
        },
        ...
      ))
      # Update timestamp
      private$timestamp(x = file.path(results.dirname, results.basename))
      # Clean up and allow chaining
      private$reset_batch()
      # Load layers and attributes
      private$load_layers()
      private$load_attributes(MARGIN = 1) # Genes (row_attrs)
      private$load_attributes(MARGIN = 2) # Cells (col_attrs)
      invisible(x = self)
    },
    map = function(
      FUN,
      MARGIN = 2,
      chunk.size = NULL,
      index.use = NULL,
      dataset.use = 'matrix',
      display.progress = TRUE,
      ...
    ) {
      if (!inherits(x = FUN, what = 'function')) {
        stop("FUN must be a function")
      }
      # Check how we store our results
      # And what the shape of our dataset is
      dataset.matrix <- any(vapply(
        X = c('matrix', 'layers'),
        FUN = grepl,
        FUN.VALUE = logical(length = 1L),
        x = dataset.use
      ))
      chunk.points <- self$chunk.points(
        chunk.size = chunk.size,
        MARGIN = MARGIN,
        dataset.use = dataset.use
      )
      # Ensure that index.use is integers within the bounds of [1, self$shape[MARGIN]]
      if (is.null(x = index.use)) {
        index.use <- 1L:self$shape[MARGIN]
      } else {
        # Filter index.use to values between 1 and self$shape[MARGIN]
        index.use <- sort(x = index.use)
        index.use <- as.integer(x = index.use)
        index.use <- index.use[index.use >= 1 & index.use <= self$shape[MARGIN]]
        index.use <- as.vector(x = na.omit(object = index.use))
        # If we still have values, figure out NAs, otherwise set index.use to NULL
        if (length(x = index.use) == 0) {
          warning("No values passed to 'index.use' fall within the data, using all values")
          index.use <- 1:self$shape[MARGIN]
        }
      }
      # Create our results holding object
      # results <- vector(mode = "list", length = length(x = batch))
      results <- vector(mode = 'list', length = ncol(x = chunk.points))
      myapply <- ifelse(test = display.progress, yes = pbsapply, no = sapply)
      unique(x = myapply(
        X = 1:ncol(x = chunk.points),
        FUN = function(i, ...) {
          # Get the indices we're iterating over
          start <- chunk.points[1, i]
          end <- chunk.points[2, i]
          chunk.indices <- start:end
          indices.use <- chunk.indices[chunk.indices %in% index.use]
          indices.use <- indices.use - chunk.indices[1] + 1
          if (length(x = indices.use) >= 1) {
            # Get the data and apply FUN
            chunk.data <- if (dataset.matrix) {
              switch(
                EXPR = MARGIN,
                '1' = {
                  # Chunk genes
                  x <- self[[dataset.use]][, chunk.indices]
                  x[, indices.use]
                },
                '2' = {
                  # Chunk cells
                  x <- self[[dataset.use]][chunk.indices, ]
                  x[indices.use, ]
                }
              )
            } else {
              x <- self[[dataset.use]][chunk.indices]
              x[indices.use]
            }
            results[[i]] <<- FUN(chunk.data, ...)
          }
          return(NULL)
        },
        ...
      ))
      # Bring result list into matrix or vector format
      if (inherits(x = results[[1]], what = c('matrix', 'Matrix', 'data.frame'))) {
        reduce.func <- switch(EXPR = MARGIN, '1' = cbind, '2' = rbind)
        results <- Reduce(f = reduce.func, x = results)
      } else {
        results <- unlist(x = results, use.names = FALSE)
      }
      private$reset_batch()
      return(results)
    },
    # Functions that modify `/matrix'
    add.cells = function(
      matrix.data,
      attributes.data = NULL,
      layers.data = NULL,
      # I added the do.transpose argument to remove the need to transpose matrix and layer data
      # in case a lot of cells are added
      do.transpose = TRUE,
      big.data = FALSE,
      display.progress = TRUE,
      ...
    ) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Check inputs
      n <- self[['matrix']]$dims[2]
      if (display.progress) {
        cate("Checking inputs...")
      }
      matrix.data <- check.matrix_data(x = matrix.data, n = n, do.transpose = do.transpose)
      layers.data <- check.layers(
        x = layers.data,
        n = n,
        layers.names = names(x = self[['layers']]),
        do.transpose = do.transpose
      )
      matrices <- vapply(attributes.data, is.matrix, FUN.VALUE = logical(1L))
      attributes.data <- check.col_attrs(
        x = attributes.data,
        obj = self,
        matrices = matrices
      )
      # Get the number of cells we're adding
      num.cells <- c(
        nCells.matrix_data(x = matrix.data, do.transpose = do.transpose),
        nCells.layers(x = layers.data, do.transpose = do.transpose),
        nCells.col_attrs(x = attributes.data)
      )
      num.cells <- max(num.cells)
      # Flesh out the input data to add
      if (display.progress) {
        cate(paste("Adding", num.cells, "to this loom file"))
      }
      matrix.data <- addCells.matrix_data(x = matrix.data, n = n, m2 = num.cells,
                                          do.transpose = do.transpose)
      layers.data <- addCells.layers(x = layers.data, n = n, m2 = num.cells,
                                     do.transpose = do.transpose)
      attributes.data <- addCells.col_attrs(x = attributes.data, m2 = num.cells, matrices = matrices)
      # Add the input to the loom file
      dims.fill <- self[['matrix']]$dims[1]
      dims.fill <- (dims.fill + 1L):(dims.fill + num.cells)
      # Matrix data
      if (display.progress) {
        cate("Adding data to /matrix")
      }
      if (do.transpose) {
        if (big.data) {
          matrix.data <- WGCNA::transposeBigData(matrix.data, ...)
        } else {
          matrix.data <- t(matrix.data)
        }
      }
      self[['matrix']][dims.fill, ] <- matrix.data
      # Layer data
      if (display.progress) {
        cate("Adding data to /layers")
        pb <- newPB()
        counter <- 0
      }
      layers.names <- names(x = self[['layers']])
      if (do.transpose) {
        if (big.data) {
          layers.data <- lapply(layers.data, WGCNA::transposeBigData, ...)
        } else {
          layers.data <- lapply(layers.data, t)
        }
      }
      for (i in layers.names) {
        self[['layers']][[i]][dims.fill, ] <- layers.data[[i]]
        if (display.progress) {
          counter <- counter + 1
          setTxtProgressBar(pb = pb, value = counter / length(x = layers.names))
        }
      }
      # Column attributes
      if (display.progress) {
        cate("Adding data to /col_attrs")
        pb <- newPB()
        counter <- 0
      }
      attrs.names <- names(x = self[['col_attrs']])
      for (i in attrs.names) {
        if (matrices[i]) {
          self[["col_attrs"]][[i]][,dims.fill] <- attributes.data[[i]]
        } else {
          self[['col_attrs']][[i]][dims.fill] <- attributes.data[[i]]
        }
        if (display.progress) {
          counter <- counter + 1
          setTxtProgressBar(pb = pb, value = counter / length(x = attrs.names))
        }
      }
      # Update self$shape
      self$shape[] <- self$matrix$dims[2:1]
      # # Load layers and attributes
      # private$load_layers()
      # private$load_attributes(MARGIN = 1)
      # private$load_attributes(MARGIN = 2)
      # self$shape <- self[['matrix']]$dims
    },
    add.loom = function(
      other,
      other.key = NULL,
      self.key = NULL,
      ...
    ) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Connect to the other loom file
      if (inherits(x = other, what = 'loom')) {
        ofile <- other
      } else if (is.character(x = other)) {
        ofile <- connect(filename = other)
      } else {
        stop("'other' must be either a loom object or a path to a loom file")
      }
      # If we have row keys to use
      if (!is.null(x = other.key) && !is.null(x = self.key)) {
        other.key <- basename(path = other.key)
        self.key <- basename(path = self.key)
        tryCatch(
          expr = other.key <- other[['row_attrs']][[other.key]][],
          error = function(e) {
            if (is.character(x = other)) {
              ofile$close_all()
            }
            stop("Failed to find the gene names dataset in the other loom file")
          }
        )
        tryCatch(
          expr = self.key <- self[['row_attrs']][[self.key]][],
          error = function(e) {
            if (is.character(x = other)) {
              ofile$close_all()
            }
            stop("Failed to find the gene names dataset in this loom file")
          }
        )
        # Match rows
        rows.use <- match(x = other.key, table = self.key, nomatch = 0)
        rows.use <- rows.use[rows.use > 0]
      } else {
        cate("Adding the loom file as-is, assuming in the correct order")
        Sys.sleep(time = 3)
        rows.use <- 1:other[['matrix']]$dims[2]
      }
      if (max(rows.use) > self[['matrix']]$dims[2]) {
        stop("More genes in the other loom file than present in this one")
      } else if (max(rows.use) < self[['matrix']]$dims[2]) {
        cate("...")
      }
      # Clean up
      if (is.character(x = other)) {
        ofile$close_all()
      }
    },
    # Datetime functions
    last.modified = function(x = 'self') {
      if (is.null(x = x)) {
        x <- c(
          'self',
          list.datasets(
            object = self,
            path = '/',
            full.names = FALSE,
            recursive = TRUE
          ),
          list.groups(
            object = self,
            path = '/',
            full.names = FALSE,
            recursive = TRUE
          )
        )
      }
      if (!all(x == 'self' | sapply(X = x, FUN = self$exists))) {
        stop("Cannot find all groups/datasets passed")
      }
      timestamps <- lapply(
        X = x,
        FUN = function(i) {
          time <- if (i == 'self') {
            tryCatch(
              expr = h5attr(x = self, which = 'last_modified'),
              error = function(e) {
                if (self$mode != 'r') {
                  private$timestamp(x = NULL)
                  return(h5attr(x = self, which = 'last_modified'))
                } else {
                  return(NULL)
                }
              }
            )
          } else {
            tryCatch(
              expr = h5attr(x = self[[i]], which = 'last_modified'),
              error = function(e) {
                if (self$mode != 'r') {
                  private$timestamp(x = i)
                  return(h5attr(x = self, which = 'last_modified'))
                } else {
                  return(NULL)
                }
              }
            )
          }
          if (!is.null(x = time)) {
            time <- unlist(x = strsplit(x = time, split = '.', fixed = TRUE))[1]
            if (substr(x = time, start = nchar(x = time), stop = nchar(x = time)) != 'Z') {
              time <- paste0(time, 'Z')
            }
            time <- as.POSIXct(x = strptime(
              x = time,
              format = '%Y%m%dT%H%M%SZ',
              tz = 'UTC'
            ))
            time <- format(x = time, tz = Sys.timezone(), usetz = TRUE)
          }
          return(time)
        }
      )
      names(x = timestamps) <- x
      return(Filter(f = Negate(f = is.null), x = timestamps))
    },
    get.changes.since = function(ts) {
      if (!inherits(x = ts, what = 'POSIXct')) {
        tryCatch(
          expr = ts <- as.POSIXct(x = ts),
          error = function(e) {
            stop("'timestamp', must be of class POSIXct")
          }
        )
      }
      all.timestamps <- as.POSIXct(x = self$last.modified(x = NULL))
      return(all.timestamps[all.timestamps >= ts])
    }
  ),
  # Private fields and methods
  # @field err_init An error message for if this object hasn't been created with loomR::create or loomR::connect
  # @field err_mode An error message for if this object is in read-only mode
  # @field it Iterator object for batch.scan and batch.next
  # @field iter.dataset Dataset for iterating on
  # @field iter.margin Margin to iterate over
  # @field iter.index Index values for iteration
  # @field skipped.validation Was validation skipped?
  # @section Methods:
  # \describe{
  #   \item{\code{load_attributes(MARGIN)}}{Load attributes of a given MARGIN into \code{self$col.attrs} or \code{self$row.attrs}}
  #   \item{\code{load_layers()}}{Load layers into \code{self$layers}}
  #   \item{\code{reset_batch()}}{Reset the batch iterator fields}
  #   \item{\code{iter_range(index.use)}}{Get the range of indices for a batch iteration}
  #   \item{\code{timestamp(x, silent = FALSE)}}{Timestamp a group or dataset, will modify timestamps recursively}
  # }
  private = list(
    # Fields
    err_init = "This loom object has not been created with either loomR::create or loomR::connect, please use these functions to create or connect to a loom file",
    err_mode = "Cannot modify a loom file in read-only mode",
    err_ambiguous = "Cannot identify the dataset provided, found too many like it; please be more specific",
    deprecated = "batch.scan and batch.next have been deprecated; please use map/apply instead. For generating chunks, please use loom$chunk.points or loom$chunk.indices",
    iter.scan = NULL,
    it = NULL,
    iter.chunksize = NULL,
    iter.dataset = NULL,
    iter.margin = NULL,
    iter.index = NULL,
    skipped.validation = FALSE,
    # Methods
    load_attributes = function(MARGIN) {
      attribute <- switch(
        EXPR = MARGIN,
        '1' = 'row_attrs',
        '2' = 'col_attrs',
        stop('Invalid attribute dimension')
      )
      group <- self[[attribute]]
      if (length(x = names(x = group)) > 0) {
        attributes <- unlist(x = lapply(
          X = names(x = group),
          FUN = function(x) {
            d <- list(group[[x]])
            names(x = d) <- x
            return(d)
          }
        ))
      } else {
        attributes <- NULL
      }
      if (MARGIN == 1) {
        self$row.attrs <- attributes
      } else if (MARGIN == 2) {
        self$col.attrs <- attributes
      }
    },
    load_layers = function() {
      if (private$skipped.validation) {
        if (getOption(x = 'verbose')) {
          warning("Not loading layers field")
        }
      } else if (length(x = names(x = self[['layers']])) > 0) {
        self$layers <- unlist(x = lapply(
          X = names(x = self[['layers']]),
          FUN = function(n) {
            d <- list(self[[paste('layers', n, sep = '/')]])
            names(x = d) <- n
            return(d)
          }
        ))
      }
    },
    reset_batch = function() {
      private$it <- NULL
      private$iter.chunksize <- NULL
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
    },
    timestamp = function(x = NULL, silent = FALSE) {
      # Check read-only mode
      if (self$mode == 'r') {
        if (silent) {
          invisible(x = NULL)
        } else {
          stop(private$err_mode)
        }
      }
      # Timestamp format, always in UTC
      # %Y -- four-digit year
      # %m -- two-digit month
      # %d -- two-digit day
      # T -- sepearte date from time
      # %H -- 24-hour time
      # %M -- two-digit minute
      # %OS6 -- seconds with 6 digits after the decimal
      # Z -- end time
      # Generate time
      time <- strftime(x = Sys.time(), format = '%Y%m%dT%H%M%OS6Z', tz = 'UTC')
      if (!is.null(x = x)) {
        # Ensure all groups and datasets exist
        if (!all(sapply(X = x, FUN = self$exists))) {
          stop("Cannot find all groups/datasets passed")
        }
        # Find all groups containing passed groups/datasets
        x <- strsplit(x = x, split = '/')
        x <- lapply(
          X = x,
          FUN = function(split.name) {
            split.name <- Filter(f = nchar, x = split.name)
            combinations <- vector(mode = 'character', length = length(x = split.name))
            for (i in 1:length(x = split.name)) {
              combinations[i] <- paste(split.name[1:i], collapse = '/')
            }
            return(combinations)
          }
        )
        x <- unique(x = unlist(x = x))
      }
      # Add timestamp
      h5attr(x = self, which = 'last_modified') <- time
      for (i in x) {
        h5attr(x = self[[i]], which = 'last_modified') <- time
      }
    }
  )
)

#' Create a loom object
#'
#' @param filename The name of the new loom file
#' @param data The data for \code{/matrix}. If cells are rows and genes are columns, set \code{do.transpose = FALSE}; otherwise, set \code{do.transpose = TRUE}
#' @param gene.attrs A named list of vectors with extra data for genes, each vector must be as long as the number of genes in \code{data}
#' @param cell.attrs A named list of vectors with extra data for cells, each vector must be as long as the number of cells in \code{data}
#' Will store in 'col_attrs/nUMI' and 'col_attrs/nGene', overwriting anything passed to \code{cel.attrs};
#' @param layers A named list of matrices to be added as layers
#' @param do.transpose Transpose the input? Should be \code{TRUE} if \code{data} has genes as rows and cells as columns
#' @param calc.numi Calculate number of UMIs and genes expressed per cell?
#' To set a custom threshold for gene expression, pass an integer value (eg. \code{calc.numi = 5} for a threshold of five counts per cell)
#' @param max.size Set maximum chunk size in terms of memory usage, unused if \code{chunk.dims} is set;
#' may pass a character string (eg. \code{3gb}, \code{1200mb}) or exact value in bytes
#' @param dtype Data type (h5type) used in matrix; auto-determined by default
#' @param chunk.dims Matrix chunk dimensions; auto-determined by default
#' @param chunk.size Maximum number of cells read/written to disk at once; auto-determined by default
#' @param overwrite Overwrite an already existing loom file?
#' @param display.progress Display a progress bar
#' @param ... Ignored for now
#'
#' @return A connection to a loom file
#'
#' @importFrom utils packageVersion setTxtProgressBar
#'
#' @seealso \code{\link{loom}}
#'
#' @export
#'
create <- function(
  filename,
  data,
  gene.attrs = NULL,
  cell.attrs = NULL,
  layers = NULL,
  do.transpose = TRUE,
  calc.numi = FALSE,
  max.size = '400mb',
  dtype = NULL,
  chunk.dims = NULL,
  chunk.size = NULL,
  overwrite = FALSE,
  display.progress = TRUE,
  ...
) {
  mode <- ifelse(test = overwrite, yes = 'w', no = 'w-')
  if (file.exists(filename) && !overwrite) {
    stop(paste('File', filename, 'already exists!'))
  }
  if (!inherits(x = data, what = c('matrix', 'Matrix'))) {
    data <- as.matrix(x = data)
  }
  new.loom <- loom$new(filename = filename, mode = mode)
  if (is.null(x = dtype)) {
    dtype <- guess_dtype(x = data[1, 1], string_len = getOption(x = "loomR.string_len"))
  } else if (!inherits(x = dtype, what = 'H5T')) {
    stop("'dtype' must be an HDF5 type, please see '?h5types'")
  }
  matrix.shape <- dim(x = data)
  if (do.transpose) {
    cate("Transposing input data: loom file will show input columns (cells) as rows and input rows (genes) as columns")
    cate("This is to maintain compatibility with other loom tools")
    matrix.shape <- rev(x = matrix.shape)
  } else {
    cate("Not tranposing data: loom file will show data exactly like input")
    cate("Please note, other loom tools will show this flipped")
  }
  matrix.space <- H5S$new(
    type = 'simple',
    dims = matrix.shape,
    maxdims = c(Inf, matrix.shape[2])
  )
  if (is.null(x = chunk.dims)) {
    mem.size <- charToBytes(x = max.size)
    if (mem.size > 4e9) {
      cate("HDF5 limits internal chunk sizes to 4 GB, setting to 4 GB")
      mem.size <- 4e9
    }
    chunk.dims <- guess_chunks(
      space_maxdims = matrix.space$maxdims,
      dtype_size = dtype$get_size(),
      chunk_size = mem.size
    )
    gc(verbose = FALSE)
  }
  chunk.dims <- as.integer(x = chunk.dims)
  chunk.dims <- pmin(chunk.dims, matrix.shape)
  new.loom$create_dataset(
    name = 'matrix',
    dtype = dtype,
    space = matrix.space,
    chunk_dims = chunk.dims,
    gzip_level = 4
  )
  if (is.null(x = chunk.size)) {
    chunk.size <- chunk.dims[1]
  }
  chunk.points <- chunkPoints(
    data.size = matrix.shape[1],
    chunk.size = chunk.size
  )
  h5attr(x = new.loom, which = 'chunks') <- paste0(
    '(',
    paste(chunk.dims, collapse = ', '),
    ')'
  )
  if (display.progress) {
    pb <- newPB()
  }
  for (col in 1:ncol(x = chunk.points)) {
    row.start <- chunk.points[1, col]
    row.end <- chunk.points[2, col]
    data.add <- if (do.transpose) {
      t(x = as.matrix(x = data[, row.start:row.end]))
    } else {
      as.matrix(x = data[row.start:row.end, ])
    }
    new.loom[['matrix']][row.start:row.end, ] <- data.add
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
    }
    gc(verbose = FALSE)
  }
  if (display.progress) {
    close(con = pb)
  }
  new.loom$matrix <- new.loom[['matrix']]
  new.loom$shape <- rev(x = new.loom[['matrix']]$dims)
  # Groups
  groups <- c('layers', 'row_attrs', 'col_attrs', 'row_graphs', 'col_graphs')
  for (group in groups) {
    new.loom$create_group(name = group)
  }
  # Check for the existance of gene or cell names
  if (!is.null(x = colnames(x = data))) {
    if (do.transpose) {
      new.loom$add.col.attribute(attribute = list('cell_names' = colnames(x = data)))
    } else {
      new.loom$add.row.attribute(attribute = list('gene_names' = colnames(x = data)))
    }
  }
  if (!is.null(x = rownames(x = data))) {
    if (do.transpose) {
      new.loom$add.row.attribute(attribute = list('gene_names' = rownames(x = data)))
    } else {
      new.loom$add.col.attribute(attribute = list('cell_names' = rownames(x = data)))
    }
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
    new.loom$add.layer(layer = layers, display.progress = display.progress)
  }
  if (!is.null(x = gene.attrs)) {
    new.loom$add.row.attribute(attribute = gene.attrs)
  }
  if (!is.null(x = cell.attrs)) {
    new.loom$add.col.attribute(attribute = cell.attrs)
  }
  # Calculate nUMI and nGene if asked
  if (calc.numi) {
    is.expr <- ifelse(test = is.numeric(x = calc.numi), yes = calc.numi, no = 0)
    calc.umi(
      object = new.loom,
      chunk.size = chunk.size,
      display.progress = display.progress,
      is.expr = is.expr
    )
  }
  # Return the connection
  return(new.loom)
}

#' Connect to a loom file
#'
#' @param filename The loom file to connect to
#' @param mode How do we connect to it? Pass 'r' for read-only or 'r+' for read/write.
#' If \code{mode} is 'r+', loomR will automatically add missing required groups during validation
#' @param skip.validate Skip the validation steps, use only for extremely large loom files
#'
#' @return A loom file connection
#'
#' @seealso \code{\link{loom}}
#'
#' @export
#'
connect <- function(filename, mode = "r", skip.validate = FALSE) {
  if (!(mode %in% c('r', 'r+'))) {
    stop("'mode' must be one of 'r' or 'r+'")
  }
  new.loom <- loom$new(filename = filename, mode = mode, skip.validate = skip.validate)
  return(new.loom)
}

#' Subset a loom file
#'
#' @param x A loom object
#' @param m Rows (cells) to subset, defaults to all rows
#' @param n Columns (genes) to subset, defaults to all columns
#' @param filename Filename for new loom object, defaults to ...
#' @param chunk.size Chunk size to iterate through \code{x}
#' @param overwrite Overwrite \code{filename} if already exists?
#' @param display.progress Display a progress bar
#' @param ... Ignored for now
#'
#' @return A loom object connected to \code{filename}
#'
#' @aliases subset
#'
#' @importFrom utils setTxtProgressBar
#'
#' @seealso \code{\link{loom}}
#'
#' @export
#' @method subset loom
#'
subset.loom <- function(
  x,
  m = NULL,
  n = NULL,
  filename = NULL,
  chunk.size = 1000,
  overwrite = FALSE,
  display.progress = TRUE,
  ...
) {
  m.all <- is.null(x = m)
  n.all <- is.null(x = n)
  if (m.all && n.all) {
    stop("Subset is set to all data, not subsetting")
  }
  # Set some defaults
  if (m.all) { # Pull all cells
    m <- 1:x[['matrix']]$dims[1]
  }
  if (n.all) { # Pull all genes
    n <- 1:x[['matrix']]$dims[2]
  }
  if (is.null(x = filename)) { # Set default filename
    filename <- paste(
      unlist(x = strsplit(x = x$filename, split = '.', fixed = TRUE)),
      collapse = '_subset.'
    )
  }
  # Ensure that m and n are within the bounds of the loom file
  if (max(m) > x[['matrix']]$dims[1] || max(n) > x[['matrix']]$dims[2]) {
    stop(paste(
      "'m' and 'n' must be less than",
      x[['matrix']]$dims[1],
      "and",
      x[['matrix']]$dims[2],
      "respectively"
    ))
  }
  extension <- rev(x = unlist(x = strsplit(x = filename, split = '.', fixed = TRUE)))[1]
  # Ensure that we call our new file a loom file
  if (extension != 'loom') {
    filename <- paste0(filename, '.loom')
  }
  # Make the loom file
  mode <- ifelse(test = overwrite, yes = 'w', no = 'w-')
  if (file.exists(filename) && !overwrite) {
    stop(paste('File', filename, 'already exists!'))
  } else if (display.progress) {
    catn("Writing new loom file to", filename)
  }
  new.loom <- loom$new(filename = filename, mode = mode)
  # Add /matrix
  matrix.dims <- c(length(x = m), length(x = n))
  new.loom$create_dataset(
    name = 'matrix',
    dtype = getDtype2(x = class(x = x[['matrix']]$get_type())[1]),
    dims = matrix.dims
  )
  new.loom$shape <- rev(x = matrix.dims)
  batch <- break.consecutive(
    x = m,
    max.length = chunk.size,
    min.length = chunk.size * 0.5
  )
  previous <- 1L
  layers <- list.datasets(object = x, path = 'layers', full.names = TRUE)
  if (display.progress) {
    msg <- "Adding data for /matrix"
    num.layers <- length(x = layers)
    if (num.layers) {
      msg <- paste(msg, 'and', num.layers, 'layers')
    }
    catn(msg)
    if (!num.layers) {
      cate("No layers found")
    }
    pb <- newPB()
  }
  if (length(x = layers) > 0) {
    for (layer in layers) {
      new.loom[['layers']]$create_dataset(
        name = basename(path = layer),
        dtype = getDtype2(x = class(x = x[[layer]]$get_type())[1]),
        dims = matrix.dims
      )
    }
  }
  for (i in 1:length(x = batch)) {
    indices.use <- batch[[i]]
    chunk.indices <- min(indices.use):max(indices.use)
    indices.return <- previous:(previous + length(x = indices.use) - 1)
    indices.use <- chunk.indices[chunk.indices %in% indices.use] - chunk.indices[1] + 1
    previous <- max(indices.return) + 1L
    # chunk.indices <- x$batch.next(return.data = FALSE)
    # indices.use <- chunk.indices[chunk.indices %in% m]
    # if (length(x = indices.use) < 1) {
    #   if (display.progress) {
    #     setTxtProgressBar(pb = pb, value = i / length(x = batch))
    #   }
    #   next
    # }
    # indices.return <- match(x = indices.use, table = m)
    # indices.use <- indices.use - chunk.indices[1] + 1
    chunk.data <- x[['matrix']][chunk.indices, ]
    chunk.data <- matrix(data = chunk.data, nrow = length(x = chunk.indices))
    chunk.data <- chunk.data[indices.use, n]
    # chunk.data <- chunk.data[indices.use, n]
    new.loom[['matrix']][indices.return, ] <- chunk.data
    gc(verbose = FALSE)
    if (length(x = layers) > 0) {
      for (layer in layers) {
        layer.data <- x[[layer]][chunk.indices, ]
        layer.data <- layer.data[indices.use, n]
        new.loom[[layer]][indices.return, ] <- layer.data
        gc(verbose = FALSE)
      }
    }
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = i / length(x = batch))
    }
  }
  if (display.progress) {
    close(con = pb)
  }
  # Add groups
  for (group in c('row_attrs', 'row_graphs', 'col_attrs', 'col_graphs', 'layers')) {
    new.loom$create_group(name = group)
  }
  # Add row attributes
  row.attrs <- list.datasets(object = x, path = 'row_attrs', full.names = TRUE)
  if (length(x = row.attrs) > 0) {
    if (display.progress) {
      catn("Adding", length(x = row.attrs), "row attributes")
      pb <- newPB()
      counter <- 0
    }
    for (row.attr in row.attrs) {
      if (n.all) {
        new.loom[['row_attrs']]$obj_copy_from(
          src_loc = x,
          src_name = row.attr,
          dst_name = basename(path = row.attr)
        )
      } else if (length(x = x[[row.attr]]$dims) == 2) {
        new.loom[[row.attr]] <- x[[row.attr]][,][, n]
      } else {
        new.loom[[row.attr]] <- x[[row.attr]][][n]
      }
      gc(verbose = FALSE)
      if (display.progress) {
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter / length(x = row.attrs))
      }
    }
    if (display.progress) {
      close(con = pb)
    }
  } else {
    cate("No row attributes found")
  }
  # Add col attributes
  col.attrs <- list.datasets(object = x, path = 'col_attrs', full.names = TRUE)
  if (length(x = col.attrs) > 0) {
    if (display.progress) {
      catn("Adding", length(x = col.attrs), "column attributes")
      pb <- newPB()
      counter <- 0
    }
    for (col.attr in col.attrs) {
      if (m.all) {
        new.loom[['col_attrs']]$obj_copy_from(
          src_loc = x,
          src_name = col.attr,
          dst_name = basename(path = col.attr)
        )
      } else if (length(x = x[[col.attr]]$dims) == 2) {
        new.loom[[col.attr]] <- x[[col.attr]][, m]
      } else {
        new.loom[[col.attr]] <- x[[col.attr]][m]
      }
      gc(verbose = FALSE)
      if (display.progress) {
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter / length(x = col.attrs))
      }
    }
    if (display.progress) {
      close(con = pb)
    }
  } else {
    cate("No column attributes found")
  }
  # # Add layers
  # layers <- list.datasets(object = x, path = 'layers', full.names = TRUE)
  # if (length(x = layers) > 0) {
  #   if (display.progress) {
  #     catn("Adding", length(x = layers), "layers")
  #     pb <- newPB()
  #   }
  #   # Initialize datasets
  #   for (layer in layers) {
  #     new.loom[['layers']]$create_dataset(
  #       name = basename(path = layer),
  #       dtype = getDtype2(x = class(x = x[[layer]]$get_type())[1]),
  #       dims = matrix.dims
  #     )
  #   }
  #   # batch <- x$batch.scan(chunk.size = chunk.size)
  #   for (i in 1:length(x = batch)) {
  #     # chunk.indices <- x$batch.next(return.data = FALSE)
  #     indices.use <- batch[[i]]
  #     chunk.indices <- min(indices.use):max(indices.use)
  #     indices.use <- chunk.indices[chunk.indices %in% m]
  #     indices.return <- match(x = indices.use, table = m)
  #     indices.use <- indices.use - chunk.indices[1] + 1
  #     # Add the data for each layer for this batch
  #     for (layer in layers) {
  #       layer.data <- x[[layer]][chunk.indices, ]
  #       layer.data <- layer.data[indices.use, n]
  #       new.loom[[layer]][indices.return, ] <- layer.data
  #       gc(verbose = FALSE)
  #     }
  #     if (display.progress) {
  #       setTxtProgressBar(pb = pb, value = i / length(x = batch))
  #     }
  #   }
  #   if (display.progress) {
  #     close(con = pb)
  #   }
  # } else {
  #   cate("No layers found")
  # }
  new.loom$flush()
  new.loom$load.fields()
  return(new.loom)
}

#' Combine loom files
#'
#' @param looms A list of \code{loom} objects or paths to loom files
#' @param filename Name for resultant loom file
#' @param chunk.size How many rows from each input loom should we stream to the merged loom file at any given time?
#' @param order.by Optional row attribute to order each input loom by, must be one dimensional
#' @param overwrite Overwrite \code{filename} if already exists?
#' @param display.progress Display a progress bar as we're copying over data
#' @param chunk.dims Chunk dimensions used for matrix; will be copied from first file by default
#' @param ... Ignored for now
#'
#' @return A loom object connected to \code{filename}
#'
#' @importFrom utils setTxtProgressBar
#'
#' @seealso \code{\link{loom}}
#'
#' @export
#'
combine <- function(
  looms,
  filename,
  chunk.size = 1000,
  order.by = NULL,
  overwrite = FALSE,
  display.progress = TRUE,
  chunk.dims = NULL,
  ...
) {
  if (is.character(x = looms)) {
    looms <- as.list(x = looms)
  }
  if (!is.list(x = looms)) {
    stop("'combine' takes a list of loom objects or paths to loom files")
  }
  # Basic checking of input arguments
  looms <- looms[sapply(X = looms, FUN = inherits, what = c('loom', 'character'))]
  if (length(x = looms) < 2) {
    stop("Need at least two loom objects or files to merge")
  }
  looms <- unlist(x = looms, recursive = TRUE, use.names = FALSE)
  # length.check <- unique(x = vapply(X = looms, FUN = length, FUN.VALUE = integer(length = 1L)))
  # if (length(x = length.check) != 1 || length.check != 1) {
  #   stop("Each entry in the list of looms to combine must have a length of one")
  # }
  # Ensure order.by comes from row_attrs
  if (!is.null(x = order.by) && dirname(path = order.by) != 'row_attrs') {
    order.by <- paste('row_attrs', order.by, sep = '/')
  }
  # Check the existance of loom files
  loom.names <- looms[vapply(X = looms, FUN = is.character, FUN.VALUE = logical(length = 1L))]
  loom.names <- unlist(x = loom.names)
  if (length(x = loom.names) > 0) {
    if (!all(file.exists(loom.names))) {
      stop(paste0(
        "Cannot find the following loom files: '",
        paste(loom.names[!file.exists(loom.names)], collapse = "', '"),
        "'"
      ))
    }
  }
  # Set mode and provide more useful error
  mode <- ifelse(test = overwrite, yes = 'w', no = 'w-')
  if (file.exists(filename) && !overwrite) {
    stop(paste('File', filename, 'already exists!'))
  }
  # Check loom contents
  # Every loom must have same number of genes (rows, MARGIN = 1)
  # and same datasets in the groups
  row.attrs <- vector(mode = 'character')
  row.types <- list()
  col.attrs <- vector(mode = 'character')
  col.dims <- list()
  col.ndims <- list()
  col.types <- list()
  layers <- vector(mode = 'character')
  layers.types <- list()
  nrows <- vector(mode = 'integer', length = length(x = looms))
  ncols <- vector(mode = 'integer', length = length(x = looms))
  matrix.type <- vector(mode = 'list', length = length(x = looms))
  loom.chunk.dims <- vector(mode = 'list', length = length(x = looms))
  if (display.progress) {
    catn("Validating", length(x = looms), "input loom files")
    pb <- newPB()
  }
  for (i in 1:length(x = looms)) {
    this <- if (is.character(x = looms[[i]])) {
      connect(filename = looms[[i]])
    } else {
      looms[[i]]
    }
    row.attrs <- unique(x = c(
      row.attrs,
      list.datasets(
        object = this,
        path = 'row_attrs',
        full.names = TRUE
      )
    ))
    # row.attrs[[i]] <- sort(x = list.datasets(
    #   object = this,
    #   path = 'row_attrs',
    #   full.names = TRUE
    # ))
    for (attr in row.attrs) {
      if (this$exists(name = attr)) {
        row.types[[attr]] <- c(
          row.types[[attr]],
          class(x = this[[attr]]$get_type())[1]
        )
      }
    }
    if (!is.null(x = order.by) && !this$exists(name = order.by)) {
      if (is.character(x = looms[[i]])) {
        this$close_all()
      }
      stop(
        "Cannot find the order.by dataset",
        order.by,
        "in the file",
        this$filename
      )
    }
    col.attrs <- unique(x = c(
      col.attrs,
      list.datasets(
        object = this,
        path = 'col_attrs',
        full.names = TRUE
      )
    ))
    # col.attrs[[i]] <- sort(x = list.datasets(
    #   object = this,
    #   path = 'col_attrs',
    #   full.names = TRUE
    # ))
    for (attr in col.attrs) {
      if (this$exists(name = attr)) {
        col.types[[attr]] <- c(
          col.types[[attr]],
          class(x = this[[attr]]$get_type())[1]
        )
        col.dims[[attr]] <- sum(col.dims[[attr]], this[[attr]]$dims[1])
        col.ndims[[attr]] <- c(
          col.ndims[[attr]],
          length(x = this[[attr]]$dims)
        )
        # I added this, since otherwise, the number of rows of 2D column attributes
        # will be doubled
        if (all(col.ndims[[attr]] == 2)) {
          col.dims[[attr]] <- this[[attr]]$dims[1]
        }
      }
    }
    layers <- unique(x = c(
      layers,
      list.datasets(
        object = this,
        path = 'layers',
        full.names = TRUE
      )
    ))
    # layers[[i]] <- sort(x = list.datasets(
    #   object = this,
    #   path = 'layers',
    #   full.names = TRUE
    # ))
    for (lay in layers) {
      if (this$exists(name = lay)) {
        layers.types[[lay]] <- c(
          layers.types[[lay]],
          class(x = this[[lay]]$get_type())[1]
        )
      }
    }
    nrows[i] <- this[['matrix']]$dims[2]
    ncols[i] <- this[['matrix']]$dims[1]
    matrix.type[[i]] <- class(x = this[['matrix']]$get_type())[1]
    # loom.chunk.dims[[i]] <- as.numeric(strsplit(h5attr(this, 'chunks'), split = '[\\(, \\)]+')[[1]][2:3])
    loom.chunk.dims[[i]] <- this[['matrix']]$chunk_dims
    if (is.character(x = looms[[i]])) {
      this$close_all()
    }
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = i / length(x = looms))
    }
  }
  if (display.progress) {
    close(con = pb)
  }
  nrows <- unique(x = nrows)
  ncells <- sum(ncols)
  if (is.null(chunk.dims)) {
    chunk.dims <- loom.chunk.dims[[1]]  # for now, just copy from first file
  }
  # if (length(x = row.attrs) != 1) {
  #   stop("Not all loom objects have the same row attributes")
  # }
  # if (length(x = col.attrs) != 1) {
  #   stop("Not all loom objects have the same column attributes")
  # }
  # if (length(x = layers) != 1) {
  #   stop("Not all loom objects have the same layers")
  # }
  if (length(x = nrows) != 1) {
    stop("Not all loom objects have the number of rows/genes (MARGIN = 1)")
  }
  # Check for the row attribute to order by
  if (!is.null(x = order.by)) { # We have something to order by
    if (!grepl(pattern = order.by, x = row.attrs)) { # Check to ensure that the order.by is in each of our row attributes
      stop(paste0("Cannot find '", order.by, "' in the row attributes for the loom files provided"))
    } else {
      # If it is, get the values to order by
      order.by <- basename(path = order.by) # Basename of the order.by attribute name
      temp <- if (is.character(x = looms[1])) { # Connect to the loom first loom file
        connect(filename = looms[1])
      } else {
        looms[1]
      }
      order.dat <- temp[['row_attrs']][[order.by]] # Ensure that our order.by attribute is one dimmensional
      if (length(x = order.dat$dims) != 1) {
        if (is.character(x = looms[1])) {
          temp$close_all()
        }
        stop("'order.by' must reference a one dimensional attribute")
      }
      # If so, pull the data
      order.use <- order.dat[]
      # If the first loom was initially closed, close it again
      if (is.character(x = looms[1])) {
        temp$close_all()
      }
    }
  }
  # Check data types:
  matrix.type <- unlist(x = unique(x = matrix.type))
  row.types <- lapply(X = row.types, FUN = unique)
  row.types.counts <- vapply(X = row.types, FUN = length, FUN.VALUE = integer(length = 1L))
  col.types <- lapply(X = col.types, FUN = unique)
  col.types.counts <- vapply(X = col.types, FUN = length, FUN.VALUE = integer(length = 1L))
  layers.types <- lapply(X = layers.types, FUN = unique)
  layers.types.counts <- vapply(X = layers.types, FUN = length, FUN.VALUE = integer(length = 1L))
  if (any(row.types.counts > 1)) {
    stop(paste0(
      "The following row attributes have multiple types across the input loom files: '",
      paste(names(x = row.types.counts[row.types.counts > 1]), collapse = "', '"),
      "'; cannot combine"
    ))
  }
  if (any(col.types.counts > 1)) {
    stop(paste0(
      "The following column attributes have multiple types across the input loom files: '",
      paste(names(x = col.types.counts[col.types.counts > 1]), collapse = "', '"),
      "'; cannot combine"
    ))
  }
  if (any(layers.types.counts > 1)) {
    stop(paste0(
      "The following layers have multiple types across the input loom files: '",
      paste(names(x = layers.types.counts[layers.types.counts > 1]), collapse = "', '"),
      "'; cannot combine"
    ))
  }
  if (length(x = matrix.type) != 1) {
    stop("Cannot combine multiple datatypes for '/matrix'")
  }
  # Check dimmensions for col_attrs
  col.ndims <- lapply(X = col.ndims, FUN = unique)
  col.ndims.counts <- vapply(X = col.ndims, FUN = length, FUN.VALUE = integer(length = 1L))
  if (any(col.ndims.counts > 1)) {
    stop(paste0(
      "The following column attributes have varying dimmensions across the input loom files: '",
      paste(names(x = col.ndims.counts[col.ndims.counts > 1]), collapse = "', '"),
      "'; cannot combine"
    ))
  } else if (any(!col.ndims %in% c(1, 2))) {
    stop("loomR only supports one- and two-dimensional attributes")
  }
  # Create the new HDF5 file and the required groups
  new.loom <- loom$new(filename = filename, mode = mode)
  groups <- c('layers', "row_attrs", "col_attrs", "row_graphs", "col_graphs")
  for (group in groups) {
    new.loom$create_group(name = group)
  }
  # Initialize the '/matrix' dataset as well as datasets in 'col_attrs' and 'layers'
  col.attrs <- unlist(col.attrs)
  row.attrs <- unlist(row.attrs)
  new.loom$create_dataset(
    name = 'matrix', # Name is '/matrix'
    dtype = getDtype2(x = matrix.type), # Use the single type that we got from above
    dims = c(ncells, nrows), # Use the number of cells from the sum of ncols above, nrows should be the same for everyone
    chunk_dims = chunk.dims,
    gzip_level = 4
  )
  for (lay in layers) {
    #if (length(x = lay) > 1) {
      new.loom$create_dataset(
        name = lay,
        dtype = getDtype2(x = layers.types[[lay]]),
        dims = c(ncells, nrows),
        chunk_dims = chunk.dims,
        gzip_level = 4
      )
    #}
  }
  for (attr in col.attrs) {
    if (length(x = attr) > 0) {
      dims.use <- switch(
        EXPR = col.ndims[[attr]],
        '1' = ncells,
        '2' = c(col.dims[[attr]], ncells)
      )
      new.loom$create_dataset(
        name = attr,
        dtype = getDtype2(x = col.types[[attr]]),
        dims = dims.use
      )
    }
  }
  # Start adding loom objects
  matrix.previous <- 0
  col.previous <- vector(mode = 'integer', length = length(x = col.attrs))
  names(x = col.previous) <- col.attrs
  for (i in 1:length(x = looms)) {
    if (display.progress) {
      catn("\nAdding loom file", i ,"of", length(x = looms))
    }
    # Open the loom file
    this <- if (is.character(x = looms[[i]])) {
      connect(filename = looms[[i]])
    } else {
      looms[[i]]
    }
    # Get the chunk points
    chunk.points <- chunkPoints(data.size = ncols[[i]], chunk.size = chunk.size)
    # Pull ordering information
    order.genes <- if (is.null(x = order.by)) {
      1:nrows # If order.by wasn't provided, just use 1:number of genes
    } else {
      # Otherwise, match the data order from our order.use
      order(match(x = this[[row.attrs]][[order.by]][], table = order.use))
    }
    # Add data to /matrix, /col_attrs, and /layers
    if (display.progress) {
      catn("Adding data to /matrix and /layers")
      pb <- newPB()
    }
    for (col in 1:ncol(x = chunk.points)) {
      cells.use <- chunk.points[1, col]:chunk.points[2, col]
      # Add /matrix for this chunk
      matrix.add <- this[['matrix']][cells.use, ]
      new.loom[['matrix']][cells.use + matrix.previous, ] <- matrix.add[, order.genes]
      # Add layers for this chunk
      for (lay in layers) {
        if (!this$exists(name = lay)) {
          if (display.progress) {
            cate("Cannot find layer", lay, "in loom file", this$filename)
          }
          next
        }
        lay.add <- this[[lay]][cells.use, ]
        new.loom[[lay]][cells.use + matrix.previous, ] <- lay.add[, order.genes]
      }
      if (display.progress) {
        setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
      }
      gc(verbose = FALSE)
    }
    if (display.progress) {
      close(con = pb)
    }
    matrix.previous <- matrix.previous + max(chunk.points)
    # Add col_attrs for this chunk
    if (display.progress) {
      catn("\nAdding data to /col_attrs")
      pb <- newPB()
    }
    for (j in 1L:length(x = col.attrs)) {
      attr <- col.attrs[j]
      start <- col.previous[attr]
      end <- start + this[['matrix']]$dims[1]
      col.previous[attr] <- end
      if (!this$exists(name = attr)) {
        if (display.progress) {
          cate("Cannot find attribute", attr, "in loom file", this$filename)
          setTxtProgressBar(pb = pb, value = j / length(x = col.attrs))
        }
        next
      } else if (col.ndims[[attr]] == 1) {
        new.loom[[attr]][(start + 1):end] <- this[[attr]][]
      } else {
        for (col in 1:ncol(x = chunk.points)) {
          col.use <- chunk.points[1, col]:chunk.points[2, col]
          new.loom[[attr]][, col.use] <- this[[attr]][, col.use]
        }
      }
      if (display.progress) {
        setTxtProgressBar(pb = pb, value = j / length(x = col.attrs))
      }
      gc(verbose = FALSE)
    }
    if (display.progress) {
      close(con = pb)
    }
    # Copy row attributes from the first loom object into the merged one
    if (i == 1) {
      if (display.progress) {
        catn("\nCopying data for /row_attrs")
        pb <- newPB()
      }
      for (attr in list.datasets(object = this, path = 'row_attrs')) {
        new.loom[['row_attrs']]$obj_copy_from(
          src_loc = this[['row_attrs']],
          src_name = attr,
          dst_name = attr
        )
        if (display.progress) {
          setTxtProgressBar(
            pb = pb,
            value = grep(
              pattern = paste0("^", attr, "$"),
              x = list.datasets(object = this[['row_attrs']])
            ) / length(x = list.datasets(object = this[['row_attrs']]))
          )
        }
      }
      if (display.progress) {
        close(con = pb)
      }
    }
    new.loom$flush()
    # Close current loom file, if not open previously
    if (is.character(x = looms[[i]])) {
      this$close_all()
    }
  }
  new.loom$load.fields()
  new.loom$update.shape()
  return(new.loom)
}

# #need to comment
# #need to add progress bar
# #but otherwise, pretty cool
# #for paul to try :
# # f <- connect("~/Downloads/10X43_1.loom")
# # mean_var = map(f,f_list = c(mean,var),chunksize = 5000)
# # nGene <- map(f, f_list = function(x) length(which(x>0)), MARGIN = 2)
# map <- function(self, f_list = list(mean, var), MARGIN=1, chunksize=1000, selection) {
#   n_func = length(f_list)
#   if (n_func == 1) {
#     f_list <- list(f_list)
#   }
#   if (MARGIN == 1) {
#     results <- list()
#     for (j in 1:n_func) {
#       results[[j]] <- numeric(0)
#     }
#     rows_per_chunk <- chunksize
#     ix <- 1
#     while (ix <= self@shape[1]) {
#       rows_per_chunk <- min(rows_per_chunk, self@shape[1] - ix + 1)
#       chunk <- self["matrix"][ix:(ix + rows_per_chunk - 1), ]
#       for (j in 1:n_func) {
#         new_results <- apply(chunk, 1, FUN = f_list[[j]])
#         results[[j]] <- c(results[[j]], new_results)
#       }
#       ix <- ix + chunksize
#     }
#   }
#   if (MARGIN == 2) {
#     results <- list()
#     for (j in 1:n_func) {
#       results[[j]] <- numeric(0)
#     }
#     cols_per_chunk <- chunksize
#     ix <- 1
#     while (ix <= self@shape[2]) {
#       cols_per_chunk <- min(cols_per_chunk, self@shape[2] - ix + 1)
#       chunk <- self["matrix"][, ix:(ix + cols_per_chunk - 1)]
#       for (j in 1:n_func) {
#         new_results <- apply(chunk, 2, FUN = f_list[[j]])
#         results[[j]] <- c(results[[j]], new_results)
#       }
#       ix <- ix + chunksize
#     }
#   }
#   if (n_func == 1) {
#     results <- results[[1]]
#   }
#   return(results)
# }
