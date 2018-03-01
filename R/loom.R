#' @include internal.R
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
#' @field shape Shape of \code{/matrix} in genes (columns) by cells (rows)
#' @field chunksize Chunks set for this dataset in columns (cells) by rows (genes)
#' @field matrix The main data matrix, stored as columns (cells) by rows (genes)
#' @field layers Additional data matricies, the same shape as \code{/matrix}
#' @field col.attrs Extra information about cells
#' @field row.attrs Extra information about genes
#'
#' @section Methods:
#' \describe{
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
#'   \item{\code{batch.scan(chunk.size, MARGIN, index.use, dataset.use, force.reset)}, \code{batch.next(return.data)}}{
#'     Scan a dataset in the loom file from \code{index.use[1]} to \code{index.use[2]}, iterating by \code{chunk.size}.
#'     \code{dataset.use} can be the name, not \code{group/name}, unless the name is present in multiple groups.
#'     Pass \code{MARGIN = 1} to iterate on genes or \code{MARGIN = 2} to iterate on cells for 'matrix' or any dataset in 'layers'.
#'     To force reevaluation of the iterator object, pass \code{force.reset = TRUE}.
#'     \code{MARGIN} does not need to be set for datasets in 'row_attrs' or 'col_attrs'.
#'     \code{chunk.size} defaults to \code{self$chunksize}, \code{MARGIN} defaults to 2,
#'     \code{index.use} defaults to \code{1:self$shape[MARGIN]}, \code{dataset.use} defaults to 'matrix'
#'   }
#'   \item{\code{apply(name, FUN, MARGIN, chunk.size, dataset.use, overwrite, display.progress, ...)}}{
#'     Apply a function over a dataset within the loom file, stores the results in the loom file.
#'     \code{name} must be the full name of the dataset ('group/name').
#'     \code{apply} will always use the entire dataset specified in \code{dataset.use}
#'     \describe{
#'       \item{\code{name}}{Name of dataset to store results to}
#'       \item{\code{FUN}}{Function to apply}
#'       \item{\code{MARGIN}}{Iterate over genes (1) or cells (2)}
#'       \item{\code{chunk.size}}{Size to chunk \code{MARGIN} by}
#'       \item{\code{dataset.use}}{Name of dataset to use}
#'       \item{\code{overwrite}}{Overite \code{name} if already exists}
#'       \item{\code{display.progress}}{Display progress}
#'       \item{\code{...}}{Extra parameters to pass to \code{FUN}}
#'     }
#'   }
#'   \item{\code{map(FUN, MARGIN, chunk.size, index.use, dataset.use, display.progress, expected, ...)}}{
#'     Map a function onto a dataset within the loom file, returns the result into R.
#'     The result will default to the shape of the dataset used; to change pass either 'vector' or 'matrix' to \code{expected}.
#'   }
#'   \item{\code{add.cells(matrix.data, attributes.data = NULL, layers.data = NULL, display.progress = TRUE)}}{
#'     Add m2 cells to a loom file.
#'     \describe{
#'       \item{\code{matrix.data}}{a list of m2 cells where each entry is a vector of length n (num genes, \code{self$shape[1]})}
#'       \item{\code{attributes.data}}{a list where each entry is named for one of the datasets in \code{self[['col_attrs']]}; each entry is a vector of length m2.}
#'       \item{\code{layers.data}}{a list where each entry is named for one of the datasets in \code{self[['layers']]}; each entry is an n-by-m2 matrix where n is the number of genes in this loom file and m2 is the number of cells being added.}
#'       \item{\code{display.progress}}{display progress}
#'     }
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
        self$shape <- rev(self[['matrix']]$dims)
        # Store the chunk size
        chunks <- tryCatch(
          expr = h5attr(x = self, which = 'chunks'),
          error = function(e) {
            hchunks <- self[['matrix']]$chunk_dims
            pchunks <- paste0('(', paste(hchunks, collapse = ', '), ')')
            if (mode != 'r') {
              h5attr(x = self, which = 'chunks') <- pchunks
            }
            return(pchunks)
          }
        )
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
    # Addding attributes and layers
    add.layer = function(layers, chunk.size = 1000, overwrite = FALSE) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Value checking
      if (!is.list(x = layers) || is.null(x = names(x = layers))) {
        stop("'layers' must be a named list")
      }
      if (is.null(x = self$shape)) {
        stop(private$err_init)
      }
      # Add layers
      for (i in 1:length(x = layers)) {
        # if (!is.matrix(x = layers[[i]])) {
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
        dtype <- getDtype(x = layers[[i]][1, 1])
        self[['layers']]$create_dataset(
          name = layer.name,
          dtype = dtype,
          dims = dim(x = layers[[i]])
        )
        chunk.points <- chunkPoints(
          data.size = dim(x = layers[[i]])[1],
          chunk.size = chunk.size
        )
        # if (display.progress) {
        #   pb <- txtProgressBar(char = '=', style = 3)
        # }
        for (col in 1:ncol(x = chunk.points)) {
          row.start <- chunk.points[1, col]
          row.end <- chunk.points[2, col]
          self[['layers']][[layer.name]][row.start:row.end, ] <- as.matrix(
            x = layers[[i]][row.start:row.end, ]
          )
          # if (display.progress) {
          #   setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
          # }
        }
        # self[['layers']]$create_dataset(
        #   name = names(x = layers)[i],
        #   robj = layers[[i]],
        #   chunk_dims = self$chunksize
        # )
      }
      self$flush()
      gc(verbose = FALSE)
      private$load_layers()
      invisible(x = self)
    },
    add.attribute = function(attribute, MARGIN, overwrite = FALSE) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Value checking
      if (is.data.frame(x = attribute)) {
        attribute <- as.list(x = attribute)
      }
      is.actual.list <- is.list(x = attribute)
      if (!is.actual.list || is.null(x = names(x = attribute))) {
        stop("Attributes must be provided as a named list")
      }
      # if (is.data.frame(x = attribute)) {
      #   attribute <- as.list(x = attribute)
      # }
      # if (!is.list(x = attribute) || is.null(x = names(x = attribute))) {
      #   stop("'attribute' must be a named list")
      # }
      if (!MARGIN %in% c(1, 2)) {
        stop("'MARGIN' must be 1 or 2")
      }
      length.use <- rev(x = self[['matrix']]$dims)[MARGIN]
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
    # Get attribute information
    get.attribute.df = function(
      attribute.layer = c("row", "col"),
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
      if (!attribute.layer %in% c("row", "col")) {
        stop("Invalid attribute.layer. Please select either 'row' or 'col'.")
      }
      attribute.layer <- paste0(attribute.layer, "_attrs")
      # check that attribute names are present
      if (!all(attribute.names %in% self[[attribute.layer]]$names)) {
        invalid.names <- attribute.names[which(!attribute.names %in% self[[attribute.layer]]$names)]
        stop(paste0("Invalid attribute.names: ", paste0(invalid.names, collapse = ", ")))
      }
      if (attribute.layer == "row_attrs") {
        combined.df <- data.frame(
          self[[paste0(attribute.layer, "/", attribute.names[1])]][],
          row.names = self[[paste0(attribute.layer, "/", row.names)]][]
        )
      } else {
        combined.df <- data.frame(
          self[[paste0(attribute.layer, "/", attribute.names[1])]][],
          row.names = self[[paste0(attribute.layer, "/", col.names)]][]
        )
      }
      if (length(x = attribute.names) > 1) {
        for (i in 2:length(x = attribute.names)) {
          combined.df[, attribute.names[i]] <- self[[paste0(attribute.layer, "/", attribute.names[i])]][]
        }
      }
      colnames(x = combined.df) <- attribute.names
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
    # Chunking functions
    batch.scan = function(
      chunk.size = NULL,
      MARGIN = 2,
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
          chunk.size <- rev(x = self$chunksize)[private$iter.margin]
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
          c('genes', 'cells')[name.check],
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
      # Ensure index.use is integers within the bounds of [1, self$shape[MARGIN]]
      if (!is.null(x = index.use)) {
        # Filter index.use to values between 1 and self$shape[MARGIN]
        index.use <- as.integer(x = index.use)
        index.use[index.use >= 1 & index.use <= self$shape[MARGIN]]
        index.use <- as.vector(x = na.omit(object = index.use))
        # If we still have values, figure out NAs, otherwise set index.use to NULL
        if (length(x = index.use) > 0) {
          # Do a trial run to figure out the class of NAs
          na.use <- NA
          if (display.progress) {
            catn("Running trial to determine class of NAs")
          }
          trial <- switch(
            EXPR = MARGIN,
            '1' = self[[dataset.use]][, 1],
            '2' = self[[dataset.use]][1, ]
          )
          trial <- FUN(trial, ...)
          if (is.list(x = trial)) {
            trial <- unlist(x = trial)
          }
          class(x = na.use) <- class(x = trial)
        } else {
          warning("No values passed to 'index.use' fall within the data, using all values")
          index.use <- NULL
        }
      }
      if (display.progress) {
        pb <- txtProgressBar(char = '=', style = 3)
      }
      # Have to initialize the dataset differently than
      # appending to it
      first <- TRUE
      for (i in 1:length(x = batch)) {
        # Get the indices we're iterating over
        these.indices <- self$batch.next(return.data = FALSE)
        if (is.null(x = index.use)) {
          chunk.indices <- these.indices
        } else {
          chunk.indices <- index.use[index.use %in% these.indices]
          chunk.na <- these.indices[!(these.indices %in% chunk.indices)]
        }
        # Get the data and apply FUN
        chunk.data <- if (dataset.matrix) {
          switch(
            EXPR = MARGIN,
            '1' = self[[dataset.use]][, chunk.indices], # Chunk genes
            '2' = self[[dataset.use]][chunk.indices, ] # Chunk cells
          )
        } else {
          self[[private$iter.datset]][chunk.indices]
        }
        chunk.data <- FUN(chunk.data, ...)
        # If this is the first iteration
        # Initialize the dataset within group, set first to FALSE
        if (first) {
          if (!is.null(x = index.use)) {
            # If we had indices to chunk on, create a holding matrix the size
            # of what the results should be, add the
            holding <- switch(
              EXPR = MARGIN,
              '1' = matrix(nrow = nrow(x = chunk.data), ncol = length(x = these.indices)),
              '2' = matrix(nrow = length(x = these.indices), ncol = ncol(x = chunk.data))
            )
            switch (
              EXPR = MARGIN,
              '1' = holding[, chunk.indices] <- chunk.data,
              '2' = holding[chunk.indices, ] <- chunk.data
            )
            chunk.data <- holding
          }
          group[[results.basename]] <- chunk.data
          first <- FALSE
        } else {
          # If we're writign to a matrix
          # Figure out which way we're writing the data
          if (results.matrix) {
            switch(
              EXPR = MARGIN,
              '1' = group[[results.basename]][, chunk.indices] <- chunk.data,
              '2' = group[[results.basename]][chunk.indices, ] <- chunk.data
            )
            if (!is.null(x = index.use)) {
              switch(
                EXPR = MARGIN,
                '1' = group[results.basename][, chunk.na] <- na.use,
                '2' = group[results.basename][chunk.na, ] <- na.use,
              )
            }
          } else {
            # Just write to the vector
            group[[results.basename]][chunk.indices] <- chunk.data
            if (!is.null(x = index.use)) {
              group[[results.basename]][chunk.na] <- na.use
            }
          }
        }
        if (display.progress) {
          setTxtProgressBar(pb = pb, value = i / length(x = batch))
        }
      }
      if (display.progress) {
        close(con = pb)
      }
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
      # Create our results holding object
      results <- vector(mode = "list", length = length(x = batch))
      if (display.progress) {
        pb <- txtProgressBar(char = '=', style = 3)
      }
      for (i in 1:length(x = batch)) {
        chunk.indices <- self$batch.next(return.data = FALSE)
        chunk.data <- if (dataset.matrix) {
          switch(
            EXPR = MARGIN,
            '1' = self[[dataset.use]][, chunk.indices], # Chunk genes
            '2' = self[[dataset.use]][chunk.indices, ] # Chunk cells
          )
        } else {
          self[[dataset.use]][chunk.indices]
        }
        results[[i]] <- FUN(chunk.data, ...)
        if (display.progress) {
          setTxtProgressBar(pb = pb, value = i / length(x = batch))
        }
      }
      # Bring result list into vector or matrix format
      if (class(results[[1]]) == 'numeric') {
        results <- unlist(x = results, use.names = FALSE)
      } else {
        if (MARGIN == 1) {
          results <- Reduce(f = cbind, x = results)
        } else if (MARGIN == 2) {
          results <- Reduce(f = rbind, x = results)
        }
      }
      if (display.progress) {
        close(con = pb)
      }
      private$reset_batch()
      return(results)
    },
    # Functions that modify `/matrix'
    add.cells = function(
      matrix.data,
      attributes.data = NULL,
      layers.data = NULL,
      display.progress = TRUE
    ) {
      if (self$mode == 'r') {
        stop(private$err_mode)
      }
      # Check inputs
      n <- self[['matrix']]$dims[2]
      if (display.progress) {
        cate("Checking inputs...")
      }
      matrix.data <- check.matrix_data(x = matrix.data, n = n)
      layers.data <- check.layers(
        x = layers.data,
        n = n,
        layers.names = names(x = self[['layers']])
      )
      attributes.data <- check.col_attrs(
        x = attributes.data,
        attrs.names = names(x = self[['col_attrs']])
      )
      # Get the number of cells we're adding
      num.cells <- c(
        nCells.matrix_data(x = matrix.data),
        nCells.layers(x = layers.data),
        nCells.col_attrs(x = attributes.data)
      )
      num.cells <- max(num.cells)
      # Flesh out the input data to add
      if (display.progress) {
        cate(paste("Adding", num.cells, "to this loom file"))
      }
      matrix.data <- addCells.matrix_data(x = matrix.data, n = n, m2 = num.cells)
      layers.data <- addCells.layers(x = layers.data, n = n, m2 = num.cells)
      attributes.data <- addCells.col_attrs(x = attributes.data, m2 = num.cells)
      # Add the input to the loom file
      dims.fill <- self[['matrix']]$dims[1]
      dims.fill <- (dims.fill + 1L):(dims.fill + num.cells)
      # Matrix data
      if (display.progress) {
        cate("Adding data to /matrix")
      }
      matrix.data <- t(x = as.matrix(x = data.frame(matrix.data)))
      self[['matrix']][dims.fill, ] <- matrix.data
      # Layer data
      if (display.progress) {
        cate("Adding data to /layers")
        pb <- new.pb()
        counter <- 0
      }
      layers.names <- names(x = self[['layers']])
      for (i in layers.names) {
        self[['layers']][[i]][dims.fill, ] <- t(x = layers.data[[i]])
        if (display.progress) {
          counter <- counter + 1
          setTxtProgressBar(pb = pb, value = counter / length(x = layers.names))
        }
      }
      # Column attributes
      if (display.progress) {
        cate("Adding data to /col_attrs")
        pb <- new.pb()
        counter <- 0
      }
      attrs.names <- names(x = self[['col_attrs']])
      for (i in attrs.names) {
        self[['col_attrs']][[i]][dims.fill] <- attributes.data[[i]]
        if (display.progress) {
          counter <- counter + 1
          setTxtProgressBar(pb = pb, value = counter / length(x = attrs.names))
        }
      }
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
              ofile$close
            }
            stop("Failed to find the gene names dataset in the other loom file")
          }
        )
        tryCatch(
          expr = self.key <- self[['row_attrs']][[self.key]][],
          error = function(e) {
            if (is.character(x = other)) {
              ofile$close
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
  # \describe{
  #   \item{\code{load_attributes(MARGIN)}}{Load attributes of a given MARGIN into \code{self$col.attrs} or \code{self$row.attrs}}
  #   \item{\code{load_layers()}}{Load layers into \code{self$layers}}
  #   \item{\code{reset_batch()}}{Reset the batch iterator fields}
  #   \item{\code{iter_range(index.use)}}{Get the range of indices for a batch iteration}
  # }
  private = list(
    # Fields
    err_init = "This loom object has not been created with either loomR::create or loomR::connect, please use these functions to create or connect to a loom file",
    err_mode = "Cannot modify a loom file in read-only mode",
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
      attributes <- unlist(x = lapply(
        X = names(x = group),
        FUN = function(x) {
          d <- list(group[[x]])
          names(x = d) <- x
          return(d)
        }
      ))
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
      } else {
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
#' @param chunk.size How many rows of \code{data} should we stream to the loom file at any given time?
#' @param overwrite Overwrite an already existing loom file?
#'
#' @return A connection to a loom file
#'
#' @importFrom utils packageVersion txtProgressBar setTxtProgressBar
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
  chunk.dims = 'auto',
  chunk.size = 1000,
  overwrite = FALSE,
  display.progress = TRUE
) {
  mode <- ifelse(test = overwrite, yes = 'w', no = 'w-')
  if (file.exists(filename) && !overwrite) {
    stop(paste('File', filename, 'already exists!'))
  }
  if (!inherits(x = data, what = c('matrix', 'Matrix'))) {
    data <- as.matrix(x = data)
  }
  if (length(x = chunk.dims) > 2 || length(x = chunk.dims) < 1) {
    stop("'chunk.dims' must be a one- or two-length integer vector or 'auto'")
  } else if (length(x = chunk.dims) == 1) {
    if (!grepl(pattern = '^auto$', x = chunk.dims, perl = TRUE)) {
      chunk.dims <- rep.int(x = as.integer(x = chunk.dims), times = 2)
    }
  } else {
    chunk.dims <- as.integer(x = chunk.dims)
  }
  new.loom <- loom$new(filename = filename, mode = mode)
  # Create the matrix
  # new.loom$create_dataset(
  #   name = 'matrix',
  #   robj = data,
  #   chunk_dims = chunk.dims
  # )
  dtype <- getDtype(x = data[1, 1])
  new.loom$create_dataset(
    name = 'matrix',
    dtype = dtype,
    dims = dim(x = data)
  )
  chunk.points <- chunkPoints(
    data.size = dim(x = data)[1],
    chunk.size = chunk.size
  )
  if (display.progress) {
    pb <- txtProgressBar(char = '=', style = 3)
  }
  for (col in 1:ncol(x = chunk.points)) {
    row.start <- chunk.points[1, col]
    row.end <- chunk.points[2, col]
    new.loom[['matrix']][row.start:row.end, ] <- as.matrix(x = data[row.start:row.end, ])
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
    }
  }
  new.loom$matrix <- new.loom[['matrix']]
  new.loom$shape <- rev(x = new.loom[['matrix']]$dims)
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

#' Connect to a loom file
#'
#' @param filename The loom file to connect to
#' @param mode How do we connect to it? Pass 'r' for read-only or 'r+' for read/write.
#' If \code{mode} is 'r+', loomR will automatically add missing required groups during validation
#' @param skip.validate Skip the validation steps, use only for extremely large loom files
#'
#' @return A loom file connection
#'
#' @seealso \code{\link{loom-class}}
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
#' @param m Rows (genes) to subset, defaults to all rows
#' @param n Columns (cells) to subset, defaults to all columns
#' @param filename Filename for new loom object, defaults to ...
#' @param overwrite Overwrite \code{filename} if already exists?
#' @param display.progress Display progress as we're copying over data
#' @param ... Ignored for now
#'
#' @return A loom object connected to \code{filename}
#'
#' @importFrom utils setTxtProgressBar
#'
#' @export subset.loom
#' @method subset loom
#'
subset.loom <- function(
  x,
  m = NULL,
  n = NULL,
  filename = NULL,
  overwrite = FALSE,
  display.progress = TRUE,
  ...
) {
  # Set some defaults
  if (is.null(x = m)) {
    m <- 1:x$shape[1]
  }
  if (is.null(x = n)) {
    n <- 1:x$shape[2]
  }
  if (is.null(x = filename)) {
    filename <- paste(
      unlist(x = strsplit(x = x$filename, split = '.', fixed = TRUE)),
      collapse = '_subset.'
    )
  }
  if (length(x = m) == 1) {
    m <- 1:m
  }
  if (length(x = n) == 1) {
    n <- 1:n
  }
  # Ensure that m and n are within the bounds of the loom file
  if (max(m) > x$shape[1] || max(n) > x$shape[2]) {
    stop(paste(
      "'m' and 'n' must be less than",
      x$shape[1],
      "and",
      x$shape[2],
      "respectively"
    ))
  }
  extension <- rev(x = unlist(x = strsplit(x = filename, split = '.', fixed = TRUE)))[1]
  # Ensure that we call our new file a loom file
  if (extension != 'loom') {
    filename <- paste0(filename, '.loom')
  }
  if (display.progress) {
    catn("Writing new loom file to", filename)
  }
  # Make the loom file
  new.loom <- create(
    filename = filename,
    # data = t(x = x[['matrix']][n, m]),
    data = x[['matrix']][n, m],
    overwrite = overwrite
  )
  # Add row attributes
  row.attrs <- list.datasets(object = x, path = 'row_attrs', full.names = TRUE)
  if (length(x = row.attrs) > 0) {
    if (display.progress) {
      catn("\nAdding", length(x = row.attrs), "row attributes")
      pb <- new.pb()
      counter <- 0
    }
    for (row.attr in row.attrs) {
      row.list <- list(x[[row.attr]][m])
      names(x = row.list) <- basename(path = row.attr)
      new.loom$add.row.attribute(attribute = row.list)
      if (display.progress) {
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter / length(x = row.attrs))
      }
    }
  } else {
    warning("No row attributes found")
  }
  # Add col attributes
  col.attrs <- list.datasets(object = x, path = 'col_attrs', full.names = TRUE)
  if (length(x = col.attrs) > 0) {
    if (display.progress) {
      catn("\nAdding", length(x = col.attrs), "column attributes")
      pb <- new.pb()
      counter <- 0
    }
    for (col.attr in col.attrs) {
      col.list <- list(x[[col.attr]][n])
      names(x = col.list) <- basename(path = col.attr)
      new.loom$add.col.attribute(attribute = col.list)
      if (display.progress) {
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter / length(x = col.attrs))
      }
    }
  } else {
    warning("No column attributes found")
  }
  # Add layers
  layers <- list.datasets(object = x, path = 'layers', full.names = TRUE)
  if (length(x = layers) > 0) {
    if (display.progress) {
      catn("\nAdding", length(x = layers), "layers")
      pb <- new.pb()
      counter <- 0
    }
    for (layer in layers) {
      layer.list <- list(x[[layer]][n, m])
      names(x = layer.list) <- basename(path = layer)
      new.loom$add.layer(layers = layer.list)
      if (display.progress) {
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter / length(x = layers))
      }
    }
  } else {
    warning("No layers found")
  }
  new.loom$flush()
  return(new.loom)
}

#' Combine loom files
#'
#' @param looms A vector of loom files or filenames
#' @param filename Name for resultant vector
#' @param chunk.size How many rows form each input loom should we stream to the merged loom file at any given time?
#' @param order.by Optional row attribute to order each input loom by, must be one dimensional
#' @param overwrite Overwrite \code{filename} if already exists?
#' @param display.progress Display progress as we're copying over data
#'
#' @return A loom object connected to \code{filename}
#'
#' @importFrom utils setTxtProgressBar
#'
#' @seealso \code{\link{loom-class}}
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
  ...
) {
  # Basic checking of input arguments
  looms <- looms[vapply(
    X = looms,
    FUN = inherits,
    FUN.VALUE = logical(length = 1L),
    what = c('loom', 'character')
  )]
  if (length(x = looms) < 2) {
    stop("Need at least two loom objects or files to merge")
  }
  # Check the existance of loom files
  loom.names <- looms[is.character(x = looms)]
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
  # Every loom must have same number of genes (rows, MARGIN = 2)
  # and same datasets in the groups
  row.attrs <- vector(mode = 'list', length = length(x = looms))
  row.types <- list()
  col.attrs <- vector(mode = 'list', length = length(x = looms))
  col.dims <- list()
  col.ndims <- list()
  col.types <- list()
  layers <- vector(mode = 'list', length = length(x = looms))
  layers.types <- list()
  nrows <- vector(mode = 'list', length = length(x = looms))
  ncols <- vector(mode = 'list', length = length(x = looms))
  matrix.type <- list()
  if (display.progress) {
    catn("Validating", length(x = looms), "input loom files")
    pb <- new.pb()
  }
  for (i in 1:length(x = looms)) {
    this <- if (is.character(x = looms[i])) {
      connect(filename = looms[i])
    } else {
      looms[i]
    }
    row.attrs[[i]] <- sort(x = list.datasets(
      object = this,
      path = 'row_attrs',
      full.names = TRUE
    ))
    for (attr in row.attrs[[i]]) {
      if (length(x = attr) > 0) {
        row.types[[attr]] <- c(
          row.types[[attr]],
          class(x = this[[attr]]$get_type())[1]
        )
      }
    }
    col.attrs[[i]] <- sort(x = list.datasets(
      object = this,
      path = 'col_attrs',
      full.names = TRUE
    ))
    for (attr in col.attrs[[i]]) {
      if (length(x = attr) > 0) {
        col.types[[attr]] <- c(
          col.types[[attr]],
          class(x = this[[attr]]$get_type())[1]
        )
        col.dims[[attr]] <- sum(col.dims[[attr]], this[[attr]]$dims[1])
        col.ndims[[attr]] <- c(
          col.ndims[[attr]],
          length(x = this[[attr]]$dims)
        )
      }
    }
    layers[[i]] <- sort(x = list.datasets(
      object = this,
      path = 'layers',
      full.names = TRUE
    ))
    for (lay in layers) {
      if (length(x = lay) > 0) {
        layers.types[[lay]] <- c(
          layers.types[[lay]],
          class(x = this[[lay]]$get_type())[1]
        )
      }
    }
    nrows[[i]] <- this[['matrix']]$dims[2]
    ncols[[i]] <- this[['matrix']]$dims[1]
    matrix.type[[i]] <- class(x = this[['matrix']]$get_type())[1]
    if (is.character(x = looms[i])) {
      this$close_all()
    }
    if (display.progress) {
      setTxtProgressBar(pb = pb, value = i / length(x = looms))
    }
  }
  row.attrs <- unique(x = row.attrs)
  col.attrs <- unique(x = col.attrs)
  layers <- unique(x = layers)
  nrows <- unlist(x = unique(x = nrows))
  ncols <- unlist(x = ncols)
  ncells <- sum(ncols)
  if (length(x = row.attrs) != 1) {
    stop("Not all loom objects have the same row attributes")
  }
  if (length(x = col.attrs) != 1) {
    stop("Not all loom objects have the same column attributes")
  }
  if (length(x = layers) != 1) {
    stop("Not all loom objects have the same layers")
  }
  if (length(x = nrows) != 1) {
    stop("Not all loom objects have the number of rows (MARGIN = 2)")
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
  }
  # Create the new HDF5 file and the required groups
  new.loom <- loom$new(filename = filename, mode = mode)
  new.loom$create_group(name = 'layers')
  new.loom$create_group(name = "row_attrs")
  new.loom$create_group(name = "col_attrs")
  new.loom$create_group(name = "row_edges")
  new.loom$create_group(name = "col_edges")
  # Initialize the '/matrix' dataset as well as datasets in 'col_attrs' and 'layers'
  new.loom$create_dataset(
    name = 'matrix', # Name is '/matrix'
    dtype = getDtype2(x = matrix.type), # Use the single type that we got from above
    dims = c(ncells, nrows) # Use the number of cells from the sum of ncols above, nrows should be the same for everyone
  )
  for (attr in col.attrs) {
    if (length(x = attr) > 0) {
      dims.use <- switch(
        EXPR = col.ndims[[attr]],
        '1' = ncells,
        '2' = c(col.dims[[attr]], ncells),
        stop("loomR supports only one- and two-dimmensional attribute datasets")
      )
      new.loom$create_dataset(
        name = attr,
        dtype = getDtype2(x = col.types[[attr]]),
        dims = dims.use
      )
    }
  }
  for (lay in layers) {
    if (length(x = lay) > 1) {
      new.loom$create_dataset(
        name = lay,
        dtype = getDtype2(x = layers.types[[lay]]),
        dims = c(ncells, nrows)
      )
    }
  }
  # Start adding loom objects
  matrix.previous <- 0
  col.previous <- vector(mode = 'integer', length = length(x = col.attrs))
  names(x = col.previous) <- unlist(x = col.attrs)
  for (i in 1:length(x = looms)) {
    if (display.progress) {
      catn("\nAdding loom file", i ,"of", length(x = looms))
    }
    # Open the loom file
    this <- if (is.character(x = looms[i])) {
      connect(filename = looms[i])
    } else {
      looms[i]
    }
    # Get the chunk points
    chunk.points <- chunkPoints(data.size = ncols[[i]], chunk.size = chunk.size)
    # print(chunk.points)
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
      pb <- new.pb()
    }
    for (col in 1:ncol(x = chunk.points)) {
      cells.use <- chunk.points[1, col]:chunk.points[2, col]
      # Add /matrix for this chunk
      matrix.add <- this[['matrix']][cells.use, ]
      new.loom[['matrix']][cells.use + matrix.previous, ] <- matrix.add[, order.genes]
      # Add layers for this chunk
      for (lay in list.datasets(object = this[['layers']])) {
        lay.add <- this[['layers']][[lay]][cells.use, ]
        new.loom[['layers']][[lay]][cells.use + matrix.previous, ] <- lay.add[, order.genes]
      }
      if (display.progress) {
        setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
      }
      gc()
    }
    matrix.previous <- matrix.previous + max(chunk.points)
    # Add col_attrs for this chunk
    if (display.progress) {
      catn("\nAdding data to /col_attrs")
      pb <- new.pb()
    }
    for (attr in list.datasets(object = this[['col_attrs']], full.names = TRUE)) {
      start <- col.previous[attr]
      end <- start + this[[attr]]$dims[length(x = this[[attr]]$dims)]
      if (col.ndims[[attr]] == 1) {
        new.loom[[attr]][(start + 1):end] <- this[[attr]][]
      } else {
        for (col in 1:ncol(x = chunk.points)) {
          col.use <- chunk.points[1, col]:chunk.points[2, col]
          new.loom[[attr]][col.use, (start + 1):end] <- this[[attr]][col.use, ]
        }
      }
      col.previous[attr] <- end
      if (display.progress) {
        setTxtProgressBar(
          pb = pb,
          value = grep(
            pattern = attr,
            x = list.datasets(object = this[['col_attrs']], full.names = TRUE)
          ) / length(x = list.datasets(object = this[['col_attrs']], full.names = TRUE))
        )
      }
    }
    # Copy row attributes from the first loom object into the merged one
    if (i == 1) {
      if (display.progress) {
        catn("\nCopying data for /row_attrs")
        pb <- new.pb()
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
              pattern = attr,
              x = list.datasets(object = this[['row_attrs']])
            ) / length(x = list.datasets(object = this[['row_attrs']]))
          )
        }
      }
    }
    new.loom$flush()
    # Close current loom file, if not open previously
    if (is.character(x = looms[i])) {
      this$close_all()
    }
  }
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
