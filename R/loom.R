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
#' @field shape Shape of \code{/matrix} in columns (cells) by rows (genes)
#' @field chunksize Chunks set for this dataset in columns (cells) by rows (genes)
#' @field matrix The main data matrix, stored as columns (cells) by rows (genes)
#' @field layers Additional data matricies, the same shape as \code{/matrix}
#' @field col.attrs Extra information about cells
#' @field row.attrs Extra information about genes
#'
#' @section Methods:
#' \describe{
#'   \item{\code{add.layer(layer, overwrite)}}{Add a data layer to this loom file, must be the same dimensions as \code{/matrix}}
#'   \item{\code{add.attribute(attribute, MARGIN, overwrite)}}{
#'     Add extra information to this loom file where
#'     \code{attribute} is a named list where each element is a vector that is as long as one dimension of \code{/matrix},
#'     \code{MARGIN} is either 1 for genes or 2 for cells, and
#'     \code{overwrite} tells us whether we can overwrite existing attributes or not
#'   }
#'   \item{\code{add.row.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 1)}}
#'   \item{\code{add.col.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 2)}}
#'   \item{\code{add.meta.data(meta.data)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 2)}}
#'   \item{\code{get.attribute.df(attribute.layer, attribute.names, row.names, col.names)}}{
#'     Extract a data.frame of \code{attribute.names} from an \code{attribute.layer} ("row" - row_attrs or "col" - col_attrs).
#'     Returns a data.frame into memory with \code{attribute.names} as the columns.
#'     Removes rows that are entirely composed of NA values.
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
    add.layer = function(layers, overwrite = FALSE) {
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
          layers[[i]] <- t(x = layers[[i]])
        }
        if (names(x = layers)[i] %in% list.datasets(object = self[['layers']])) {
          if (overwrite) {
            self[['layers']]$link_delete(name = names(x = layers)[i])
          } else {
            stop(paste(
              "A layer with the name",
              names(x = layers)[i],
              "already!"
            ))
          }
        }
        self[['layers']]$create_dataset(
          name = names(x = layers)[i],
          robj = layers[[i]],
          chunk_dims = self$chunksize
        )
      }
      self$flush()
      gc(verbose = FALSE)
      private$load_layers()
      invisible(x = self)
    },
    add.attribute = function(attribute, MARGIN, overwrite = FALSE) {
      if (self$mode == 'r') {
        stop("Cannot add attributes in read-only mode")
      }
      # Value checking
      if (is.data.frame(x = attribute)) {
        attribute <- as.list(x = attribute)
      }
      if (!is.list(x = attribute) || is.null(x = names(x = attribute))) {
        stop("'attribute' must be a named list")
      }
      for (i in 1:length(x = attribute)) {
        if (!is.vector(x = attribute[[i]]) && !is.factor(x = attribute[[i]])) {
          if (length(x = dim(x = attribute[[i]])) > 1) {
            print(length(x = attribute[[i]]))
            stop("All attributes must be one-dimensional vectors")
          } else {
            attribute[[i]] <- as.vector(x = attribute[[i]])
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
      grp.name <- c('row_attrs', 'col_attrs')[MARGIN]
      grp <- self[[grp.name]]
      for (i in 1:length(x = attribute)) {
        if (length(attribute[[i]]) != rev(x = self$shape)[MARGIN])
          stop(paste(
            "All",
            switch(EXPR = MARGIN, '1' = 'gene', '2' = 'cell'),
            "attributes must be of length",
            self$shape[MARGIN]
          ))
        if (names(x = attribute)[i] %in% list.datasets(object = grp)) {
          if (overwrite) {
            grp$link_delete(name = names(x = attribute)[i])
          } else {
            stop(paste(
              "An attribute with the name",
              names(x = attribute)[i],
              "already exists!"
            ))
          }
        }
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
      if(attribute.layer == "row_attrs") {
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
        stop("Cannot write to disk in read-only mode")
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
            cat("Running trial to determine class of NAs\n")
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
      # Determine the shape of our results
      index.use <- private$iter_range(index.use = index.use)
      # Create our results holding object
      if (results.matrix) {
        switch(
          EXPR = private$iter.margin,
          '1' = results <- matrix( # Genes, make matrix with nCells rows and range(index.use) columns
            nrow = self$shape[1],
            ncol = length(x = index.use[1]:index.use[2])
          ),
          '2' = results <- matrix( # Cells, make matrix with range(index.use) rows and nGenes columns
            nrow = length(x = index.use[1]:index.use[2]),
            ncol = self$shape[2]
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
            '1' = self[[dataset.use]][, chunk.indices], # Chunk genes
            '2' = self[[dataset.use]][chunk.indices, ] # Chunk cells
          )
        } else {
          self[[dataset.use]][chunk.indices]
        }
        chunk.data <- FUN(chunk.data, ...)
        if (results.matrix) {
          if (MARGIN == 1) {
            results[, chunk.indices] <- chunk.data
          } else if (MARGIN == 2) {
            results[chunk.indices, ] <- chunk.data
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
      # matrix.data is a vector of data for one cell or a list of data for several cells
      # each entry in matrix.data must be the same length as number of genes
      # attributes.data is an optional list or vector (with optional names) for col_attrs
      # each entry in col_attrs must be the same length as the number of cells being added (NAs added for those that aren't)
      # layers.data is an optional list (with optional names) for layers
      # each entry in layers.data must be an N by M matrix where N is the number of genes and M is the number of cells
      if (is.vector(x = matrix.data) && !is.list(x = matrix.data)) {
        matrix.data <- list(matrix.data)
      }
      list.check <- vapply(
        X = list(matrix.data, attributes.data, layers.data),
        FUN = is.list,
        FUN.VALUE = logical(length = 1L)
      )
      if (!all(list.check)) {
        stop("'matrix.data', 'attributes.data', and 'layers.data' must be lists")
      }
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
      num.cells.added <- max(lengths)
      return(num.cells.added)
      # private$load_layers()
      # private$load_attributes(MARGIN = 1)
      # private$load_attributes(MARGIN = 2)
      return(lengths)
    }
    # add.loom = function() {}
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
        index.use <- c(1, rev(x = self$shape)[private$iter.margin])
      } else if (length(x = index.use) == 1) {
        # If one index was provided, start at one and go to index
        index.use <- c(1, index.use)
      } else {
        # Use index.use[1] and index.use[2]
        index.use <- c(index.use[1], index.use[2])
      }
      # Ensure the indices provided fit within the range of the dataset
      index.use[1] <- max(1, index.use[1])
      index.use[2] <- min(index.use[2], rev(x = self$shape)[private$iter.margin])
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
#' @param overwrite Overwrite an already existing loom file?
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
  chunk.dims = 'auto',
  overwrite = FALSE
) {
  mode <- ifelse(test = overwrite, yes = 'w', no = 'w-')
  if (file.exists(filename) && !overwrite) {
    stop(paste('File', filename, 'already exists!'))
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
  new.loom <- loom$new(filename = filename, mode = mode)
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

#' Connect to a loom file
#'
#' @param filename The loom file to connect to
#' @param mode How do we connect to it? Pass 'r' for read-only or 'r+' for read/write
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
#'
#' @return A loom object connected to \code{filename}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
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
  new.pb <- function() {return(txtProgressBar(style = 3, char = '='))}
  # Set some defaults
  if (is.null(x = m)) {
    m <- 1:x$shape[1]
  }
  if (is.null(x = n)) {
    n <- 1:x$shape[2]
  }
  if (is.null(x = filename)) {
    filename <- paste(
      unlist(x = strsplit(x = lfile$filename, split = '.', fixed = TRUE)),
      collapse = '_subset.'
    )
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
  extension <- rev(x = unlist(x = strsplit(x = filename, split = '.', fixed = TRUE)))
  # Ensure that we call our new file a loom file
  if (extension != 'loom') {
    filename <- paste0(filename, '.loom')
  }
  if (display.progress) {
    cat("Writing new loom file to", filename, '\n')
  }
  # Make the loom file
  new.loom <- create(
    filename = filename,
    data = t(x = x[['matrix']][n, m]),
    overwrite = overwrite
  )
  # Add row attributes
  row.attrs <- list.datasets(object = x, path = 'row_attrs', full.names = TRUE)
  if (length(x = row.attrs) > 0) {
    if (display.progress) {
      cat("\nAdding", length(x = row.attrs), "row attributes\n")
      pb <- new.pb()
      counter <- 0
    }
    for (row.attr in row.attrs) {
      base.row <- basename(path = row.attr)
      new.loom$add.row.attribute(attribute = base.row = x[[row.attr]][m]))
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
      cat("\nAdding", length(x = col.attrs), "row attributes\n")
      pb <- new.pb()
      counter <- 0
    }
    for (col.attr in col.attrs) {
      base.col <- basename(path = col.attr)
      new.loom$add.col.attribute(attribute = list(base.col = x[[col.attr]][n]))
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
      cat("\nAdding", length(x = layers), "row attributes\n")
      pb <- new.pb()
      counter <- 0
    }
    for (layer in layers) {
      base.layer <- basename(path = layer)
      new.loom$add.layer(layers = list(base.layer = x[[layer]][n, m]))
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
