#' @import hdf5r
#' @importFrom R6 R6Class
NULL

#' A class for loom files
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @return Object of \code{\link{R6Class}} to generate \code{loom} objects
#' @format An \code{\link{R6Class}} object
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
#'   \item{\code{add.attribute(attribute, MARGIN)}}{Add extra information to this loom file; \code{attribute} is a named list where each element is a vector that is as long as one dimension of \code{/matrix}, \code{MARGIN} is either 1 for cells or 2 for genes}
#'   \item{\code{add.row.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 2)}}
#'   \item{\code{add.col.attribute(attribute)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 1)}}
#'   \item{\code{add.meta.data(meta.data)}}{A wrapper for \code{add.attribute(attribute, MARGIN = 1)}}
#' }
#'
#' @importFrom utils packageVersion
#'
#' @export
#'
loom <- R6Class(
  classname = 'loom',
  inherit = hdf5r::H5File,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
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
    initialize = function(filename = NULL, mode = c('a', 'r', 'r+', 'w', 'w-'), ...) {
      # If the file exists, run validation steps
      do.validate <- file.exists(filename) && !(mode %in% c('w', 'w+'))
      super$initialize(filename = filename, mode = mode, ...)
      if (do.validate) {
        # Run the validation steps
        validateLoom(object = self)
        # Store the shape of /matrix
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
        self$version <- packageVersion(pkg = 'loomR')
      }
    },
    add.layer = function(layer, name) {
      # Layers have to be matricies
      if (!is.matrix(x = layer)) {
        layer <- as.matrix(x = layer)
      }
      if (is.null(x = self$shape)) {
        stop(private$err_msg)
      }
      do.transpose <- FALSE
      if (any(dim(x = layer) != self$shape)) {
        if (all(rev(x = dim(x = layer)) == self$shape)) {
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
      # Transpose the matrix since hdf5r uses column x row
      if (do.transpose) {
        self[['layers', name]] <- t(x = layer)
      } else {
        self[[layers, name]] <- layer
      }
      private$load_layers()
    },
    add.attribute = function(attribute, MARGIN) {
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
        if (length(attribute[i]) != self$shape[MARGIN])
          stop(paste(
            "All",
            switch(EXPR = MARGIN, '1' = 'row', '2' = 'column'),
            "attributes must be of length",
            self$shape[MARGIN]
          ))
        grp[[names(x = attribute)[i]]] <- attribute[[i]]
      }
      gc(verbose = FALSE)
      # Load the attributes for this margin
      private$load_attributes(MARGIN = MARGIN)
    },
    add.row.attribute = function(attribute) {
      self$add.attribute(attribute = attribute, MARGIN = 2)
    },
    add.col.attribute = function(attribute) {
      self$add.attribute(attribute = attribute, MARGIN = 1)
    },
    add.meta.data = function(meta.data) {
      self$add.col.attribute(attribute = meta.data)
    }
  ),
  private = list(
    err_msg = "This loom object has not been created with either loomR::create or loomR::connect, please use these function to create or connect to a loom file",
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
  if (length(x = chunk.dims) > 2 || length(x = chunk.dims < 1)) {
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
    robj = t(x = data),
    chunk_dims = chunk.dims
  )
  if (!is.null(x = colnames(x = data))) {
    new.loom$add.row.attribute(attribute = list('gene_names' = colnames(x = data)))
  }
  if (!is.null(x = rownames(x = data))) {
    new.loom$add.col.attribute(attribute = list('cell_names' = colnames(x = data)))
  }
  # Store some constants as HDF5 attributes
  h5attr(x = new.loom, which = 'version') <- new.loom$version
  h5attr(x = new.loom, which = 'chunks') <- paste0(
    '(',
    paste(new.loom[['matrix']]$chunk_dims, collapse = ', '),
    ')'
  )
  # Groups
  new.loom$create_group(name = 'layers')
  new.loom$create_group(name = 'row_attrs')
  new.loom$create_group(name = 'col_attrs')
  # Add layers
  for (ly in layers) {
    new.loom$add.layer(layer = ly)
  }
  new.loom$add.row.attribute(attribute = gene.attrs)
  new.loom$add.col.attribute(attribute = cell.attrs)
  # Set last bit of information
  new.loom$shape <- new.loom[['matrix']]$dims
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
#' @return None, errors if object is an invalid loom connection
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
  dim.matrix <- object[[root.datasets]]$dims # Rows x Columns
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
#' @export
#'
connect <- function(filename, mode = "r") {
  if (!(mode %in% c('r', 'r+'))) {
    stop("'mode' must be one of 'r' or 'r+'")
  }
  new.loom <- loom$new(filename = filename, mode = mode)
  return(new.loom)
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
