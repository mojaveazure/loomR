#' @import hdf5r
#' @importFrom R6 R6Class
NULL

#' A class for loom
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @return Object of class \code{\link{loom}}
#' @seealso \code{\link{hdf5r::H5File}}
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
    # Methods
    initialize = function(
      filename = NULL,
      mode = c('a', 'r', 'r+'),
      ...
    ) {
      do.validate <- file.exists(filename)
      super$initialize(filename = filename, mode = mode, ...)
      if (do.validate) {
        validateLoom(object = self)
        self$shape <- self[['matrix']]$dims
        chunks <- h5attr(x = self, which = 'chunks')
        chunks <- gsub(pattern = '(', replacement = '', x = chunks, fixed = TRUE)
        chunks <- gsub(pattern = ')', replacement = '', x = chunks, fixed = TRUE)
        chunks <- unlist(x = strsplit(x = chunks, split = ','))
        self$chunks <- as.integer(x = chunks)
        self$version <- as.character(x = tryCatch(
          expr = h5attr(x = self, which = 'version'),
          error = function(e) packageVersion(pkg = 'loomR')
        ))
      } else {
        # self$version <- packageVersion(pkg = 'loomR')
        print()
      }
    }
  )
)

#' Create a loom object
#'
#' @param filename ...
#' @param data ...
#' @param row.attrs ...
#' @param col.attrs ...
#'
#' @return A connection to a loom file
#'
#' @seealso \code{\link{loom-class}}
#'
create <- function(filename, data, row.attrs, col.attrs) {
  if (file.exists(filename)) {
    stop(paste('File', file, 'already exists!'))
  }
  new.loom <- loom$new(filename = filename, mode = 'r')
}

# #' @importFrom utils packageVersion
# #'
# setMethod(
#   f = 'initialize',
#   signature = 'loom',
#   definition = function(.Object, name, mode = 'a') {
#     .Object <- callNextMethod(
#       .Object,
#       name = name,
#       mode = mode
#     )
#     validateLoom(object = .Object)
#     #.Object@version <- packageVersion(pkg = 'loom')
#     # .Object@filename <- name
#     .Object@shape <- dim(.Object['/matrix'])
#     return(.Object)
#   }
# )


#' Validate a loom object
#'
#' @param object A loom object
#'
#' @return None, errors if object is an invalid loom object
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
  names(dim.matrix) <- required.groups[c(2, 1)]
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
connect <- function(filename, mode = "r+") {
  self <- new("loom", filename, mode)
  # self@filename <- filename
  self@shape <- self["matrix"]@dim
  return(self)
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
  if (n_func == 1) f_list=list(f_list)
  if (MARGIN == 1) {
    results=list();
    for (j in 1:n_func) {
      results[[j]] <- numeric(0)
    }
    rows_per_chunk <- chunksize
    ix <- 1
    while (ix <= self@shape[1]) {
      rows_per_chunk <- min(rows_per_chunk, self@shape[1]-ix+1)
      chunk <- self["matrix"][ix:(ix + rows_per_chunk -1), ]
      for(j in 1:n_func) {
        new_results <- apply(chunk, 1, FUN = f_list[[j]])
        results[[j]] <- c(results[[j]], new_results)
      }
      ix <- ix + chunksize
    }
  }
  if (MARGIN == 2) {
    results=list();
    for (j in 1:n_func) {
      results[[j]] <- numeric(0)
    }
    cols_per_chunk <- chunksize
    ix <- 1
    while (ix <= self@shape[2]) {
      cols_per_chunk <- min(cols_per_chunk, self@shape[2]-ix+1)
      chunk <- self["matrix"][,ix:(ix + cols_per_chunk -1)]
      for(j in 1:n_func) {
        new_results <- apply(chunk, 2, FUN = f_list[[j]])
        results[[j]] <- c(results[[j]], new_results)
      }
      ix <- ix + chunksize
    }
  }
  if (n_func == 1) return(results[[1]])
  return(results)
}
