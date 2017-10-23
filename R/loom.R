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
  #i'm not sure what we should store as slots, and what we should store as attributes or groups
  slots = c(
    version = 'ANY',
    filename = 'ANY',
    shape = "vector"
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
    #.Object@version <- packageVersion(pkg = 'loom')
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

connect = function(filename, mode = "r+") {
  self <- new("loom", filename, mode)
  self@filename <- filename
  self@shape = dim(self["matrix"])
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
  if (n_func==1) f_list=list(f_list)
  if (MARGIN==1) {
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
  
  if (MARGIN==2) {
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
