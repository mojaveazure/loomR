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
