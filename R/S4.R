#' @import R6
#' @importFrom methods setOldClass
NULL

.onLoad <- function(libname, pkgname) {
  setOldClass(Classes = 'loom')
}

.onAttach <- function(libname, pkgname) {
  .onLoad(libname = libname, pkgname = pkgname)
}

# # Examples setting S4 methods for R6 classes
# # Requires the setOldClass from above
# setGeneric(
#   name = 'foo',
#   signature = 'x',
#   def = function(x) {
#     return(standardGeneric(f = 'foo'))
#   }
# )
#
# setMethod(
#   f = 'foo',
#   signature = c('x' = 'loom'),
#   definition = function(x) {
#     print("Hello loom!")
#     print(' -From foo')
#   }
# )
