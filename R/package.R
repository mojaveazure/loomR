#' @importFrom methods setOldClass
NULL

#' An R interface for loom files
#'
#' @docType package
#' @name loomR
#'
NULL

# Hooks to set loom as an S4 class upon
# loadNamespace or library/require
.onLoad <- function(libname, pkgname) {
  setOldClass(Classes = 'loom')
}

# # Examples setting S4 methods for R6 classes
# # Requires the setOldClass from above
# #' @export foo
# methods::setGeneric(
#   name = 'foo',
#   signature = 'x',
#   def = function(x) {
#     return(standardGeneric(f = 'foo'))
#   }
# )
#
# #' @exportMethod foo
# methods::setMethod(
#   f = 'foo',
#   signature = c('x' = 'loom'),
#   definition = function(x) {
#     print("Hello loom!")
#     print(' -From foo')
#   }
# )
