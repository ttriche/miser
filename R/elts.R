#' get elements from a consistently separated character string
#' 
#' Another function that I keep rewriting; no more!
#' 
#' @param x       the string(s) to split
#' @param sep     the separator (default is "_")
#' @param elt     the element[s] to get (default is element 1)
#' @param naify   return NA_character_ if y is not present? (FALSE) 
#' 
#' @return        a string, joined by `sep` if `elt` is a vector
#' 
#' @examples
#' 
#' string <- "foo_123_bar.baz.qux"
#'
#' elts(string)                     # [1] "foo" 
#' elts(string, elts=2:3)           # [1] "123_bar.baz.qux"
#' elts(string, sep=".", elts=2:3)  # [1] "baz.qux"
#' 
#' strings <- paste0(letters[1:3], "_string")
#'
#' elts(strings)              # a_string b_string c_string 
#'                            #      "a"      "b"      "c" 
#' 
#' strings2 <- c(strings, "nothing") 
#' 
#' elts(strings2)             #  a_string  b_string  c_string   nothing 
#'                            #       "a"       "b"       "c" "nothing" 
#' elts(strings2, naify=TRUE) #  a_string b_string c_string  nothing 
#'                            #       "a"      "b"      "c"       NA 
#' 
#' @export 
elts <- function(x, sep="_", elt=1, naify=FALSE) {
  if (length(x) > 1) return(vapply(x, elts, sep=sep, elt=elt, naify=naify, ""))
  if (naify & !grepl(sep, x)) return(NA_character_) 
  paste(sapply(strsplit(x, sep, fixed=TRUE), `[`, elt), collapse=sep)
}
