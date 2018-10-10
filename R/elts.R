#' get elements from a consistently separated character string
#' 
#' Another function that I keep rewriting; no more!
#' 
#' @param x       the string
#' @param y       the separator
#' @param z       the element[s] to get
#' @param naify   return NA_character_ if y is not present? (FALSE) 
#' 
#' @return        a string, joined by the separator if z is a vector
#' 
#' @examples
#' 
#' string <- "foo_123_bar.baz.qux"
#'
#' elts(string)               # [1] "foo" 
#' elts(string, z=2:3)        # [1] "123_bar.baz.qux"
#' elts(string, y=".", z=2:3) # [1] "baz.qux"
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
elts <- function(x, y="_", z=1, naify=FALSE) {
  if (length(x) > 1) return(vapply(x, elts, y=y, z=z, naify=naify, ""))
  if (naify & !grepl(y, x)) return(NA_character_) 
  paste(sapply(strsplit(x, y, fixed=TRUE), `[`, z), collapse=y)
}
