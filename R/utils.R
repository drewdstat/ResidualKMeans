getmode <- function(v){
  unique(v)[which.max(tabulate(match(v, unique(v))))]
}

complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y))
  rho * sd(x) * y + x * sd(y) * sqrt(1 - rho^2)
}

binpred<-function(x, b0, b1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  missprobpred <- 1/(1 + exp(-1 * (b0 + x * b1)))
  rbinom(length(missprobpred), 1, missprobpred)
}

decimalplaces  <-  function(x, beforedecimal = F) {
  if(beforedecimal){
    if (abs(x - round(x))  >  .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0 + $', '', as.character(x)), ".",
                     fixed = TRUE)[[1]][[1]])
    } else {
      nchar(x)
    }
  } else {
    if (abs(x - round(x))  >  .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0 + $', '', as.character(x)), ".",
                     fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
}

reorder_clusters <- function(x){
  fx <- as.factor(x)
  summfx <- summary(fx)
  newlabs <- summfx[order(-summfx)]
  fx <- factor(fx, levels = names(newlabs))
  return(as.numeric(fx))
}
