# setwd("put/your/local/repo/location/here") ## comment out of submitted
a <- scan("https://www.gutenberg.org/files/4300/4300-0.txt", what="character", skip=73, nlines=32858-73, fileEncoding="UTF-8")
a <- gsub("_(","", a, fixed=TRUE) ## remove "_("

split_punct <- function(vec, char) {
  vec
}