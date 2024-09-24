# setwd("put/your/local/repo/location/here") ## comment out of submitted
a <- scan("https://www.gutenberg.org/files/4300/4300-0.txt", what="character", skip=73, nlines=32858-73, fileEncoding="UTF-8")
a <- gsub("_(","", a, fixed=TRUE) ## remove "_("

split_punct <- function(vec, char) {
  esc_char <- paste0('\\', char)
  
  ii <- grep(esc_char, vec)
  
  out <- vector(mode = 'character', length = length(ii) + length(vec))
  
  out[ii + 1:length(ii)] <- char
  
  out[-c(ii + 1:length(ii))] <- gsub(esc_char, '', vec)
  
  return(out)
}

punct_to_remove <- c(",", ".", ";", "!", ":", "?")

for(punct in punct_to_remove) {
  a <- split_punct(a, punct)
} 
