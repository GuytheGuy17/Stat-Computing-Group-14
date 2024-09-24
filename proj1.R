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

# Function generating the m (=1000) most commonly occurring words for a vector
most_common_words <- function(vec, m = 1000) {
  # Convert words to lowercase 
  vec <- tolower(vec)
  
  # Find the vector of unique words
  unique <- unique(vec)
  
  # Determine the index vector
  index <- match(vec, unique)
  
  # Count the number of occurrences
  occurrences <- tabulate(index)
  
  # Determine the threshold frequency required to be included in the top 1000 words
  sorted_occurrences <- occurrences[order(occurrences, decreasing = TRUE)]
  threshold_frequency <- sorted_occurrences[m]
  
  # Store the indices that beat the threshold frequency 
  most_common_indices <- c()
  for(ind in c(1:length(occurrences))){
    if(occurrences[ind] >= threshold_frequency){
      most_common_indices = append(most_common_indices, ind)
    }
  }
  
  # Store the words corresponding to these indices (top 1000 words)
  b <- c()
  
  for(ind in most_common_indices){
    b = append(b, unique[ind])
  }
  
  return(b)
}

length(most_common_words(a))

