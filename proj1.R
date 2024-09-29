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

# Function generating the m(=1000) most commonly occurring words for a vector

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
      most_common_indices <- append(most_common_indices, ind)

    }
  }
  
  # Store the words corresponding to these indices (top 1000 words)
  b <- c()
  
  for(ind in most_common_indices){
    b <- append(b, unique[ind])

  }
  
  return(b)
}

length(most_common_words(a))

b <- most_common_words(a)

create_M <- function(a, b, mlag = 4) {
  match_ind <- match(a, b)
  
  n <- length(a)
  
  M <- matrix(NA, nrow = n - mlag, ncol = mlag + 1)
  
  for(i in 1:(mlag + 1)) {
    M[, i] <- match_ind[1:(n - mlag) + i - 1]
  }
  
  return(M)
}

M <- create_M(tolower(a), b, 4)

# 8.
simulate_text <- function(token_matrix, words, nw = 50, mlag = 4) {
  output <- integer(nw)
  
  first_column <- token_matrix[,1]
  valid_tokens <- first_column[!is.na(first_column)]
  output[1] <- sample(valid_tokens, 1)
  
  for (i in 2:nw) {
    for (j in (mlag):1) if (i>j) {
      context <- output[(i-j):(i-1)]
      
      matched_rows <- which(apply(token_matrix[,1:j, drop=FALSE], 1, function(x) all(x == context)))
      
      if (length(matched_rows) > 0) {
        next_tokens <- token_matrix[matched_rows, j+1]
        token_table <- table(next_tokens)
        
        token_table <- token_table[!is.na(names(token_table))]
        probs <- as.vector(token_table) / sum(token_table)
        
        if (length(token_table) > 0) {
          output[i] <- as.integer(sample(names(token_table), 1, prob = probs))
          break 
        }
      }
    }
    
    if (is.na(output[i])) {
      output[i] <- sample(valid_tokens, 1)
    }
  }
  
  return(words[output])
}

generated_text <- simulate_text(M, b, mlag = 4)
cat(generated_text, sep=" ")

# Function generating a section of words of a given size(=50), by independently
# drawing each word from a fixed collection of words with their associated weights

frequency_simulation <- function(words, weights, size = 50){
  n <- length(weights)
  
  j <- sample(1:n, size, replace = TRUE, prob = weights)
  
  section <- words[j]
  
  return(paste(section, collapse = " "))
}

# calculate the weights of each word
freq <- tabulate(match(tolower(a), b))

frequency_simulation(tolower(a), freq, size = 50)

# find the modified b vector
low <- tolower(a)

# loop over b
b_mod <- lapply(b, function(x) {
  # find all the matches in the lowercase text and get the original occurances of them
  matches <- a[which(low == x)]
  
  # sort this and find the most common occurance
  most_common <- sort(table(matches), decreasing = TRUE)[1]
  
  names(most_common)
}) |>
  unlist()

# add whitespace in front of non punctuation
b_mod[!b_mod %in% punct_to_remove] <- paste0(' ', b_mod[!b_mod %in% punct_to_remove])

# example of generating text with modified b
generated_text <- simulate_text(M, b_mod, mlag = 4)

# need to remove the leading whitespace from first letter
generated_text[1] <- gsub(' ', '', generated_text[1])

cat(generated_text, sep = '')
