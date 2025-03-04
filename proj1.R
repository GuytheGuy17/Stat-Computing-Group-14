#' Group Members: 
#' Guy McClennan - s2036567
#' Alexandru Girban - s2148980
#' Louis Bennett - s2744241
#' 
#' This project was completed largely collaboratively with each member contributing 
#' ideas for a range of answers. Then it would be left to one person to implement 
#' this in code. Louis created the split_punct, create_M functions as well as the code 
#' for  Q10, Alex added the most_common_words and frequency simulation functions
#' as well as helping with how to solve Q8. Guy set up the Github repo, added the
#' text for parts 2 and 3 and created the base simulate_text function code for Q8 with 
#' the group as a whole working on ways to optimise this function.

# setwd("~/Stat-Computing-Group-14")
a <- scan("4300-0.txt", what="character", skip=73, nlines=32858-73, fileEncoding="UTF-8")
a <- gsub("_(","", a, fixed=TRUE) ## remove "_("

#' split_punct
#' 
#' @description A function to search for punctuation marks in a word vector, remove them 
#' from their associated word and input them as a new unique entry in the word vector.
#' @param vec The word vector to search for punctuation marks.
#' @param char The punctuation marks that are being searched for in the word vector.
#' @return A word vector with the punctuation marks of the original word vector
#' marked as unique inputs after the word they were originally in.

split_punct <- function(vec, char) {
  # Escape the character as grep interprets some punctuation as regex
  esc_char <- paste0('\\', char)
  
  ii <- grep(esc_char, vec)
  
  # Initialise output vector
  out <- vector(mode = 'character', length = length(ii) + length(vec))
  
  # Store the (un-escaped) character in the right place
  out[ii + 1:length(ii)] <- char
  
  # Add the rest of the tokens back in, with the character removed
  out[-c(ii + 1:length(ii))] <- gsub(esc_char, '', vec)
  
  return(out)
}

punct_to_remove <- c(",", ".", ";", "!", ":", "?")

for(punct in punct_to_remove) {
  a <- split_punct(a, punct)
} 

#' most_common_words
#' 
#' @description Function generating the m(=1000) most commonly occurring words for a vector.
#' @param vec A word vector the function is searching over.
#' @param m The number of most common words being identified.
#' @return A vector containing the m most common words in the input word vector.

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
  
  # Print the threshold frequency
  print(paste('The threshold frequency found was:', threshold_frequency))
  
  # Store the indices that beat the threshold frequency 
  most_common_indices <- which(occurrences >= threshold_frequency)
  
  # Store the words corresponding to these indices (top 1000 words)
  b <- unique[most_common_indices]
  
  return(b)
}

b <- most_common_words(a)

length(b)

#' create_M
#' 
#' @description Function creating the matrix M
#' @param full_text A word vector the function is searching over.
#' @param words A vector containing the m most common words, generated by most_common_words
#' @param mlag The order of the Markov model to use
#' @return The matrix M as defined in the question sheet

create_M <- function(full_text, words, mlag = 4) {
  # Tokenise the full text vector
  match_ind <- match(full_text, words)
  
  # Initialise the matrix M
  n <- length(full_text)
  M <- matrix(NA, nrow = n - mlag, ncol = mlag + 1)
  
  # Loop over and insert each column of the matrix, shifting the text along by one each time
  for(i in 1:(mlag + 1)) {
    M[, i] <- match_ind[1:(n - mlag) + i - 1]
  }
  
  return(M)
}

M <- create_M(tolower(a), b, 4)

head(M)

#' simulate_text
#' 
#' @description Function to simulate text of length nw(=50)
#' @param full_text A word vector the function is searching over.
#' @param words A vector containing the m most common words, generated by most_common_words
#' @param mlag The order of the Markov model to use
#' @param nw The number of tokens/words to generate in the text
#' @return Simulated text from full_text using a mlag order Markov model

simulate_text <- function(full_text, words, nw = 50, mlag = 4) {
  # Create the matrix M 
  # The 'tolower' and 'gsub' ensure the text and words are in the same form and will be matched
  token_matrix <- create_M(tolower(full_text), tolower(gsub(' ', '',words)), mlag)
  
  # Initialise output
  output <- integer(nw)
  
  # Sample the first token randomly from all non NAs
  first_column <- token_matrix[, 1]
  valid_tokens <- first_column[!is.na(first_column)]
  output[1] <- sample(valid_tokens, 1)
  
  # Loop up to nw
  for (i in 2:nw) {
    # Possible that i <= j for the first couple of words (before we've generated the first mlag words)
    for (j in (mlag):1) if (i>j) { 
      # Relevant sequence of words
      context <- output[(i - j):(i - 1)]
      
      # Initialise the test matrix to find matching sequences in M
      test <- matrix(NA, nrow = NROW(token_matrix), ncol = j + 1)
      test[, 1] <- TRUE
      
      # Loop over the word section
      for(k in 1:j) {
        # Check if the last column is TRUE and if the current word also matches
        test[, k + 1] <- !is.na(token_matrix[, k]) & token_matrix[, k] == context[k] & test[, k]
      }
      
      # Find rows which matched the whole sequence
      matched_rows <- which(test[, j + 1])
      
      # If no matched rows try lower order model
      if (length(matched_rows) > 0) {
        # Get the next tokens
        next_tokens <- token_matrix[matched_rows, j + 1]
        token_table <- table(next_tokens)
        
        token_table <- token_table[!is.na(names(token_table))]
        
        # Calculate the probabilities of each token
        probs <- as.vector(token_table) / sum(token_table)
        
        # If there's only one match or all are NA try lower order model
        if (length(token_table) > 1) {
          # Sample the next word then break from the loop
          output[i] <- as.integer(sample(names(token_table), 1, prob = probs))
          break 
        }
      }
    }
    
    # If the first order model failed, use the 0th order (sampling randomly from all words)
    if (is.na(output[i]) | output[i] == 0) {
      output[i] <- sample(valid_tokens, 1)
    }
  }
  
  return(words[output])
}

## Generating a 50 word sequence with order 4 from the given text
generated_text_1 <- simulate_text(a, b, mlag = 4)
cat(generated_text_1, sep=" ")

#' frequency_simulation
#' 
#' @description # Function generating a section of words of a given size(=50), 
#' by independently drawing each word from a fixed collection of words 
#' with their associated weights/probability of being sampled.
#' @param words The word space to generate the simulation from.
#' @param weights The associated weights of the words we are simulating from.
#' @param size The size of the text we are constructing with this simulation
#' @return A text simulated according to the procedure in the description.
 
frequency_simulation <- function(words, weights, size = 50){
  n <- length(weights)
  
  j <- sample(1:n, size, replace = TRUE, prob = weights)
  
  section <- words[j]
  
  return(section)
}

# Convert a to lowercase
low <- tolower(a)

# Calculate the weights of each word
freq <- tabulate(match(low, b))

freq_generated_text <- frequency_simulation(low, freq, size = 50)

cat(freq_generated_text, collapse = " ")

## Find the modified b vector
# Loop over b
b_mod <- sapply(b, function(x) {
  # Find all the matches in the lowercase text and get the original occurrences of them
  matches <- a[which(low == x)]
  
  # Sort this and find the most common occurrence 
  # find the most common form of the word rather than just comparing capital/non-capital cases
  most_common <- sort(table(matches), decreasing = TRUE)[1]
  
  names(most_common)
})

# Add whitespace in front of non punctuation
b_mod[!b_mod %in% punct_to_remove] <- paste0(' ', b_mod[!b_mod %in% punct_to_remove])

# Example of generating text with modified b
generated_text_2 <- simulate_text(a, b_mod, mlag = 4)

# Need to remove the leading whitespace from first letter
generated_text_2[1] <- gsub(' ', '', generated_text_2[1])

cat(generated_text_2, sep = '')