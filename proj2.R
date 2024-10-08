## File for second project
setwd("~/Stat-Computing-Group-14")

# Only using the first 150 rows of the data for this practical.
initial_data <- read.table("engcov.txt", nrows = 150)

# Define the modified Pearson statistic to use as a metric for the goodness of fit
# observed_data, simulated_data represent the vectors we are computing the Pearson statistic for 
pearson_stat <- function(observed_data, simulated_data) {
  sum <- 0
  
  for (i in 1:length(observed_data)){
    sum <- sum + ((observed_data[i] - simulated_data[i]) ^ 2) / (max(1, simulated_data[i]))
  }
  
  return(sum)
}

basic_algorithm <- function(data, num_rows = 150, meanlog = 3.152, sdlog = 0.451, max_duration = 80) {
  
  # Set up infection-to-death distribution
  probabilities <- dlnorm(1:max_duration, meanlog, sdlog)
  probabilities <- probabilities / sum(probabilities)  # Normalize
  
  # Calculate total deaths and create death days vector
  total_deaths <- sum(data$nhs)
  death_days <- rep(data$julian, data$nhs)
  
  # Generate infection-to-death durations
  durations <- sample(1:max_duration, total_deaths, replace = TRUE, prob = probabilities)
  
  # Calculate initial guesses for infection days
  t0 <- death_days - durations
  
  return(list(
    data = data,
    probabilities = probabilities,
    total_deaths = total_deaths,
    death_days = death_days,
    t0 = t0
  ))
}

basic_algorithm(data = initial_data)

