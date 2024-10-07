## File for second project
setwd("~/Stat-Computing-Group-14")

# Only using the first 150 rows of the data for this practical.
initial_data <- read.table("engcov.txt", nrows = 150)

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

