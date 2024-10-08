
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

pearson_stat_vec <- function(observed_data, simulated_data) {
  denom <- simulated_data
  denom[simulated_data <= 1] <- 1
  
  sum((observed_data - simulated_data) ^ 2 / denom)
}

calc_days_to_death_probs <- function(max_duration = 80, meanlog = 3.152, sdlog = 0.451) {
  # Set up infection-to-death distribution
  probabilities <- dlnorm(1:max_duration, meanlog, sdlog)
  probabilities / sum(probabilities)  # Normalize
}

create_t0 <- function(days, deaths, max_duration = 80) {
  probabilities <- calc_days_to_death_probs(max_duration)
  
  # Calculate total deaths and create death days vector
  total_deaths <- sum(deaths)
  death_days <- rep(days, deaths)
  
  # Generate infection-to-death durations
  durations <- sample(1:max_duration, total_deaths, replace = TRUE, prob = probabilities)
  
  # Calculate initial guesses for infection days
  t0 <- death_days - durations
  
  # set negative values to 0
  t0[t0 <= 0] <- 0
  
  return(t0)
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  if(is.null(t0)) {
    t0 <- create_t0(t, deaths)
  }
  
  # extend deaths to 1 to 350 days
  ext_deaths <- numeric(length = 350)
  ext_deaths[t] <- deaths
  
  # initial vectors/matrices to store history of P and infections
  P_hist <- numeric(length = n.rep)
  inft <- matrix(nrow = length(ext_deaths), ncol = n.rep)
  
  # calculate the probs for drawing from later
  probs <- calc_days_to_death_probs()
  
  for(j in 1:n.rep) {
    # initialise options for days
    if(j <= 50) {
      days_opt <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    } else if(j <= 75) {
      days_opt <- c(-4, -2, -1, 1, 2, 4)
    } else {
      days_opt <- c(-2, -1, 1, 2)
    }
    
    days_to_death <- sample(1:80, length(t0), replace = TRUE, prob = probs)
    
    # create pred_deaths - nbins 350 to match the length
    pred_deaths <- tabulate(t0 + days_to_death, nbins = 350)
    
    # calculate the P
    P <- pearson_stat_vec(pred_deaths, ext_deaths)
    
    # sample without replacement to find the random ordering
    ordering <- sample(1:length(t0), length(t0), replace = FALSE)
    
    # sort both vectors in same ordering
    t0_sorted <- t0[ordering]
    
    days_to_death_sorted <- days_to_death[ordering]
    
    # sample random changes
    rand_change <- sample(days_opt, length(t0), replace = TRUE)
    
    for(i in 1:length(t0)) {
      # Add 2 new vectors to store the proposed changes
      t0_proposed <- t0_sorted
      pred_deaths_proposed <- pred_deaths
      
      # add random change - don't want to propose a negative number of days so replace with 0
      t0_proposed[i] <- max(t0_sorted[i] + rand_change[i], 0)
      
      # recalcuate day of death and the previous day of death
      proposed_death_day <- t0_proposed[i] + days_to_death_sorted[i]
      
      old_death_day <- t0_sorted[i] + days_to_death_sorted[i]
      
      # inject these into the predicted deaths
      pred_deaths_proposed[proposed_death_day] <- pred_deaths_proposed[proposed_death_day] + 1
      pred_deaths_proposed[old_death_day] <- pred_deaths_proposed[old_death_day] - 1
      
      P_proposed <- pearson_stat_vec(pred_deaths_proposed, ext_deaths)
      
      # if P_proposed < P then update all our vectors (strictly less than - in case of equality don't update for efficiency)
      if(P_proposed < P) {
        pred_deaths <- pred_deaths_proposed
        t0_sorted <- t0_proposed
        P <- P_proposed
      }
    }
    
    # undo the sorting of t0 and store it
    reordering <- numeric(length = length(ordering))
    reordering[ordering] <- 1:length(ordering)
    
    t0 <- t0_sorted[reordering]
    
    # create plot
    plot(x = 1:350, y = ext_deaths, type = 'l', xlab = 'Days', ylab = 'Count')
    lines(x = 1:350, y = tabulate(t0, nbins = 350), col = 'blue')
    lines(x = 1:350, y = pred_deaths, col = 'red')
    legend(x = 175, y = 900, legend = c('Estimated Incidence', 'Simulated Deaths', 'Real Deaths'), fill = c('blue', 'red', 'black'))
    
    # save current state of variables
    P_hist[j] <- P
    inft[, j] <- pred_deaths
  }
  
  # produce the plot here
  
  list(
    P = P_hist,
    inft = inft,
    t0 = t0
  )
}

# just using n.rep = 50 as example
run <- deconv(initial_data$julian, initial_data$nhs, n.rep = 200)

# create plot - need to add confidence intervals into this once bootstrapping
cases <- tabulate(run$t0, nbins = 350)
plot(x = 1:350, y = cases, type = 'l', xlab = 'Days', ylab = 'Count')
lines(x = initial_data$julian, y = initial_data$nhs, col = 'blue')
abline(v = 84, lty = 'dashed')
text(x = 84, y = max(cases), labels = 'UK Lockdown', pos = 4)
legend(x = 350, y = max(cases), legend = c('Estimated Incidences', 'Deaths'), col = c('black', 'blue'), lty = 'solid', xjust = 1, yjust = 1, cex = 0.8)

## cases coming down before lockdown

## comparison of 2 pearson stat functions 
days_to_death <- sample(1:80, length(t0), replace = TRUE, prob = calc_days_to_death_probs())
t0 <- create_t0(initial_data$julian, initial_data$nhs)
pred_deaths <- tabulate(t0 + days_to_death, nbins = 350)

ext_deaths <- numeric(length = 350)
ext_deaths[initial_data$julian] <- initial_data$nhs

t <- Sys.time()
test <- lapply(1:100000, \(x) pearson_stat_vec(pred_deaths, ext_deaths))
print(Sys.time() - t)

t <- Sys.time()
test2 <- lapply(1:100000, \(x) pearson_stat(pred_deaths, ext_deaths))
print(Sys.time() - t)
