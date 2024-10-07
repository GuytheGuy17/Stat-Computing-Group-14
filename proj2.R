## File for second project
setwd("~/Stat-Computing-Group-14")

# Only using the first 150 rows of the data for this practical.
data <- read.table("engcov.txt", nrows = 150)

x <- rlnorm(80, meanlog=3.152, sdlog=0.451)
norm_x <- x/sum(x)

y <- sample(1:80, 100, replace = TRUE, prob = norm_x) # Sampled length of time from infection to death

total_deaths <- sum(data$nhs)

death_days <- rep(data$julian, data$nhs) # Day of death vector

durations <- sample(1:80, total_deaths, replace = TRUE, prob = norm_x)
                    
t0 <- death_days - durations
