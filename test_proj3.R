library(lme4)

## test script to compare multiple different model structures with the output from lme4::lmer

# generate data w random effects
generate_random_data <- function(num_rows, num_fixed, num_random) {
  fixed_effects <- as.data.frame(matrix(rnorm(num_rows * num_fixed), ncol = num_fixed))
  names(fixed_effects) <- paste0("X", 1:num_fixed)
  
  random_effects <- lapply(1:num_random, \(i) factor(sample(1:3, num_rows, replace = TRUE)))
  names(random_effects) <- paste0("R", 1:num_random)
  
  data <- cbind(fixed_effects, random_effects)
  
  data$y <- rnorm(num_rows)
  
  return(data)
}

# compare the 2 models
test_lmm <- function(form, dat, ref) {
  cust_lmm <- lmm(form, dat, ref)
  
  if(length(ref) > 0) {
    rand_effects <- sapply(ref, \(x) paste(x, collapse = ':'))
    
    lmer_form <- as.formula(paste0(
      deparse(form), " + (1|", paste(rand_effects, collapse = ") + (1|"), ")"
    ))
    
    lmer_lmm <- lme4::lmer(lmer_form, dat, REML = FALSE)
  } else {
    lmer_lmm <- lm(form, dat)
  }
  
  return(list(cust = cust_lmm, lmer = lmer_lmm))
}

# simple model, no random effects
data1 <- generate_random_data(100, 2, 0)
formula1 <- y ~ X1 + X2
ref_list1 <- list()
test_lmm(formula1, data1, ref_list1)

# simple model with one random effect
data2 <- generate_random_data(100, 2, 1)
formula2 <- y ~ X1 + X2
ref_list2 <- list("R1")
test_lmm(formula2, data2, ref_list2)

# model with two random effects (nested structure)
data3 <- generate_random_data(100, 2, 2)
formula3 <- y ~ X1 + X2
ref_list3 <- list("R1", "R2")
test_lmm(formula3, data3, ref_list3)

# model with interaction term and one random effect
data4 <- generate_random_data(100, 3, 1)
formula4 <- y ~ X1 + X2 * X3
ref_list4 <- list("R1")
test_lmm(formula4, data4, ref_list4)

# invalid random effect w variable missing from data
data5 <- generate_random_data(100, 2, 1)
formula5 <- y ~ X1 + X2
ref_list5 <- list("NonExistentVar")
tryCatch(test_lmm(formula5, data5, ref_list5), error = function(e) cat(e$message))

# empty data frame
data6 <- data.frame()
formula6 <- y ~ X1
ref_list6 <- list()
tryCatch(test_lmm(formula6, data6, ref_list6), error = function(e) cat(e$message))

# missing values in the data
data7 <- generate_random_data(100, 2, 1)
data7$X1[1:10] <- NA  # Introduce NAs
formula7 <- y ~ X1 + X2
ref_list7 <- list("R1")
tryCatch(test_lmm(formula7, data7, ref_list7), error = function(e) cat(e$message))

## more complex random effect structures ---

# one interaction - this one isn't working as Z gets too wide (if too many levels of R1/R2 are included in data)
data8 <- generate_random_data(100, 2, 2)
formula8 <- y ~ X1 + X2
ref_list8 <- list('R1', c("R1", "R2"))
test_lmm(formula8, data8 , ref_list8)

# one interaction & 3 effects - this one isn't working & not sure why
data9 <- generate_random_data(100, 2, 3)
formula9 <- y ~ X1
ref_list9 <- list('R1', c("R1", "R2"), 'R3')
test_lmm(formula9, data9 , ref_list9)


## -----

library(haven) 

popular2data <- read_sav(file ="https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")
popular2data <- popular2data[c('pupil', 'class', 'extrav', 'sex', 'texp', 'popular')]

test_lmm(popular ~ sex + extrav + texp, dat = popular2data, ref = list('class'))


## ---

lmm.data <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

test_lmm(extro ~ open + agree + social, lmm.data, ref = list('school', c('school', 'class')))
