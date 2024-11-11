library(nlme);library(lme4)
lmm(score ~ Machine, Machines, list("Worker",c("Worker","Machine")))

lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), data = Machines, REML=FALSE)

