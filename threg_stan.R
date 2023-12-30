library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(threg)
library(survminer)

setwd("~/R/threg")

## load the data "lkr"
data("lkr", package = "threg")

## Transform the "treatment2" variable into factor variable "f.treatment2" .
lkr$f.treatment2 <- factor(lkr$treatment2)
table(lkr$treatment2)

## generate the Kaplan-Meier survival curves for the drug B group and 
## the standard group.
fit.km <- survfit(Surv(weeks, relapse) ~ f.treatment2, data = lkr) 

ggsurvplot(fit.km, risk.table = F)

#it is important to stack the data so that all the completed event rows are listed first
data_list <- list(
  Nuc = length(which(lkr$relapse==1)),
  Nc = length(which(lkr$relapse==0)),
  N = dim(lkr)[1],
  y = c(unlist(lkr$weeks[which(lkr$relapse==1)]),unlist(lkr$weeks[which(lkr$relapse==0)])),
  T=40,
  K=2,
  Xvar = as.matrix(cbind(c(rep(1,dim(lkr)[1])),
                         c(unlist(lkr$treatment2[which(lkr$relapse==1)]),unlist(lkr$treatment2[which(lkr$relapse==0)])))) 
)


###########################################################################

file <- file.path("threg.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


fit_threg <- mod$sample(
  data = data_list,
  seed = 123,
  init = 2,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
#  output_dir = "~/temp"
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = .8
  
)

fit_threg$summary(variables = c("b_lny0","b_mu"))

## fit the threshold regression model on the factor variable "f.treatment2"  
standard.threg <- threg(Surv(weeks, relapse) ~ f.treatment2 | f.treatment2, data = lkr)
standard.threg
## coefficients are similar to the posterior means from Stan

treatment0 <- fit_threg$summary(c("Survival0"))

treatment0_ribbons <-
  treatment0 %>%
  select(mean, q5, q95) %>%
  mutate(month = c(1:41))


treatment1 <- fit_threg$summary(c("Survival1"))

treatment1_ribbons <-
  treatment1 %>%
  select(mean, q5, q95) %>%
  mutate(month = c(1:41))


ggsurvplot(fit.km, risk.table = F)$plot + 
  geom_line(data = treatment0_ribbons, aes( x=month, y=mean, color='model0'), linetype='dashed') +
  geom_line(data = treatment1_ribbons, aes( x=month, y=mean, color='model1'), linetype='dashed') +
  scale_color_manual(values=c('red','blue','red','blue'))  
  
###############################################################


file <- file.path("threg_normal.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


fit_threg2 <- mod$sample(
  data = data_list,
  seed = 123,
  init = 2,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  #  output_dir = "~/temp"
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = .8
  
)

fit_threg2$summary(variables = c("a_mu","a_sig","b_mu","b_lny0"))

## fit the threshold regression model on the factor variable "f.treatment2"  
standard.threg <- threg(Surv(weeks, relapse) ~ f.treatment2 | f.treatment2, data = lkr)
standard.threg
## coefficients are similar to the posterior means from Stan

treatment0 <- fit_threg2$summary(c("Survival0"))

treatment0_ribbons <-
  treatment0 %>%
  select(mean, q5, q95) %>%
  mutate(month = c(1:41))


treatment1 <- fit_threg2$summary(c("Survival1"))

treatment1_ribbons <-
  treatment1 %>%
  select(mean, q5, q95) %>%
  mutate(month = c(1:41))


ggsurvplot(fit.km, risk.table = F)$plot + 
  geom_line(data = treatment0_ribbons, aes( x=month, y=mean, color='model0'), linetype='dashed') +
  geom_line(data = treatment1_ribbons, aes( x=month, y=mean, color='model1'), linetype='dashed') +
  scale_color_manual(values=c('red','blue','red','blue'))
