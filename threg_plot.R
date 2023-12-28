library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(ggfortify)
library(threg)
library(survminer)
library(survival)




table(lung$sex)
lung$sex2 <- as.numeric(factor(lung$sex))-1
table(lung$sex2 )


lung$status2 <- as.numeric(factor(lung$status))-1
table(lung$status2)

lung$time2 <- lung$time/7
summary(lung$time2)

###Kaplan-Meier###
kmfit = survfit(Surv(time2,status2)~sex2, data=lung)


####Standard Threg####

std.threg <- threg(Surv(time2, status2) ~ sex2 | sex2, data = lung)
std.threg


####Bayesian Threg####
data_list <- list(
  Nuc = length(which(lung$status2==1)),
  Nc = length(which(lung$status2==0)),
  N = dim(lung)[1],
  y = c(unlist(lung$time2[which(lung$status2==1)]),unlist(lung$time2[which(lung$status2==0)])),
  T=150,
  K=2,
  Xvar = as.matrix(cbind(c(rep(1,dim(lung)[1])),lung$sex2)) 
)


file <- file.path("~/threg/threg.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


bayes.threg <- mod$sample(
  data = data_list,
  seed = 123,
  init = 2,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  output_dir = "~/temp"
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = .8
  
)

bayes.threg$summary(variables = c("b_lny0","b_mu"))
std.threg

posterior1 <- bayes.threg$draws(format="df")

mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars =c("b_lny0[1]","b_lny0[2]") )
mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars =c("b_mu[1]","b_mu[2]") )



model_draws_Surv <-
  posterior1 %>%
  select(.chain, .iteration, .draw, starts_with("Survival"))
#names(model_draws_Surv)
colnames(model_draws_Surv) = gsub("(\\[)|(\\])", "", colnames(model_draws_Surv))
#names(model_draws_Surv)



#summarize the posterior survival curves
model_draws_Surv0 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival0"),
    names_to = "time",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 0)

model_draws_Surv1 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival1"),
    names_to = "time",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 1)


# ggplot() + 
#   geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
#   geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
#   scale_color_manual(name = "", 
#                      values = c("treatment = 0"="blue", "treatment = 1" = "red"))

autoplot(kmfit, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=time, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=time, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red")) +
  ggtitle("with Covar")


##############################################################################


file <- file.path("~/threg/threg_normal.stan")
mod_n <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


fit_threg_n <- mod_n$sample(
  data = data_list,
  seed = 123,
  init = 2,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  output_dir = "~/temp"
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = .8
  
)

bayes.threg$summary(variables = c("b_lny0","b_mu"))
fit_threg_n$summary(variables = c("b_lny0","b_mu"))
fit_threg_n$summary(variables = c("a_mu","a_sig"))

posterior_n <- fit_threg_n$draws(format="df")
mcmc_areas(posterior_n, prob = 0.90, prob_outer = 1, pars =c("a_mu","a_sig") )
mcmc_areas(posterior_n, prob = 0.90, prob_outer = 1, pars =c("b_mu[1]","b_mu[2]") )


model_draws_Surv <-
  posterior_n %>%
  select(.chain, .iteration, .draw, starts_with("Survival"))
#names(model_draws_Surv)
colnames(model_draws_Surv) = gsub("(\\[)|(\\])", "", colnames(model_draws_Surv))
#names(model_draws_Surv)


#summarize the posterior survival curves
model_draws_Surv0 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival0"),
    names_to = "time",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 0)

model_draws_Surv1 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival1"),
    names_to = "time",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 1)


# ggplot() + 
#   geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
#   geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
#   scale_color_manual(name = "", 
#                      values = c("treatment = 0"="blue", "treatment = 1" = "red"))

autoplot(kmfit, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=time, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=time, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red")) +
  ggtitle("With unobs-heter")


#####################################################################################

file <- file.path("~/threg/threg_normal0.stan")
mod_n0 <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


fit_threg_n0 <- mod_n0$sample(
  data = data_list,
  seed = 123,
  init = 2,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  output_dir = "~/temp"
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = .8
  
)

fit_threg$summary(variables = c("b_lny0","b_mu"))
fit_threg_n$summary(variables = c("b_lny0","b_mu"))
fit_threg_n$summary(variables = c("a_mu","a_sig"))

fit_threg_n0$summary(variables = c("a_mu","a_sig","b_mu"))


posterior_n0 <- fit_threg_n0$draws(format="df")

mcmc_areas(posterior_n0, prob = 0.90, prob_outer = 1, pars =c("a_mu","a_sig") )
mcmc_areas(posterior_n0, prob = 0.90, prob_outer = 1, pars =c("b_mu[1]","b_mu[2]") )


model_draws_Surv <-
  posterior_n0 %>%
  select(.chain, .iteration, .draw, starts_with("Survival"))
names(model_draws_Surv)
colnames(model_draws_Surv) = gsub("(\\[)|(\\])", "", colnames(model_draws_Surv))
names(model_draws_Surv)


#summarize the posterior survival curves
model_draws_Surv0 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival0"),
    names_to = "time",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 0)

model_draws_Surv1 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("Survival1"),
    names_to = "time",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (time = as.numeric(time) ) %>%
  group_by(time) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 1)


# ggplot() + 
#   geom_line(data = model_draws_Surv0, aes( x=time, y=surv_median, color = "treatment = 0")) +
#   geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
#   scale_color_manual(name = "", 
#                      values = c("treatment = 0"="blue", "treatment = 1" = "red"))

autoplot(kmfit, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=time, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=time, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red"))
