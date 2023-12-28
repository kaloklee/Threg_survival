library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(ggfortify)
library(threg)

## load the data "lkr"
data("lkr", package = "threg")

#use threg to check the coefficients with just intercept
fit0 <- threg(Surv(weeks, relapse) ~ 1 | 1, data = lkr)
fit0

#prepare Stan data
data_list <- list(
  Nuc = length(which(lkr$relapse==1)), 
  yuc = unlist(lkr$weeks[which(lkr$relapse==1)]),
  Nc = length(which(lkr$relapse==0)),
  yc = unlist(lkr$weeks[which(lkr$relapse==0)]),
  T=40
)

#compile Stan model.  The 'stanc_options =' can be omitted.
file <- file.path("~/threg/threg_nocov.stan")
mod0 <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)

fit_threg0 <- mod0$sample(
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


fit_threg0$summary(variables = c("y0","mu","lny0"))

posterior0 <- fit_threg0$draws(format="df")

#check parameters convergence and diagnostics
mcmc_trace(posterior0, pars = c("y0","mu"))

mcmc_intervals(posterior0, pars =c("y0","mu"))

mcmc_areas(posterior0, prob = 0.90, prob_outer = 1, pars =c("y0") )
mcmc_areas(posterior0, prob = 0.90, prob_outer = 1, pars =c("mu") )
 
mcmc_dens_overlay(posterior0, pars=c("y0","mu")) 

mcmc_pairs(posterior0,pars=c("y0","mu"))

mcmc_hex(posterior0, pars=c("y0","mu"))

mcmc_scatter(posterior0,pars=c("y0","mu"))

####################################################################


model_draws_Surv <-
  posterior0 %>%
  select(.chain, .iteration, .draw, starts_with("expected["))
names(model_draws_Surv)
colnames(model_draws_Surv) = gsub("(\\[)|(\\])", "", colnames(model_draws_Surv))
names(model_draws_Surv)

#summarize the posterior survival curves
model_draws_Surv2 <-
  model_draws_Surv %>%
  pivot_longer(
    cols = starts_with("expected"),
    names_to = "week",
    names_prefix = "expected",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() 
#  %>%
#  arrange(week)

View(model_draws_Surv2)


fit.km<-survfit(Surv(weeks, relapse) ~ 1, data = lkr, conf.int=T)  

km_df <- data.frame(fit.km$time,fit.km$surv)
 
names(km_df)
 
km_df <- km_df %>%
rename(time = fit.km.time, surv = fit.km.surv) 

model_draws_Surv2$week=as.numeric(model_draws_Surv2$week)


autoplot(fit.km, conf.int = F, surv.colour ='red') + 
  geom_line(data = model_draws_Surv2, aes( x=week , y=surv_median), inherit.aes = F, color='blue' ) +
  geom_ribbon(data = model_draws_Surv2, aes(x=week, ymin=surv_p5, ymax=surv_p95), alpha=.2, inherit.aes = F)

+ scale_color_manual(name = "", 
                   labels = c("K-M", "posterior prediction") ,
                   values = c("red","blue"))


autoplot(fit.km, conf.int = F) + 
  geom_line(data = model_draws_Surv2, aes( x=week, y=surv_median, color = "posterior prediction")) +
  geom_ribbon(data = model_draws_Surv2, aes(x=week, ymin=surv_p5, ymax=surv_p95), alpha=.2) 

  labs(x = "\n Escalation Time (Year) ", y = "% Not Escalated \n", 
       title = "Time to Escalate from Initial Dx") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .1), limits=c(0,1)) +
  geom_vline(xintercept = 1, linetype="dashed", 
             color = "red", size=.5)

#Plot them side by side.  Output graph is omitted.
ggplot() + 
  geom_line(data = model_draws_Surv, aes( x=timing, y=surv_median, color = "posterior prediction")) +
  geom_ribbon(data = model_draws_Surv, aes(x=timing, ymin=surv_p5, ymax=surv_p95), alpha=.2) +
  geom_line(data = km_df, aes( x=time, y=surv, color = "K-M")) + 
  scale_color_manual(name = "", 
                     values = c("posterior prediction"="blue", "K-M" = "red"))


######################################################################################

## Transform the "treatment2" variable into factor variable "f.treatment2" .
lkr$f.treatment2 <- factor(lkr$treatment2)

## generate the Kaplan-Meier survival curves for the drug B group and 
## the standard group.
fit.km <- survfit(Surv(weeks, relapse) ~ f.treatment2, data = lkr) 
plot(fit.km, lty = 1:2, xlab="time",ylab="Estimated S(t)") 
legend(20, 1, c("Standard","Drug B"), lty = 1:2) 


## fit the threshold regression model on the factor variable "f.treatment2", 
fit.threg <- threg(Surv(weeks, relapse) ~ f.treatment2 | f.treatment2, data = lkr)
fit.threg

## generate the threshold regression predicted survival curves for the drug B group and 
## the standard group.
plot(fit, var = f.treatment2, graph = "sv", nolegend = 1, nocolor = 1)
legend(20, 1, c("Standard", "Drug B"), lty = 1:2) 


lkr$f.treatment2num <- as.numeric(factor(lkr$treatment2))-1
table(lkr$f.treatment2num)

data_list <- list(
  Nuc = length(which(lkr$relapse==1)),
  Nc = length(which(lkr$relapse==0)),
  N = dim(lkr)[1],
  y = c(unlist(lkr$weeks[which(lkr$relapse==1)]),unlist(lkr$weeks[which(lkr$relapse==0)])),
  T=40,
  K=2,
  Xvar = as.matrix(cbind(c(rep(1,dim(lkr)[1])),lkr$f.treatment2num)) 
)


file <- file.path("~/threg/threg.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)


fit_threg <- mod$sample(
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
fit

posterior1 <- fit_threg$draws(format="df")


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
    names_to = "week",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
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
    names_to = "week",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
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

autoplot(fit.km, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red")) +
  ggtitle("with Covar")


##############################################################################


data_list <- list(
  Nuc = length(which(lkr$relapse==1)),
  Nc = length(which(lkr$relapse==0)),
  N = dim(lkr)[1],
  y = c(unlist(lkr$weeks[which(lkr$relapse==1)]),unlist(lkr$weeks[which(lkr$relapse==0)])),
  T=40,
  K=2,
  Xvar = as.matrix(cbind(c(rep(1,dim(lkr)[1])),lkr$f.treatment2num)) 
)


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

fit_threg$summary(variables = c("b_lny0","b_mu"))
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
    names_to = "week",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
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
    names_to = "week",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
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

autoplot(fit.km, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red")) +
  ggtitle("With un-heter")


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
    names_to = "week",
    names_prefix = "Survival0",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
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
    names_to = "week",
    names_prefix = "Survival1",
    values_to = "Surv",
    values_drop_na = TRUE
  ) %>%
  mutate (week = as.numeric(week) ) %>%
  group_by(week) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate ( treatment = 1)


ggplot() + 
  geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red"))

autoplot(fit.km, conf.int = F) + 
  geom_line(data = model_draws_Surv0, aes( x=week, y=surv_median, color = "treatment = 0")) +
  geom_line(data = model_draws_Surv1, aes( x=week, y=surv_median, color = "treatment = 1")) +
  scale_color_manual(name = "", 
                     values = c("treatment = 0"="blue", "treatment = 1" = "red"))
