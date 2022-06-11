## estimate missing abundance in 2020,
## fit Bayesian regression model, 
## and forecast 2022 value

library(tidyverse)
library(mice)
library(rstan)
library(brms)
library(bayesplot)

source("./scripts/stan_utils.R")

theme_set(theme_bw())


# load borealization DFA trend

trend <- read.csv("./output/dfa_trend.csv")

trend <- trend %>%
  select(t, estimate) %>%
  rename(year = t, trend = estimate)

# load snow crab abundance

abundance <- read.csv("./Data/imm_abun.csv", row.names = 1)

# clean up and combine

abundance <- abundance %>%
  rename(year = AKFIN_SURVEY_YEAR) %>%
  mutate(log_abundance = log(ABUNDANCE_MIL), .keep = "unused")

dat <- left_join(trend, abundance)

plot.dat <- dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)


# estimate 2020 with mice ---------------------------------------
library(mice)

# load prediction time series
pred_dat <- read.csv("./Data/imputation_data.csv")

# clean up and remove survey data - can consider these in a lagged application if needed
pred_dat <- pred_dat[1:42,1:8] %>%
  select(-survey_female_mat_biomass, -survey_male_mat_biomass)

# and log transform

pred_dat[,2:6] <- log(pred_dat[,2:6])

dat <- abundance %>%
  rbind(.,
        data.frame(year = 2020,
                   log_abundance = NA)) %>%
  arrange(year) %>%
  left_join(., pred_dat)

# examine correlations!
cor(dat[,-1], use = "p")

imp <- mice(data = dat, method = "norm", m = 100, seed = 957)

# pull out 2020 estimates
estimated_2020 <- NA


for(i in 1:100){

estimated_2020[i] <- complete(imp, i)$log_abundance[41]

}


# plot estimate value compared with observed time series
estimated <- data.frame(year = 2020,
                        log_abundance = mean(estimated_2020),
                        LCI = quantile(estimated_2020, 0.025),
                        UCI = quantile(estimated_2020, 0.975))


abundance.plot <- abundance %>%
  mutate(LCI = NA,
         UCI = NA)

abundance.plot <- rbind(abundance.plot, estimated)

ggplot(abundance.plot, aes(year, log_abundance)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
  ylab("log(immature snow crab abundance)") +
  theme(axis.title.x = element_blank())

## process imp and pass mice imputations to brms ------------------

imputed.data <- list()

# get a df to plot y values for diagnostics
y <- data.frame()

for(i in 1:100){

  temp <- complete(imp, i)[,1:2] %>%
    mutate(log_abundance_lead1 = lead(log_abundance))

  y <- rbind(y, data.frame(
    year = temp$year[1:41],
    log_abundance_lead1 = temp$log_abundance_lead1[1:41])
  )

  # add borealization trend
  temp <- left_join(temp, trend)

  imputed.data[[i]] <- temp

  }

# set up and run brms

form <- bf(log_abundance_lead1 ~ log_abundance + s(trend))

## fit
mice_brm <- brm_multiple(form,
              data = imputed.data,
              seed = 992,
              cores = 4, chains = 4, iter = 3000,
              save_pars = save_pars(all = TRUE),
              control = list(adapt_delta = 0.99999999, max_treedepth = 12))

saveRDS(mice_brm, file = "output/mice_brm.rds") # too big to push to github!

summary(mice_brm)

# diagnostics
mice_brm <- readRDS("./output/mice_brm.rds")
check_hmc_diagnostics(mice_brm$fit)
neff_lowest(mice_brm$fit)


bayes_R2(mice_brm) # this only fits for the first model

plot(conditional_effects(mice_brm), ask = FALSE)

y <- y %>%
  group_by(year) %>%
  summarise(log_abundance_lead1 = mean(log_abundance_lead1))

y <- y$log_abundance_lead1

yrep_mice_brm  <- fitted(mice_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_mice_brm[sample(nrow(yrep_mice_brm), 25), ]) +
  ggtitle("mice_brm")

ggsave("./Figs/mice_brms_ppc.png", width = 6, height = 4, units = 'in')

trace_plot(mice_brm$fit)

# plot mice_brm

## 95% CI
ce1s_1 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mice_brm, effect = "trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(trend$trend[trend$year %in% 1980:2020]), amount = 0.01),
                          rep(NA, 100-length(trend$trend[trend$year %in% 1980:2020])))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Borealization index", y = "Log abundance") +
  geom_rug(aes(x=rug.anom, y=NULL))

print(g1)

ggsave("./Figs/mice_borealization_abundance_regression.png", width = 4, height = 3, units = 'in')


# now predict for 2022 survey
new.dat <- data.frame(log_abundance = abundance$log_abundance[abundance$year == 2021],
                      trend = trend$trend[trend$year == 2021])


pred.2022 <- posterior_epred(mice_brm, newdata = new.dat)


mean(pred.2022)
LCI.80 <- exp(quantile(pred.2022, 0.1))
UCI.80 <- exp(quantile(pred.2022, 0.9))

LCI.95 <- exp(quantile(pred.2022, 0.025))
UCI.95 <- exp(quantile(pred.2022, 0.975))

overall.mean <- exp(mean(abundance$log_abundance[abundance$year %in% 1980:2019]))

LCI.80 / overall.mean
UCI.80 / overall.mean

LCI.95 / overall.mean
UCI.95 / overall.mean

# probability of greater value in 2022 than 2021
sum(pred.2022 > abundance$log_abundance[abundance$year == 2021]) / length(pred.2022)

# add blank year to abundance.plot
xtra.80 <- data.frame(year = 2022,
                   log_abundance = 5.5, # dummy value to make plot work
                   LCI = quantile(pred.2022, 0.1),
                   UCI = quantile(pred.2022, 0.9))

# abundance.plot <- rbind(abundance.plot, xtra)

# create data frame for 95% CI
xtra.95 <- data.frame(year = 2022,
                      log_abundance = 5.5, # dummy value to make plot work
                      LCI = quantile(pred.2022, 0.025),
                      UCI = quantile(pred.2022, 0.975))

ggplot(abundance.plot, aes(year, log_abundance)) +
  geom_line() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), color = "dark grey") +
  geom_point(size = 2) +
  geom_errorbar(data = xtra.95, aes(x = year, ymin = LCI, ymax = UCI), color = "gold") +
  geom_errorbar(data = xtra.80, aes(x = year, ymin = LCI, ymax = UCI), color = "firebrick") +

  ylab("Log abundance") +
  theme(axis.title.x = element_blank())

# save plot
ggsave("./Figs/imputed_and_predicted_survey_abundance.png", width = 5, height = 3, units = 'in')
