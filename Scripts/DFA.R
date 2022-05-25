# combine multiple indicators of borealization with 
# a Dynamic Factor Analysis model

library(tidyverse)
library(MARSS)

theme_set(theme_bw())

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## data processing -------------------------------------

d1 <- read.csv("./Data/bloom_timing.csv")

d1 <- d1 %>%
  filter(sub_6domain == "south_middle") %>%
  rename(value = peak_mean) %>%
  mutate(name = "Bloom timing") %>%
  select(year, name, value)
  
d2 <- read.csv("./Data/bloom_type.csv")

d2 <- d2 %>%
  filter(!is.na(bloom_type),
         sub_6domain == "south_middle")

ggplot(d2, aes(year, count, color = bloom_type)) +
  geom_line() +
  geom_point() 

d2 <- d2 %>%
  pivot_wider(names_from = bloom_type,
              values_from = count) %>%
  rename(`Ice-edge bloom` = ice_full,
         `Open water bloom` = ice_free) %>%
  pivot_longer(cols = c(-year, -sub_6domain)) %>%
  select(year, name, value)

d3 <- read.csv("./Data/ice.csv")

d3 <- d3 %>%
  rename(`Jan-Feb ice cover` = JanFeb_ice,
         `Mar-Apr ice cover` = MarApr_ice) %>%
  filter(year >= 1972) %>%
  pivot_longer(cols = -year)

d4 <- read.csv("./Data/chl_a.csv")

d4 <- d4 %>%
  select(-north_int_chla, -north_fract_chla) %>%
  rename(`Bloom size` = south_int_chla,
         `Phytoplankton size` = south_fract_chla) %>%
  pivot_longer(cols = -year)

d5 <- read.csv("./Data/bcs_prev.csv", row.names = 1)

plot <- d5 %>%
  pivot_longer(cols = -year)

ggplot(plot, aes(year, value, color = name)) +
  geom_line() +
  geom_point()

d5 <- d5 %>% 
  select(year, Population) %>%
  rename(`Hematodinium` = Population) %>%
  pivot_longer(cols = -year)
  
d6 <- read.csv("./Data/date_corrected_bottom_temp.csv")

d6 <- d6 %>%
  rename(value = bottom.temp) %>%
  mutate(name = 'Bottom temperature')

d7 <- read.csv("./Data/groundfish_med_cpue.csv", row.names = 1)

d7 <- d7 %>%
  rename(year = YEAR,
         `Pacific cod` = med_cod_CPUE,
         `Arctic groundfish` = med_arctic_CPUE) %>%
  pivot_longer(cols = -year)

d8 <- read.csv("./Data/summarized_zooplankton.csv")

d8 <- d8 %>%
  select(-LCI, -UCI) %>%
  rename(name = group, 
         value = log_abundance)

dat <- rbind(d1, d2, d3, d4, d5, d6, d7, d8)

# reorder for plot
plot.order <- data.frame(name = unique(dat$name),
                         order = c(6,5,4,1,2,7,8,11,3,12,13,9,10))

dat <- left_join(dat, plot.order)

dat$name <- reorder(dat$name, dat$order)

ggplot(dat, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 5) +
  theme(axis.title.x = element_blank())

ggsave("./Figs/borealization_time_series.png", width = 12, height = 5, units = 'in')

dfa.dat <- dat %>%
  select(-order) %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(year) %>%
  select(-year) %>%
  t()
  
colnames(dfa.dat) <- unique(d3$year)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models & store results
for(R in levels.R) {
  for(m in 1:2) {  # find best single-trend model
    
    dfa.model = list(A="zero", R=R, m=m)
    
    kemz = MARSS(dfa.dat, model=dfa.model,
                 form="dfa", z.score=TRUE, control=cntl.list)
    
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data # diagonal and unequal is the best model but doesn't converge

# save model selection table
write.csv(model.data, "./output/dfa_model_selection_table.csv",
          row.names = F)

## fit the best model --------------------------------------------------
model.list = list(A="zero", m=2, R="diagonal and unequal") # best model 

# not sure that these changes to control list are needed for this best model, but using them again!
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# rotate
# get the inverse of the rotation matrix
Z.est <- coef(mod, type = "matrix")$Z

H.inv <- varimax(coef(mod, type = "matrix")$Z)$rotmat

# rotate factor loadings
Z.rot <- Z.est %*% H.inv

# rotate trends
trends.rot <- solve(H.inv) %*% mod$states

# Add CIs to marssMLE object
mod <- MARSSparamCIs(mod)

# Use coef() to get the upper and lower CIs
Z.low <- coef(mod, type = "Z", what = "par.lowCI")
Z.up <- coef(mod, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv

plot.CI <- data.frame(names=rownames(dfa.dat),
  mean = as.vector(Z.rot),
  upCI = as.vector(Z.rot.up),
  lowCI = as.vector(Z.rot.low)
)
plot.CI

# none of the loadings can be distinguished from 0! 

# fit the second-best model (1 trend diagonal and unequal)

model.list = list(A="zero", m=1, R="diagonal and unequal") # second-best model 

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                          mean=CI$par$Z[1:13],
                          upCI=CI$par.upCI$Z[1:13],
                          lowCI=CI$par.lowCI$Z[1:13])

dodge <- position_dodge(width=0.9)


plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:13])

loadings.plot <- ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=9), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

# plot trend
trend <- data.frame(t=1972:2021,
                        estimate=as.vector(mod$states),
                        conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                        conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))


trend.plot <- ggplot(trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[2]) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[2]) + xlab("") + ylab("Trend")


# save
png("./Figs/borealization_DFA_loadings_trend.png", width = 9, height = 3.5, units = 'in', res = 300)

ggpubr::ggarrange(loadings.plot,
                  trend.plot,
                  ncol = 2,
                  widths = c(0.45, 0.55),
                  labels = "auto")

dev.off()

# and save loadings and trend
write.csv(plot.CI, "./Output/dfa_loadings.csv", row.names = F)
write.csv(trend, "./Output/dfa_trend.csv", row.names = F)
