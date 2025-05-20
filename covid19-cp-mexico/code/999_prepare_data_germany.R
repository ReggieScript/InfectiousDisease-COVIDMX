library(tidyverse)
library(readr)
library(lubridate)
library(gamlss)
library(lubridate)

file_path <- file.path("../data/raw/RKI-2021-01-05.csv")

#Read and convert Date cols to dates
covid <- read_csv(file_path) %>%  
  mutate(across(c(Refdatum, Meldedatum), as.Date)) %>% dplyr::select(-Altersgruppe2)


#' Create age-group and select disease onset (if reference date is disease onset date), 
#' reporting date, age-group, federal state, and case counts
covid %>% distinct(Altersgruppe)

covid = covid %>% mutate(age_group = factor(Altersgruppe, levels = c("A00-A04",
                                                             "A05-A14",
                                                             "A15-A34",
                                                             "A35-A59",
                                                             "A60-A79",
                                                             "A80+",
                                                             "unbekannt"),
                                    labels = c("0-14", "0-14", "15-59", "15-59", "60-79", "80+", "unknown"))) %>%
  mutate(onset=if_else(IstErkrankungsbeginn==1, as.character(Refdatum), NA_character_)) %>%
  dplyr::select(age_group, onset, rep_date = Meldedatum, federal_state = Bundesland, num_case = AnzahlFall) %>%
  filter(onset>=ymd("2020-02-15")|is.na(onset),
         onset<=ymd("2020-09-28")|is.na(onset),
         rep_date>=ymd("2020-03-01"),
         rep_date<=ymd("2020-06-01"),
         num_case>=1,
         onset<=rep_date|is.na(onset),
         !is.na(age_group)) %>%
  mutate(onset=ymd(onset))

# Time-series disease onset if available
# Germany
covid %>% filter(!is.na(onset)) %>% group_by(date=onset) %>% summarise(n=sum(num_case)) %>%
  ggplot() + geom_line(aes(date, n)) + 
  geom_line(aes(date, n), covid %>% filter(!is.na(onset)) %>% 
              group_by(date=rep_date) %>% 
              summarise(n=sum(num_case)), col = "red")

# Impute missing disease onsets
imputation_dat = covid %>% slice(rep(1:n(), times = num_case)) %>%
  dplyr::select(-num_case) %>%
  mutate(rep_date_weekday = weekdays(rep_date),
         rep_week = as.numeric(strftime(rep_date, format = "%V")),
         delay = as.numeric(rep_date - onset),
         # Add age-group unknown to biggest age-group
         age_group_model =factor(age_group, 
                                 levels = c("0-14", "15-59","60-79", "80+","unknown"),
                                 labels = c("0-14", "15-59","60-79", "80+","15-59")))

# Training data with available onset
train_dat = imputation_dat %>% dplyr::filter(!is.na(onset))
# Data for imputation
imp_dat = imputation_dat %>% dplyr::filter(is.na(onset))

fed_states = imputation_dat %>% pull(federal_state) %>% unique()

imputed_dat_fed_state = list()
plot_dat_fed_state = list()
#' Estimate imputation model per state and impute missing onsets, 
#' safe imputed data in list and list for plot
#' The estimation of the model is currently not performed,
#' previously estimated models are loaded from the folder
#' ./results/imputation_ger
#' uncomment the corresponding rows for the estimation of the models

set.seed(124124)
for(state in fed_states) {
  # Restrict training data to specific state
  train_dat_state = train_dat %>%
    dplyr::filter(federal_state==state) %>%
    dplyr::select(onset, rep_date, rep_week, delay, rep_date_weekday, age_group_model, federal_state) %>%
    mutate(delay_prep = if_else(delay == 0, 1e-2, delay))
  
  # Estimate imputation model
  # imp_mod = gamlss::gamlss(
  #   delay_prep  ~ gamlss::cs(rep_week) + age_group_model + rep_date_weekday,
  #   sigma.formula = ~ gamlss::cs(rep_week) + age_group_model + rep_date_weekday,
  #   family = gamlss.dist::WEI2,
  #   data = train_dat_state,
  #   method = mixed(400, 20))
  # 
  # save(imp_mod, file = paste0("../results/imputation_ger/imp_model", state, ".RData")) 
  load(paste0("../results/imputation_ger/imp_model", state, ".RData"))
  # Create dataset to plot imputation results
  plot_imp_mod_dat <- data.frame(expand.grid(rep_week = unique(train_dat_state$rep_week), 
                                             age_group_model=unique(train_dat_state$age_group_model),
                                             rep_date_weekday = unique(train_dat_state$rep_date_weekday)))
  plot_dat_fed_state[[state]] = plot_imp_mod_dat %>% 
    mutate(mu =exp(predict(imp_mod, newdata = plot_imp_mod_dat, 
                           what = "mu",
                           data = train_dat_state)),
           sigma = exp(predict(imp_mod, newdata = plot_imp_mod_dat, 
                               what="sigma",
                               data = train_dat_state)),
           med = gamlss.dist:::qWEI2(0.5, mu=mu, sigma=sigma),
           q25 = gamlss.dist:::qWEI2(0.25, mu=mu, sigma=sigma),
           q75 = gamlss.dist:::qWEI2(0.75, mu=mu, sigma=sigma),
           rep_week_local = factor(rep_week, levels = sort(unique(rep_week)), 
                                   labels = sort(unique(rep_week))),
           rep_date_local_weekday = factor(rep_date_weekday,
                                           levels = c("Monday", "Tuesday", "Wednesday",
                                                      "Thursday", "Friday", "Saturday",
                                                      "Sunday")),
           age = as.numeric(as.character(factor(age_group_model,
                                                levels = c("0-14", 
                                                           "15-59", 
                                                           "60-79", 
                                                           "80+"),
                                                labels = c("7", "37",
                                                           "70", "90")))),
           federal_state=state)
  
  # Impute missing onsets
  imp_dat_state = imp_dat %>%
    dplyr::filter(federal_state==state) 
  
  imputed_dat_fed_state[[state]] = imp_dat_state %>%
    mutate(
      # Get estimated person-specific location and scale
      mu = exp(predict(imp_mod, newdata = imp_dat_state %>% dplyr::select(rep_week, age_group_model, 
                                                                          rep_date_weekday, federal_state),
                       what = "mu", data = train_dat_state)),
      sigma = exp(predict(imp_mod, newdata = imp_dat_state %>% dplyr::select(rep_week, age_group_model, 
                                                                             rep_date_weekday, federal_state),
                          what = "sigma", data = train_dat_state)),
      # Randomly sample from person-specific delay distribution
      delay_imputed = gamlss.dist::rWEI2(nrow(imp_dat_state), mu = mu, sigma = sigma),
      # Get day of symptom onset assuming all cases reported at noon (i.e. delay <0.5 -> reporting at same 
      # day as symptom onset)
      onset = ymd(rep_date - delay_imputed + 0.5))
  
}

#
plot_imp = do.call(rbind, plot_dat_fed_state)
imp_median_weibull = ggplot(aes(age, med, col = rep_week_local), data = plot_imp) +
  geom_line() + facet_grid(federal_state~rep_date_local_weekday) +
  scale_color_discrete(name="Week") +
  theme_bw() + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 4)) +
  xlab(label = "Age") + ylab("Delay")

ggsave(imp_median_weibull, filename = "../results/imputation_ger/figure_median.png", width = 10, height = 25)

# Derive overall dataset
imp_dat = do.call(rbind, imputed_dat_fed_state)
# Observed and imputed
onset_age_group_ger = rbind(train_dat %>% dplyr::select(age_group, date=onset, federal_state),
                imp_dat %>% dplyr::select(age_group, date=onset, federal_state)) %>%
  group_by(age_group, federal_state, date) %>%
  summarise(onsets=n()) %>%
  dplyr::filter(date>=ymd("2020-02-15"))
  

write.table(onset_age_group_ger, file = "../data/onset_age_group_ger.csv", row.names = FALSE,
            sep="\t", quote = FALSE)

onset_rep_age_group_ger = onset_age_group_ger %>% full_join(covid %>% 
                                                              group_by(date=rep_date, 
                                                                       age_group,
                                                                       federal_state) %>% 
                                                              summarise(reported=sum(num_case)))
write.table(onset_rep_age_group_ger, file = "../data/onset_reported_age_group_ger.csv", row.names = FALSE,
            sep = "\t", quote = FALSE)
# Only observed (sensitivity analysis)
onset_age_group_ger_obs = train_dat %>% dplyr::select(age_group, date=onset, federal_state) %>%
  group_by(age_group, federal_state, date) %>%
  summarise(onsets=n()) %>%
  dplyr::filter(date>=ymd("2020-02-15"))


write.table(onset_age_group_ger_obs, file = "../data/onset_age_group_ger_obs.csv", row.names = FALSE,
            sep="\t", quote = FALSE)

# Observed and non-observed reporting date
onset_age_group_ger_obsrep2 = rbind(train_dat %>% dplyr::select(age_group, date=onset, federal_state),
                                imp_dat %>% dplyr::mutate(date = rep_date)  %>%
                                  dplyr::select(age_group, date, federal_state)) %>%
  group_by(age_group, federal_state, date) %>%
  summarise(onsets=n()) %>%
  dplyr::filter(date>=ymd("2020-02-15"))
write.table(onset_age_group_ger_obsrep2, file = "../data/onset_age_group_ger_obsrep2.csv", row.names = FALSE,
            sep="\t", quote = FALSE)

# Plot Imputed disease onset, reported disease onsets and reported cases per day
onset_age_group_ger %>% group_by(date) %>%
  summarise(n=sum(onsets)) %>% ggplot() + 
  geom_line(aes(date, n)) + 
  geom_line(aes(date, n), 
            data = covid %>% filter(!is.na(onset)) %>%
              group_by(date=onset) %>%
              summarise(n=n()), lty=2) +
  geom_line(aes(date, n), 
            data = covid %>%
              group_by(date=rep_date) %>%
              summarise(n=n()), col = "red")

