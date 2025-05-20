# library(tidyverse)
# theme_set(theme_bw())
# library(lubridate)
# Load prepared data: timeseries number of disease onsets per day
bav = read_tsv("../data/onset_age_group_bav.csv")
ger = read_tsv("../data/onset_age_group_ger.csv")

# Comparison epidemic curve in Bavaria based on impuation of German and Bavarian data
bav %>% group_by(date) %>% summarise(n=sum(onsets)) %>%
  filter(date<=ymd("2020-05-01"),
         date>=ymd("2020-02-15")) %>%
  mutate(Datasource="Bavaria (LGL)") %>%
  rbind(ger %>% filter(federal_state=="Bayern") %>%
          group_by(date) %>%
          summarise(n=sum(onsets)) %>%
          filter(date<=ymd("2020-05-01"),
                 date>=ymd("2020-02-15")) %>%
          mutate(Datasource="Germany (RKI)")) %>%
ggplot() +
  geom_line(aes(date, n, col = Datasource)) +
  ylab("Number") +
  ggtitle("Bavaria (disease onsets obs + imputed)")+
  ggsave("../results/imputation_ger/comp_onset_imputed_ger_bav_dat.png", width=8, height = 4)


bav_full = bav %>% group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))

ger_full = ger %>% group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))

## Estimate changepoints
# library(future)
# library(future.apply)
# library(segmented)
# library(surveillance)
# library(MASS)
# library(nlme)
# plan("multisession")
# options(mc.cores = 6)

source("cp_fun.R")

cp_res_bav_full = perform_cp_analysis(data = bav_full,
                                      type = "both",
                                      cp_max_onset = 6,
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T,
                                      use_disc_optim_results = T,
                                      name_disc =  "bav_full")


cp_res_ger_full = perform_cp_analysis(data = ger_full,
                                      type = "both",
                                      cp_max_onset = 6,
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T,
                                      use_disc_optim_results = T,
                                      name_disc =  "ger_full")

save(cp_res_bav_full, file = "../results/changepoint/cp_results/cp_res_bav.RData")
save(cp_res_ger_full, file = "../results/changepoint/cp_results/cp_res_ger.RData")

## Save plots
lapply(cp_res_bav_full$cp_segmented_list_backpro, function(x){
  ggsave(x$plot, file = paste0("../results/changepoint/cp_results/figures/bav_backpro/cp_", nrow(x$breakpoints), ".png"),
         width = 10, height = 6)
})
lapply(cp_res_bav_full$cp_segmented_list_onset, function(x){
  ggsave(x$plot, file = paste0("../results/changepoint/cp_results/figures/bav_onset/cp_", nrow(x$breakpoints), ".png"),
         width = 10, height = 6)
})

lapply(cp_res_ger_full$cp_segmented_list_backpro, function(x){
  tryCatch(ggsave(x$plot, file = paste0("../results/changepoint/cp_results/figures/ger_backpro/cp_", nrow(x$breakpoints), ".png"),
         width = 10, height = 6), error=function(e) print("missing"))
})
lapply(cp_res_ger_full$cp_segmented_list_onset, function(x){
  tryCatch(ggsave(x$plot, file = paste0("../results/changepoint/cp_results/figures/ger_onset/cp_", nrow(x$breakpoints), ".png"),
         width = 10, height = 6), error=function(e) print("missing"))
})

## Selected models based on BIC
# Backprojection
# Bavaria
xtable::xtable(cp_res_bav_full$cp_segmented_list_backpro[[which.min(cp_res_bav_full$bic_backpro)]]$breakpoints)
xtable::xtable(cp_res_bav_full$cp_segmented_list_backpro[[which.min(cp_res_bav_full$bic_backpro)]]$coef)

# Germany
xtable::xtable(cp_res_ger_full$cp_segmented_list_backpro[[which.min(cp_res_ger_full$bic_backpro)]]$breakpoints)
xtable::xtable(cp_res_ger_full$cp_segmented_list_backpro[[which.min(cp_res_ger_full$bic_backpro)]]$coef)

## Selected models based on BIC
# Onsets
# Bavaria
xtable::xtable(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$breakpoints)
xtable::xtable(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$coef)

# Germany
xtable::xtable(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_backpro)]]$breakpoints)
xtable::xtable(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_backpro)]]$coef)

# Figure S1. Results of change point models for onset data in Bavaria/Germany
theme = theme(
  axis.text=element_text(size = rel(1.3)),
  axis.title=element_text(size = rel(1.3)),
  legend.text = element_text(size = rel(1.3)),
  legend.title =element_text(size = rel(1.3)))

fig_s1_onset_bav_ger = ggpubr::ggarrange(
  cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot +
    ylab("Number Onsets") + ggtitle("Bavaria") + 
    scale_x_continuous(breaks = seq(max(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$t), 
                                    min(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$t), by = -14),
                       labels = paste0(strftime(
                         seq(max(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$date), 
                             min(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$date), by = -14), 
                         format = "%d.%m."), "\n(",
                       strftime(
                         seq(max(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$date), 
                             min(cp_res_bav_full$cp_segmented_list_onset[[which.min(cp_res_bav_full$bic_onset)]]$plot$plot_env$dat$date), by = -14)-5, 
                         format = "%d.%m."), ")")) + theme,  
  cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot +
    ylab("Number Onsets") + ggtitle("Germany") +
    scale_x_continuous(breaks = seq(max(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$t), 
                                    min(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$t), by = -14),
                       labels = paste0(strftime(
                         seq(max(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$date), 
                             min(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$date), by = -14), 
                         format = "%d.%m."), "\n(",
                         strftime(
                           seq(max(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$date), 
                               min(cp_res_ger_full$cp_segmented_list_onset[[which.min(cp_res_ger_full$bic_onset)]]$plot$plot_env$dat$date), by = -14)-5, 
                           format = "%d.%m."), ")")) + theme, 
  labels = "AUTO")

ggsave(fig_s1_onset_bav_ger, filename = "../results/changepoint/fig_s1_bponsets_ger_bav.png",
       width=12, height = 4)
