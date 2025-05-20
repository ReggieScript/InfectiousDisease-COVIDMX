# library(tidyverse)
# theme_set(theme_bw())
# library(lubridate)
# Load prepared data number of disease onsets per day
bav = read_tsv("../data/onset_age_group_bav_obs.csv")
ger = read_tsv("../data/onset_age_group_ger_obs.csv")
bav2 = read_tsv("../data/onset_age_group_bav_obsrep2.csv")
ger2 = read_tsv("../data/onset_age_group_ger_obsrep2.csv")

# Data editing
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


bav_full2 = bav2 %>% group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))

ger_full2 = ger2 %>% group_by(date) %>%
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
# plan("multisession")
# options(mc.cores = 6)
source("cp_fun.R")

cp_res_bav_full = perform_cp_analysis(data = bav_full, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "bav_full_sens")

cp_res_ger_full = perform_cp_analysis(data = ger_full, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "ger_full_sens")

save(cp_res_bav_full, file = "../results/changepoint/cp_results/cp_res_bav_sens.RData")
save(cp_res_ger_full, file = "../results/changepoint/cp_results/cp_res_ger_sens.RData")

# Sensitivity 2
cp_res_bav_full2 = perform_cp_analysis(data = bav_full2, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "bav_full_sens2")

cp_res_ger_full2 = perform_cp_analysis(data = ger_full2, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "ger_full_sens2")

save(cp_res_bav_full2, file = "../results/changepoint/cp_results/cp_res_bav_sens2.RData")
save(cp_res_ger_full2, file = "../results/changepoint/cp_results/cp_res_ger_sens2.RData")

# Results
load("../results/changepoint/cp_results/cp_res_bav_sens.RData")
cp_res_bav_full1 = cp_res_bav_full
load("../results/changepoint/cp_results/cp_res_bav.RData")
load("../results/changepoint/cp_results/cp_res_bav_sens2.RData")

plot_bav = tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
         pred = exp(predict(cp_res_bav_full$cp_segmented_list_backpro[[which.min(cp_res_bav_full$bic_backpro)]]$segmented_model)),
       type = "orig") %>%
  rbind(tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
        pred = exp(predict(cp_res_bav_full1$cp_segmented_list_backpro[[which.min(cp_res_bav_full1$bic_backpro)]]$segmented_model)),
        type = "sens 1")) %>%
  rbind(tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
               pred = exp(predict(cp_res_bav_full2$cp_segmented_list_backpro[[which.min(cp_res_bav_full2$bic_backpro)]]$segmented_model)),
               type = "sens 2")) %>%
  ggplot() + geom_line(aes(date, pred, lty = type)) +
  theme(legend.title = element_blank()) + 
  ylab("Number") + xlab("Date") + ggsave(filename = "../results/changepoint/cp_results/fig_sens_bav.png", width = 6, height = 6)

load("../results/changepoint/cp_results/cp_res_ger_sens.RData")
cp_res_ger_full1 = cp_res_ger_full
load("../results/changepoint/cp_results/cp_res_ger.RData")
load("../results/changepoint/cp_results/cp_res_ger_sens2.RData")

plot_ger = tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
       pred = exp(predict(cp_res_ger_full$cp_segmented_list_backpro[[which.min(cp_res_ger_full$bic_backpro)]]$segmented_model)),
       type = "orig") %>%
  rbind(tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
               pred = exp(predict(cp_res_ger_full1$cp_segmented_list_backpro[[which.min(cp_res_ger_full1$bic_backpro)]]$segmented_model)),
               type = "sens 1")) %>%
  rbind(tibble(date= seq(ymd("2020-02-28"), ymd("2020-05-01"), by="1 day"),
               pred = exp(predict(cp_res_ger_full2$cp_segmented_list_backpro[[which.min(cp_res_ger_full2$bic_backpro)]]$segmented_model)),
               type = "sens 2")) %>%
  ggplot() + geom_line(aes(date, pred, lty = type)) +
  theme(legend.title = element_blank()) + 
  ylab("Number") + xlab("Date") + ggsave(filename = "../results/changepoint/cp_results/fig_sens_ger.png", width = 6, height = 6)

theme = theme(
  axis.text=element_text(size = rel(1.3)),
  axis.title=element_text(size = rel(1.3)),
  legend.text = element_text(size = rel(1.3)),
  legend.title =element_blank())
fig_s2 = ggpubr::ggarrange(plot_bav+theme, plot_ger+theme, common.legend = T, legend = "bottom", labels = "AUTO")

ggsave(fig_s2, filename = "../results/changepoint/fig_s2_sens_impu.png", width = 12, height = 4)
