# library(tidyverse)
# theme_set(theme_bw())
# library(lubridate)
# Load prepared data number of disease onsets per day
bav = read_tsv("../data/onset_age_group_bav.csv")
ger = read_tsv("../data/onset_age_group_ger.csv")


bav_014 = bav %>% dplyr::filter(age_group=="0-14") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
bav_1559 = bav %>% dplyr::filter(age_group=="15-59") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
bav_6079 = bav %>% dplyr::filter(age_group=="60-79") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
bav_80 = bav %>% dplyr::filter(age_group=="80+") %>%
  group_by(date) %>%
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

cp_res_bav_014 = perform_cp_analysis(data = bav_014, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "bav_014")

cp_res_bav_1559 = perform_cp_analysis(data = bav_1559, 
                                     type = "backpro",
                                     cp_max_onset = 6, 
                                     cp_max_backpro = 6,
                                     save_disc_optim_results = T, 
                                     use_disc_optim_results = T, 
                                     name_disc =  "bav_1559")
cp_res_bav_6079 = perform_cp_analysis(data = bav_6079, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "bav_6079")
cp_res_bav_80 = perform_cp_analysis(data = bav_80, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "bav_80")

res_bav014_selected_backpro = cp_res_bav_014$cp_segmented_list_backpro[[which.min(cp_res_bav_014$bic_backpro)]]
res_bav1559_selected_backpro = cp_res_bav_1559$cp_segmented_list_backpro[[which.min(cp_res_bav_1559$bic_backpro)]]
res_bav6079_selected_backpro = cp_res_bav_6079$cp_segmented_list_backpro[[which.min(cp_res_bav_6079$bic_backpro)]]
res_bav80_selected_backpro = cp_res_bav_80$cp_segmented_list_backpro[[which.min(cp_res_bav_80$bic_backpro)]]

# Table change points
res_bav014_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                   fac_014=paste0(format(round(mult_factor,2), digits = 2), 
                                                                 " [", 
                                                                 format(round(CI_lwr,2), digits = 2), 
                                                                 ", ", 
                                                                 format(round(CI_upr,2), digits = 2), 
                                                                 "]")) %>%
  dplyr::select(ID, fac_014) %>% 
  full_join(res_bav1559_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                        fac_1559=paste0(format(round(mult_factor,2), digits = 2), 
                                                                       " [", 
                                                                       format(round(CI_lwr,2), digits = 2), 
                                                                       ", ", 
                                                                       format(round(CI_upr,2), digits = 2), 
                                                                       "]")) %>%
              dplyr::select(ID, fac_1559)) %>%
  full_join(res_bav6079_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                         fac_6079=paste0(format(round(mult_factor,2), digits = 2), 
                                                                         " [", 
                                                                         format(round(CI_lwr,2), digits = 2), 
                                                                         ", ", 
                                                                         format(round(CI_upr,2), digits = 2), 
                                                                         "]")) %>%
              dplyr::select(ID, fac_6079)) %>%
  full_join(res_bav80_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                         fac_80=paste0(format(round(mult_factor,2), digits = 2), 
                                                                         " [", 
                                                                         format(round(CI_lwr,2), digits = 2), 
                                                                         ", ", 
                                                                         format(round(CI_upr,2), digits = 2), 
                                                                         "]")) %>%
              dplyr::select(ID, fac_80)) %>%
dplyr::select(-ID) %>% xtable()

res_bav014_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                   BP_014=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                 " [", 
                                                                 gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                 ", ", 
                                                                 gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                 "]")) %>%
  dplyr::select(ID, BP_014) %>%
  full_join(res_bav1559_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                                BP_1559=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                               " [", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                               ", ", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                               "]")) %>%
              dplyr::select(ID, BP_1559)) %>%
  full_join(res_bav6079_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                                BP_6079=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                               " [", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                               ", ", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                               "]")) %>%
              dplyr::select(ID, BP_6079)) %>%
  full_join(res_bav80_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                              BP_80=paste0(gsub("2020-", "", gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T))), 
                                                                           " [", 
                                                                           gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                           ", ", 
                                                                           gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                           "]")) %>%
              dplyr::select(ID, BP_80)) %>% dplyr::select(-ID) %>% xtable()

# Age-Distribution Bavaria
age_bav = read.csv("../data/age_dist_bav.csv")
age_bav_group = age_bav %>% mutate(age_group=cut(Age, c(-1,14,59,79,120), labels = c("0-14", "15-59", "60-79", "80+"))) %>%
  group_by(age_group) %>%
  summarise(n=sum(Num_191231))

plot_bav_age_group = res_bav014_selected_backpro$segmented_model$model %>% 
  dplyr::select(t, logbackpro) %>%
  mutate(age_group="0-14") %>%
  cbind(pred=res_bav014_selected_backpro$segmented_model$fitted) %>%
  rbind(res_bav1559_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="15-59") %>%
          cbind(pred=res_bav1559_selected_backpro$segmented_model$fitted)) %>%
  rbind(res_bav6079_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="60-79") %>%
          cbind(pred=res_bav6079_selected_backpro$segmented_model$fitted)) %>%
  rbind(res_bav80_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="80+") %>%
          cbind(pred=res_bav80_selected_backpro$segmented_model$fitted)) %>%
  mutate(t=ymd("2020-02-27")+t) %>%
  right_join(age_bav_group) %>%
  mutate(pred_per_100k = (exp(pred)/n)*100000) %>%
  ggplot() +
  #geom_line(aes(t, exp(logbackpro), col=age_group), lty=2) +
  geom_line(aes(t, pred_per_100k, lty=age_group)) +
  ylab("Number infections\n per 100k") + xlab("Date") +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.3)),
    legend.title =element_text(size = rel(1.3)))+
  guides(lty=guide_legend(title="Age group"))

# Figure 2
load("../results/changepoint/cp_results/cp_res_bav.RData")
theme = theme(
  axis.text=element_text(size = rel(1.3)),
  axis.title=element_text(size = rel(1.3)),
  legend.text = element_text(size = rel(1.3)),
  legend.title =element_text(size = rel(1.3)))
fig_2_bav_ov_age_grp = ggpubr::ggarrange(cp_res_bav_full$cp_segmented_list_backpro[[which.min(cp_res_bav_full$bic_backpro)]]$plot +
                    ggtitle("Overall") + theme + ylab("Number infections"),
                  plot_bav_age_group + ggtitle("Age group") + theme, labels = "AUTO")
ggsave(fig_2_bav_ov_age_grp, filename = "../results/changepoint/fig2_results_bav_backpro.png", width=12, height = 4)

## Germany
ger_014 = ger %>% dplyr::filter(age_group=="0-14") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
ger_1559 = ger %>% dplyr::filter(age_group=="15-59") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
ger_6079 = ger %>% dplyr::filter(age_group=="60-79") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))
ger_80 = ger %>% dplyr::filter(age_group=="80+") %>%
  group_by(date) %>%
  summarise(onsets=sum(onsets)) %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  arrange(date) %>%
  mutate(onsets=replace_na(onsets,0))


cp_res_ger_014 = perform_cp_analysis(data = ger_014, 
                                     type = "backpro",
                                     cp_max_onset = 6, 
                                     cp_max_backpro = 6,
                                     save_disc_optim_results = T, 
                                     use_disc_optim_results = T, 
                                     name_disc =  "ger_014")

cp_res_ger_1559 = perform_cp_analysis(data = ger_1559, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "ger_1559")
cp_res_ger_6079 = perform_cp_analysis(data = ger_6079, 
                                      type = "backpro",
                                      cp_max_onset = 6, 
                                      cp_max_backpro = 6,
                                      save_disc_optim_results = T, 
                                      use_disc_optim_results = T, 
                                      name_disc =  "ger_6079")
cp_res_ger_80 = perform_cp_analysis(data = ger_80, 
                                    type = "backpro",
                                    cp_max_onset = 6, 
                                    cp_max_backpro = 6,
                                    save_disc_optim_results = T, 
                                    use_disc_optim_results = T, 
                                    name_disc =  "ger_80")

res_ger014_selected_backpro = cp_res_ger_014$cp_segmented_list_backpro[[which.min(cp_res_ger_014$bic_backpro)]]
res_ger1559_selected_backpro = cp_res_ger_1559$cp_segmented_list_backpro[[which.min(cp_res_ger_1559$bic_backpro)]]
res_ger6079_selected_backpro = cp_res_ger_6079$cp_segmented_list_backpro[[which.min(cp_res_ger_6079$bic_backpro)]]
res_ger80_selected_backpro = cp_res_ger_80$cp_segmented_list_backpro[[which.min(cp_res_ger_80$bic_backpro)]]

# Table change points
res_ger014_selected_backpro$coef %>% mutate(ID=1:n(), 
                                            fac_014=paste0(format(round(mult_factor,2), digits = 2), 
                                                           " [", 
                                                           format(round(CI_lwr,2), digits = 2), 
                                                           ", ", 
                                                           format(round(CI_upr,2), digits = 2), 
                                                           "]")) %>%
  dplyr::select(ID, fac_014) %>% 
  full_join(res_ger1559_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                         fac_1559=paste0(format(round(mult_factor,2), digits = 2), 
                                                                         " [", 
                                                                         format(round(CI_lwr,2), digits = 2), 
                                                                         ", ", 
                                                                         format(round(CI_upr,2), digits = 2), 
                                                                         "]")) %>%
              dplyr::select(ID, fac_1559)) %>%
  full_join(res_ger6079_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                         fac_6079=paste0(format(round(mult_factor,2), digits = 2), 
                                                                         " [", 
                                                                         format(round(CI_lwr,2), digits = 2), 
                                                                         ", ", 
                                                                         format(round(CI_upr,2), digits = 2), 
                                                                         "]")) %>%
              dplyr::select(ID, fac_6079)) %>%
  full_join(res_ger80_selected_backpro$coef %>% mutate(ID=1:n(), 
                                                       fac_80=paste0(format(round(mult_factor,2), digits = 2), 
                                                                     " [", 
                                                                     format(round(CI_lwr,2), digits = 2), 
                                                                     ", ", 
                                                                     format(round(CI_upr,2), digits = 2), 
                                                                     "]")) %>%
              dplyr::select(ID, fac_80)) %>%
  dplyr::select(-ID) %>% xtable()

res_ger014_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                   BP_014=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                 " [", 
                                                                 gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                 ", ", 
                                                                 gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                 "]")) %>%
  dplyr::select(ID, BP_014) %>%
  full_join(res_ger1559_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                                BP_1559=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                               " [", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                               ", ", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                               "]")) %>%
              dplyr::select(ID, BP_1559)) %>%
  full_join(res_ger6079_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                                BP_6079=paste0(gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T)), 
                                                                               " [", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                               ", ", 
                                                                               gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                               "]")) %>%
              dplyr::select(ID, BP_6079)) %>%
  full_join(res_ger80_selected_backpro$breakpoints %>% mutate(ID=1:n(), 
                                                              BP_80=paste0(gsub("2020-", "", gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP, perl=T))), 
                                                                           " [", 
                                                                           gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_lwr, perl=T)), 
                                                                           ", ", 
                                                                           gsub("2020-", "", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", BP_CI_upr, perl=T)), 
                                                                           "]")) %>%
              dplyr::select(ID, BP_80)) %>% dplyr::select(-ID) %>% xtable()

# Age-Distribution Germany
age_ger = read.csv("../data/age_dist_ger.csv")
age_ger_group = age_ger %>% mutate(age_group=cut(Age, c(-1,14,59,79,120), labels = c("0-14", "15-59", "60-79", "80+"))) %>%
  group_by(age_group) %>%
  summarise(n=sum(Num_191231))

plot_ger_age_group = res_ger014_selected_backpro$segmented_model$model %>% 
  dplyr::select(t, logbackpro) %>%
  mutate(age_group="0-14") %>%
  cbind(pred=res_ger014_selected_backpro$segmented_model$fitted) %>%
  rbind(res_ger1559_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="15-59") %>%
          cbind(pred=res_ger1559_selected_backpro$segmented_model$fitted)) %>%
  rbind(res_ger6079_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="60-79") %>%
          cbind(pred=res_ger6079_selected_backpro$segmented_model$fitted)) %>%
  rbind(res_ger80_selected_backpro$segmented_model$model %>% 
          dplyr::select(t, logbackpro) %>%
          mutate(age_group="80+") %>%
          cbind(pred=res_ger80_selected_backpro$segmented_model$fitted)) %>%
  mutate(t=ymd("2020-02-27")+t) %>%
  right_join(age_ger_group) %>%
  mutate(pred_per_100k = (exp(pred)/n)*100000) %>%
  ggplot() +
  #geom_line(aes(t, exp(logbackpro), col=age_group), lty=2) +
  geom_line(aes(t, pred_per_100k, lty=age_group)) +
  ylab("Number infections\n per 100k") + xlab("Date") +
  theme(
    axis.text=element_text(size = rel(1.5)),
    axis.title=element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))) +
  guides(lty=guide_legend(title="Age group"))

# Figure 3
load("../results/changepoint/cp_results/cp_res_ger.RData")
theme = theme(
  axis.text=element_text(size = rel(1.3)),
  axis.title=element_text(size = rel(1.3)),
  legend.text = element_text(size = rel(1.3)),
  legend.title =element_text(size = rel(1.3)))
fig_3_ger_ov_age_grp = ggpubr::ggarrange(cp_res_ger_full$cp_segmented_list_backpro[[which.min(cp_res_ger_full$bic_backpro)]]$plot +
                                           theme + ylab("Number infections") +
                                           ggtitle("Overall") ,
                                         plot_ger_age_group + theme + ggtitle("Age group") , labels = "AUTO")
ggsave(fig_3_ger_ov_age_grp, filename = "../results/changepoint/fig3_results_ger_backpro.png", width=12, height = 4)
