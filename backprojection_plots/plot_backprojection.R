# plot back-projection

rm(list = ls())
library(tidyverse)
theme_set(theme_bw())
library(lubridate)
library(segmented)
library(surveillance)
library(MASS)
library(nlme)
#library(ggpubr)
library("nleqslv")
##################
### Full year#####
##################
# Read dataset containing diseae onsets and number of reported cases per day


onset_rep_age_group_bav = read_csv("data/age_group_data.csv", show_col_types = FALSE)
onset_bav = onset_rep_age_group_bav %>%
 # right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  group_by(date) %>% summarise(onsets=sum(onsets, na.rm = T),
                               reported=sum(reported, na.rm = T))
#' Backprojection of data 
#' We take a literature based approach to deduce an incubation time distribution.
#' Lauer et al. (2020) - log normal distribution - same as in Dehning et al. (2020)
#' Source: [Lauer et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32150748).
quantiles_incu = data.frame(q=c(0.5, 0.975), value=c(5.1, 11.5))
# Fit distribution
f_target2 = function(theta) {
  qlnorm(quantiles_incu$q[c(1,2)], theta[1], theta[2]) - quantiles_incu$value[c(1,2)]
}
incu_lnorm = nleqslv::nleqslv(c(1,1), f_target2)$x
# Compare observed and fitted value
data.frame(
  q        = quantiles_incu$q,
  observed = quantiles_incu$value,
  fitted   = qlnorm( quantiles_incu$q, incu_lnorm[1], incu_lnorm[2])
)

inc_pdf = data.frame(t = seq(0,15, length=1000)) %>%
  mutate(pdf = dlnorm(t, incu_lnorm[1], incu_lnorm[2]))

# Discretize incubation period distribution
cdf = plnorm(seq(0,14,by=1), incu_lnorm[1], incu_lnorm[2])
pmf = structure(c(0,diff(cdf)), names=seq(0,14,by=1))
# Normalize the discrete incubation period distribution
pmf = pmf/sum(pmf)
df = data.frame(days=as.numeric(names(pmf)), pmf=pmf)

#' The backprojection can be done using the function `surveillance::backprojNP`
#' The backprojected curve shows the expected number of infections per day and can be
#' compared to interventions similar to Werber et al. (2013) (https://doi.org/10.1093/aje/kwt069)

perform_backpro = function(date, onsets, pmf) {
  # Onset data
  sts_symp = sts(
    epoch       = date,
    observed    = matrix(onsets, ncol = 1, nrow = length(onsets)),
    epochAsDate = TRUE)
  
  # Perform back projection with smoothing to adjust for weekday effects (k=6)
  bp = backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
  dat_backpro = bp %>% as.data.frame() %>%
    dplyr::select(date=epoch,
                  backpro = upperbound) %>%
    dplyr::filter(date>=ymd("2020-02-28"),
                  date<=max(date)-(length(pmf)-1))
  dat_backpro
}

backpro_bav = perform_backpro(onset_bav$date, onset_bav$onsets, pmf = pmf)


plot_dat_bav = backpro_bav %>% left_join(onset_bav) %>% 
  mutate(reported_7 = as.numeric(stats::filter(reported, sides=2, rep(1/7,7)))) %>%
  pivot_longer(cols = c("backpro", "onsets", "reported", "reported_7"), 
               values_to = "count", names_to = "type") %>%
  mutate(type = factor(type, levels = c("backpro", "onsets", "reported", "reported_7"),
                       labels = c("Backprojection", "Onsets", "Reported Raw", "Reported")))


plot_bav = ggplot() + 
  geom_line(aes(date, count, col=type, lty=type), plot_dat_bav %>%
              filter(type %in% c("Backprojection", "Onsets", "Reported")),
            cex=1.1) +
  # geom_point(aes(date, count, col = type), plot_dat_bav %>%
  #              filter(type %in% c("Reported Raw")) %>%
  #              mutate(type="Reported"), show.legend = FALSE, pch = 3) +
  # geom_line(aes(date, count, col = type), plot_dat_bav %>%
  #             filter(type %in% c("Reported Raw")) %>%
  #             mutate(type="Reported"), cex=.2, lty = 2) +
  scale_color_manual(values = c("Backprojection"=rgb(0.1,.1,.1,1),
                                "Onsets"=rgb(.4,.4,.4,1),
                                "Reported" = rgb(.8,.8,.8,1))) +
  scale_linetype_manual(values = c("Backprojection"=1,
                                   "Onsets"=1,
                                   "Reported" = 1)) +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.1)),
    legend.title = element_blank(),
    legend.position = "bottom") + 
  ylab("Number") + 
  xlab("Date") +
  ggtitle("Mexico")

ggsave(
  "Plot_data_full.png",
  plot = plot_bav
)

##################
#March and April##
##################
onset_bav = onset_rep_age_group_bav %>%
   right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  group_by(date) %>% summarise(onsets=sum(onsets, na.rm = T),
                               reported=sum(reported, na.rm = T))
#' Backprojection of data 
#' We take a literature based approach to deduce an incubation time distribution.
#' Lauer et al. (2020) - log normal distribution - same as in Dehning et al. (2020)
#' Source: [Lauer et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32150748).
quantiles_incu = data.frame(q=c(0.5, 0.975), value=c(5.1, 11.5))
# Fit distribution
f_target2 = function(theta) {
  qlnorm(quantiles_incu$q[c(1,2)], theta[1], theta[2]) - quantiles_incu$value[c(1,2)]
}
incu_lnorm = nleqslv::nleqslv(c(1,1), f_target2)$x
# Compare observed and fitted value
data.frame(
  q        = quantiles_incu$q,
  observed = quantiles_incu$value,
  fitted   = qlnorm( quantiles_incu$q, incu_lnorm[1], incu_lnorm[2])
)

inc_pdf = data.frame(t = seq(0,15, length=1000)) %>%
  mutate(pdf = dlnorm(t, incu_lnorm[1], incu_lnorm[2]))

# Discretize incubation period distribution
cdf = plnorm(seq(0,14,by=1), incu_lnorm[1], incu_lnorm[2])
pmf = structure(c(0,diff(cdf)), names=seq(0,14,by=1))
# Normalize the discrete incubation period distribution
pmf = pmf/sum(pmf)
df = data.frame(days=as.numeric(names(pmf)), pmf=pmf)

#' The backprojection can be done using the function `surveillance::backprojNP`
#' The backprojected curve shows the expected number of infections per day and can be
#' compared to interventions similar to Werber et al. (2013) (https://doi.org/10.1093/aje/kwt069)

perform_backpro = function(date, onsets, pmf) {
  # Onset data
  sts_symp = sts(
    epoch       = date,
    observed    = matrix(onsets, ncol = 1, nrow = length(onsets)),
    epochAsDate = TRUE)
  
  # Perform back projection with smoothing to adjust for weekday effects (k=6)
  bp = backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
  dat_backpro = bp %>% as.data.frame() %>%
    dplyr::select(date=epoch,
                  backpro = upperbound) %>%
    dplyr::filter(date>=ymd("2020-02-28"),
                  date<=max(date)-(length(pmf)-1))
  dat_backpro
}

backpro_bav = perform_backpro(onset_bav$date, onset_bav$onsets, pmf = pmf)


plot_dat_bav = backpro_bav %>% left_join(onset_bav) %>% 
  mutate(reported_7 = as.numeric(stats::filter(reported, sides=2, rep(1/7,7)))) %>%
  pivot_longer(cols = c("backpro", "onsets", "reported", "reported_7"), 
               values_to = "count", names_to = "type") %>%
  mutate(type = factor(type, levels = c("backpro", "onsets", "reported", "reported_7"),
                       labels = c("Backprojection", "Onsets", "Reported Raw", "Reported")))

plot_bav2 = ggplot() + 
  geom_line(aes(date, count, col=type, lty=type), plot_dat_bav %>%
              filter(type %in% c("Backprojection", "Onsets", "Reported")),
            cex=1.1) +
  scale_color_manual(values = c("Backprojection"=rgb(0.1,.1,.1,1),
                                "Onsets"=rgb(.4,.4,.4,1),
                                "Reported" = rgb(.8,.8,.8,1))) +
  scale_linetype_manual(values = c("Backprojection"=1,
                                   "Onsets"=1,
                                   "Reported" = 1)) +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.1)),
    legend.title = element_blank(),
    legend.position = "bottom") + 
  ylab("Number") + 
  xlab("Date") +
  ggtitle("Mexico")

ggsave("Plot_2_Month.png",  plot = plot_bav2)
