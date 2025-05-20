# library(tidyverse)
# theme_set(theme_bw())
# library(lubridate)
# library(segmented)
# library(surveillance)
# library(MASS)
# library(nlme)
# Read dataset containing diseae onsets and number of reported cases per day
onset_rep_age_group_ger = read_tsv("../data/onset_reported_age_group_ger.csv")
onset_ger = onset_rep_age_group_ger %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  group_by(date) %>% summarise(onsets=sum(onsets, na.rm = T),
                               reported=sum(reported, na.rm = T))
  
onset_rep_age_group_bav = read_tsv("../data/onset_reported_age_group_bav.csv")
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

backpro_ger = perform_backpro(onset_ger$date, onset_ger$onsets, pmf = pmf)
backpro_bav = perform_backpro(onset_bav$date, onset_bav$onsets, pmf = pmf)

# Create dataset for comparison of time-series
# Download JHU data
# jhu <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv") 
# # Restrict to Germany
# jhu_de <- jhu %>% filter(`Country/Region` == "Germany") %>% dplyr::select(-`Province/State`, -Lat, -Long) 
# # Convert to long format
# jhu_de_ts <- jhu_de %>% 
#   pivot_longer(cols=-`Country/Region`, values_to="Cum_Counts") %>% 
#   mutate(Date = mdy(name), 
#          jhu_Counts = Cum_Counts - lag(Cum_Counts),
#          # Smooth
#          jhu_Counts7 = as.numeric(stats::filter(jhu_Counts, sides=2, rep(1/7,7)))) %>%
#   dplyr::select(date=Date, jhu_smooth=jhu_Counts7)
# write.table(jhu_de_ts, file = "../data/raw/jhu_germany.csv", sep="\t", row.names = FALSE, quote = FALSE)
jhu = read_tsv("../data/raw/jhu_germany.csv")

plot_dat_ger = backpro_ger %>% left_join(onset_ger) %>% left_join(jhu) %>%
  mutate(reported_7 = as.numeric(stats::filter(reported, sides=2, rep(1/7,7)))) %>%
  pivot_longer(cols = c("backpro", "onsets", "reported", "reported_7", "jhu_smooth"), 
               values_to = "count", names_to = "type") %>%
  mutate(type = factor(type, levels = c("backpro", "onsets", "reported", "reported_7", "jhu_smooth"),
                       labels = c("Backprojection", "Onsets", "Reported Raw", "Reported (RKI)", "Reported (JHU)")))

plot_dat_bav = backpro_bav %>% left_join(onset_bav) %>% 
  mutate(reported_7 = as.numeric(stats::filter(reported, sides=2, rep(1/7,7)))) %>%
  pivot_longer(cols = c("backpro", "onsets", "reported", "reported_7"), 
               values_to = "count", names_to = "type") %>%
  mutate(type = factor(type, levels = c("backpro", "onsets", "reported", "reported_7"),
                       labels = c("Backprojection", "Onsets", "Reported Raw", "Reported")))
plot_ger = ggplot() + 
  geom_line(aes(date, count, col=type, lty=type), plot_dat_ger %>%
            filter(type %in% c("Backprojection", "Onsets", "Reported (RKI)")),
            cex=1.1) +
  geom_point(aes(date, count, col = type), plot_dat_ger %>%
               filter(type %in% c("Reported Raw")) %>%
               mutate(type="Reported (RKI)"), show.legend = FALSE, pch=3) +
  geom_line(aes(date, count, col = type), plot_dat_ger %>%
                filter(type %in% c("Reported Raw")) %>%
                mutate(type="Reported (RKI)"), cex=.2, lty = 2) +
  scale_color_manual(values = c("Backprojection"=rgb(0.1,.1,.1,1),
                                "Onsets"=rgb(.4,.4,.4,1),
                                "Reported (RKI)" = rgb(.8,.8,.8,1),
                                "Reported (JHU)" = rgb(.8,.8,.8,1))) +
  scale_linetype_manual(values = c("Backprojection"=1,
                                   "Onsets"=1,
                                   "Reported (RKI)" = 1,
                                   "Reported (JHU)" =2)) +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_blank(),
    legend.position = "bottom") + 
  ylab("Number") + 
  xlab("Date") +
  ggtitle("Germany")

# Alternative plot including also reporting data by JHU
plot_ger_jhu = ggplot() + 
  geom_line(aes(date, count, col=type, lty=type), plot_dat_ger %>%
              filter(type %in% c("Backprojection", "Onsets", "Reported (RKI)", "Reported (JHU)")),
            cex=1.1) +
  geom_point(aes(date, count, col = type), plot_dat_ger %>%
               filter(type %in% c("Reported Raw")) %>%
               mutate(type="Reported (RKI)"), show.legend = FALSE, pch=3) +
  # geom_line(aes(date, count, col = type), plot_dat_ger %>%
  #              filter(type %in% c("Reported Raw")) %>%
  #              mutate(type="Reported (RKI)"), cex=.2, lty = 2) +
  scale_color_manual(values = c("Backprojection"=rgb(0.1,.1,.1,1),
                                "Onsets"=rgb(.4,.4,.4,1),
                                "Reported (RKI)" = rgb(.8,.8,.8,1),
                                "Reported (JHU)" = rgb(.8,.8,.8,1))) +
  scale_linetype_manual(values = c("Backprojection"=1,
                                   "Onsets"=1,
                                   "Reported (RKI)" = 1,
                                   "Reported (JHU)" =2)) +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_blank(),
    legend.position = "bottom") + 
  ylab("Number") + 
  xlab("Date") +
  ggtitle("Germany")

plot_bav = ggplot() + 
  geom_line(aes(date, count, col=type, lty=type), plot_dat_bav %>%
              filter(type %in% c("Backprojection", "Onsets", "Reported")),
            cex=1.1) +
  geom_point(aes(date, count, col = type), plot_dat_bav %>%
               filter(type %in% c("Reported Raw")) %>%
               mutate(type="Reported"), show.legend = FALSE, pch = 3) +
  geom_line(aes(date, count, col = type), plot_dat_bav %>%
               filter(type %in% c("Reported Raw")) %>%
               mutate(type="Reported"), cex=.2, lty = 2) +
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
  ggtitle("Bavaria")

fig_1 = ggpubr::ggarrange(plot_bav, plot_ger, nrow=1, common.legend = F, legend = "bottom", 
                          labels = "AUTO")
ggsave(fig_1, filename = "../results/changepoint/fig1_backpro.png", width=12, height = 4)

## Sensitivity analysis Old Age
# Alternative incubation period distribution
quantiles_incu_sens = data.frame(q=c(0.5, 0.95), value=c(7.7, 14.1))
# Fit distribution
f_target_sens = function(theta) {
  qlnorm(quantiles_incu_sens$q[c(1,2)], theta[1], theta[2]) - quantiles_incu_sens$value[c(1,2)]
}
incu_lnorm_sens = nleqslv::nleqslv(c(1,1), f_target_sens)$x
# Compare observed and fitted value

data.frame(
  q        = quantiles_incu_sens$q,
  observed = quantiles_incu_sens$value,
  fitted   = qlnorm( quantiles_incu_sens$q, incu_lnorm_sens[1], 
                     incu_lnorm_sens[2])
)

inc_pdf_sens = data.frame(t = seq(0,19, length=1000)) %>%
  mutate(pdf = dlnorm(t, incu_lnorm_sens[1], incu_lnorm_sens[2]))

# Discretize incubation period distribution
cdf_sens = plnorm(seq(0,18,by=1), incu_lnorm_sens[1], incu_lnorm_sens[2])
pmf_sens = structure(c(0,diff(cdf_sens)), names=seq(0,18,by=1))
# Normalize the discrete incubation period distribution
pmf_sens = pmf_sens/sum(pmf_sens)
df_sens = data.frame(days=as.numeric(names(pmf_sens)), pmf=pmf_sens)

onset_ger_80 = onset_rep_age_group_ger %>% 
  filter(age_group=="80+") %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  group_by(date) %>% summarise(onsets=sum(onsets, na.rm = T),
                               reported=sum(reported, na.rm = T))

onset_bav_80 = onset_rep_age_group_bav %>%
  filter(age_group=="80+") %>%
  right_join(tibble(date=seq(ymd("2020-02-24"), ymd("2020-05-15"), by = "1 day"))) %>%
  group_by(date) %>% summarise(onsets=sum(onsets, na.rm = T),
                               reported=sum(reported, na.rm = T))


backpro_80_bav = perform_backpro(onset_bav_80$date, onset_bav_80$onsets, pmf = pmf_sens) %>%
  mutate(type="Sensitivity") %>% 
  rbind(perform_backpro(onset_bav_80$date, onset_bav_80$onsets, pmf = pmf) %>%
          mutate(type="Original"))

backpro_80_ger = perform_backpro(onset_ger_80$date, onset_ger_80$onsets, pmf = pmf_sens) %>%
  mutate(type="Sensitivity") %>% 
  rbind(perform_backpro(onset_ger_80$date, onset_ger_80$onsets, pmf = pmf) %>%
          mutate(type="Original"))

sens_bav_plot = ggplot(backpro_80_bav) + 
 geom_line(aes(date, backpro, lty = type)) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_blank(),
    legend.position = "bottom") + ylab("Number") +xlab("Date") +
  ggtitle("Bavaria")

sens_ger_plot = ggplot(backpro_80_ger) + 
  geom_line(aes(date, backpro, lty = type)) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size = rel(1.3)),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_blank(),
    legend.position = "bottom") + ylab("Number") +xlab("Date") +
  ggtitle("Germany")

fig_sens_backpro = ggpubr::ggarrange(sens_bav_plot, sens_ger_plot, nrow=1, common.legend = T, legend = "bottom", 
                          labels = "AUTO")
ggsave(fig_sens_backpro, filename = "../results/changepoint/fig_s3_sens_backpro.png", width=12, height = 4)
