#' Function to estimate a piecewise-linear model for the analysis of backprojected
#' new COVID-19 infections per day. We model the logarithmic number of expected infections
#' after the backprojection based on a normal distribution with AR-1 error structure.
#' Breakpoints are found based on discrete optimization
#' over all potential combinations of breakpoints with at least 3 days in between 2 breakpoints
#' @param data dataset of subsequent days including a column 'date' and a column 'backpro' which contains
#' a backprojection for the estimated number of new infections per day
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @return list with entries: 'mod_min_dev': the fitted model with minimum deviance; 'bps': the breakpoints
#' of the model with minimum deviance (in days), 'deviance': the deviance of the model with minimum deviance;
#' 'overdispersion': overdispersion of model with minimum deviance; 'res_disc_opt': tibble with deviance and overdispersion
#' for each combination of breakpoints
estimate_bp_disc_optim_backpro = function(data, bp=3) {
  dat = data %>% mutate(t=1:n())
  day_vec = lapply(1:bp, function(x) 3:(nrow(data)-3))
  bp_scen = tibble(bp_1 = day_vec[[1]])
  for(i in 2:bp) {
    bp_scen = expand_grid(bp_scen, !!paste0("bp_", i) := day_vec[[i]]) %>%
      dplyr::filter(!!as.symbol(paste0("bp_", i)) > !!as.symbol(paste0("bp_", i-1)),
                    (!!as.symbol(paste0("bp_", i)) - !!as.symbol(paste0("bp_", i-1))) >3)
  }
  create_data_est_model = function(bp_vec, data = dat) {
    # Create data matrix, linear trends between changepoints
    dat_mod = cbind(data %>% dplyr::select(backpro, t),
                    do.call(cbind, lapply(1:length(bp_vec),
                                          function(x) pmax(0, data$t - bp_vec[x] + 1))))
    colnames(dat_mod) = c("backpro", "t", paste0("t", 1:length(bp_vec)))

    # Estimate model
    mod = tryCatch(gls(log(backpro)~., data = dat_mod, correlation=corAR1(form= ~1),method="ML"), error = NULL)
    # Out: changepoints and corresponding deviance
    if (is.null(mod) || is.null(mod$deviance)) {
      manual_deviance = NA
      } else {
        manual_deviance = mod$deviance
        }
    if (is.null(mod) || is.null(mod$dispersion)) {
      manual_dispersion = NA
      } else {
        manual_dispersion = summary(mod)$dispersion
        }
    c(bp_vec, deviance = manual_deviance,
      dispersion = manual_dispersion)
  }
  res = future_apply(bp_scen, 1, create_data_est_model)
  # res = apply(bp_scen, 1, create_data_est_model)
  res = as_tibble(t(res))
  # Get Breakpoints with lowest deviance
  bp_ts = res %>% dplyr::filter(rank(deviance) == 1)
  dat_mod = cbind(dat %>% dplyr::select(backpro, t),
                  do.call(cbind, lapply(1:bp,
                                        function(x) pmax(0, dat$t - unlist(bp_ts[1,1:bp])[x] + 1))))
  colnames(dat_mod) = c("backpro", "t", paste0("t", 1:bp))

  # Re-estimate model
  mod = gls(log(backpro)~., data = dat_mod, correlation=corAR1(form= ~1),method="ML")

  list(mod_min_dev = mod,
       bps = unlist(bp_ts[1:bp]),
       deviance = mod$deviance,
       overdispersion = summary(mod)$dispersion,
       res_disc_opt = res)
}

#' Function to estimate a piecewise-linear model for the analysis of backprojected COVID-19 infections per
#' day using the 'segmented' function of the 'segmented' package. We model the logarithmic number of expected infections
#' after the backprojection based on a normal distribution with AR-1 error structure.
#' @param data dataset of subsequent days including a column 'date' and a column 'backpro' which contains
#' a backprojection for the estimated expected number of new infections per day
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @param start_bp starting values for the estimated breakpoints (might be derived from discrete optimization), length
#' of the vector has to correspond with the 'bp' parameter
#' @param plot_main desired title of the computed ggplot object
#' @param n_boot number of bootstrap samples in the bootstrap restarting algorithm of the 'segmented' function. Default to 50.
#' Bootstrapping can help to escape local minima of the objective function. Depending on the data situation, optimization can
#' fail to converge for few or several of the bootstrap samples. Please check resulting warning messages and results
#' carefully and potentially adjust the number of breakpoints. Changing starting values and seed can help to check stability
#' of results.
#' @param segmented_seed seed used in the 'segmented' function for reproducibility of (bootstrap) results.
#' @return list with entries: 'segmented_model': the result object of the 'segmented' funtion,
#' 'coef': tibble with information on multiplicative effects on daily new cases in the segments of the
#' estimated epidemic curve, 'breakpoints': tibble of the estimated breakpoints, 'plot': ggplot
#' of the estimated epidemic curve


estimate_bp_segmented_backpro = function(
  data,
  bp = 3,
  start_bp = c(18,23,54),
  n_boot = 50,
  segmented_seed = 1135235) {

  stopifnot(length(start_bp)==bp)
  dat = data %>% mutate(t=1:n(),
                        logbackpro = log(backpro))

  mod1 = gls(logbackpro~t, data=dat, correlation=corAR1(form= ~1), method="ML")
  segmented = segmented(
    mod1, seg.Z=~t, psi = start_bp,
    control = seg.control(n.boot = n_boot, it.max = 1000, seed = segmented_seed))

  # Results
  # get regression coefficients
  gamma = summary(segmented)$tTable[2:(2+bp),1]
  # linear effect on log scale per interval
  st = cumsum(gamma)
  # get varcov of gammas
  m = vcov(segmented)[2:(2+bp),2:(2+bp)]
  # derive variance of linear combinations / linear effect per interval
  lin_comb = lower.tri(matrix(rep(0,length(st)^2), ncol=length(st)), diag = T)
  var = diag(lin_comb %*% m %*% t(lin_comb))
  # multiplicative effects on cases per day and 95%-CI
  mod_res = tibble(mult_factor = exp(st),
                   CI_lwr = exp(st - 2*sqrt(var)),
                   CI_upr = exp(st+2*sqrt(var)))
  # breakpoints
  break_points = as_tibble(confint.segmented(segmented)) %>%
    rename(bp=Est., bp_lwr=`CI(95%).low`, bp_upr = `CI(95%).up`) %>%
    mutate(bp_date = min(dat$date) -1 + round(bp),
           bp_lwr_date = min(dat$date)-1 + floor(bp_lwr),
           bp_upr_date = min(dat$date)-1 + ceiling(bp_upr))

  bp_dates_res = break_points %>%
    mutate(BP =
             paste0(as.character(round(bp,1)), " (", as.character(bp_date), ")"),
           BP_CI_lwr =
             paste0(as.character(round(bp_lwr,1)), " (", as.character(bp_lwr_date), ")"),
           BP_CI_upr =
             paste0(as.character(round(bp_upr,1)), " (", as.character(bp_upr_date), ")")) %>%
    dplyr::select(BP, BP_CI_lwr, BP_CI_upr)

  dat = dat %>% mutate(pred_seg = segmented$fitted.values)
  bp_plot = ggplot() +
    geom_col(data = dat, aes(x=t, y=backpro),
             col = "grey", fill = "lightgrey") +
    geom_rect(aes(xmin=floor(break_points$bp_lwr) - 0.49,
                  xmax=ceiling(break_points$bp_upr) + 0.49, ymin=0, ymax=Inf), fill = "steelblue", alpha = 0.25)+
    geom_line(aes(x=seq(1,max(dat$t), 0.1),
                  y=exp(predict_segmented_gls(seq(1, 64, by=.1), bpts=break_points$bp,
                    beta0=coef(segmented)[1], gamma=gamma))),
              col = "black", lwd = 1.1) +
    geom_vline(xintercept=break_points$bp, lty=2, col ="steelblue", size=0.8) +
    theme_bw()+
    scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -14),
                       labels = strftime(seq(max(dat$date), min(dat$date), by = "-2 weeks"),
                                         format = "%d.%m.")) +
    scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    labs(x = "Date", y = "Number of infections (expected)") +
    theme(
      axis.text=element_text(size = rel(1.5)),
      axis.title=element_text(size = rel(1.9)),
      legend.text = element_text(size = rel(1.5)),
      legend.title =element_text(size = rel(1.5)),
      legend.position = "bottom")
  list(segmented_model = segmented, coef=mod_res, breakpoints = bp_dates_res, plot = bp_plot)
}


#' Function to estimate a piecewise-linear negative-binomial model for the analysis
#' of new COVID-19 disease onsets per day. Breakpoints are found based on discrete optimization
#' over all potential combinations of breakpoints with at least 3 days in between 2 breakpoints
#' @param data dataset of subsequent days including a column 'date' and a column 'onset' which contains
#' the number of reported disease onsets per day (epidemic curve)
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @return list with entries: 'mod_min_dev': the fitted negative-binomial GLM with minimum deviance; 'bps': the breakpoints
#' of the model with minimum deviance (in days), 'deviance': the deviance of the model with minimum deviance;
#' 'overdispersion': overdispersion of model with minimum deviance; 'res_disc_opt': tibble with deviance and overdispersion
#' for each combination of breakpoints
estimate_bp_disc_optim_onset = function(data, bp=3) {
  dat = data %>% mutate(t=1:n())
  day_vec = lapply(1:bp, function(x) 3:(nrow(data)-3))
  bp_scen = tibble(bp_1 = day_vec[[1]])
  for(i in 2:bp) {
    bp_scen = expand_grid(bp_scen, !!paste0("bp_", i) := day_vec[[i]]) %>%
      dplyr::filter(!!as.symbol(paste0("bp_", i)) > !!as.symbol(paste0("bp_", i-1)),
                    (!!as.symbol(paste0("bp_", i)) - !!as.symbol(paste0("bp_", i-1))) >3)
  }
  create_data_est_model = function(bp_vec, data = dat) {tryCatch({
    dat_mod = cbind(data %>% dplyr::select(backpro, t),
                    do.call(cbind, lapply(1:length(bp_vec),
                                          function(x) pmax(0, data$t - bp_vec[x] + 1))))
    colnames(dat_mod) = c("backpro", "t", paste0("t", 1:length(bp_vec)))

    mod = gls(log(backpro)~., data = dat_mod, correlation=corAR1(form= ~1),method="ML")

    c(bp_vec, deviance = mod$deviance, dispersion = summary(mod)$dispersion)
  }, error = function(e) {
    warning(sprintf("Model failed for bp_vec = %s: %s", paste(bp_vec, collapse = ","), e$message))
    c(rep(NA, length(bp_vec)), deviance = NA, dispersion = NA)
  })
}
  res = future_apply(bp_scen, 1, create_data_est_model)
  res = as_tibble(t(res))
  # Get Breakpoints with lowest deviance
  bp_ts = res %>% dplyr::filter(rank(deviance) == 1)
  dat_mod = cbind(dat %>% dplyr::select(onsets, t),
                  do.call(cbind, lapply(1:bp,
                                        function(x) pmax(0, dat$t - unlist(bp_ts[1,1:bp])[x] + 1))))
  colnames(dat_mod) = c("onsets", "t", paste0("t", 1:bp))
  # Re-estimate model
  mod = glm.nb(onsets~., data = dat_mod)

  list(mod_min_dev = mod, bps = unlist(bp_ts[1:bp]),
       deviance = mod$deviance, overdispersion = mod$theta, res_disc_opt = res)
}

#' Function to estimate a piecewise-linear model negative-binomial model for the analysis
#' of new COVID-19 disease onsets per day using the 'segmented' function of the 'segmented' package.
#' @param data dataset of subsequent days including a column 'date' and a column 'onset' which contains
#' the number of disease onsets per day
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @param start_bp starting values for the estimated breakpoints (might be derived from discrete optimization), length
#' of the vector has to correspond with the 'bp' parameter
#' @param n_boot number of bootstrap samples in the bootstrap restarting algorithm of the 'segmented' function. Default to 50.
#' Bootstrapping can help to escape local minima of the objective function. Depending on the data situation, optimization can
#' fail to converge for few or several of the bootstrap samples. Please check resulting warning messages and results
#' carefully and potentially adjust the number of breakpoints. Changing starting values and seed can help to check stability
#' of results.
#' @param segmented_seed seed used in the 'segmented' function for reproducibility of (bootstrap) results.
#' @return list with entries: 'segmented_model': the result object of the 'segmented' funtion,
#' 'coef': tibble with information on multiplicative effects on daily new cases in the segments of the
#' estimated epidemic curve, 'breakpoints': tibble of the estimated breakpoints, 'plot': ggplot
#' of the estimated epidemic curve


estimate_bp_segmented_onset = function(data, bp = 3, start_bp = c(18,23,54),
                                 n_boot = 50,
                                 segmented_seed = 1135235) {
  stopifnot(length(start_bp)==bp)
  dat = data %>% mutate(t=1:n())
  mod1 = glm.nb(onsets~t, data=dat)
  segmented = segmented(mod1, seg.Z=~t,
                        psi = start_bp,
                        control = seg.control(n.boot = n_boot, it.max = 1000, seed = segmented_seed))

  # Results
  # get regression coefficients
  gamma = summary(segmented)$coefficients[2:(2+bp),1]
  # linear effect on log scale per interval
  st = cumsum(gamma)
  # get varcov of gammas
  m = vcov(segmented)[2:(2+bp),2:(2+bp)]
  # derive variance of linear combinations / linear effect per interval
  lin_comb = lower.tri(matrix(rep(0,length(st)^2), ncol=length(st)), diag = T)
  var = diag(lin_comb %*% m %*% t(lin_comb))
  # multiplicative effects on cases per day and 95%-CI
  mod_res = tibble(mult_factor = exp(st), CI_lwr = exp(st - 2*sqrt(var)),
                   CI_upr = exp(st+2*sqrt(var)))
  # breakpoints
  break_points <- as_tibble(confint(segmented)) %>%
    rename(bp=Est., bp_lwr=`CI(95%).low`, bp_upr = `CI(95%).up`) %>%
    mutate(bp_date = min(dat$date) -1 + round(bp),
           bp_lwr_date = min(dat$date)-1 + floor(bp_lwr),
           bp_upr_date = min(dat$date)-1 + ceiling(bp_upr))

  bp_dates_infect = break_points %>% mutate(bp_inf = bp_date - 5,
                                            bp_lwr_inf = bp_lwr_date-5,
                                            bp_upr_inf = bp_upr_date-5)
  bp_dates_res = bp_dates_infect %>%
    mutate(BP =
             paste0(as.character(round(bp,1)), " (", as.character(bp_date), ")"),
           BP_CI_lwr =
             paste0(as.character(round(bp_lwr,1)), " (", as.character(bp_lwr_date), ")"),
           BP_CI_upr =
             paste0(as.character(round(bp_upr,1)), " (", as.character(bp_upr_date), ")")) %>%
    dplyr::select(BP, BP_CI_lwr, BP_CI_upr)

  dat = dat %>% mutate(pred_seg = segmented$fitted.values)
  bp_plot = ggplot() +
    geom_col(data = dat, aes(x=t, y=onsets),
             col = "grey", fill = "lightgrey") +
    geom_rect(aes(xmin=floor(break_points$bp_lwr) - 0.49,
                  xmax=ceiling(break_points$bp_upr) + 0.49, ymin=0, ymax=Inf), fill = "steelblue", alpha = 0.25)+
    geom_line(aes(x=seq(1,max(dat$t), 0.1),
                  y=predict(segmented, newdata = data.frame(t=seq(1,max(dat$t), 0.1)),
                            type = "response")),
              col = "black", lwd = 1.1) +
    geom_vline(xintercept=break_points$bp, lty=2, col ="steelblue", size=0.8) +
    theme_bw()+
    scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -14),
                       labels = strftime(seq(max(dat$date), min(dat$date), by = "-2 weeks"),
                                         format = "%d.%m.")) +
    scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    labs(x = "Date", y = "Number of disease onsets") +
    theme(
      axis.text=element_text(size = rel(1.5)),
      axis.title=element_text(size = rel(1.9)),
      legend.text = element_text(size = rel(1.5)),
      legend.title =element_text(size = rel(1.5)),
      legend.position = "bottom")

  list(segmented_model = segmented, coef=mod_res, breakpoints = bp_dates_res, plot = bp_plot)
}

#' Function that predicts the segmented regression based on a gls object.
#' @param timepoints The fitted gls object returned by the segmented package, which is not an object of class gls.
#' @param bp The estimated break points
#' @param beta0 the intercept of the model
#' @param st Vector of estimated slope coefficients in each segment (on the log scale)

predict_segmented_gls <- function(
  timepoints,
  bpts,
  beta0,
  gamma
) {

  bpts <- c(round(bpts, 1), max(timepoints))
  res <- beta0 + gamma[1] * timepoints
  for(i in seq_along(bpts[-1])) {
    ind_segment <- timepoints >= bpts[i]
    res[ind_segment] <- res[ind_segment] + (timepoints[ind_segment]-bpts[i])*gamma[i+1]
  }

  res

}

#' Function to perform change point analysis on backprojected daily expected numbers of infections and
#' daily disease onsets based on piecewise linear models. The function is based on three analysis steps:
#' 1) backprojection of number of disease onsets, 2) estimation of change points for backprojected number of infections,
#' 3) estimation of change points for disease onsets. Estimation of change-points is performed based on discrete optimisation
#' for two and three change points and subsequent estimation via the segmented::segmented function using results of the discrete 
#' optimisation a starting values. For change points >3 the additional starting values of the change points are selected based 
#' on a heuristical approach at the end of the observations window (with minimum distance >3 days to other CPs)
#' @param data dataset of subsequent days including a column 'date' and a column 'onsets' which contains
#' the number of disase onsets per day
#' @param type what kind of change point analysis to perform: either 'both' for the analysis of the backprojected
#' number of infections and the number of diseases onsets per day, or "backpro" or "onset" for the respective counts
#' @param cp_max_onset maximum number of change points for disease onset model, integer between 2 and 6 (larger numbers than 3 require long computation time and available memory)
#' @param cp_max_backpro number of change points for backprojected infections model, integer between 2 and 6 (larger numbers than 3 require long computation time and available memory)
#' @param save_disc_optim_results should the results of the discrete optimization be safed for later use?
#' @param name_disc identifier for saved results of the discrete optimization (string)
#' @param use_disc_optim_results should the results of previous discrete optimization (if available)
#' be loaded and used instead of being newly estimated?
#'
#' @return list with entries: 
#' aic_backpro: Vector of length cp_max_backpro-1, the AIC of the change point models on backprojected infections
#' bic_backpro: Vector of length cp_max_backpro-1, the BIC of the change point models on backprojected infections
#' cp_segmented_list_backpro: List contaning results of the estimated change point models on backprojected infections
#' aic_onset: as above but for onset models
#' bic_onset: as above but for onset models
#' cp_segmented_list_onset: as above but for onset models
#'
perform_cp_analysis = function(data,
                               type=c("both", "backpro", "onset"),
                               cp_max_onset=2,
                               cp_max_backpro = 2,
                               save_disc_optim_results=TRUE,
                               name_disc="",
                               use_disc_optim_results=TRUE) {
  stopifnot(cp_max_backpro>=2 & cp_max_backpro<=6)
  stopifnot(cp_max_onset>=2 & cp_max_onset<=6)

  if (type %in% c("both", "backpro")) {
    print("perform analysis of backprojected infections")
  #' ## Backprojection of the epidemic curve
  #'
  #' Non-parametric back-projection as in Becker et al. (1991). The exposure time is
  #' the relevant time scale to assess interventions.
  #'
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
    fitted   = qlnorm( quantiles_incu$q, incu_lnorm[1], incu_lnorm[2]))

  inc_pdf = data.frame(t = seq(0,15, length=1000)) %>%
    mutate(pdf = dlnorm(t, incu_lnorm[1], incu_lnorm[2]))

  # Discretize incubation period distribution
  cdf = plnorm(seq(0,14,by=1), incu_lnorm[1], incu_lnorm[2])
  pmf = structure(c(0,diff(cdf)), names=seq(0,14,by=1))
  # Normalize the discrete incubation period distribution
  pmf = pmf/sum(pmf)
  df = data.frame(days=as.numeric(names(pmf)), pmf=pmf)

  #' The backprojection can be done using the function `surveillance::backprojNP`
  #' The backprojected curve shows the number of infections per day and can be
  #' compared to interventions similar to Werber et al. (2013) (https://doi.org/10.1093/aje/kwt069)

  # Extract data from nowcast
  sts_symp = sts(
    epoch       = data$date,
    observed    = matrix(data$onsets, ncol = 1, nrow = nrow(data)),
    epochAsDate = TRUE)

  # Perform back projection with smoothing to adjust for weekday effects (k=6)
  bp = backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
  dat_backpro = bp %>% as.data.frame() %>%
    dplyr::select(date=epoch,
                  backpro = upperbound) %>%
    dplyr::filter(date>=ymd("2020-02-28"),
                  date<=max(date)-14)

  # Perform discrete optimization for backprojected numbers of infection
  if (use_disc_optim_results &
      file.exists(paste0("results/changepoint/discrete_optimization/backpro_disc_optim_",
                         name_disc,
                         ".RData"))) {
    load(paste0("results/changepoint/discrete_optimization/backpro_disc_optim_", name_disc,".RData"))
  } else {
    backpro_discrete = list()
  }
  if (is.null(backpro_discrete[["two_bp"]])) {
    print("perform discrete optimization for 2 change points infections")
    backpro_discrete[["two_bp"]] = estimate_bp_disc_optim_backpro(data=dat_backpro, bp=2)
  }
  if (cp_max_backpro>2 & is.null(backpro_discrete[["three_bp"]])) {
    print("perform discrete optimization for 3 change points infections")
    backpro_discrete[["three_bp"]] = estimate_bp_disc_optim_backpro(data=dat_backpro, bp=3)
  }
  if(save_disc_optim_results)
    save(backpro_discrete, file = paste0("results/changepoint/discrete_optimization/backpro_disc_optim_", name_disc,".RData"))
  # Estimate changepoint models based on segmented package
  backpro_seg = list()

  print("estimate change point models based on segmented package infections")
  backpro_seg[["two_bp"]] = tryCatch(estimate_bp_segmented_backpro(data = dat_backpro, bp = 2,
                                                  start_bp = backpro_discrete$two_bp$bps),
                                     error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))

  if (cp_max_backpro>2) {
    message(paste0("n bp:", 3))
    backpro_seg[["three_bp"]] = tryCatch(estimate_bp_segmented_backpro(data=dat_backpro, bp=3,
                                                      start_bp = backpro_discrete$three_bp$bps),
                                         error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
    
  }
  if (cp_max_backpro>3) {
    message(paste0("n bp:", 4))
    start_bp_4 = sort(c(backpro_discrete$three_bp$bps,
                        max(setdiff(3:(nrow(dat_backpro)-3),
                                    c(backpro_discrete$three_bp$bps[1]+(-3:3),
                                      backpro_discrete$three_bp$bps[2]+(-3:3),
                                      backpro_discrete$three_bp$bps[3]+(-3:3))))))

    backpro_seg[["four_bp"]] = tryCatch(estimate_bp_segmented_backpro(data=dat_backpro, bp=4,
                                                             start_bp = start_bp_4),
                                        error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
    }
  if (cp_max_backpro>4) {
    message(paste0("n bp:", 5))
    start_bp_5 = sort(c(start_bp_4,
                        max(setdiff(3:(nrow(dat_backpro)-3),
                                    c(start_bp_4[1]+(-3:3),
                                      start_bp_4[2]+(-3:3),
                                      start_bp_4[3]+(-3:3),
                                      start_bp_4[4]+(-3:3))))))

    backpro_seg[["five_bp"]] = tryCatch(estimate_bp_segmented_backpro(data=dat_backpro, bp=5,
                                                             start_bp = start_bp_5),
                                        error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
    
  }
  if (cp_max_backpro>5) {
    message(paste0("n bp:", 6))
    start_bp_6 = sort(c(start_bp_5,
                        max(setdiff(3:(nrow(dat_backpro)-3),
                                    c(start_bp_5[1]+(-3:3),
                                      start_bp_5[2]+(-3:3),
                                      start_bp_5[3]+(-3:3),
                                      start_bp_5[4]+(-3:3),
                                      start_bp_5[5]+(-3:3))))))
    backpro_seg[["six_bp"]] = tryCatch(estimate_bp_segmented_backpro(data=dat_backpro, bp=6,
                                                    start_bp = start_bp_6),
                                       error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
  }
  aic_backpro = sapply(backpro_seg, function(x) tryCatch(AIC(x$segmented_model), error=function(e) NA))
  bic_backpro = sapply(backpro_seg, function(x) tryCatch(BIC(x$segmented_model), error=function(e) NA))
  }

  if (type %in% c("both", "onset")) {
    print("perform analysis of onsets")

    # Perform discrete optimization for onsets
    dat_onset=data %>%
      dplyr::filter(date>=ymd("2020-02-28"),
                    date<=max(date)-14)
    if (use_disc_optim_results &
        file.exists(paste0("results/changepoint/discrete_optimization/onset_disc_optim_",
                           name_disc,
                           ".RData"))) {
      load(paste0("results/changepoint/discrete_optimization/onset_disc_optim_", name_disc,".RData"))
    } else {
      onset_discrete = list()
    }
    if (is.null(onset_discrete[["two_bp"]])) {
      print("perform discrete optimization for 2 change points onsets")
      onset_discrete[["two_bp"]] = estimate_bp_disc_optim_onset(data=dat_onset, bp=2)
    }
    if (cp_max_onset>2 & is.null(onset_discrete[["three_bp"]])) {
      print("perform discrete optimization for 3 change points onsets")
      onset_discrete[["three_bp"]] = estimate_bp_disc_optim_onset(data=dat_onset, bp=3)
    }
    if(save_disc_optim_results)
      save(onset_discrete, file = paste0("results/changepoint/discrete_optimization/onset_disc_optim_", name_disc,".RData"))
    # Estimate changepoint models based on segmented package
    onset_seg = list()
    print("estimate change point models based on segmented package onset")
    onset_seg[["two_bp"]] = tryCatch(estimate_bp_segmented_onset(data=dat_onset, bp = 2,
                                                        start_bp = onset_discrete$two_bp$bps),
                                     error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
    if (cp_max_onset>2)
      onset_seg[["three_bp"]] = tryCatch(estimate_bp_segmented_onset(data=dat_onset, bp=3,
                                                            start_bp = onset_discrete$three_bp$bps),
                                         error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
    
    if (cp_max_onset>3) {
      start_bp_4 = sort(c(onset_discrete$three_bp$bps,
                          max(setdiff(3:(nrow(dat_onset)-3),
                                      c(onset_discrete$three_bp$bps[1]+(-3:3),
                                        onset_discrete$three_bp$bps[2]+(-3:3),
                                        onset_discrete$three_bp$bps[3]+(-3:3))))))

      onset_seg[["four_bp"]] = tryCatch(estimate_bp_segmented_onset(data=dat_onset, bp=4,
                                                               start_bp = start_bp_4),
                                        error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
      
    }
    if (cp_max_onset>4) {
      start_bp_5 = sort(c(start_bp_4,
                          max(setdiff(3:(nrow(dat_onset)-3),
                                      c(start_bp_4[1]+(-3:3),
                                        start_bp_4[2]+(-3:3),
                                        start_bp_4[3]+(-3:3),
                                        start_bp_4[4]+(-3:3))))))

      onset_seg[["five_bp"]] = tryCatch(estimate_bp_segmented_onset(data=dat_onset, bp=5,
                                                               start_bp = start_bp_5),
                                        error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
      
    }
    if (cp_max_onset>5) {
      start_bp_6 = sort(c(start_bp_5,
                          max(setdiff(3:(nrow(dat_onset)-3),
                                      c(start_bp_5[1]+(-3:3),
                                        start_bp_5[2]+(-3:3),
                                        start_bp_5[3]+(-3:3),
                                        start_bp_5[4]+(-3:3),
                                        start_bp_5[5]+(-3:3))))))
      onset_seg[["six_bp"]] = tryCatch(estimate_bp_segmented_onset(data=dat_onset, bp=6,
                                                              start_bp = start_bp_6),
                                       error=function(e) list(segmented_model = NA, coef=NA, breakpoints = NA, plot = NA))
      
    }
    aic_onset = sapply(onset_seg, function(x) tryCatch(AIC(x$segmented_model), error=function(e) NA))
    bic_onset = sapply(onset_seg, function(x) tryCatch(BIC(x$segmented_model), error=function(e) NA))
  }

  # Changepoint models onset
  if (type == "both") {
    res = list(aic_backpro=aic_backpro, bic_backpro=bic_backpro, cp_segmented_list_backpro = backpro_seg,
               aic_onset=aic_onset, bic_onset=bic_onset, cp_segmented_list_onset=onset_seg)
  } else if (type=="backpro") {
    res = list(aic_backpro=aic_backpro, bic_backpro=bic_backpro, cp_segmented_list_backpro = backpro_seg)
  } else if (type=="onset") {
    res = list(aic_onset=aic_onset, bic_onset=bic_onset, cp_segmented_list_onset=onset_seg)
  }
  res
}
