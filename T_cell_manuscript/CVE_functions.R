# This is a helper function that restricts the estimates
# and marker range of ests$cr to between low and high

truncate_ce <- function(ests, low, high){
  ind_exclude = which(ests$cr$s < low | ests$cr$s > high)
  
  if(length(ind_exclude) == 0)
    return(ests)
  
  else {
    ests$cr$s = ests$cr$s[-ind_exclude]
    ests$cr$est = ests$cr$est[-ind_exclude]
    ests$cr$se = ests$cr$se[-ind_exclude]
    ests$cr$ci_lower = ests$cr$ci_lower[-ind_exclude]
    ests$cr$ci_upper = ests$cr$ci_upper[-ind_exclude]
    return(ests)
  }
}


# This is a helper function that restricts the estimates
# and marker range of ests$cve to between low and high 

truncate_cve <- function(ests, low, high){
  ind_exclude = which(ests$cve$s < low | ests$cve$s > high)
  
  if(length(ind_exclude) == 0)
    return(ests)
  
  else {
    ests$cve$s = ests$cve$s[-ind_exclude]
    ests$cve$est = ests$cve$est[-ind_exclude]
    ests$cve$se = ests$cve$se[-ind_exclude]
    ests$cve$ci_lower = ests$cve$ci_lower[-ind_exclude]
    ests$cve$ci_upper = ests$cve$ci_upper[-ind_exclude]
    return(ests)
  }
}

# Conduct a basic CoR analysis using the vaccine package
CoR_D1 <- function(dt, density,
                   arm_name, time, event, 
                   marker_name, adj_vars, t0,
                   title_name, xlab_name,
                   dir_name, save_pdf = TRUE){
  
  dat <- load_data(
    time = time,
    event = event,
    vacc = arm_name,
    marker = marker_name,
    covariates = adj_vars,
    weights = "wt.D15",
    ph2 = "casecontrol",
    data = dt)
  
  ests_cox <- est_ce(dat = dat, type="Cox", t_0 = t0, return_p_value = TRUE)
  ests_np <- est_ce(dat = dat, type="NP", t_0 = t0, return_p_value = TRUE)
  
  low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  
  Cox = truncate_ce(ests_cox, low, high)
  Nonparametric = truncate_ce(ests_np, low, high)
  
  p = plot_ce(Cox, Nonparametric,
              zoom_y = c(0, 1.1), zoom_x = c(1, 6)) +
    geom_density(aes(x = x, y = after_stat(scaled)), 
                 data = data.frame(x = density$s), 
                 inherit.aes = F, bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted risk ', t0, ' days post D15')) +
    labs(title = title_name) +
    theme_bw(base_size = 24) + 
    theme(legend.position = 'bottom')
  
  if (save_pdf)
    ggsave(dir_name, plot = p, width = 12, height = 8)
  
  return(list(plot = p,
              Cox_CoR = ests_cox,
              NP_CoR = ests_np))
}


# Helper function that plots both Cox-based and NP-based CoR curves
# ests_cox and ests_np are returned objects from 
# the vaccine package.
# mk_dist is ph2 marker measurements with weights

plot_cr <- function(ests_cox, ests_np, mk_dist,
                     xlab_name = 'log10 CD4+ IFNg OR IL-2 (%)',
                     title_name = 'Hybird Immunity Per-Protocol Group (HIPP)',
                     t0 = 165, x_range = c(-2, 1),
                     x_breaks = c(-2, -1, 0, 1),
                    x_labels = c('-2', '-1', '0', '1')){
  
  density = list(s = as.vector(unlist(mk_dist[,1])),
                 w = as.vector(unlist(mk_dist[,2])))
  low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  
  Cox = truncate_ce(ests_cox, low, high)
  Nonparametric = truncate_ce(ests_np, low, high)
  
  p = plot_ce(Cox, Nonparametric,
              zoom_y = c(0, 1.1), zoom_x = x_range) +
    geom_density(aes(x = x, y = after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w), 
                 inherit.aes = F, 
                 bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted risk \n', t0, ' days post M1')) +
    labs(title = title_name) +
    theme_bw(base_size = 16) + 
    theme(legend.position = 'bottom')
  
  print(p)
}

# Helper function that plots the NP-based CoR curve

plot_cr2 <- function(ests_cox, ests_np, mk_dist, overall_incidence,
                    xlab_name = 'log10 CD4+ IFNg OR IL-2 (%)',
                    title_name = 'Hybird Immunity Per-Protocol Group (HIPP)',
                    t0 = 165, x_range = c(-2, 1),
                    x_breaks = c(-2, -1, 0, 1),
                    x_labels = c('-2', '-1', '0', '1')){
  
  density = list(s = as.vector(unlist(mk_dist[,1])),
                 w = as.vector(unlist(mk_dist[,2])))
 low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  
  Cox = truncate_ce(ests_cox, low, high)
  Nonparametric = truncate_ce(ests_np, low, high)
  
  p = plot_ce(Nonparametric,
              zoom_y = c(0, 1.1), zoom_x = x_range) +
    geom_density(aes(x = x, y = after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w), 
                 inherit.aes = F, 
                 bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    geom_hline(yintercept = overall_incidence, 
               color = 'gray41', linetype = 'dashed') +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted risk \n', t0, ' days post M1')) +
    labs(title = title_name) +
    theme_bw(base_size = 18) + 
    theme(legend.position = 'bottom')
  
  print(p)
}


# Helper function that plots the NP-based CoR curve for fold-rise markers

plot_cr_fr <- function(ests_np, mk_dist, overall_incidence,
                       xlab_name = 'log10 CD4+ IFNg OR IL-2 (%)',
                       title_name = 'Hybird Immunity Per-Protocol Group (HIPP)',
                       t0 = 165, x_range = c(-2, 1),
                       x_breaks = c(-2, -1, 0, 1)){
  
  density = list(s = as.vector(unlist(mk_dist[,1])),
                 w = as.vector(unlist(mk_dist[,2])))
  low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  
  Cox = truncate_ce(ests_cox, low, high)
  Nonparametric = truncate_ce(ests_np, low, high)
  
  p = plot_ce(Nonparametric,
              zoom_y = c(0, 1.1), zoom_x = x_range) +
    geom_density(aes(x = x, y = after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w), 
                 inherit.aes = F, 
                 bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    geom_hline(yintercept = overall_incidence, 
               color = 'gray41', linetype = 'dashed') +
    scale_x_continuous(breaks = x_breaks ) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted risk \n', t0, ' days post M1')) +
    labs(title = title_name) +
    theme_bw(base_size = 16) + 
    theme(legend.position = 'bottom')
  
  print(p)
}

# Helper function that plots the Cox-based and NP-based CVE or RR curves

plot_cve <- function(ests_cox, ests_np, mk_dist, mk_VIPP_med,
                     xlab_name = 'log10 CD4+ IFNg OR IL-2 (%)',
                     title_name = 'Hybird Immunity Per-Protocol Group (HIPP)',
                     t0 = 165, x_range = c(-2, 1), y_range = c(-1.5, 2),
                     x_breaks = c(-2, -1, 0, 1), 
                     x_labels = c('-2', '-1', '0', '1'),
                     y_breaks = c(-1, -0.5, 0, 0.5, 1)){
  
  density = list(s = as.vector(unlist(mk_dist[,1])),
                 w = as.vector(unlist(mk_dist[,2])))
  low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  Cox = truncate_cve(ests_cox, mk_VIPP_med, high)
  Nonparametric = truncate_cve(ests_np, mk_VIPP_med, high)
  
  p = plot_ce(Cox, Nonparametric, which = 'CVE',
              zoom_y = y_range, zoom_x = x_range) +
    geom_density(aes(x = x, y = after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w), 
                 inherit.aes = F, 
                 bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    geom_vline(xintercept = mk_VIPP_med, color = 'dark red',
               linetype = 'dashed', linewidth = 1.1) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(breaks = y_breaks, labels = function(x) paste0(x*100, "%")) +
    labs(title = 'Controlled Relative Risk',
         subtitle = 'Marker-Specific HIPP Risk vs. Overall VIPP Risk') +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted controlled \n relative risk', t0, ' days post M1')) +
    theme_bw(base_size = 16) + 
    theme(legend.position = 'bottom')
  
  print(p)
}

# Helper function that plots the NP-based CVE or RR curves

plot_cve2 <- function(ests_cox, ests_np, mk_dist, mk_VIPP_med,
                     xlab_name = 'log10 CD4+ IFNg OR IL-2 (%)',
                     title_name = 'Hybird Immunity Per-Protocol Group (HIPP)',
                     t0 = 165, x_range = c(-2, 1), y_range = c(-0.2, 3.2),
                     x_breaks = c(-2, -1, 0, 1), 
                     x_labels = c('-2', '-1', '0', '1'),
                     y_breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)){
  
  density = list(s = as.vector(unlist(mk_dist[,1])),
                 w = as.vector(unlist(mk_dist[,2])))
  low = wtd.quantile(density$s, density$w, c(0.025, 0.975))[1]
  high = wtd.quantile(density$s, density$w, c(0.025, 0.975))[2]
  
  Cox = truncate_cve(ests_cox, mk_VIPP_med, high)
  Nonparametric = truncate_cve(ests_np, mk_VIPP_med, high)
  
  p = plot_ce(Nonparametric, which = 'CVE',
              zoom_y = y_range, zoom_x = x_range) +
    geom_density(aes(x = x, y = 3*after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w), 
                 inherit.aes = F, 
                 bounds = c(min(density$s), max(density$s)),
                 fill = "forestgreen", alpha = 0.3, color = NA) +
    geom_vline(xintercept = mk_VIPP_med, color = 'dark red',
               linetype = 'dashed', linewidth = 1.1) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = 'Controlled Relative Risk',
         subtitle = 'Marker-Specific HIPP Risk vs. Overall VIPP Risk') +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted controlled \n relative risk', t0, ' days post M1')) +
    theme_bw(base_size = 18) + 
    theme(legend.position = 'bottom')
  
  print(p)
}
