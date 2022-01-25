library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)


#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/SMA-SAE/"))

source("analysis/models.R")
res_dir <- "results/"
dir.create(file.path(res_dir), showWarnings = FALSE)
cluster_res_dir <- "results/cluster/"
dir.create(file.path(cluster_res_dir), showWarnings = FALSE)
sim_res_dir <- "results/cluster/sims/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)

#### GENERATE POPULATION ####
gen_pop_X <- function(N = 20000, 
                      C = 80) {
  pop_dat <-
    data.frame(id_unit = 1:N,
               id_area = rep(1:C, each = N / C)) %>%
    mutate(x1 = rbinom(N, size = 1, prob = 0.3 + 0.5 * id_area / C),
           x2 = rbinom(N, size = 1, prob = 0.2))
  pop_dat
}

gen_pop_Y <- function(pop_dat, eta_c_sd = .15,
                      e_ch_sd = .5,
                      beta = c(3, .03, .04)) {
  N <- nrow(pop_dat)
  C <- length(unique(pop_dat$id_area))
  pop_dat %>%
    mutate(eta_c = rep(rnorm(C, sd = eta_c_sd), each = N / C),
           e_ch = rnorm(N, sd = e_ch_sd),
           y = exp(beta[1] + beta[2] * x1 + beta[3] * x2 + eta_c + e_ch))
}


gen_sample_ind <- 
  function(
    pop_dat, strata,
    prop_sample_per_strata = .02
  ) {
    strata_vec <- as.vector(model.matrix(strata, pop_dat)[, 2])
    n_strata <- length(unique(strata_vec))
    pop_dat$strata <- strata_vec
    pop_dat$idx <- 1:nrow(pop_dat)
    ns_per_area_table <- pop_dat %>%
      group_by(strata) %>%
      summarize(n = round(n() * prop_sample_per_strata))
    ns_per_area <- ns_per_area_table$n
    
    pop_dat <- pop_dat %>%
      left_join(ns_per_area_table, by = "strata") %>%
      group_by(strata) %>%  
      mutate(wt = n() / n) %>%
      ungroup()
    
    sample_dat <- pop_dat %>%
      group_by(strata) %>%
      nest() %>%             
      ungroup() %>% 
      mutate(samp = 
               map2(data, ns_per_area,
                    function(x, y) 
                      slice_sample(x, n = y,
                                   weight_by = 1/x$wt))) %>%
      dplyr::select(-data) %>%
      unnest(samp)
    return(list(ind = sample_dat$idx, pop_dat = pop_dat))
  }

gen_sample_ind_inf <- 
  function(
    pop_dat, strata,
    prop_sample_per_strata = .02
  ) {
    strata_vec <- as.vector(model.matrix(strata, pop_dat)[, 2])
    n_strata <- length(unique(strata_vec))
    pop_dat$strata <- strata_vec
    pop_dat$idx <- 1:nrow(pop_dat)
    ns_per_area_table <- pop_dat %>%
      group_by(strata) %>%
      summarize(n = round(n() * prop_sample_per_strata))
    ns_per_area <- ns_per_area_table$n
    pop_dat <- pop_dat %>%
      mutate(a = 1 + 1 * (x2 == 1)) %>%
      left_join(ns_per_area_table, by = "strata") %>%
      group_by(strata) %>%  
      mutate(wt = sum(a) / a / n) %>%
      ungroup()
    sample_dat <- pop_dat %>%
      group_by(strata) %>%  
      nest() %>%             
      ungroup() %>% 
      mutate(samp = 
               map2(data, ns_per_area,
                    function(x, y) 
                      slice_sample(x, n = y,
                                   weight_by = 1/x$wt))) %>%
      dplyr::select(-data) %>%
      unnest(samp)
    return(list(ind = sample_dat$idx, pop_dat = pop_dat))
  }
#### RUN SIMULATIONS ####


run_one_sim <- function(sample_ind, pop_dat_X, pop_gen_fn, i, ...) {
  pop_dat <- pop_gen_fn(pop_dat_X, ...) %>%
    mutate(region = as.character(id_area)) %>%
    mutate(z = 1 * (y < 12))
  sample_dat <- pop_dat[sample_ind, ]
  pop_means <- pop_dat %>%
    group_by(region) %>%
    summarize(pop_mean = mean(z)) %>%
    as.data.frame()
  
  sample_dat$id <- 1:nrow(sample_dat)
  sample_des <- svydesign(ids = ~id, strata = ~region,
                          weights = ~wt, data = sample_dat)
  direct_est <- get_direct(~z, ~region, sample_des)
  iid_sdir_logit_est <- get_iid_sdir_logit(direct_est)
  iid_sdir_est <- get_iid_sdir(direct_est)
  working_fit <- svyglm(z ~ x1, sample_des, family = quasibinomial())
  greg_est <- get_greg(working_fit, z ~ x1, ~region,
                       pop_dat, sample_des)
  iid_sgreg_logit_est <- get_iid_sdir_logit(greg_est)
  iid_sgreg_est <- get_iid_sdir(greg_est)
  
  working_full_fit <- svyglm(z ~ x1 + x2, 
                             sample_des, family = quasibinomial())
  greg_full_est <- get_greg(working_full_fit, z ~ x1 + x2,
                            ~region, pop_dat, sample_des) %>%
    mutate(method = "GREGFull")
  iid_sgreg_full_logit_est <- get_iid_sdir_logit(greg_full_est) %>%
    mutate(method = "iidSGREGFull")
  iid_sgreg_full_est <- get_iid_sdir(greg_full_est) %>%
    mutate(method = "iidSLogitGREGFull")
  
  iid_y_ulm_est <- get_iid_gau_ulm(log(y) ~ x1, 
                                   ~region, pop_dat, sample_dat,
                                   h = function(x) exp(x) < 12)
  iid_z_ulm_est <- get_iid_bin_ulm(z ~ x1,
                                   ~region, 
                                   sample_dat = sample_dat,
                                   pop_dat = pop_dat)
  
  iid_y_full_ulm_est <- get_iid_gau_ulm(log(y) ~ x1 + x2, 
                                        ~region, pop_dat, sample_dat,
                                        h = function(x) exp(x) < 12) %>%
    mutate(method = "iidGauFullULM")
  iid_z_full_ulm_est <- get_iid_bin_ulm(z ~ x1 + x2, 
                                        ~region, 
                                        sample_dat = sample_dat,
                                        pop_dat = pop_dat) %>%
    mutate(method = "iidBinFullULM")
  res <- bind_rows(direct_est, 
                   iid_sdir_est, 
                   iid_sdir_logit_est, 
                   greg_est,
                   iid_sgreg_est,
                   iid_sgreg_logit_est,
                   greg_full_est,
                   iid_sgreg_full_est,
                   iid_sgreg_full_logit_est,
                   iid_y_ulm_est,
                   iid_z_ulm_est,
                   iid_y_full_ulm_est,
                   iid_z_full_ulm_est)
  res$id_sim = i
  res <- left_join(res, pop_means, by = "region")
  return(res)
}

set.seed(1125)
pop_dat <- gen_pop_X()
sample_ind_res <- gen_sample_ind(pop_dat, ~id_area, .2)
inf_sample_ind_res <- gen_sample_ind_inf(pop_dat, ~id_area, .2)
sample_ind <- sample_ind_res$ind
pop_dat <- sample_ind_res$pop_dat
inf_sample_ind <- inf_sample_ind_res$ind
inf_pop_dat <- inf_sample_ind_res$pop_dat

args = commandArgs(TRUE)
# supplied at the command line
k = as.numeric(args[1])
set.seed(k)
res <- do.call(
  rbind,
  lapply(1:10, function(i) run_one_sim(sample_ind,  pop_dat, gen_pop_Y, i + (k - 1) * 10))
  )
saveRDS(res, paste0(sim_res_dir, "CMN_res_", k, ".rds"))

k = as.numeric(args[1])
set.seed(k)
res <- do.call(
  rbind,
  lapply(1:10, function(i) run_one_sim(inf_sample_ind,  inf_pop_dat, gen_pop_Y, i + (k - 1) * 10))
)
saveRDS(res, paste0(sim_res_dir, "CMN-inf-smp_res_", k, ".rds"))
