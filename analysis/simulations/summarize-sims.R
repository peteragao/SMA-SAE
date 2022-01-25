#### 0 LIBRARIES AND PATHS #####################################################
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

summarize_res <- function(res) {
  res %>% 
    group_by(method) %>%
    summarize(rmse = sqrt(mean((median - pop_mean)^2)),
              bias = mean(median - pop_mean),
              cov90 = mean(lower < pop_mean & upper > pop_mean))
}
#### 1 CMN 2020 WEAK MODEL #####################################################
CMN_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) readRDS(paste0(sim_res_dir, "CMN_res_", x, ".rds"))
          )) 
saveRDS(summarize_res(CMN_res), "results/simulations/CMN_res.rds")
CMN_inf_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) readRDS(paste0(sim_res_dir, "CMN-inf-smp_res_", x, ".rds"))
          )) 
saveRDS(summarize_res(CMN_inf_res), "results/simulations/CMN_inf_res.rds")


#### 2 CMN 2020 IMPROVED MODEL #################################################
CMN_imp_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "CMN-improved-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "CMN-improved-model_res_", x, ".rds"))
              }
            }
          )) 
saveRDS(summarize_res(CMN_imp_res), "results/simulations/CMN_imp_res.rds")
CMN_imp_inf_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "CMN-improved-model-inf-smp_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "CMN-improved-model-inf-smp_res_", x, ".rds"))
              }
            }
          )) 
saveRDS(summarize_res(CMN_imp_inf_res), "results/simulations/CMN_imp_inf_res.rds")

#### 3 SPATIAL MODEL ###########################################################
spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "CMN-cts-spatial-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "CMN-cts-spatial-model_res_", x, ".rds"))
              }
            })) 
saveRDS(summarize_res(spa_res), "results/simulations/cts_spatial_model_res.rds")
spa_inf_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "CMN-cts-spatial-model-inf-smp_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "CMN-cts-spatial-model-inf-smp_res_", x, ".rds"))
              }
            })) 
saveRDS(summarize_res(spa_inf_res), "results/simulations/cts_spatial_model_inf_res.rds")

#### 4 FIGURES #################################################################
CMN_res <- readRDS("results/simulations/CMN_res.rds")
CMN_inf_res <- readRDS("results/simulations/CMN_inf_res.rds")
CMN_imp_res <- readRDS("results/simulations/CMN_imp_res.rds")
CMN_imp_inf_res <- readRDS("results/simulations/CMN_imp_inf_res.rds")
spa_res <- readRDS("results/simulations/cts_spatial_model_res.rds")
spa_inf_res <- readRDS("results/simulations/cts_spatial_model_inf_res.rds")

pub_names <- rbind(
  c("GREG", "GREG (reduced)"),
  c("GREGFull", "GREG (full)"),
  c("HT", "Hájek"),
  c("iidBinFullULM", "Binomial ULM (full)"),
  c("iidBinULM", "Binomial ULM (reduced)"),
  c("iidGauFullULM", "Gaussian ULM (full)"),
  c("iidGauULM", "Gaussian ULM (reduced)"),
  c("iidSGREG", "SGREG (reduced) alt"),
  c("iidSGREGFull", "SGREG (full) alt"),
  c("iidSLogitGREG", "SGREG (reduced)"),
  c("iidSLogitGREGFull", "SGREG (full)"),
  c("iidSHT", "SH alt"),
  c("iidSLogitHT", "SH"),
  c("bymBinFullULM", "Spatial Binomial ULM (full)"),
  c("bymBinULM", "Spatial Binomial ULM (reduced)"),
  c("bymGauFullULM", "Spatial Gaussian ULM (full)"),
  c("bymGauULM", "Spatial Gaussian ULM (reduced)"),
  c("bymSGREG", "Spatial SGREG (reduced) alt"),
  c("bymSGREGFull", "Spatial SGREG (full) alt"),
  c("bymSLogitGREG", "Spatial SGREG (reduced)"),
  c("bymSLogitGREGFull", "Spatial SGREG (full)"),
  c("bymSHT", "Spatial SH alt"),
  c("bymSLogitHT", "Spatial SH")
) %>%
  as.data.frame() %>%
  setNames(c("internal", "publication")) 

selected_methods = 
  c("GREG", "HT", 
    "GREGFull", "iidSLogitHT",
    "iidSLogitGREG", "iidSLogitGREGFull",
    "iidBinULM", "iidBinFullULM")

pub_order <- c("Population Mean",
               "Hájek","SH", "Spatial SH",  "GREG (reduced)", 
               "Spatial SGREG (reduced)",
               "SGREG (reduced)", 
               "Spatial Binomial ULM (reduced)",
               "Binomial ULM (reduced)",
               "Spatial Gaussian ULM (reduced)",
               "Gaussian ULM (reduced)",
               "Spatial SH alt", "Spatial SGREG (reduced) alt",
               "SH alt", "SGREG (reduced) alt",
               "GREG (full)",
               "Spatial SGREG (full)",
               "SGREG (full)",
               "Spatial Binomial ULM (full)",
               "Binomial ULM (full)",
               "Spatial Gaussian ULM (full)",
               "Gaussian ULM (full)",
               "Spatial SGREG (full) alt",
               "SGREG (full) alt")



pub_names_spa <- rbind(
  c("GREG", "GREG (reduced)"),
  c("GREGFull", "GREG (full)"),
  c("HT", "Hájek"),
  c("iidGauFullULM", "Gaussian ULM (full)"),
  c("iidGauULM", "Gaussian ULM (reduced)"),
  c("iidSGREG", "SGREG (reduced)"),
  c("iidSGREGFull", "SGREG (full)"),
  c("iidSHT", "SH"),
  c("bymGauFullULM", "Spatial Gaussian ULM (full)"),
  c("bymGauULM", "Spatial Gaussian ULM (reduced)"),
  c("bymSGREG", "Spatial SGREG (reduced)"),
  c("bymSGREGFull", "Spatial SGREG (full)"),
  c("bymSHT", "Spatial SH")
) %>%
  as.data.frame() %>%
  setNames(c("internal", "publication")) 
pub_order_spa<- c("Population Mean",
                  "Hájek","SH","Spatial SH",  "GREG (reduced)", 
                   "SGREG (reduced)",
                  "Binomial ULM (reduced)",
                  "Gaussian ULM (reduced)",
                  "Spatial SGREG (reduced)",
                  "Spatial Binomial ULM (reduced)",
                  "Spatial Gaussian ULM (reduced)",
                  "Spatial SH alt", "Spatial SGREG (reduced) alt",
                  "SH alt", "SGREG (reduced) alt",
                  "GREG (full)",
                  "SGREG (full)",
                  "Binomial ULM (full)",
                  "Gaussian ULM (full)",
                  "Spatial SGREG (full)",
                  "Spatial Binomial ULM (full)",
                  "Spatial Gaussian ULM (full)",
                  "Spatial SGREG (full) alt",
                  "SGREG (full) alt")
selected_spatial_methods = 
  c("GREG", "HT", 
    "GREGFull", "iidSHT",
    "iidSGREG", "iidSGREGFull",
    "iidGauULM", "iidGauFullULM", "bymSHT",
    "bymSGREG", "bymSGREGFull",
    "bymGauULM", "bymGauFullULM")
fmt_tbl <- function(res, file, methods = unique(res$method), 
                    pub_names = pub_names, pub_order = pub_order,
                    rows = NULL) {
  res <- res %>%
    filter(method %in% methods) %>%
    mutate(method = 
             pub_names$publication[match(method, pub_names$internal)]) %>%
    arrange(match(method, pub_order)) %>%
    mutate(rmse = rmse * 100, bias = bias * 1000, cov90 = cov90 * 100) %>%
    setNames(c("Method", 
               paste0("RMSE (x 100)"),
               paste0("Bias (x 1000)"),
               "90% Interval Cov.")) %>%
    knitr::kable(digits = c(0, 2, 2, 0), format = "latex", booktabs = T,
                 linesep = "")
  if (!is.null(rows)) {
    res <- res %>%  kableExtra::row_spec(rows, hline_after = T) 
  }
  res %>% writeLines(file)
  res
}
fmt_tbl(CMN_res, 'paper/figures/CMN_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(CMN_res, 'paper/figures/sel_CMN_res_summary.tex', selected_methods,
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(CMN_inf_res, 'paper/figures/CMN_inf_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(CMN_inf_res,
        'paper/figures/sel_CMN_inf_res_summary.tex', selected_methods,
        pub_names = pub_names, pub_order = pub_order)

fmt_tbl(CMN_imp_res, 'paper/figures/CMN_imp_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(CMN_imp_res, 'paper/figures/sel_CMN_imp_res_summary.tex', selected_methods,
        pub_names = pub_names, pub_order = pub_order,
        rows = c(2, 5))
fmt_tbl(CMN_imp_inf_res, 'paper/figures/CMN_imp_inf_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(CMN_imp_inf_res,
        'paper/figures/sel_CMN_imp_inf_res_summary.tex', selected_methods,
        pub_names = pub_names, pub_order = pub_order,
        rows = c(2, 5))

fmt_tbl(spa_res, 'paper/figures/spatial_sim_res_summary.tex',
        pub_names = pub_names_spa, pub_order = pub_order_spa)
fmt_tbl(spa_res,
        'paper/figures/sel_spatial_sim_res_summary.tex', selected_spatial_methods,
        pub_names = pub_names_spa, pub_order = pub_order_spa,
        rows = c(3, 8))
fmt_tbl(spa_inf_res, 'paper/figures/spatial_inf_sim_res_summary.tex',
        pub_names = pub_names_spa, pub_order = pub_order_spa)
fmt_tbl(spa_inf_res,
        'paper/figures/sel_spatial_inf_sim_res_summary.tex', 
        selected_spatial_methods, pub_names = pub_names_spa, pub_order = pub_order_spa,
        rows = c(3, 8))

#### 5 FIGURES #################################################################
i <- 2
x <- round(i / 10) * 10 + 1
CMN_imp_inf_res_x <- readRDS(paste0(sim_res_dir, "CMN-improved-model-inf-smp_res_", x, ".rds"))
res <- CMN_imp_inf_res

selected_methods = 
  c("GREG", "HT", 
    "iidSLogitGREG", 
    "iidBinFullULM", "iidGauULM")
est_plot <- function(res, methods = selected_methods) {
  pop_means <- res %>% filter(method == "HT" & id_sim == i) %>%
    mutate(est = pop_mean, method = "Population Mean") 
  res <- res %>%
    filter(method %in% methods) %>%
    mutate(method = 
             pub_names$publication[match(method, pub_names$internal)]) %>%
    arrange(match(method, pub_order))
  ranks <- pop_means %>%
    dplyr::select(region, est) %>%
    mutate(rank = rank(est, ties.method = "first")) %>%
    dplyr::select(region, rank)
  plot_dat <- bind_rows(pop_means, res %>% filter(id_sim == i)) %>%
    left_join(ranks, by = "region")  %>%
    arrange(match(method, pub_order)) %>%
    mutate(method = factor(method, levels = pub_order))
  ggplot(plot_dat, aes(x = rank, y = est, color = method, shape = method)) + geom_point()
  
}

metric_plot <- function(res, methods = selected_methods) {
  pop_means <- res %>% filter(method == "HT" & id_sim == i) %>%
    mutate(est = pop_mean, method = "Population Mean") 
  ranks <- pop_means %>%
    dplyr::select(region, est) %>%
    mutate(rank = rank(est, ties.method = "first")) %>%
    dplyr::select(region, rank)
  metrics <- res %>% 
    filter(method %in% methods) %>%
    mutate(method = 
             pub_names$publication[match(method, pub_names$internal)]) %>%
    arrange(match(method, pub_order)) %>%
    group_by(method, region) %>%
    summarize(rmse = sqrt(mean((median - pop_mean)^2)),
              bias = mean(median - pop_mean),
              cov90 = mean(lower < pop_mean & upper > pop_mean)) %>%
    ungroup() %>%
    left_join(ranks)
  ggplot(metrics, aes(x = rank, y = cov90, color = method, shape = method)) + geom_point()
  
}

