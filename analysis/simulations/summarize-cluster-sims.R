#### 0 LIBRARIES AND PATHS #####################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)
library(kableExtra)
#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/SMA-SAE/"))

source("analysis/models.R")
res_dir <- "results/"
cluster_res_dir <- "results/cluster/"
sim_res_dir <- "results/cluster/spatial-cluster-sims/"

summarize_res <- function(res) {
  res %>% 
    group_by(method) %>%
    summarize(rmse = sqrt(mean((median - pop_mean)^2)),
              abs_bias = mean(abs(median - pop_mean)),
              cov90 = mean(lower < pop_mean & upper > pop_mean),
              int_len = mean(upper - lower))
}
#### 3 SPATIAL MODEL ###########################################################
spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "FULL-bin-spatial-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "FULL-bin-spatial-model_res_", x, ".rds"))
              }
            })) 
saveRDS(summarize_res(spa_res), "results/simulations/bin_spatial_model_res.rds")
spa_inf_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "FULL-bin-spatial-model-inf-smp_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "FULL-bin-spatial-model-inf-smp_res_", x, ".rds"))
              }
            })) 
saveRDS(summarize_res(spa_inf_res), "results/simulations/bin_spatial_model_inf_res.rds")

#### 4 FIGURES #################################################################
# set names for manuscript tables
pub_names <- rbind(
  c("GREG", "MA  (red.)"),
  c("GREGFull", "MA (full)"),
  c("HT", "Direct (Hájek)"),
  c("iidBinFullULM", "Bin. ULM (full)"),
  c("iidBinULM", "Bin. ULM  (red.)"),
  c("iidLNBinFullULM", "LoNo-Bin. ULM (full)"),
  c("iidLNBinULM", "LoNo-Bin. ULM  (red.)"),
  c("iidBBinFullULM", "Betabin. ULM (full)"),
  c("iidBBinULM", "Betabin. ULM  (red.)"),
  c("iidSLogitGREG", "SMA  (red.)"),
  c("iidSLogitGREGFull", "SMA (full)"),
  c("iidSLogitHT", "SH"),
  c("bymBinFullULM", "Spatial Bin. ULM (full)"),
  c("bymBinULM", "Spatial Bin. ULM  (red.)"),
  c("bymLNBinFullULM", "Spatial LoNo-Bin. ULM (full)"),
  c("bymLNBinULM", "Spatial LoNo-Bin. ULM  (red.)"),
  c("bymBBinFullULM", "Spatial Betabin. ULM (full)"),
  c("bymBBinULM", "Spatial Betabin. ULM  (red.)"),
  c("bymSLogitGREG", "Spatial SMA  (red.)"),
  c("bymSLogitGREGFull", "Spatial SMA (full)"),
  c("bymSLogitHT", "Spatial SH")
) %>%
  as.data.frame() %>%
  setNames(c("internal", "publication")) 


pub_order <- c("Population Mean",
               "Direct (Hájek)","SH", "Spatial SH",  
               "MA  (red.)", 
               "SMA  (red.)", 
               "Bin. ULM  (red.)",
               "Betabin. ULM  (red.)",
               "LoNo-Bin. ULM  (red.)",
               "Spatial SMA  (red.)",
               "Spatial Bin. ULM  (red.)",
               "Spatial Betabin. ULM  (red.)",
               "Spatial LoNo-Bin. ULM  (red.)",
               "MA (full)",
               "SMA (full)",
               "Bin. ULM (full)",
               "Betabin. ULM (full)",
               "LoNo-Bin. ULM (full)",
               "Spatial SMA (full)",
               "Spatial Bin. ULM (full)",
               "Spatial Betabin. ULM (full)",
               "Spatial LoNo-Bin. ULM (full)")

fmt_tbl <- function(res, file, methods = unique(res$method), 
                    pub_names = pub_names, pub_order = pub_order,
                    rows = NULL) {
  res <- res %>%
    filter(method %in% methods) %>%
    mutate(method = 
             pub_names$publication[match(method, pub_names$internal)]) %>%
    arrange(match(method, pub_order)) %>%
    mutate(rmse = rmse * 100, abs_bias = abs_bias * 100, cov90 = cov90 * 100, int_len = int_len * 100) %>%
    setNames(c("Method", 
               paste0("RMSE (x 100)"),
               paste0("MAE (x 100)"),
               "90% Cov.",
               "Int. Len. (x 100)")) %>%
    knitr::kable(digits = c(0, 2, 2, 0, 2), format = "latex", booktabs = T,
                 linesep = "")
  if (!is.null(rows)) {
    #res <- res %>%  kableExtra::row_spec(rows, hline_after = T) 
  }
  res %>% writeLines(file)
  res
}

fmt_tbl(summarize_res(spa_res), 'paper/figures/bin_spatial_sim_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
fmt_tbl(summarize_res(spa_inf_res), 'paper/figures/bin_spatial_inf_sim_res_summary.tex',
        pub_names = pub_names, pub_order = pub_order)
