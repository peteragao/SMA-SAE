#### 0 LIBRARIES AND PATHS ####################################################
library(rgdal)
#library(readstata13)
library(geosphere)
library(spdep)
library(survey)
library(INLA)
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
library(sf)
library(patchwork)
library(tiff)
# Set working directory (not necessary if using .Rproj)
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/SMA-SAE/"))
res_dir <- "results/"
dir.create(file.path(res_dir), showWarnings = FALSE)
cluster_res_dir <- "results/cluster/"
dir.create(file.path(cluster_res_dir), showWarnings = FALSE)
mcv_res_dir <- "results/cluster/mcv/"
dir.create(file.path(mcv_res_dir), showWarnings = FALSE)

#### 0.1 COUNTRY INFO ##########################################################
country <- "Nigeria"
survey_year <- 2018
gadm_abbrev <- "NGA"
pop_abbrev <- 'nga'

# dhsStata, which contains survey data
dhs_file <- "dhsStata/NGKR7BDT/NGKR7BFL.DTA"
# survey GPS
dhsFlat_file <- "NGGE7BFL"

country_dir <- paste0("../", country, "/")
poly_path <- paste0(country_dir, "shapeFiles_gadm")
res_dir <- "./results"
proj_data_dir <-  paste0("./data/" , country, "/", sep="")
dir.create(file.path(proj_data_dir), showWarnings = FALSE)
# link to the main dropbox folder for population rasters 
pop_raw_dir <- paste0(country_dir, 'Population/')
cov_raw_dir <- paste0(country_dir, 'covariates/')

#### 1 LOAD DATA ##############################################################
 
#### 1.1 GEO DATA ##############################################################

#### 1.1.1 Load GADM polygons ####
poly_layer_adm0 <- paste('gadm36', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm36', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm36', gadm_abbrev,
                         '2', sep = "_")

poly_adm0 <- readOGR(dsn = poly_path, encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly_layer_adm0)) 
# use encoding to read special characters
poly_adm1 <- readOGR(dsn = poly_path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly_layer_adm1))

if(sum(grepl(paste('gadm36', gadm_abbrev,
                   '2', sep = "_"), list.files(poly_path))) != 0){
  poly_adm2 <- readOGR(dsn = poly_path,encoding = "UTF-8", use_iconv = TRUE,
                       layer = as.character(poly_layer_adm2))}

if (exists("poly_adm2")) {
  proj4string(poly_adm0) <- proj4string(poly_adm1)  <- proj4string(poly_adm2)
}else {
  proj4string(poly_adm0) <- proj4string(poly_adm1)
}
poly_adm2$NAME_1 <- as.character(poly_adm2$NAME_1)
poly_adm1$NAME_1 <- as.character(poly_adm1$NAME_1)
poly_adm2$NAME_1[poly_adm2$NAME_1 == "Federal Capital Territory"] <- "Abuja"
poly_adm1$NAME_1[poly_adm1$NAME_1 == "Federal Capital Territory"] <- "Abuja"

#### 1.1.2 Create adjacency matrices and tables ####
if(exists("poly_adm1")){
  admin1_mat <- poly2nb(SpatialPolygons(poly_adm1@polygons))
  admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)
  colnames(admin1_mat) <- 
    rownames(admin1_mat) <-
    paste0("admin1_", 1:dim(admin1_mat)[1])
  admin1_names <- data.frame(GADM = poly_adm1@data$NAME_1,
                             Internal = rownames(admin1_mat))
}else{
  message("There is no Admin1 polygon file.")
}
if(exists("poly_adm2")){
  admin2_mat <- poly2nb(SpatialPolygons(poly_adm2@polygons))
  admin2_mat <- nb2mat(admin2_mat, zero.policy = TRUE)
  colnames(admin2_mat) <- 
    rownames(admin2_mat) <- 
    paste0("admin2_", 1:dim(admin2_mat)[1])
  admin2_names <- data.frame(GADM = poly_adm2@data$NAME_2,
                             Internal = rownames(admin2_mat))
}else{
  message("There is no Admin2 polygon file.")
}

#### 1.2 DHS DATA ##############################################################
if (TRUE) {
  svy_dat <- readRDS(paste0(proj_data_dir, "clean_DHS_data_DELETE-ME.rds"))
} else {
  #### 1.2.1 Load EA locations ####
  ea_locs_path <- paste0(country_dir, "dhsFlat/", dhsFlat_file)
  ea_locs <- readOGR(dsn = path.expand(ea_locs_path),
                     layer = as.character(dhsFlat_file))
  
  # Remove EAs with missing geo information
  missing_idx <- ea_locs$LATNUM == 0
  ea_locs <- ea_locs[!missing_idx, ]
  
  #### 1.2.2 Load Recode data ####
  in_dat = read.dta13(paste0(country_dir, dhs_file))

  svy_dat = data.frame(
    cluster = in_dat$v001, 
    hshold = in_dat$v002,
    stratum = in_dat$v023, 
    h9 = in_dat$h9,
    mcv_yes = 1 * (in_dat$h9 == "vaccination date on card" | 
                     in_dat$h9 == "reported by mother" | 
                     in_dat$h9 == "vaccination marked on card"),
    alive = in_dat$b5 == "yes",
    doi = in_dat$v008,
    dob = in_dat$b3,
    wt = in_dat$v005/1000000
    )
  svy_dat <- subset(svy_dat, doi-dob <= 23 & doi-dob >= 12)
  svy_dat <- subset(svy_dat, !is.na(mcv_yes))
  svy_dat = subset(svy_dat, alive)
  
  #### 1.2.3 Merge geographic info ####
  ea_dat <- data.frame(cluster = ea_locs$DHSCLUST,
                       urban = ea_locs$URBAN_RURA,
                       lon = ea_locs$LONGNUM,
                       lat = ea_locs$LATNUM)
  svy_dat = merge(svy_dat, ea_dat, by = "cluster")

  points_frame <- as.data.frame(svy_dat[,c("lon", "lat")])
  points_frame <- SpatialPoints(points_frame)
  
  # assign points to admin 2
  poly_over_adm2 <- SpatialPolygons(poly_adm2@polygons)
  proj4string(points_frame) <-
    proj4string(poly_over_adm2) <- 
    proj4string(poly_adm2)  <- 
    proj4string(poly_adm1)  
  admin2_key <- over(points_frame, poly_over_adm2)
  miss_frame_adm2 <- 
    matrix(unique(points_frame@coords[which(is.na(admin2_key)),]), ncol = 2)
  
  if(dim(miss_frame_adm2)[1] != 0){
    miss_poly_adm2 <- dist2Line( miss_frame_adm2, poly_over_adm2)
    for(i in 1:dim(miss_poly_adm2)[1]){
      long_ids <- 
        which(points_frame@coords[,c("lon")] %in% miss_frame_adm2[i,1])
      lat_ids <-
        which(points_frame@coords[,c("lat")] %in% miss_frame_adm2[i,2])
      ids <- intersect(long_ids, lat_ids)
      admin2_key[ids] <- rep(miss_poly_adm2[i, 'ID'], length(ids))
    }
  }
  svy_dat$admin2 <- admin2_key
  svy_dat$admin2_char <- paste0("admin2_", admin2_key)
  svy_dat$admin2_name <- as.character(poly_adm2@data$NAME_2)[admin2_key]
  
  # assign points to admin 2
  poly_over_adm1 <- SpatialPolygons(poly_adm1@polygons)
  proj4string(points_frame) <-
    proj4string(poly_over_adm1) <- 
    proj4string(poly_adm1) 
  admin1_key <- over(points_frame, poly_over_adm1)
  miss_frame_adm1 <- 
    matrix(unique(points_frame@coords[which(is.na(admin1_key)),]), ncol = 2)
  if(dim(miss_frame_adm1)[1] != 0){
    miss_poly_adm1 <- dist2Line(miss_frame_adm1, poly_over_adm1)
    for(i in 1:dim(miss_poly_adm1)[1]){
      long_ids <- 
        which(points_frame@coords[,c("lon")] %in% miss_frame_adm1[i,1])
      lat_ids <- 
        which(points_frame@coords[,c("lat")] %in% miss_frame_adm1[i,2])
      ids <- intersect(long_ids, lat_ids)
      admin1_key[ids] <- rep(miss_poly_adm1[i, 'ID'], length(ids))
    }
  }
  svy_dat$admin1 <- admin1_key
  svy_dat$admin1_char <- paste0("admin1_", admin1_key)
  svy_dat$admin1_name <- as.character(poly_adm1@data$NAME_1)[admin1_key]
  
  # add number of trials
  svy_dat$n_trials <- 1
  saveRDS(svy_dat,
          file = paste0(proj_data_dir, "clean_DHS_data_DELETE-ME.rds"))
}
#### 1.2 POPULATION DATA #######################################################
if (F) {
  source("./analysis/Nigeria/process-pop-cov.R")
  load_pop_data(poly_adm1, poly_adm2)
} 
# adm1_pix_tbl <- readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl.rds"))
adm1_pop_totals <- readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pop_totals.rds"))

#adm1_pix_tbl <- 
#  readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_1.rds"))
#### 1.3 COVARIATE DATA ########################################################

# access raster
access <- 
  raster("../data-SMA-SAE/Nigeria/2015_accessibility_to_cities_v1.0.tif")
# poverty raster
poverty <- 
  raster("../data-SMA-SAE/Nigeria/nga10povcons200.tif")

#### 1.3.1 cluster covariates
svy_dat_coords <- st_as_sf(svy_dat[, c("lon", "lat")], coords = c(1:2))
svy_dat$access <- extract(access, svy_dat_coords)
svy_dat$l1a <- log(1 + svy_dat$access)
svy_dat$poverty <- extract(poverty, svy_dat_coords)
svy_dat$poverty <- ifelse(is.na(svy_dat$poverty), mean(svy_dat$poverty, na.rm = T), svy_dat$poverty)
svy_dat$urban <- as.factor(svy_dat$urban)
#### 2 DESIGN-BASED METHODS ####################################################
source("analysis/models.R")

if (F) {
  #### 2.1.1 Hajek ####
  sample_des <- svydesign(id = ~cluster + hshold,
                          strata = ~stratum, nest=T, 
                          weights = ~wt, data=svy_dat)
  direct_est <- get_direct(~mcv_yes, ~admin1_char, sample_des)
  
  #### 2.1.1 Model-Assisted ####
  working_fit <- svyglm(mcv_yes ~ urban + l1a + poverty, 
                        sample_des, family = quasibinomial())
  greg_est_list <- list()
  for (i in 1:37) {
    pop_dat_i <- 
      readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
    pop_dat_i$l1a = log(1 + pop_dat_i$access)
    pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
    pop_dat_i$admin1_char <-
      admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
    pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                              c("U", "R"))
    greg_est_list[[i]] <- get_greg(working_fit, mcv_yes ~ urban + l1a + poverty,
                                   ~admin1_char,
                                   pop_dat_i, sample_des, ~svy_1_5_pop)
  }
  greg_est <- do.call(rbind, greg_est_list)
  int_working_fit <- svyglm(mcv_yes ~ (urban + l1a + poverty)*admin1_char, 
                        sample_des, family = quasibinomial())
  igreg_est_list <- list()
  for (i in 1:37) {
    pop_dat_i <- 
      readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
    pop_dat_i$l1a = log(1 + pop_dat_i$access)
    pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
    pop_dat_i$admin1_char <-
      admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
    pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                              c("U", "R"))
    igreg_est_list[[i]] <- 
      get_greg(int_working_fit, mcv_yes ~ (urban + l1a + poverty)*admin1_char,
               ~admin1_char,
               pop_dat_i, sample_des, ~svy_1_5_pop)
  }
  igreg_est <- do.call(rbind, igreg_est_list)
}


#### 3 MODEL-BASED METHODS #####################################################

#### 3.1 BINOMIAL MODELS #######################################################
if (!file.exists(paste0(mcv_res_dir, "iid_bin_ulm_fit.rds"))) {
  set.seed(1201)
  iid_bin_ulm_fit <- get_iid_bin_ulm_fit(mcv_yes ~ urban + l1a + poverty,
                                         ~admin1_char, svy_dat)
  saveRDS(iid_bin_ulm_fit, file = paste0(mcv_res_dir, "iid_bin_ulm_fit.rds"))
} else {
  iid_bin_ulm_fit <- readRDS(paste0(mcv_res_dir, "iid_bin_ulm_fit.rds"))
}
if (!file.exists(paste0(mcv_res_dir, "bym_bin_ulm_fit.rds"))) {
  set.seed(1201)
  bym_bin_ulm_fit <- get_bym_bin_ulm_fit(mcv_yes ~ urban + l1a + poverty,
                                         ~admin1_char, admin1_mat, svy_dat)
  saveRDS(bym_bin_ulm_fit, file = paste0(mcv_res_dir, "bym_bin_ulm_fit.rds"))
} else {
  bym_bin_ulm_fit <- readRDS(paste0(mcv_res_dir, "bym_bin_ulm_fit.rds"))
}


args = commandArgs(TRUE)
# supplied at the command line
i = as.numeric(args[1])
if (!file.exists(paste0(mcv_res_dir, "iid_bin_ulm_est_admin1_", i, ".rds"))) {
  pop_dat_i <- 
    readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
  pop_dat_i$l1a = log(1 + pop_dat_i$access)
  pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
  pop_dat_i$admin1_char <-
    admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
  
  pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                            c("U", "R"))
  iid_bin_ulm_est <- get_iid_bin_ulm_est(mcv_yes ~ urban + l1a + poverty,
                                         ~admin1_char,
                                         iid_bin_ulm_fit$sample_dat,
                                         pop_dat_i, 
                                         iid_bin_ulm_fit$ulm_fit, 
                                         iid_bin_ulm_fit$ulm_fit_sample, 
                                         ~svy_1_5_pop)
  saveRDS(iid_bin_ulm_est, 
          file = paste0(mcv_res_dir, "iid_bin_ulm_est_admin1_", i, ".rds"))
}
if (!file.exists(paste0(mcv_res_dir, "bym_bin_ulm_est_admin1_", i, ".rds"))) {
  pop_dat_i <- 
    readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
  pop_dat_i$l1a = log(1 + pop_dat_i$access)
  pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
  pop_dat_i$admin1_char <-
    admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
  pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                            c("U", "R"))
  bym_bin_ulm_est <- get_bym_bin_ulm_est(mcv_yes ~ urban + l1a + poverty,
                                         ~admin1_char, admin1_mat,
                                         bym_bin_ulm_fit$sample_dat,
                                         pop_dat_i, 
                                         bym_bin_ulm_fit$ulm_fit, 
                                         bym_bin_ulm_fit$ulm_fit_sample, 
                                         ~svy_1_5_pop)
  saveRDS(bym_bin_ulm_est, 
          file = paste0(mcv_res_dir, "bym_bin_ulm_est_admin1_", i, ".rds"))
}
iid_bin_ulm_est_list <- list()
for (i in 1:37) {
  iid_bin_ulm_est_list[[i]] <- readRDS(paste0(mcv_res_dir, "iid_bin_ulm_est_admin1_", i, ".rds"))
}
iid_bin_ulm_est <- do.call(rbind, iid_bin_ulm_est_list)

bym_bin_ulm_est_list <- list()
for (i in 1:37) {
  bym_bin_ulm_est_list[[i]] <- readRDS(paste0(mcv_res_dir, "bym_bin_ulm_est_admin1_", i, ".rds"))
}
bym_bin_ulm_est <- do.call(rbind, bym_bin_ulm_est_list)
#### 3.2 BETABINOMIAL MODELS ###################################################
if (F) {
  svy_cluster_dat <- svy_dat %>%
    dplyr::select(-hshold, -h9, -doi, -dob, -alive) %>%
    group_by(cluster) %>%
    mutate(mcv_yes = sum(mcv_yes), 
           n_trials = sum(n_trials)) %>%
    ungroup() %>%
    unique()
  if (!file.exists(paste0(mcv_res_dir, "iid_bbin_ulm_fit.rds"))) {
    set.seed(1201)
    iid_bbin_ulm_fit <- get_iid_bbin_ulm_fit(mcv_yes ~ urban + l1a + poverty,
                                             ~admin1_char, svy_dat)
    saveRDS(iid_bbin_ulm_fit, file = paste0(mcv_res_dir, "iid_bbin_ulm_fit.rds"))
  } else {
    iid_bbin_ulm_fit <- readRDS(paste0(mcv_res_dir, "iid_bbin_ulm_fit.rds"))
  }
  if (!file.exists(paste0(mcv_res_dir, "bym_bbin_ulm_fit.rds"))) {
    set.seed(1201)
    bym_bbin_ulm_fit <- get_bym_bbin_ulm_fit(mcv_yes ~ urban + l1a + poverty,
                                             ~admin1_char, admin1_mat, svy_dat)
    saveRDS(bym_bbin_ulm_fit, file = paste0(mcv_res_dir, "bym_bbin_ulm_fit.rds"))
  } else {
    bym_bbin_ulm_fit <- readRDS(paste0(mcv_res_dir, "bym_bbin_ulm_fit.rds"))
  }
  
  args = commandArgs(TRUE)
  # supplied at the command line
  i = as.numeric(args[1])
  if (!file.exists(paste0(mcv_res_dir, "iid_bbin_ulm_est_admin1_", i, ".rds"))) {
    pop_dat_i <- 
      readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
    pop_dat_i$l1a = log(1 + pop_dat_i$access)
    pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
    pop_dat_i$admin1_char <-
      admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
    
    pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                              c("U", "R"))
    iid_bbin_ulm_est <- get_iid_bbin_ulm_est(mcv_yes ~ urban + l1a + poverty,
                                             ~admin1_char,
                                             iid_bbin_ulm_fit$sample_dat,
                                             pop_dat_i, 
                                             iid_bbin_ulm_fit$ulm_fit, 
                                             iid_bbin_ulm_fit$ulm_fit_sample, 
                                             ~svy_1_5_pop)
    saveRDS(iid_bbin_ulm_est, 
            file = paste0(mcv_res_dir, "iid_bbin_ulm_est_admin1_", i, ".rds"))
  }
  if (!file.exists(paste0(mcv_res_dir, "bym_bbin_ulm_est_admin1_", i, ".rds"))) {
    pop_dat_i <- 
      readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
    pop_dat_i$l1a = log(1 + pop_dat_i$access)
    pop_dat_i <- filter(pop_dat_i, !is.na(l1a) & !is.na(poverty))
    pop_dat_i$admin1_char <-
      admin1_names$Internal[match(pop_dat_i$admin1_name, admin1_names$GADM)]
    pop_dat_i$urban <- factor(ifelse(pop_dat_i$urban, "U", "R"),
                              c("U", "R"))
    bym_bbin_ulm_est <- get_bym_bbin_ulm_est(mcv_yes ~ urban + l1a + poverty,
                                             ~admin1_char, admin1_mat,
                                             bym_bbin_ulm_fit$sample_dat,
                                             pop_dat_i, 
                                             bym_bbin_ulm_fit$ulm_fit, 
                                             bym_bbin_ulm_fit$ulm_fit_sample, 
                                             ~svy_1_5_pop)
    saveRDS(bym_bbin_ulm_est, 
            file = paste0(mcv_res_dir, "bym_bbin_ulm_est_admin1_", i, ".rds"))
  }
  iid_bbin_ulm_est_list <- list()
  for (i in 1:37) {
    iid_bbin_ulm_est_list[[i]] <- readRDS(paste0(mcv_res_dir, "iid_bbin_ulm_est_admin1_", i, ".rds"))
  }
  iid_bbin_ulm_est <- do.call(rbind, iid_bbin_ulm_est_list)
  
  bym_bbin_ulm_est_list <- list()
  for (i in 1:37) {
    bym_bbin_ulm_est_list[[i]] <- readRDS(paste0(mcv_res_dir, "bym_bbin_ulm_est_admin1_", i, ".rds"))
  }
  bym_bbin_ulm_est <- do.call(rbind, bym_bbin_ulm_est_list)
}
#### 3.3 FAY-HERRIOT MODELS ####################################################

if (F) {
  #### 3.3.1 Smoothed Hajek ####
  iid_sdir_logit_est <- get_iid_sdir_logit(direct_est)
  iid_sdir_est <- get_iid_sdir(direct_est)
  
  bym2_sdir_est <- get_bym2_sdir(direct_est, admin1_mat)
  bym2_sdir_logit_est <- get_bym2_sdir_logit(direct_est, admin1_mat)
  
  #### 3.3.2 Smoothed Model ####
  iid_sgreg_logit_est <- get_iid_sdir_logit(greg_est)
  iid_sgreg_est <- get_iid_sdir(greg_est)
  
  bym2_sgreg_est <- get_bym2_sdir(greg_est, admin1_mat)
  bym2_sgreg_logit_est <- get_bym2_sdir_logit(greg_est, admin1_mat)
  iid_sigreg_logit_est <- get_iid_sdir_logit(igreg_est)
  iid_sigreg_est <- get_iid_sdir(igreg_est)
  
  bym2_sigreg_est <- get_bym2_sdir(igreg_est, admin1_mat)
  bym2_sigreg_logit_est <- get_bym2_sdir_logit(igreg_est, admin1_mat)
  save(direct_est,
       greg_est, 
       iid_sdir_est,
       iid_sdir_logit_est,
       bym2_sdir_est,
       bym2_sdir_logit_est,
       iid_sgreg_est,
       iid_sgreg_logit_est,
       bym2_sgreg_est,
       bym2_sgreg_logit_est,
       iid_sigreg_est,
       iid_sigreg_logit_est,
       bym2_sigreg_est,
       bym2_sigreg_logit_est, 
       iid_bin_ulm_est,
       bym_bin_ulm_est,
       iid_bbin_ulm_est,
       bym_bbin_ulm_est,
       file = paste0(mcv_res_dir, "adm1_est.RData"))
}


#### 4 FIGURES #################################################################
#### 4.1 MAP OF NIGERIA ########################################################
set.seed(1204)
nga_map <- ggplot(data = st_as_sf(poly_adm2)) +
  geom_sf(lwd = .08, fill = NA) + 
  geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
  geom_point(data = svy_dat %>%
               mutate(urban = as.factor(ifelse(urban == "U", "Urban", "Rural"))),
             aes(x = lon, y = lat, color = urban),
             shape = 3, alpha = 1, size = .55) +
  #scale_fill_manual(values = sample(colorRampPalette(c("#008753", "#FFFFFF"))(37))) + 
  scale_color_manual(values = c("mediumblue", "gold"), name = NULL) + 
  guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
  theme_bw() + guides(fill="none") +
  theme(legend.position="bottom",
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/nga_map.pdf"), nga_map,
       width = 7, height = 7)
# ggsave(paste0("paper/figures/nga_map.tiff"),
#        nga_map, dpi = 200, width = 7, height = 7)
#### 4.2 MAPS OF ADMIN-1 ESTIMATES ##############################################
load(paste0(mcv_res_dir, "adm1_est.RData"))
print("loaded estimates")
adm1_est <- bind_rows(direct_est %>%
                        mutate(method = "Hájek"),
                      greg_est %>%
                        mutate(method = "GREG"),
                      iid_bbin_ulm_est %>%
                        mutate(method = "ULM"),
                      bym_bbin_ulm_est %>%
                        mutate(method = "Spatial ULM"),
                      iid_sdir_est %>%
                        mutate(method = "SH alt"),
                      iid_sdir_logit_est %>%
                        mutate(method = "SH"),
                      bym2_sdir_est %>%
                        mutate(method = "Spatial SH alt."),
                      bym2_sdir_logit_est  %>%
                        mutate(method = "Spatial SH"),
                      # iid_sigreg_est %>%
                      #   mutate(method = "SGREG int alt"),
                      # iid_sigreg_logit_est %>%
                      #   mutate(method = "SGREG int"),
                      # bym2_sigreg_est %>%
                      #   mutate(method = "Spatial SGREG int alt"),
                      # bym2_sigreg_logit_est  %>%
                      #   mutate(method = "Spatial SGREG int"),
                      iid_sgreg_est %>%
                        mutate(method = "SGREG alt."),
                      iid_sgreg_logit_est  %>%
                        mutate(method = "SGREG"),
                      bym2_sgreg_est %>%
                        mutate(method = "Spatial SGREG alt"),
                      bym2_sgreg_logit_est  %>%
                        mutate(method = "Spatial SGREG")) %>%
    left_join(admin1_names, by = c("region" = "Internal"))
selected_methods = c("Hájek", "GREG", 
                     "Spatial ULM",
                     "Spatial SH", "Spatial SGREG")

mcv_lims <- range(adm1_est$median)
length_lims <- range(adm1_est$upper - adm1_est$lower)
adm1_maps <- st_as_sf(poly_adm1) %>% 
  dplyr::select(NAME_1) %>%
  left_join(adm1_est, by = c("NAME_1" = "GADM"))
sel_maps <- adm1_maps %>%
  filter(method %in% selected_methods) %>%
  mutate(method = factor(method, 
                         levels = c("Hájek", "Spatial SH", 
                                    "GREG", "Spatial SGREG",
                                    "Spatial ULM")))
all_ests <- ggplot(adm1_maps, aes(fill = median)) + geom_sf(lwd = 0) + 
  facet_wrap(~method) + scale_fill_viridis_c(direction = -1, name = "MCV")  +
  theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/all_adm1_mcv_est_maps.pdf"),
       all_ests, width = 8, height = 7)
# ggsave(paste0("paper/figures/all_adm1_mcv_est_maps.tiff"),
#        all_ests, dpi = 200, width = 8, height = 7)

sel_ests <- ggplot(sel_maps,
                   aes(fill = median)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, name = "MCV", limits = mcv_lims)  +
  facet_wrap(~method, nrow = 3) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/sel_adm1_mcv_est_maps.pdf"),
       sel_ests, width = 6, height = 7)
# ggsave(paste0("paper/figures/sel_adm1_mcv_est_maps.tiff"),
#        sel_ests, dpi = 200, width = 6, height = 7)

#### 4.3 MAPS OF ADMIN-1 CI LENGTHS ############################################

all_lengths <- ggplot(adm1_maps, aes(fill = upper - lower)) + 
  geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, option = "magma",
                       name = "90% CI length") +
  facet_wrap(~method) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())  + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/all_adm1_mcv_len_maps.pdf"), 
       all_lengths, width = 8, height = 7)
# ggsave(paste0("paper/figures/all_adm1_mcv_len_maps.tiff"),
#        all_lengths, dpi = 200, width = 8, height = 7)
sel_lengths <- ggplot(sel_maps,
                   aes(fill = upper - lower)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, option = "magma",
                       name = "90% CI length", limits = length_lims)  +
  facet_wrap(~method, nrow = 3) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("paper/figures/sel_adm1_mcv_len_maps.pdf"), 
       sel_lengths, width = 6, height = 7)
# ggsave(paste0("paper/figures/sel_adm1_mcv_len_maps.tiff"),
#        sel_lengths, dpi = 200, width = 6, height = 7)

comb_maps <- sel_ests + labs(caption = "(a)") + 
  theme(plot.caption = element_text(hjust = .5 , size = 15)) +
  sel_lengths + labs(caption = "(b)") + 
  theme(plot.caption = element_text(hjust = .5 , size = 15))
ggsave(paste0("paper/figures/sel_adm1_mcv_maps.pdf"), 
       comb_maps, width = 12, height = 7)
# ggsave(paste0("paper/figures/sel_adm1_mcv_maps.tiff"),
#        comb_maps, dpi = 200, width = 12, height = 7)

#### 4.4 MAPS OF ADMIN-1 CI LENGTHS ############################################

ht_comp <- adm1_est %>%
  filter(method != "HT") %>%
  left_join(adm1_est %>% 
              filter(method == "HT") %>% 
              dplyr::select(region, est, var) %>% 
              rename(HT = est,
                     HT_var = var),
            by = "region")

gg <- ggplot(ht_comp, aes(x = HT, y = est, color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) 
ggsave(plot = gg, height = 6, width = 7,
       filename = paste0("paper/figures/all_adm1_mcv_est_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_est_scatter.tiff"))
gg <- ggplot(ht_comp, aes(x = sqrt(HT_var), y = sqrt(var), color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) 
ggsave(plot = gg, height = 6, width = 7,
       filename = paste0("paper/figures/all_adm1_mcv_se_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_se_scatter.tiff"))


