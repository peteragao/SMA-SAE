library(raster)
library(dplyr)
library(tidyr)
library(raster)

get_thresh <- function(this_area_pop, adm1, poppa) {
  # based on John Paige's getRegionThresh from SUMMER package
  prop_urb <- poppa$pctUrb[poppa[["GADMarea"]] == adm1] / 100
  prop_rural <- 1 - prop_urb
  if (round(prop_rural, 10) == 1) {
    return(-Inf)
  } else if (round(prop_rural, 10) == 0) {
    return(Inf)
  }
  area_pop <- this_area_pop
  total_area_pop <- sum(area_pop)
  
  area_pop_sorted = sort(area_pop)
  cumsum_area_pop = cumsum(area_pop_sorted)
  
  thresh_idx = match(1, cumsum_area_pop >= total_area_pop * prop_rural)
  if ((thresh_idx != 1) & (thresh_idx != length(area_pop))) {
    thresh = area_pop_sorted[thresh_idx]
  } else {
    # make sure not all pixels are urban or all are rural
    if(thresh_idx == 1) {
      thresh = mean(c(area_pop_sorted[1], area_pop_sorted [2]))
    } else {
      thresh = mean(c(area_pop_sorted [length(area_pop)], 
                      area_pop_sorted[length(area_pop) - 1]))
    }
  }
  return(thresh)
}

# return both population summary as well as urban/rural raster (?)
load_pop_data <- function(poly_adm1, poly_adm2) {
  
  # Age 1-5 population rasters by sex 
  # worldpop estimates for year of census
  cen_f1_pop <- raster("data/Nigeria/nga_f_1_2006.tif")
  cen_m1_pop <- raster("data/Nigeria/nga_m_1_2006.tif")
  cen_1_5_pop <- calc(stack(cen_f1_pop, cen_m1_pop), sum)
  
  # worldpop estimates for year of survey
  svy_f1_pop <- raster("data/Nigeria/nga_f_1_2018.tif")
  svy_m1_pop <- raster("data/Nigeria/nga_m_1_2018.tif")
  svy_1_5_pop <- calc(stack(svy_f1_pop, svy_m1_pop), sum)
  
  # Total population rasters
  # worldpop estimates for year of census
  cen_pop <- raster("data/Nigeria/nga_ppp_2006_UNadj.tif")
  
  # worldpop estimates for year of survey
  svy_pop <- raster("data/Nigeria/nga_ppp_2018_UNadj.tif")
  
  # access raster
  access <- raster("data/Nigeria/2015_accessibility_to_cities_v1.0.tif")
  # poverty raster
  poverty <- raster("data/Nigeria/nga10povcons200.tif")
  
  svy_poppa <- read.csv("data/Nigeria/nga_2018_poppa.csv")
  
  pop_stk <- stack(cen_1_5_pop, svy_1_5_pop, cen_pop, svy_pop)
  names(pop_stk) <- c("cen_1_5_pop", "svy_1_5_pop", "cen_pop", "svy_pop")
  
  adm1_pix_list <- list()
  adm1_pop_totals_list <- list()
  # make pixel table
  adm1_names <- as.character(poly_adm1$NAME_1)
  for (i in 1:length(adm1_names)) {
    adm1 <- adm1_names[i]
    print(adm1)
    adm1_pix <-
      mask(crop(pop_stk, subset(poly_adm2, NAME_1 == adm1), snap='out'),
           subset(poly_adm2, NAME_1 == adm1))
    adm1_pix <- st_join(st_as_sf(rasterToPoints(adm1_pix, spatial = T)),
                        st_as_sf(poly_adm2)[, c("NAME_1", "NAME_2")])
    adm1_pix <- rename(adm1_pix, "admin1_name" = "NAME_1")
    adm1_pix <- rename(adm1_pix, "admin2_name" = "NAME_2")
    
    # USE SURVEY POPULATION TO GET URBAN/RURAL THRESHOLD?
    thresh <- get_thresh(adm1_pix$cen_pop, adm1, svy_poppa)
    adm1_pix$urban <- adm1_pix$cen_pop >= thresh
    
    
    adm1_pix$access <- extract(access, adm1_pix)
    adm1_pix$poverty <- extract(poverty, adm1_pix)
    
    adm1_pix_coords <- st_coordinates(adm1_pix)
    adm1_pix$lon <- adm1_pix_coords[, 1] 
    adm1_pix$lat <- adm1_pix_coords[, 2]
    adm1_pix <- st_set_geometry(adm1_pix, NULL)
    adm1_pop_totals <-
      aggregate(.~admin1_name + admin2_name,
                dplyr::select(adm1_pix, -lon, -lat, -urban, -access, -poverty), sum)
    adm1_urban_pop_totals <-
      aggregate(.~admin1_name + admin2_name + urban,
                dplyr::select(adm1_pix, -lon, -lat, -access, -poverty), sum)
    adm1_urban_pop_totals <- subset(adm1_urban_pop_totals, urban == T)
    adm1_urban_pop_totals <- dplyr::select(adm1_urban_pop_totals, -urban)
    colnames(adm1_urban_pop_totals)[3:6] <- 
      paste0(colnames(adm1_urban_pop_totals)[3:6], "_urb")
    adm1_pop_totals <- merge(adm1_pop_totals, 
                             adm1_urban_pop_totals[, ], all.x = T)
    adm1_pop_totals[is.na(adm1_pop_totals)] <- 0
    saveRDS(adm1_pix, 
            paste0("data/Nigeria/nga_pix_tbl_admin1_", i, ".rds"))
    adm1_pix_list[[i]] <- adm1_pix
    adm1_pop_totals_list[[i]] <- adm1_pop_totals
    
    
  }
  saveRDS(do.call(rbind, adm1_pix_list), 
          paste0("data/Nigeria/nga_pix_tbl.rds"))
  adm1_pop_totals <- do.call(rbind, adm1_pop_totals_list)
  adm1_pop_totals <- aggregate(.~admin1_name + admin2_name, data = adm1_pop_totals, sum)
  saveRDS(adm1_pop_totals, 
          paste0("data/Nigeria/nga_pop_totals.rds"))
}
