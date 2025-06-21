library(tidyverse)
library(sf)
library(PresenceAbsence)
library(dismo)
library(raster)
library(jaccard)

####Import Sightings and Predictions####
best_mod_list <- readRDS("best_mod_list_smr.RDS")

#####Extract prediction layers#####
pa_list <- list()
thresh_list <- list()
pa_preds_list <- list()
for(sp in c('humpback whale', 'fin whale', 'harbour porpoise', 'dalls porpoise')){
  pa_list[[sp]] <- list()
  thresh_list[[sp]] <- list()
  pa_preds_list[[sp]] <- list()
  for(i in 1:length(best_mod_list[[sp]])){
    e.mx <- readRDS(paste0("eval_obj_smr_list-", sp, "_", i, ".RDS"))
    best_mod <- best_mod_list[[sp]][[i]]
    tune_args <- as.character(best_mod$tune.args)
    occs <- e.mx@occs%>%
      rownames_to_column(., var="rows")%>%
      mutate(pres = 1,
             index = paste(sp, i, "pres", rows, sep="_"))
    bg <- e.mx@bg%>%
      rownames_to_column(., var="rows")%>%
      mutate(pres = 0,
             index = paste(sp, i, "abs", rows, sep="_"))
    points <- rbind(occs, bg)%>%
      st_as_sf(coords = c("Longitude", "Latitude"),
               crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    
    points$preds <- raster::extract(e.mx@predictions[[tune_args]], points)
    pa_list[[sp]][[i]] <- points %>% dplyr::select(index, pres, preds)
    thresh_list[[sp]][[i]] <- error.threshold.plot(DATA=as.data.frame(pa_list[[sp]][[i]] %>%                                           dplyr::select(-geometry)),
                                                   which.model=1,
                                                   #main="SPlow",
                                                   na.rm=TRUE,
                                                   color=TRUE,
                                                   opt.thresholds=TRUE,
                                                   opt.methods="MaxSens+Spec")
    criteria <- c(0,thresh_list[[sp]][[i]]$threshold, 0,
                  thresh_list[[sp]][[i]]$threshold, 1, 1)
    crit_matrix <- matrix(criteria, ncol=3, byrow=TRUE)
    pa_preds_list[[sp]][[i]] <- reclassify(e.mx@predictions[[tune_args]], crit_matrix)
  }
}

saveRDS(thresh_list, "sdm_thresh_list.RDS")

pa_preds_stacks <- lapply(pa_preds_list, stack)

pa_preds_ensembles <- lapply(pa_preds_stacks, function(x){calc(x, sum)})
plot(pa_preds_ensembles[[1]])
plot(pa_preds_ensembles[[2]])
plot(pa_preds_ensembles[[3]])
plot(pa_preds_ensembles[[4]])

for(sp in names(pa_preds_ensembles)){
  rast <- pa_preds_ensembles[[sp]]
  writeRaster(rast, paste0(sp, "_pa_preds_ensembles.tif"))
  rast[rast[] < 5] = NA
  writeRaster(rast, paste0(sp, "_pa_preds_ens5.tif"))
}

saveRDS(pa_preds_ensembles, "pa_preds_ensembles.RDS")

saveRDS(pa_preds_stacks, "pa_preds_stacks.RDS")

####Summer 2018 Predictions#####
#Import rasters for 2018 summer means of each predictor
envs_2018 <- raster::stack(c(raster("Sea Surface Temperature_smr_2018_mean.tif"),
                        raster("Chlorophyll-A_smr_2018_mean.tif"),
                        raster("gebco_2021_resampled_4km.tif")))

names(envs_2018) <- names(envs)

pa_list <- list()
thresh_list <- list()
pa_preds_list <- list()
for(sp in c('humpback whale', 'fin whale', 'harbour porpoise', 'dalls porpoise')){
  pa_list[[sp]] <- list()
  thresh_list[[sp]] <- list()
  pa_preds_list[[sp]] <- list()
  for(i in 1:length(best_mod_list[[sp]])){
    e.mx <- readRDS(paste0("eval_obj_smr_list-", sp, "_", i, ".RDS"))
    best_mod <- best_mod_list[[sp]][[i]]
    tune_args <- as.character(best_mod$tune.args)
    occs <- e.mx@occs%>%
      rownames_to_column(., var="rows")%>%
      mutate(pres = 1,
             index = paste(sp, i, "pres", rows, sep="_"))
    bg <- e.mx@bg%>%
      rownames_to_column(., var="rows")%>%
      mutate(pres = 0,
             index = paste(sp, i, "abs", rows, sep="_"))
    points <- rbind(occs, bg)%>%
      st_as_sf(coords = c("Longitude", "Latitude"),
               crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    sp_mod <- maxent(envs, 
                     occs %>% dplyr::select(c(Longitude, Latitude)), 
                     bg %>% dplyr::select(c(Longitude, Latitude)), 
                     tune_args, 
                     factors=NULL, 
                     removeDuplicates=TRUE)
    preds <- predict(sp_mod, envs_2018)
    #writeRaster(ens_map, paste0(getwd(), "/SDM Null Ensembles/", sp, "_null_pa_preds_ens_", i, ".tif"), overwrite=TRUE)
    points$preds <- raster::extract(preds, points)
    pa_list[[sp]][[i]] <- points %>% dplyr::select(index, pres, preds)
    thresh_list[[sp]][[i]] <- error.threshold.plot(DATA=as.data.frame(pa_list[[sp]][[i]] %>%                                           dplyr::select(-geometry)),
                                                   which.model=1,
                                                   #main="SPlow",
                                                   na.rm=TRUE,
                                                   color=TRUE,
                                                   opt.thresholds=TRUE,
                                                   opt.methods="MaxSens+Spec")
    criteria <- c(0,thresh_list[[sp]][[i]]$threshold, 0,
                  thresh_list[[sp]][[i]]$threshold, 1, 1)
    crit_matrix <- matrix(criteria, ncol=3, byrow=TRUE)
    pa_preds_list[[sp]][[i]] <- reclassify(e.mx@predictions[[tune_args]], crit_matrix)
  }
}

saveRDS(thresh_list, "sdm_thresh_list_2018.RDS")

pa_preds_stacks <- lapply(pa_preds_list, stack)

pa_preds_ensembles <- lapply(pa_preds_stacks, function(x){calc(x, sum)})
#Fix HB ensemble - only one subset performed better than null models
pa_preds_ensembles[["humpback whale"]] <- pa_preds_stacks[["humpback whale"]][[3]]*5

plot(pa_preds_ensembles[[1]])
plot(pa_preds_ensembles[[2]])
plot(pa_preds_ensembles[[3]])
plot(pa_preds_ensembles[[4]])

pa_preds_ens5_list <- list()
for(sp in names(pa_preds_ensembles)){
  rast <- pa_preds_ensembles[[sp]]
  writeRaster(rast, paste0(sp, "_pa_preds_ensembles_2018.tif"), overwrite=TRUE)
  rast[rast[] < 5] = NA
  writeRaster(rast, paste0(sp, "_pa_preds_ens5_2018.tif"), overwrite=TRUE)
  pa_preds_ens5_list[[sp]] <- rast
}

###Import Wright et al Rasters#####
wright_files <- list.files("Spatial Model Rasters", pattern = "rescaled.*\\.tif$", full.names=TRUE)

wright_list <- lapply(wright_files, raster)

#Convert Wright et al densities to presence-absence using 10% threshold
wright_list_10p <- lapply(wright_list, function(x){
  val_matrix <- matrix(c(0,0.1, 0, 0.1, 1, 1), ncol=3, byrow=TRUE)
  reclassify(x, val_matrix)
})

lapply(wright_list_10p, plot)

names(wright_list_10p) <- c("dalls porpoise", "fin whale", "harbour porpoise", "humpback whale")

for(sp in names(wright_list_10p)){
  lyr <- wright_list_10p[[sp]]
  writeRaster(lyr, paste0("./Spatial Model Density Threshold Rasters/", sp, "_DSM_pa_10p.tif", overwrite=TRUE))
}

#####Import Wright et al centroids####
#Grid cell centroids derived in QGIS
centroids_files <- list.files("./Spatial Model Data csvs", pattern=".shp$", full.names=TRUE)
centroids_list <- lapply(centroids_files, read_sf)
names(centroids_list) <- c("dalls porpoise", "fin whale", "harbour porpoise", "humpback whale")


pa_preds_ens5_list1 <- lapply(pa_preds_ens5_list, function(x){
  #y <- calc(x, function(a){a/5})
  z <- projectRaster(x, crs = crs(wright_list_10p[["dalls porpoise"]]))
  return(z)
})
plot(pa_preds_ens5_list1[[1]])
plot(pa_preds_ens5_list1[[2]])
plot(pa_preds_ens5_list1[[3]])
plot(pa_preds_ens5_list1[[4]])


###Sample Both rasters####
sampled_list <- list()
for (sp in names(wright_list_10p)){
  pts <- centroids_list[[sp]]
  sdm_sampled <- raster::extract(pa_preds_ens5_list1[[sp]], pts)
  wright_sampled <- raster::extract(wright_list_10p[[sp]], pts)
  df <- cbind(pts %>% dplyr::select(c(X, Y)), wright_sampled, sdm_sampled)%>%
    mutate(sdm_sampled_pa = case_when(sdm_sampled == 5~1, 
                                      TRUE~0),
           wright_sampled_pa = ifelse(wright_sampled==1, 1, 0))
  sampled_list[[sp]] <- df
}

#Jaccard Test####
jaccard_list <- lapply(sampled_list, function(x){
  x <- x %>%
    dplyr::filter(!is.na(sdm_sampled_pa), !is.na(wright_sampled_pa))
  jaccard.test(x$sdm_sampled_pa, x$wright_sampled_pa, method="bootstrap")
})

saveRDS(jaccard_list, "jaccard_tests.RDS")