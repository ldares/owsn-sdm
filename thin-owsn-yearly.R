library(tidyverse)
library(dismo)
library(corrplot)
library(lubridate)
library(ENMeval)
library(ecospat)
library(spdep)
library(gstat)
library(spThin)
library(MuMIn) #For AICc
library(pROC)
library(kuenm)

set.seed(123)

#####Import Data####
sight <- read.csv("sightings_env_sampled_monthly.csv")

sight1 <- sight%>%
  mutate(date = as_date(sightingdate, format='%m/%d/%Y')) %>%
  mutate(month = month(date)) %>%
  mutate(year = year(date))%>%
  mutate(season = ifelse(month %in% seq(5,9), 'summer', 'winter')) %>%
  filter_at(vars(sst_mean, 
                 cha_mean, 
                 dac_mean, 
                 depth), 
            all_vars(!is.na(.)))

sight2 <- sight1 %>% filter(confidence == 'certain')

View(sight2 %>%
       group_by(species, year)%>%
       summarize(count=n())%>%
       pivot_wider(names_from=year, values_from=count))

sp_of_interest <- c("dalls porpoise", 'fin whale', "grey whale", "harbour porpoise", "humpback whale", "killer whale", "minke whale", "pacific white-sided dolphin", "rissos dolphin", "sperm whale")



win_sight <- sight2 %>% filter(season == 'winter')
smr_sight <- sight2 %>% filter(season == 'summer')

rescale_env <- function(df){
  sst_mst <-  scale(df$sst_mean, center=TRUE, scale=TRUE)
  cha_mst <- scale(df$cha_mean, center = TRUE, scale=TRUE)
  depth_mst <- scale(df$depth, center = TRUE, scale=TRUE)
  df2= cbind(df, sst_mst, cha_mst, depth_mst)
  return(df2)
}

win_sight <- rescale_env(win_sight)
smr_sight <- rescale_env(smr_sight)

smr_sp_list <- split(smr_sight, smr_sight$species)
smr_sp_list <- lapply(smr_sp_list, function(x){
  x %>% mutate(pres = 1)
})

for(i in 1:length(smr_sp_list)){
  spp <- unique(smr_sp_list[[i]]$species)
  df <- smr_sp_list[[i]]
  write.csv(df, paste0(spp, "_summer.csv"), row.names=FALSE)
}

smr_sp_list_yr <- lapply(smr_sp_list, function(x){split(x, x$year)})

win_sp_list <- split(win_sight, win_sight$species)
win_sp_list <- lapply(win_sp_list, function(x){
  x %>% mutate(pres = 1)
})

for(i in 1:length(win_sp_list)){
  spp <- unique(win_sp_list[[i]]$species)
  df <- win_sp_list[[i]]
  write.csv(df, paste0(spp, "_winter.csv"), row.names=FALSE)
}

win_sp_list_yr <- lapply(win_sp_list, function(x){split(x, x$year)})

smr_list_yr <- split(sight2, sight2$year)

#####Thin datasets by year####
thinned_pres_smr <- list()
for (i in 1:length(smr_sp_list_yr)){
  thinned_yr <- list()
  smr_list <- smr_sp_list_yr[[i]]
  for (j in 1:length(smr_list)){
    df <- smr_list[[j]]
    thinned_yr[[j]] <- thin(df, lat.col="latitudedec", long.col="longitudedec", spec.col="species",
                                thin.par=10, reps=100, locs.thinned.list.return = TRUE, write.files=FALSE)
   #plotThin(thinned_yr[[j]])
  }
  thinned_pres_smr[[i]] <- thinned_yr
}
#Give names to each level of thinned_pres_smr
names(thinned_pres_smr) <- sp_of_interest
for (i in 1:length(thinned_pres_smr)){
  yrs <- names(smr_sp_list_yr[[i]])
  names(thinned_pres_smr[[i]]) <- yrs
}

#Save thinned list
saveRDS(thinned_pres_smr,'thinned_pres_smr_yr_10km.RDS')

#Repeat for winter presences
thinned_pres_win <- list()
for (i in 1:length(win_sp_list_yr)){
  thinned_yr <- list()
  win_list <- win_sp_list_yr[[i]]
  for (j in 1:length(win_list)){
    df <- win_list[[j]]
    thinned_yr[[j]] <- thin(df, lat.col="latitudedec", long.col="longitudedec", spec.col="species",
                            thin.par=10, reps=100, locs.thinned.list.return = TRUE, write.files=FALSE)
    #plotThin(thinned_yr[[j]])
  }
  thinned_pres_win[[i]] <- thinned_yr
}

#Give names to each level of thinned_pres_win
names(thinned_pres_win) <- sp_of_interest
for (i in 1:length(thinned_pres_win)){
  yrs <- names(win_sp_list_yr[[i]])
  names(thinned_pres_win[[i]]) <- yrs
}

saveRDS(thinned_pres_win,'thinned_pres_win_yr_10km.RDS')

#####Paste together absences and Re-thin#####
smr_sp_list_thinned <- list()
for (i in 1:length(thinned_pres_smr)){
  smr_sp_list_thinned_yr <- list()
  thinned_pres_smr_sp <- thinned_pres_smr[[i]]
  thinned_abs_smr_sp <- thinned_pres_smr[-i]
  abs_full_smr_sp <- smr_sp_list_yr[-i]
  abs_list_sp_yr <- lapply(thinned_abs_smr_sp, function(x){lapply(x, '[[', sample(seq(1, 100), 1))}) #Randomly selects one thinned df per year for each sp from 100 generated in spThin procedure
  for (yr in names(thinned_pres_smr_sp)){
    pres_loc <- thinned_pres_smr_sp[[yr]][[sample(seq(1, 100), 1)]]
    pres_full <- smr_sp_list_yr[[i]][[yr]]
    pres <- left_join(pres_loc, pres_full, by=c("Longitude"="longitudedec", "Latitude"="latitudedec"))%>%
      distinct(Latitude, Longitude, .keep_all=TRUE)%>%
      mutate(pres=1)
    #bind together 100th thinned dataset for each abs species for each year
    abs_loc_yr_df <- bind_rows(lapply(abs_list_sp_yr, '[[', yr)) %>%
      mutate(species = 'abs')
    #thin absences so all at least 10 km away
    abs_loc_yr_thinned <- thin(abs_loc_yr_df, lat.col='Latitude', long.col = 'Longitude', spec.col = 'species',
                          thin.par=10, reps=100, locs.thinned.list.return = TRUE, write.files=FALSE)
    #join 100th df with abs df for year to get rest of sightings info
    abs_full_yr_df <- bind_rows(lapply(abs_full_smr_sp, '[[', yr))
    abs <- left_join(abs_loc_yr_thinned[[100]], abs_full_yr_df, by=c("Longitude"="longitudedec", "Latitude"="latitudedec"))%>%
      distinct(Latitude, Longitude, .keep_all=TRUE)%>%
      mutate(pres = 0)
    smr_sp_list_thinned_yr[[yr]] <- rbind(pres, abs)
  }
  smr_sp_list_thinned[[i]] <- smr_sp_list_thinned_yr
}
names(smr_sp_list_thinned) <- sp_of_interest

saveRDS(smr_sp_list_thinned,'smr_sp_list_thinned.RDS')

#Thin winter absences
win_sp_list_thinned <- list()
for (i in 1:length(thinned_pres_win)){
  win_sp_list_thinned_yr <- list()
  thinned_pres_win_sp <- thinned_pres_win[[i]]
  thinned_abs_win_sp <- thinned_pres_win[-i]
  abs_full_win_sp <- win_sp_list_yr[-i]
  abs_list_sp_yr <- lapply(thinned_abs_win_sp, function(x){lapply(x, '[[', sample(seq(1, 100), 1))})
  for (yr in names(thinned_pres_win_sp)){
    pres_loc <- thinned_pres_win_sp[[yr]][[sample(seq(1, 100), 1)]]
    pres_full <- win_sp_list_yr[[i]][[yr]]
    pres <- left_join(pres_loc, pres_full, by=c("Longitude"="longitudedec", "Latitude"="latitudedec"))%>%
      distinct(Latitude, Longitude, .keep_all=TRUE)%>%
      mutate(pres=1)
    #bind together 100th thinned dataset for each abs species for each year
    abs_loc_yr_df <- bind_rows(lapply(abs_list_sp_yr, '[[', yr)) %>%
      mutate(species = 'abs')
    #thin absences so all at least 10 km away
    abs_loc_yr_thinned <- thin(abs_loc_yr_df, lat.col='Latitude', long.col = 'Longitude', spec.col = 'species',
                               thin.par=10, reps=100, locs.thinned.list.return = TRUE, write.files=FALSE)
    #join 100th df with abs df for year to get rest of sightings info
    abs_full_yr_df <- bind_rows(lapply(abs_full_win_sp, '[[', yr))
    abs <- left_join(abs_loc_yr_thinned[[100]], abs_full_yr_df, by=c("Longitude"="longitudedec", "Latitude"="latitudedec"))%>%
      distinct(Latitude, Longitude, .keep_all=TRUE)%>%
      mutate(pres = 0)
    win_sp_list_thinned_yr[[yr]] <- rbind(pres, abs)
  }
  win_sp_list_thinned[[i]] <- win_sp_list_thinned_yr
}
names(win_sp_list_thinned) <- sp_of_interest

saveRDS(win_sp_list_thinned,'win_sp_list_thinned.RDS')