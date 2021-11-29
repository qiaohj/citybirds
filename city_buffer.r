setwd("/media/huijieqiao/WD12T/eBird/Script/ebird_migration")
froms = seq(1, 86000, by=1000)
fff<-NA
for (from in froms){
  to = ifelse(from==85001, 85305, from+999)
  fff<-sprintf("../../Data_eBird_2020/city_buffer_data/city_buffer_data_%d_%d.rda", from, to)
  if (file.exists(fff)){
    next()
  }
  saveRDS(NULL, fff)
  break()
}
print(fff)
print(tempdir())

library(raster)
library(rgdal)
library(ggplot2)
library(data.table)
library(Rmisc)
library(landscapemetrics)
library(landscapetools)


mask<-raster("../../Data_eBird_2020/mask_2.5km.tif")
#Pop density
if (F){
  rr<-raster("../../Data_eBird_2020/TIF/GHS_POP_E2015_GLOBE_R2019A_54009_1K_V1_0/GHS_POP_E2015_GLOBE_R2019A_54009_1K_V1_0.tif")
  rr_sinu<-resample(rr, mask)
  writeRaster(rr_sinu, "../../Data_eBird_2020/POP_2.5km.tif")
  
  rr<-raster("../../Data_eBird_2020/TIF/ndvi/mean_ndvi1.tif")
  rr<-setMinMax(rr)
  rr_sinu<-projectRaster(rr, crs=crs(mask), res=res(mask))
  #rr_sinu<-resample(rr, mask, method="bilinear")
  
  writeRaster(rr_sinu, "../../Data_eBird_2020/ndvi_2.5km.tif")
  
  rr<-raster("../../Data_eBird_2020/TIF/HFP2009_footprint.tif")
  rr_sinu<-resample(rr, mask)
  writeRaster(rr_sinu, "../../Data_eBird_2020/HFP2009_2.5km.tif")
  
  city_buffer<-readOGR(dsn="../../Data_eBird_2020/Shape/named_cities", layer="b_cities")
  city_buffer_prj<-spTransform(city_buffer, CRS(proj4string(mask)))
  
  writeOGR(city_buffer_prj, dsn="../../Data_eBird_2020/Shape/named_cities", layer="b_cities_sinu",
           driver="ESRI Shapefile")
}

city_buffer<-readOGR(dsn="../../Data_eBird_2020/Shape/named_cities", layer="b_cities_sinu")
#city_buffer$OBJECTID_1

dim(city_buffer)
length(unique(city_buffer$OBJECTID_1))
length(unique(city_buffer$OBJECTID))
length(unique(city_buffer$OBJECTID_2))
length(unique(city_buffer$NAME_1))

modis<-raster(sprintf("/media/huijieqiao/Speciation_Extin/land_sparing_sharing/Data/MCD12Q1/TIF/%d_LC_Type5.tif",
                      2019))

i=1
handle_v<-function(rr, polygon, type){
  v<-crop(rr, extent(polygon))
  v<-mask(v, polygon)
  v<-values(v)
  v<-v[!is.na(v)]
  if (length(v)>0){
    
    data.table(type=type, id=polygon$OBJECTID, 
               city_name=polygon$NAME_1, 
               country_name=polygon$NAME,
               iso2=polygon$ISO2, 
               iso3=polygon$ISO3,
               mean_v=mean(v), 
               sd=sd(v), CI=CI(v)[1]-CI(v)[2],
               class=1)
  }else{
    NULL
  }
  
}
bind<-function(df, item_df){
  if (is.null(df)){
    item_df
  }else{
    rbindlist(list(df, item_df))
  }
}
city_buffer$NAME<-as.character(city_buffer$NAME)
city_buffer$NAME_1<-as.character(city_buffer$NAME_1)
city_buffer$ISO2<-as.character(city_buffer$ISO2)
city_buffer$ISO3<-as.character(city_buffer$ISO3)



pop<-raster("../../Data_eBird_2020/POP_2.5km.tif")
ndvi<-raster("../../Data_eBird_2020/ndvi_2.5km.tif")
HFP<-raster("../../Data_eBird_2020/HFP2009_2.5km.tif")
endanger_gain_n_events<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_gain_2.5km_N_EVENTS.tif")
endanger_gain_n_observations<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_gain_2.5km_N_OBSERVATIONS.tif")
endanger_gain_n_species<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_gain_2.5km_N_SPECIES.tif")
endanger_loss_n_events<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_loss_2.5km_N_EVENTS.tif")
endanger_loss_n_observations<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_loss_2.5km_N_OBSERVATIONS.tif")
endanger_loss_n_species<-raster("../../Data_eBird_2020/TIF/Turnover/endanger_loss_2.5km_N_SPECIES.tif")

gain_n_events<-raster("../../Data_eBird_2020/TIF/Turnover/gain_2.5km_N_EVENTS.tif")
gain_n_observations<-raster("../../Data_eBird_2020/TIF/Turnover/gain_2.5km_N_OBSERVATIONS.tif")
gain_n_species<-raster("../../Data_eBird_2020/TIF/Turnover/gain_2.5km_N_SPECIES.tif")
loss_n_events<-raster("../../Data_eBird_2020/TIF/Turnover/loss_2.5km_N_EVENTS.tif")
loss_n_observations<-raster("../../Data_eBird_2020/TIF/Turnover/loss_2.5km_N_OBSERVATIONS.tif")
loss_n_species<-raster("../../Data_eBird_2020/TIF/Turnover/loss_2.5km_N_SPECIES.tif")

objs<-list("pop"=pop,
           "ndvi"=ndvi,
           "HFP"=HFP,
           "endanger_gain_n_events"=endanger_gain_n_events,
           "endanger_gain_n_observations"=endanger_gain_n_observations,
           "endanger_gain_n_species"=endanger_gain_n_species,
           "endanger_loss_n_events"=endanger_loss_n_events,
           "endanger_loss_n_observations"=endanger_loss_n_observations,
           "endanger_loss_n_species"=endanger_loss_n_species,
           "gain_n_events"=gain_n_events,
           "gain_n_observations"=gain_n_observations,
           "gain_n_species"=gain_n_species,
           "loss_n_events"=loss_n_events,
           "loss_n_observations"=loss_n_observations,
           "loss_n_species"=loss_n_species)
i=2

all_df<-NULL
for (i in c(from:to)){
#for (i in c(1:nrow(city_buffer))){
  print(paste(i, nrow(city_buffer)))
  item<-city_buffer[i,]
  
  modis_r<-crop(modis, extent(item))
  modis_r<-mask(modis_r, item)
  
  c_pland <- lsm_c_pland(modis_r)
  ## landscape configuration 1: clumpiness index
  c_clumpy <- lsm_c_clumpy(modis_r)
  ## landscape configuration 2: landscape division index
  c_division <- lsm_c_division(modis_r)
  
  c_pland<-data.table(type=c_pland$metric, id=item$OBJECTID, city_name=item$NAME_1,
                      country_name=item$NAME, iso2=item$ISO2,
                      iso3=item$ISO3,mean_v=c_pland$value, sd=0, CI=0, class=c_pland$class)
  if (!is.null(c_pland)){
    all_df<-bind(all_df, c_pland)
  }
  
  c_clumpy<-data.table(type=c_clumpy$metric, id=item$OBJECTID, city_name=item$NAME_1,
                      country_name=item$NAME, iso2=item$ISO2,
                      iso3=item$ISO3,mean_v=c_clumpy$value, sd=0, CI=0, class=c_clumpy$class)
  if (!is.null(c_clumpy)){
    all_df<-bind(all_df, c_clumpy)
  }
  
  c_division<-data.table(type=c_division$metric, id=item$OBJECTID, city_name=item$NAME_1,
                      country_name=item$NAME, iso2=item$ISO2,
                      iso3=item$ISO3,mean_v=c_division$value, sd=0, CI=0, class=c_division$class)
  if (!is.null(c_division)){
    all_df<-bind(all_df, c_division)
  }
  ii<-names(objs)[1]
  for (ii in names(objs)){
    layer_df<-handle_v(objs[[ii]], item, ii)
    if (!is.null(layer_df)){
      all_df<-bind(all_df, layer_df)
    }
  }
  
  
}
  
saveRDS(all_df, fff)
#fix the ndvi error
if (F){
  froms = seq(1, 86000, by=1000)
  fff<-NA
  alldf<-NULL
  from=1
  for (from in froms){
    print(from)
    to = ifelse(from==85001, 85305, from+999)
    fff<-sprintf("../../Data_eBird_2020/city_buffer_data_no_ndvi/city_buffer_data_%d_%d.rda", from, to)
    dfitem<-readRDS(fff)
    i=to
    for (i in c(from:to)){
      print(i)
      item<-city_buffer[i,]
      layer_df<-handle_v(objs[[2]], item, "ndvi")
      dfitem<-bind(dfitem, layer_df)
    }
    saveRDS(dfitem, sprintf("../../Data_eBird_2020/city_buffer_data/city_buffer_data_%d_%d.rda", from, to))
    
  }
}
#merge the data tables
if (F){
  froms = seq(1, 86000, by=1000)
  fff<-NA
  alldf<-NULL
  for (from in froms){
    print(from)
    to = ifelse(from==85001, 85305, from+999)
    fff<-sprintf("../../Data_eBird_2020/city_buffer_data/city_buffer_data_%d_%d.rda", from, to)
    dfitem<-readRDS(fff)
    alldf<-bind(alldf, dfitem)
  }
  
  alldf<-alldf[!is.na(class)]
  dt<-NULL
  types<-unique(alldf[, c("type", "class")])
  i=1
  for (i in c(1:nrow(types))){
    obj<-types[i,]
    item<-alldf[(type==obj$type)&class==obj$class]
    names(item)[c(1, 7:10)]<-paste(obj$type, obj$class, names(item)[c(1, 7:10)], sep="_")
    if (is.null(dt)){
      dt<-item
    }else{
      dt<-merge(dt, item, by=c("id", "city_name", "country_name", "iso2", "iso3"), all=T)
    }
  }
  
  removed_column<-expand.grid(a=c("pland", "clumpy", "division"), b=c(0:11), c=c("sd", "CI", "class", "type"))
  removed_column$x<-paste(removed_column$a, removed_column$b, removed_column$c, sep="_")
  removed_columns<-removed_column$x
  dt_removed<-dt
  dt_removed[, c(removed_columns):=NULL]
  removed_columns<-c(paste(names(objs), 1, "type", sep = "_"), paste(names(objs), 1, "class", sep = "_"))
  dt_removed[, c(removed_columns):=NULL]
  dt_removed<-dt_removed[(!is.na(gain_n_species_1_mean_v))|(!is.na(loss_n_species_1_mean_v))]
  
  
  View(dt)
  dim(dt)
  write.table(dt, "../../Data_eBird_2020/city_buffer_data.csv", row.names = F, sep=",")
  write.table(names(dt), "../../Data_eBird_2020/city_buffer_columns.csv", row.names=F)
  names(dt)
}

