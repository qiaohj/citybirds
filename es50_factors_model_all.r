library(rgdal)
library(raster)
library(data.table)
library(randomForest)
library(ggplot2)
library(caret)


setwd("/media/huijieqiao/WD12T/eBird/Script/ebird_migration")
source("functions.r")

if (F){
  keys<-c("nbiomes", "broad", "urbpop_b7", "US", "UK", "Canada", "AUSTRALIA", 
          "SAfrica", "Thailand")
  
  key<-"nbiomes"
  df_list<-list()
  for (key in keys){
    print(key)
    if (key %in% c("broad", "urbpop_b7")){
      shp<-readOGR(dsn="../../Data_eBird_2020/es50/20210630",
                   layer=key)
    }
    if (key %in% c("nbiomes")){
      shp<-readOGR(dsn="../../Data_eBird_2020/es50/20210616/nbiom101",
                   layer="nbiomes")
    }
    if (key %in% c("US", "UK", "Canada", "AUSTRALIA", "SAfrica", "Thailand")){
      shp<-readOGR(dsn="../../eBird_Pendemic_2021/Shape/inatbirds_countries",
                   layer=key)
    }
    
    colnames(shp@data)
    df<-data.table(shp@data)
    df$AREA<-area(shp)
    
    
    area_threshold<-2e7
    df_sub<-df[AREA>=area_threshold]
    if (key=="urbpop_b7"){
      df_es50<-df_sub[es50>0] 
      hist(df_es50$es50)
      colnames(df_es50)[which(colnames(df_es50)=="es50")]<-"D_es50"
      colnames(df_es50)[which(colnames(df_es50)=="Mean_coast")]<-"coast_mean"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_riv")]<-"riv_mean"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_pop")]<-"pop_mean"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_stdnd")]<-"std_ndvi_m"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_mndvi")]<-"mean_ndvi_"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_fua")]<-"mean_fua"
      colnames(df_es50)[which(colnames(df_es50)=="MEAN_light")]<-"mean_light"
      colnames(df_es50)[which(colnames(df_es50)=="perfo")]<-"perfor"
      colnames(df_es50)[which(colnames(df_es50)=="perim")]<-"per_im"
      
      
    }else{
      df_es50<-df_sub[D_es50>0]  
      #hist(df_es50$D_es50)
    }
    #cols<-c("D_es50", "MEAN_can", "MEAN_dens", "coast_mean", "riv_mean",
    #        "pop_mean", "std_ndvi_m", "mean_ndvi_", "mean_fua",
    #        "mean_light", "perfor", "per_im")
    #df_es50<-df_es50[, ..cols]
    df_list[[key]]<-df_es50
  }
  saveRDS(df_list, "../../eBird_Pendemic_2021/Objects/es50.rda")
}
if (F){
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  names(df_list)
  name<-"US"
  sample_size<-1000
  for (name in names(df_list)){
    print(name)
    df_es50<-df_list[[name]]
    tunegrid <- expand.grid(.mtry=c(1:5))
    train_control <- trainControl(method="repeatedcv", number=5, repeats=10)
    if (nrow(df_es50)<sample_size){
      df_es50_sub<-df_es50
    }else{
      df_es50_sub<-df_es50[sample(nrow(df_es50), sample_size),]
    }
    rf_no_rank <- train(as.formula("D_es50 ~ MEAN_can+MEAN_dens+coast_mean+riv_mean+
                       pop_mean+std_ndvi_m+mean_ndvi_+mean_fua+
                       mean_light+perfor+per_im"), data=df_es50_sub, 
                        trControl=train_control, method="rf", 
                        tuneGrid=tunegrid, ntree=1000)
    saveRDS(rf_no_rank, sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s.rda", name))
  }
}

if (F){
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  names(df_list)
  name<-"US"
  
  for (name in names(df_list)){
    print(name)
    df_es50<-df_list[[name]]
    rf_no_rank<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s.rda", name))
    best_model<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
    model_file<-sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best.rda", name)
    if (file.exists(model_file)){
      rf<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best.rda", name))
    }else{
      rf <- randomForest(D_es50 ~ MEAN_can+MEAN_dens+coast_mean+riv_mean+
                           pop_mean+std_ndvi_m+mean_ndvi_+mean_fua+
                           mean_light+perfor+per_im, 
                         data=df_es50, importance=TRUE, ntree=1000, mtry=best_model$mtry,
                         proximity=TRUE)
      
      saveRDS(rf, sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best.rda", name))  
    }
    
    df_es50$predicted<-rf$predicted
    cor<-cor(df_es50$D_es50, df_es50$predicted)
    
    importances<-data.table(importance(rf))
    colnames(importances)[1]<-"IncMSE"
    importances$mtry<-best_model$mtry
    importances$RMSE<-best_model$RMSE
    importances$Label<-name
  }
}
if (F){
  cols<-c("MEAN_can", "MEAN_dens", "coast_mean", "riv_mean",
          "pop_mean", "std_ndvi_m", "mean_ndvi_", "mean_fua",
          "mean_light", "perfor", "per_im")
  aus_cols<-c(cols, "protected", "resource_p", "mini_use", "grazing", "crop", "plantation", "horticultu",
              "animals", "transition", "INTENSIVE", "MANUFACTUR", "residentia", "SERVICES", "TRANSPORT_",
              "MINING", "WASTE_TREA", "MARSH_WETL", "open_water")
  canada_cols<-c(cols, "tem_ne_f", "sub_ne_f", "tem_bd_f", "mix_f", "tem_sp_sh", "tem_sp_gr", "spp_shl",
                 "spp_gl", "spp_bl", "wetland", "crop", "barren", "urb", "water", "snow")
  thailand_cols<-c(cols, "Agricultur", "Forest", "Mosaic", "Urban", "Water")
  south_africa_cols<-c(cols, "FOREST", "SHRUB", "GRASS", "WATER", "WETLANDS", "BARREN", "CROPS", "URBAN", "MINES")
  usa_cols<-c(cols, "forest", "scrub", "desert", "flooded", "savanna", "grassland", "plantation", "water", "mining", 
              "developed", "tree_crops")
  uk_cols<-c(cols, "grass_1", "marsh_1", "heather_1", "rocky_1", "water_1", "urban_1", "suburban_1", "arable_1", 
             "woodland_1")
  
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  keys<-c("US", "UK", "Canada", "AUSTRALIA", 
          "SAfrica", "Thailand")
  formula_template<-"D_es50 ~ %s"
  key<-"Thailand"
  sample_size<-1000
  for (key in keys){
    print(key)
    df_es50<-df_list[[key]]
    df_es50<-data.frame(df_es50)
    if (key=="US"){
      formula_str<-as.formula(sprintf(formula_template, paste(usa_cols, collapse="+")))  
      col<-"tem_bd_f"
      for (col in usa_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    if (key=="UK"){
      formula_str<-as.formula(sprintf(formula_template, paste(uk_cols, collapse="+")))  
      for (col in uk_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    if (key=="Canada"){
      formula_str<-as.formula(sprintf(formula_template, paste(canada_cols, collapse="+")))  
      for (col in canada_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    if (key=="AUSTRALIA"){
      formula_str<-as.formula(sprintf(formula_template, paste(aus_cols, collapse="+")))  
      for (col in aus_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    if (key=="SAfrica"){
      formula_str<-as.formula(sprintf(formula_template, paste(south_africa_cols, collapse="+")))  
      for (col in south_africa_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    if (key=="Thailand"){
      formula_str<-as.formula(sprintf(formula_template, paste(thailand_cols, collapse="+")))  
      for (col in thailand_cols){
        df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
      }
    }
    tunegrid <- expand.grid(.mtry=c(1:10))
    train_control <- trainControl(method="repeatedcv", number=5, repeats=10)
    if (nrow(df_es50)<sample_size){
      df_es50_sub<-df_es50
    }else{
      df_es50_sub<-df_es50[sample(nrow(df_es50), sample_size),]
    }
    rf_no_rank <- train(formula_str, data=df_es50_sub, 
                        trControl=train_control, method="rf", 
                        tuneGrid=tunegrid, ntree=1000)
    saveRDS(rf_no_rank, sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_all_field.rda", key))
    
  }
  
}


if (F){
  cols<-c("MEAN_can", "MEAN_dens", "coast_mean", "riv_mean",
          "pop_mean", "std_ndvi_m", "mean_ndvi_", "mean_fua",
          "mean_light", "perfor", "per_im")
  aus_cols<-c(cols, "protected", "resource_p", "mini_use", "grazing", "crop", "plantation", "horticultu",
              "animals", "transition", "INTENSIVE", "MANUFACTUR", "residentia", "SERVICES", "TRANSPORT_",
              "MINING", "WASTE_TREA", "MARSH_WETL", "open_water")
  canada_cols<-c(cols, "tem_ne_f", "sub_ne_f", "tem_bd_f", "mix_f", "tem_sp_sh", "tem_sp_gr", "spp_shl",
                 "spp_gl", "spp_bl", "wetland", "crop", "barren", "urb", "water", "snow")
  thailand_cols<-c(cols, "Agricultur", "Forest", "Mosaic", "Urban", "Water")
  south_africa_cols<-c(cols, "FOREST", "SHRUB", "GRASS", "WATER", "WETLANDS", "BARREN", "CROPS", "URBAN", "MINES")
  usa_cols<-c(cols, "forest", "scrub", "desert", "flooded", "savanna", "grassland", "plantation", "water", "mining", 
              "developed", "tree_crops")
  uk_cols<-c(cols, "grass_1", "marsh_1", "heather_1", "rocky_1", "water_1", "urban_1", "suburban_1", "arable_1", 
             "woodland_1")
  
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  keys<-c("US", "UK", "Canada", "AUSTRALIA", 
          "SAfrica", "Thailand")
  formula_template<-"D_es50 ~ %s"
  key<-"Thailand"
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  names(df_list)
  name<-"US"
  
  for (key in keys){
    print(key)
    df_es50<-df_list[[key]]
    df_es50<-data.frame(df_es50)
    rf_no_rank<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_all_field.rda", key))
    best_model<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
    model_file<-sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best_all_field.rda", key)
    
    if (file.exists(model_file)){
      rf<-readRDS(model_file)
    }else{
      if (key=="US"){
        formula_str<-as.formula(sprintf(formula_template, paste(usa_cols, collapse="+")))  
        col<-"tem_bd_f"
        for (col in usa_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      if (key=="UK"){
        formula_str<-as.formula(sprintf(formula_template, paste(uk_cols, collapse="+")))  
        for (col in uk_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      if (key=="Canada"){
        formula_str<-as.formula(sprintf(formula_template, paste(canada_cols, collapse="+")))  
        for (col in canada_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      if (key=="AUSTRALIA"){
        formula_str<-as.formula(sprintf(formula_template, paste(aus_cols, collapse="+")))  
        for (col in aus_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      if (key=="SAfrica"){
        formula_str<-as.formula(sprintf(formula_template, paste(south_africa_cols, collapse="+")))  
        for (col in south_africa_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      if (key=="Thailand"){
        formula_str<-as.formula(sprintf(formula_template, paste(thailand_cols, collapse="+")))  
        for (col in thailand_cols){
          df_es50[, col]<-as.numeric(as.character(df_es50[, col]))
        }
      }
      rf <- randomForest(formula_str, 
                         data=df_es50, importance=TRUE, ntree=1000, mtry=best_model$mtry,
                         proximity=TRUE)
      
      saveRDS(rf, model_file)  
    }
  }
}






