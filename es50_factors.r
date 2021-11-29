library(rgdal)
library(raster)
library(data.table)
library(randomForest)
library(ggplot2)

setwd("/media/huijieqiao/WD12T/eBird/Script/ebird_migration")
source("functions.r")
biomes<-list.files("../../Data_eBird_2020/Shape/newbiomes", pattern = "\\.shp")
biome<-biomes[1]

if (F){
  df<-NULL
  for (biome in biomes){
    print(biome)
    biome<-gsub("\\.shp", "", biome)
    shp<-readOGR(dsn="../../Data_eBird_2020/Shape/newbiomes", layer=biome)
    item<-shp@data
    #item$AREA2<-area(shp)
    if (!("AREA" %in% names(item))){
      print("Error Layers")
      next()
      item$AREA<-area(shp)
    }
    #plot(item$AREA, item$AREA2)
    #hist(item$AREA)
    item$biome<-biome
    item<-data.table(item)
    df<-bind(df, item, fill=TRUE)
    
  }
  saveRDS(df, "../../Data_eBird_2020/es50/es50_factor.rda")
}


df<-readRDS("../../Data_eBird_2020/es50/es50_factor.rda")
area_threshold<-2e7
df_sub<-df[AREA>=area_threshold]
df_es50<-df_sub[D_es50>0]
hist(df_es50$D_es50)

df$D_es50
library(Metrics)

rf <- randomForest(D_es50 ~ mean_coast+max_coast+MAX_riv+MEAN_riv+
                     MAX_dens+MEAN_dens+MAX_can+MEAN_can+
                     MAX+MEAN, 
                   data=df_es50, importance=TRUE,
                        proximity=TRUE)
saveRDS(rf, "../../Data_eBird_2020/es50/es50_model/es50_rf.rda")
rf<-readRDS("../../Data_eBird_2020/es50/es50_model/es50_rf.rda")

print(rf)
## Look at variable importance:
round(importance(rf), 2)
predicted<-predict(rf, df_es50)
df_es50$predicted<-predicted
rmse(df_es50$D_es50, predicted)
Metrics::auc(df_es50$D_es50, df_es50$predicted)
Metrics::mae(df_es50$D_es50, df_es50$predicted)
summary(rf)
saveRDS(predicted, "../../Data_eBird_2020/es50/es50_model/es50_rf_predict.rda")
predicted<-readRDS("../../Data_eBird_2020/es50/es50_model/es50_rf_predict.rda")
plot(df_es50$D_es50, predicted)
cor(df_es50$D_es50, predicted, method = "pearson")
cor(df_es50$D_es50, predicted, method = "kendall")
cor(df_es50$D_es50, predicted, method = "spearman")

biome_label<-unique(df_es50$biome)[1]
for (biome_label in unique(df_es50$biome)){
  print(biome_label)
  df_es50_item<-df_es50[biome==biome_label]
  rf <- randomForest(D_es50 ~ mean_coast+max_coast+MAX_riv+MEAN_riv+
                       MAX_dens+MEAN_dens+MAX_can+MEAN_can+
                       MAX+MEAN, 
                     data=df_es50_item, importance=TRUE,
                     proximity=TRUE)
  saveRDS(rf, 
          sprintf("../../Data_eBird_2020/es50/es50_model/es50_rf_%s.rda", biome_label))
  print(rf)
  ## Look at variable importance:
  round(importance(rf), 2)
  predicted<-predict(rf, df_es50_item)
  #plot(df_es50_item$D_es50, predicted)
  saveRDS(predicted, 
          sprintf("../../Data_eBird_2020/es50/es50_model/es50_rf_predict_%s.rda", biome_label))
}

plot(df_es50$D_es50, rf$predicted)
cor(df_es50$D_es50, predicted)


df_es50
vars<-c("mean_coast", "max_coast", "MAX_riv", "MEAN_riv",
          "MAX_dens", "MEAN_dens", "MAX_can", "MEAN_can",
          "MAX", "MEAN")
v<-vars[1]

df_ggplot<-NULL
for (v in vars){
  myVector<-c("D_es50", "biome", v)
  item<-df_es50[, ..myVector]
  colnames(item)<-c("D_es50", "biome", "v")
  item$v_type<-v
  df_ggplot<-bind(df_ggplot, item)
}
p<-ggplot(df_ggplot[biome!="tundra"], aes(x=D_es50, y=v, colour=factor(biome)))+
  geom_point(size=0.3, alpha=0.3)+
  geom_smooth(method=loess, formula=(y~x))+
  theme_bw()+
  facet_wrap(~v_type, scale="free", ncol=4)
ggsave(p, filename = "../../Data_eBird_2020/Figures/es50_factors/curves_loess.png", width=15, height=10)

rf<-readRDS("../../Data_eBird_2020/es50/es50_model/es50_rf.rda")
rf$importance


lmm<-glm(D_es50 ~ mean_coast+max_coast+MAX_riv+MEAN_riv+
           MAX_dens+MEAN_dens+MAX_can+MEAN_can+
           MAX+MEAN+biome, data=df_es50)
#df_es50$biome
summary(lmm)
lmm$aic
df_es50_item<-df_es50
df_es50_item$pred<-predict(lmm, df_es50)



ggplot(df_es50_item[biome!="tundra"], aes(x=D_es50, y=pred))+
  geom_point(size=0.3, alpha=0.3)+
  theme_bw()+
  facet_wrap(~biome, scale="free", ncol=4)


biome_label<-unique(df_es50$biome)[8]
rf_biome<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_model/es50_rf_%s.rda", biome_label))

item<-df_es50_item[biome==biome_label]
item$pred<-rf_biome$predicted

ggplot(item, aes(x=D_es50, y=pred))+
  geom_point(size=0.3, alpha=0.3)+
  theme_bw()+
  facet_wrap(~biome, scale="free", ncol=4)

lm_result<-NULL
for (biome_label in unique(df_es50$biome)){
  if (biome_label=="tundra"){
    next()
  }
  item<-df_es50[biome==biome_label]
  item$D_es50<-scale(item$D_es50)
  item<-data.frame(item)
  var<-vars[1]
  for (var in vars){
    item[,var]<-scale(item[,var])

    f<-as.formula(sprintf("D_es50 ~ %s", var))
    lmm<-lm(f, data=item)
    s<-summary(lmm)
    intercept<-s$coefficients[1,1]
    slope<-s$coefficients[2,1]
    
    #pred<-predict(lmm, item)
    #plot(item$D_es50, pred)
    r2<-s$r.squared
    f <- s$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    p.type<-""
    if (p<=0.05){
      p.type<-"*"
    }
    if (p<=0.01){
      p.type<-"**"
    }
    if (p<=0.001){
      p.type<-"***"
    }
    
    item_df<-data.table(biome=biome_label, variable=var,
                        intercept=intercept, slope=slope,
                        r2=r2, p.value=p, p.type=p.type)
    if (is.null(lm_result)){
      lm_result<-item_df
    }else{
      lm_result<-rbindlist(list(lm_result, item_df))
    }
  }
}

write.csv(lm_result, "../../Data_eBird_2020/es50/lm_result.csv", row.names = F)
