library(ggplot2)
library(data.table)
library(randomForest)

if (F){
  keys<-c("nbiomes", "broad", "urbpop_b7", "US", "UK", "Canada", "AUSTRALIA", 
          "SAfrica", "Thailand")
  name<-"US"
  df_list<-readRDS("../../eBird_Pendemic_2021/Objects/es50.rda")
  importances<-list()
  for (name in names(df_list)){
    print(name)
    df_es50<-df_list[[name]]
    rf_no_rank<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s.rda", name))
    best_model<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
    rf<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best.rda", name))
    df_es50$predicted<-rf$predicted
    cor<-cor(df_es50$D_es50, df_es50$predicted)
    
    importances_part<-importance(rf)
    importances_part<-data.table(var_name=rownames(importances_part), IncMSE=importances_part[,1],
                                 IncNodePurity=importances_part[,2])
    
    importances_part$mtry<-best_model$mtry
    importances_part$RMSE<-best_model$RMSE
    importances_part$Label<-name
    importances_part$cor<-cor
    importances_part$field<-"part"
    if (file.exists(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_all_field.rda", name))){
      rf_no_rank<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_all_field.rda", name))
      best_model<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
      rf<-readRDS(sprintf("../../eBird_Pendemic_2021/Objects/es50/rf/%s_best_all_field.rda", name))
      df_es50$predicted_all_field<-rf$predicted
      cor<-cor(df_es50$D_es50, df_es50$predicted_all_field)
      importances_full<-importance(rf)
      importances_full<-data.table(var_name=rownames(importances_full), IncMSE=importances_full[,1], 
                                   IncNodePurity=importances_full[,2])
      importances_full$mtry<-best_model$mtry
      importances_full$RMSE<-best_model$RMSE
      importances_full$Label<-name
      importances_full$cor<-cor
      importances_full$field<-"all_field"
      
      importances_item<-rbind(importances_part, importances_full)
    }else{
      importances_item<-importances_part
    }
    importances[[name]]<-importances_item
  }
  
  saveRDS(importances, "../../eBird_Pendemic_2021/Objects/es50/importances.rda")
}

importances<-readRDS("../../eBird_Pendemic_2021/Objects/es50/importances.rda")
keys<-c("nbiomes", "broad", "urbpop_b7", "US", "UK", "Canada", "AUSTRALIA", 
        "SAfrica", "Thailand")
importances<-rbindlist(importances)
importances_part<-importances[field=="part"]
importances_part<-importances_part[Label!="nbiomes"]
p<-ggplot(importances_part)+geom_bar(aes(x=var_name, y=IncMSE), stat = "identity")+
  facet_wrap(~Label, nrow=4)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(p, filename="../../eBird_Pendemic_2021/Figures/importance_part.png", width=8, height=12)


importances_full<-importances[field=="all_field"]
importances_full[var_name %in% c("crop", "Agricultur", "CROPS", "tree_crops", "arable_1")]$var_name<-"crop"
importances_full[var_name %in% c("open_water", "water", "Water", "WATER", "water_1")]$var_name<-"water"
importances_full[var_name %in% c("plantation", "plantations")]$var_name<-"plantation"
importances_full[var_name %in% c("forest", "Forest", "FOREST", "woodland_1")]$var_name<-"forest"
importances_full[var_name %in% c("wetland", "MARSH_WETL", "WETLANDS", "marsh_1", "flooded")]$var_name<-"wetland"
importances_full[var_name %in% c("grass_1", "grassland", "GRASS")]$var_name<-"grassland"
importances_full[var_name %in% c("urb", "Urban", "urban_1", "URBAN")]$var_name<-"urban"
importances_full[var_name %in% c("scrub", "SHRUB")]$var_name<-"scrub"
df<-unique(rbind(importances_full, importances_part)[, ..cols])

importances_full[order(Label, -IncMSE), indx := seq_len(.N), "Label"]
importances_full$top_5<-F
importances_full[indx<=5]$top_5<-T
p<-ggplot(importances_full)+geom_bar(aes(x=var_name, y=IncMSE, fill=top_5), stat = "identity")+
  facet_wrap(~Label, nrow=4)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
ggsave(p, filename="../../eBird_Pendemic_2021/Figures/importance_full.png", width=15, height=8)
write.csv(importances_full, "../../eBird_Pendemic_2021/Figures/importance_full.csv")
write.csv(importances_part, "../../eBird_Pendemic_2021/Figures/importances_part.csv")
cols<-c("Label", "cor", "field")
df<-unique(rbind(importances_full, importances_part)[, ..cols])
ggplot()+geom_p
