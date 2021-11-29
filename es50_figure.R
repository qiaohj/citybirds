library(tibble)
library(dggridR)
library(dplyr)
library(purrr)
library(entropart)
detach("package:gbifrasters", unload=TRUE)
library(gbifrasters)
library(sf)


setwd("/media/huijieqiao/WD12T/eBird/Script/ebird_migration")
spacing<-20
grid<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", spacing))
grid_no_na<-grid%>%filter(!is.na(es50))
grid_no_na$var<-grid_no_na$es50
breaks<-round(quantile(grid_no_na$D_es50, seq(0, 1, by=0.1)))
Start = 2:length(breaks) %>% map_dbl(~breaks[.x-1])
Finish = 1:(length(breaks)-1) %>% map_dbl(~breaks[.x+1])

labels<-paste0(Start,"-",Finish)
breaks[11]<-50.1
grid_no_na$xx<-cut(grid_no_na$var,breaks,labels)
grid_no_na[which(is.na(grid_no_na$xx)),]
p = gbifrasters::plotPolyMap(grid_no_na,
                             variable="ES50",
                             legend.position = "bottom",
                             breaks=breaks,
                             labels=labels,
                             polygon_text_size=3,
                             polygon_alpha = 1,
                             labelType="",
                             keywidth=0.01,
                             keyheight=0.3,
                             legend_text_size=15,
                             path_alpha=0,
                             mappath="/media/huijieqiao/WD12T/GISLayers/continents/continent.shp"
)
#p
ggsave(p, filename=sprintf("../../Data_eBird_2020/Figures/es50_maps/spacing_%d.png", spacing),
       width=12, height=8)

#"US", "UK", "Canada", "AUSTRALIA", "SAfrica", "Thailand"
target_countries<-c("USA", "GBR", "ZAF", "THA", "AUS")
countries<-st_read("/media/huijieqiao/WD12T/GISLayers/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")


grid_3<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", 3))
grid_no_na_3 = grid_3 %>% filter(!is.na(D_es50))
grid_no_na_3$var<-grid_no_na_3$D_es50
country<-"USA"
pp<-p
for (country in target_countries){
  print(country)
  layer<-countries[which(countries$ISO3==country),]
  box<-matrix(st_bbox(layer), ncol=2)
  if (country=="USA"){
    box[2, 2]<-60
    box[1, 1]<--141
    box[1, 2]<--52
    abline<-data.frame(x=c(box[1, 1], box[1, 2]), y=c(49, 49))
  }
  if (country=="AUS"){
    box[2, 1]<--43.7
    box[1, 2]<-154.5
  }
  if (country=="ZAF"){
    box[2, 1]<--35.7
    box[1, 2]<-33.4
  }
  
  boundary<-data.frame(x=c(box[1, 1], box[1, 1], box[1, 2], box[1, 2], box[1, 1]),
                       y=c(box[2, 1], box[2, 2], box[2, 2], box[2, 1], box[2, 1]))
  if (F){
    countriesx = ggplot2::map_data("world")
    ggplot()+
      geom_polygon(data=countriesx,aes(x=long, y=lat, group=group), 
                   fill=NA, color="#D8DACF",alpha=0.8)+
      
      geom_path(data=boundary, aes(x=x, y=y))+
      geom_path(data=abline, aes(x=x, y=y), linetype=2)
    
  }
  pp<-pp+geom_path(data=boundary, aes(x=x, y=y), color="black", alpha=0.6)
  if (country=="USA"){
    pp<-pp+geom_path(data=abline, aes(x=x, y=y), linetype=2, alpha=0.6)
  }
  
  #next()
  
  
  width<-10
  height<-width * ((box[2, 2] - box[2, 1])/(box[1, 2] - box[1, 1]))
  grid_item<-grid_no_na_3[which((grid_no_na_3$long>=box[1, 1])&
                                  (grid_no_na_3$long<=box[1, 2])&
                                  (grid_no_na_3$lat>=box[2, 1])&
                                  (grid_no_na_3$lat<=box[2, 2])),]
  p_item = gbifrasters::plotPolyMap(grid_item,
                                    variable="D_es50",
                                    legend.position = "none",
                                    breaks=breaks,
                                    labels=labels,
                                    zoom_x=c(box[1, 1], box[1, 2]),
                                    zoom_y=c(box[2, 1], box[2, 2]),
                                    polygon_text_size=3,
                                    polygon_alpha = 1,
                                    labelType="",
                                    keywidth=0.01,
                                    keyheight=0.3,
                                    legend_text_size=15,
                                    path_alpha=0,
                                    mappath="/media/huijieqiao/WD12T/GISLayers/continents/continent.shp"
  )
  ggsave(p_item, filename = sprintf("../../Data_eBird_2020/Figures/es50_maps/%s.png", country),
         width=width, height=height)
}
#pp
ggsave(pp, filename=sprintf("../../Data_eBird_2020/Figures/es50_maps/spacing_%d_with_box.png", spacing),
       width=12, height=8)
