library(gbifrasters)
library(tibble)
library(dggridR)
library(dplyr)
library(purrr)
library('entropart')

setwd("/media/huijieqiao/WD12T/eBird/Script/ebird_migration")


# NOT RUN {

if (F){
  spacing<-10
  dggs          <- dgconstruct(spacing=spacing, metric=FALSE) 
  
  fff<-list.files("../../Data_eBird_2020/Tables/Species_With_Urban", 
                  pattern = "\\.rda", full.names = T)
  i=1
  D_df<-NULL
  for (i in c(1:length(fff))){
    print(paste(i, length(fff)))
    f<-fff[i]
    sp_df<-readRDS(f)
    if (is.null(sp_df)){
      next()
    }
    if (nrow(sp_df)==0){
      next()
    }
    sp_df$cell <- dgGEO_to_SEQNUM(dggs,sp_df$LONGITUDE, sp_df$LATITUDE)$seqnum
    saveRDS(sp_df, gsub("Species_With_Urban", 
                        sprintf("Species_with_Urban_cell_index_%d", spacing), f))
    D<-sp_df%>%dplyr::group_by(cell)%>%
      dplyr::summarise(N=n(), sp_id=i, sp=sp_df[1, "SCIENTIFIC_NAME"])
    if (is.null(D_df)){
      D_df<-D
    }else{
      D_df<-bind_rows(D_df, D)
    }
  }
  
  saveRDS(D_df, sprintf("../../Data_eBird_2020/es50/es50_df_%d.rda", spacing))

}
if (F){
  dggs_100 <- dggridR::dgconstruct(spacing=100, metric=FALSE)
  polygongrid_100<-gbifrasters::getPolygonGrid(dggs_100, 100,landOnly=FALSE)
  D_df_100<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_df_%d.rda", 100))
  D_df_100$taxonkey<-D_df_100$sp_id
  D_df_100$count<-D_df_100$N
  D_df_100$cell_raw<-D_df_100$cell
  grid_100 = gbifrasters::mergeToGrid(polygongrid_100, D_df_100)
  saveRDS(grid_100, sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", 100))
}
if (F){
  spacing<-50
  dggs <- dggridR::dgconstruct(spacing=spacing, metric=FALSE)
  polygongrid<-gbifrasters::getPolygonGrid(dggs, spacing,landOnly=FALSE)
  D_df<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_df_%d.rda", spacing))
  D_df$taxonkey<-D_df$sp_id
  D_df$count<-D_df$N
  D_df$cell_raw<-D_df$cell
  grid = gbifrasters::mergeToGrid(polygongrid, D_df)
  
  
  D<-D_df
  esNum<-50
  
  D_group<-D %>%
    group_split(cell) %>%
    map(~ .x %>%
          group_by(taxonkey) %>%
          summarize(occCount = sum(count))
    ) %>%
    map(~ deframe(.x))
  
  cell = D %>%
    group_split(cell) %>%
    map(~ .x %>% mutate(cell = as.character(cell))) %>%
    map_chr(~ unique(.x$cell))
  
  cellSpCounts <- D %>%
    group_by(cell) %>%
    summarise(spCount=n()) %>%
    mutate(cell = as.character(cell))
  cellOccCounts <- D %>%
    group_by(cell) %>%
    summarise(occCounts = sum(count))
  
  
  D_es50<-D_group%>%
    modify_if(~ length(.x) <= esNum, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, esNum)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()
  
  es50Table = tibble(cell,D_es50)
  
  D_es20<-D_group%>%
    modify_if(~ length(.x) <= 20, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 20)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()
  
  es20Table = tibble(cell,D_es20)
  
  D_es100<-D_group%>%
    modify_if(~ length(.x) <= 100, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 100)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()
  
  es100Table = tibble(cell,D_es100)
  
  D_es200<-D_group%>%
    modify_if(~ length(.x) <= 200, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 200)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()
  
  es200Table = tibble(cell,D_es200)
  
  D_es1000<-D_group%>%
    modify_if(~ length(.x) <= 1000, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 1000)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()
  
  es1000Table = tibble(cell,D_es1000)
  
  
  grid <- dggridR::dgcellstogrid(dggs,D$cell,frame=TRUE,wrapcells=TRUE)
  
  grid = merge(grid,cellSpCounts,id="cell",all.x=TRUE) %>%
    merge(es50Table,id="cell",all.x=TRUE) %>%
    merge(es20Table,id="cell",all.x=TRUE) %>%
    merge(es100Table,id="cell",all.x=TRUE) %>%
    merge(es200Table,id="cell",all.x=TRUE) %>%
    merge(es1000Table,id="cell",all.x=TRUE) %>%
    merge(cellOccCounts,id="cell",all.x=TRUE) %>%
    mutate(spCountLog10 = log10(spCount)) %>%
    mutate(occCountsLog10 = log10(occCounts)) %>%
    arrange(cell, order)
  
  
  saveRDS(grid, sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", spacing))
  
  grid<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", spacing))
  grid_no_na<-grid%>%filter(!is.na(D_es50))
  grid_no_na$var<-grid_no_na$D_es50
  
  detach("package:gbifrasters", unload=TRUE)
  library(gbifrasters)
  
  p = gbifrasters::plotPolyMap(grid_no_na,
                               variable="ES50",
                               legend.position = "bottom",
                               breaks=NULL,
                               labels=NULL,
                               polygon_text_size=3,
                               polygon_alpha = 1,
                               labelType="",
                               keywidth=0.01,
                               keyheight=0.3,
                               legend_text_size=15,
                               path_alpha=0,
                               mappath="/media/huijieqiao/WD12T/GISLayers/continents/continent.shp"
  )
  p
  
  ggsave(p, filename=sprintf("../../Data_eBird_2020/Figures/es50_maps/spacing_%d.png", spacing),
         width=12, height=8)
}

grid_100<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", 100))
spacing<-3
dggs <- dggridR::dgconstruct(spacing=spacing, metric=FALSE)

D_df<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_df_%d.rda", spacing))

D_df$taxonkey<-D_df$sp_id
D_df$count<-D_df$N
D_df$cell_raw<-D_df$cell
#min(abs(polygongrid$cell- D_df[1, "cell"]))
#polygongrid%>%dplyr::filter(!(cell %in% D_df$cell))
#D_df%>%dplyr::filter(!(cell %in% polygongrid$cell))


D<-D_df
esNum<-50

D_group<-D %>%
  group_split(cell) %>%
  map(~ .x %>%
        group_by(taxonkey) %>%
        summarize(occCount = sum(count))
  ) %>%
  map(~ deframe(.x))

cell = D %>%
  group_split(cell) %>%
  map(~ .x %>% mutate(cell = as.character(cell))) %>%
  map_chr(~ unique(.x$cell))

cellSpCounts <- D %>%
  group_by(cell) %>%
  summarise(spCount=n()) %>%
  mutate(cell = as.character(cell))
cellOccCounts <- D %>%
  group_by(cell) %>%
  summarise(occCounts = sum(count))


D_es50<-D_group%>%
  modify_if(~ length(.x) <= esNum, ~ NA) %>% # run only if more than 50 species
  modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, esNum)) %>%
  map(~ unname(.x)) %>%
  flatten_dbl()

es50Table = tibble(cell,D_es50)

D_es20<-D_group%>%
  modify_if(~ length(.x) <= 20, ~ NA) %>% # run only if more than 50 species
  modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 20)) %>%
  map(~ unname(.x)) %>%
  flatten_dbl()

es20Table = tibble(cell,D_es20)

D_es100<-D_group%>%
  modify_if(~ length(.x) <= 100, ~ NA) %>% # run only if more than 50 species
  modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 100)) %>%
  map(~ unname(.x)) %>%
  flatten_dbl()

es100Table = tibble(cell,D_es100)

D_es200<-D_group%>%
  modify_if(~ length(.x) <= 200, ~ NA) %>% # run only if more than 50 species
  modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 200)) %>%
  map(~ unname(.x)) %>%
  flatten_dbl()

es200Table = tibble(cell,D_es200)

D_es1000<-D_group%>%
  modify_if(~ length(.x) <= 1000, ~ NA) %>% # run only if more than 50 species
  modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, 1000)) %>%
  map(~ unname(.x)) %>%
  flatten_dbl()

es1000Table = tibble(cell,D_es1000)


grid <- dggridR::dgcellstogrid(dggs,D$cell,frame=TRUE,wrapcells=TRUE)

grid = merge(grid,cellSpCounts,id="cell",all.x=TRUE) %>%
  merge(es50Table,id="cell",all.x=TRUE) %>%
  merge(es20Table,id="cell",all.x=TRUE) %>%
  merge(es100Table,id="cell",all.x=TRUE) %>%
  merge(es200Table,id="cell",all.x=TRUE) %>%
  merge(es1000Table,id="cell",all.x=TRUE) %>%
  merge(cellOccCounts,id="cell",all.x=TRUE) %>%
  mutate(spCountLog10 = log10(spCount)) %>%
  mutate(occCountsLog10 = log10(occCounts)) %>%
  arrange(cell, order)

saveRDS(grid, sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", spacing))
grid<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", spacing))
grid_no_na = grid %>% filter(!is.na(D_es50))

grid_no_na$var<-grid_no_na$D_es50

p = gbifrasters::plotPolyMap(grid_no_na,
                             variable="D_es50",
                             legend.position = "bottom",
                             breaks=NULL,
                             labels=NULL,
                             polygon_text_size=3,
                             polygon_alpha = 1,
                             labelType="",
                             keywidth=0.01,
                             keyheight=0.3,
                             legend_text_size=15,
                             path_alpha=0,
                             mappath="/Users/huijieqiao/media/huijieqiao/WD12T/GISLayers/continents/continent.shp",
                             var_lab="ES50"
)
if (F){
  detach("package:gbifrasters", unload=TRUE)
  library(gbifrasters)
}
p
ggsave(p, filename="../../Data_eBird_2020/es50/es50_3.png", width=12, height=9)
grid_100<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", 100))
grid_100_no_na<-grid_100%>%filter(!is.na(es50))
grid_100_no_na$var<-grid_100_no_na$es50
p = gbifrasters::plotPolyMap(grid_100_no_na,
                             variable="es50",
                             legend.position = "bottom",
                             breaks=NULL,
                             labels=NULL,
                             polygon_text_size=3,
                             polygon_alpha = 1,
                             labelType="",
                             keywidth=0.01,
                             keyheight=0.3,
                             legend_text_size=15,
                             path_alpha=0,
                             mappath="/Users/huijieqiao/media/huijieqiao/WD12T/GISLayers/continents/continent.shp",
                             var_lab="ES50"
)
p

library(sf)
library(sfheaders)

sf <- sfheaders::sf_polygon(
  obj = grid_100_no_na
  , x = "long"
  , y = "lat"
  , polygon_id = "cell"
)
sf::st_crs( sf ) <- 4326

sp_df<-as(sf, 'Spatial')

data_df<-unique(grid_100_no_na[, c("cell", "lonCenter", "latCenter",     
                                   "spCount", "es50", "es20", "es100", "es200",         
                                   "es1000", "occCounts", "spCountLog10", "occCountsLog10")])
sp_df@data<-data_df
plot(sp_df, col=sp_df@data$es50)
writeOGR(sp_df, dsn="../../Data_eBird_2020/es50/shape", layer=sprintf("res_%d", 100),
         driver="ESRI Shapefile")


grid_3<-readRDS(sprintf("../../Data_eBird_2020/es50/es50_grid_%d.rda", 3))
grid_3_no_na<-grid%>%filter(!is.na(D_es50))
library(sf)
library(sfheaders)

sf <- sfheaders::sf_polygon(
  obj = grid_3_no_na
  , x = "long"
  , y = "lat"
  , polygon_id = "cell"
)
sf::st_crs( sf ) <- 4326

sp_df<-as(sf, 'Spatial')

data_df<-unique(grid_3_no_na[, c("cell",    
                                   "spCount", "D_es50", "D_es20", "D_es100", "D_es200",         
                                   "D_es1000", "occCounts", "spCountLog10", 
                                 "occCountsLog10")])
sp_df@data<-data_df
#plot(sp_df, col=sp_df@data$es50)
writeOGR(sp_df, dsn="../../Data_eBird_2020/es50/shape", layer=sprintf("res_%d", 3),
         driver="ESRI Shapefile")


