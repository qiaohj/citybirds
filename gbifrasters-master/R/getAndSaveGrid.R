getAndSaveGrid = function(spacing) {

  dggs = dggridR::dgconstruct(spacing=spacing, metric=FALSE, resround='down')

  grid = gbifrasters::getPolygonGrid(dggs,spacing)

  grid %>% saveRDS(file="C:/Users/ftw712/Desktop/es50/data/gridframe_" %+% spacing %+% ".rda")

  return(grid)
}
