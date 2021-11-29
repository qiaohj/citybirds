getAndSavePolyRaster = function(grid,spacing,taxonKey,facet="speciesKey") {

  gridData = grid %>% select(cell,geometry) %>% unique()

  polyraster = gbifrasters::getData(taxonKey,gridData,facet)

  polyraster %>% saveRDS(file="C:/Users/ftw712/Desktop/es50/data/polygon_" %+% spacing %+% "_" %+% taxonKey %+%".rda")

  return(polyraster)
}
