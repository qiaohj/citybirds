getPolygonGrid = function(dggs,spacing=300,landOnly=TRUE) {

  D = gbifrasters::getGridData(1) # no really small polygons

  if(landOnly) {
  D =  D %>%
    CoordinateCleaner::cc_sea(lon="longitudeMid",lat="latitudeMid",value = "clean")
  }

  D = D %>% mutate(cell = dggridR::dgGEO_to_SEQNUM(dggs,longitudeMid,latitudeMid)$seqnum)

  grid = dggridR::dgcellstogrid(dggs,D$cell,frame=TRUE,wrapcells=TRUE)


  cellList = grid %>%
    mutate(long = format(long, scientific = FALSE)) %>%
    mutate(lat = format(lat, scientific = FALSE)) %>%
    split(.$cell)

  # print(cellList)
  cell = cellList %>% map_chr(~ .x$cell %>% unique() %>% unname()) %>% as.numeric()

  # this could probably be better
  geometry = cellList %>%
    modify_if(~ nrow(.x) == 6, ~
                "geometry=POLYGON((" %+%
                .x$long[1] %+% "%20" %+% .x$lat[1] %+% "," %+%
                .x$long[2] %+% "%20" %+% .x$lat[2] %+% "," %+%
                .x$long[3] %+% "%20" %+% .x$lat[3] %+% "," %+%
                .x$long[4] %+% "%20" %+% .x$lat[4] %+% "," %+%
                .x$long[5] %+% "%20" %+% .x$lat[5] %+% "," %+%
                .x$long[6] %+% "%20" %+% .x$lat[6] %+% "," %+%
                .x$long[1] %+% "%20" %+% .x$lat[1] %+% "))&"
    ) %>%
    modify_if(~ class(.x) == "data.frame", ~
                "geometry=POLYGON((" %+%
                .x$long[1] %+% "%20" %+% .x$lat[1] %+% "," %+%
                .x$long[2] %+% "%20" %+% .x$lat[2] %+% "," %+%
                .x$long[3] %+% "%20" %+% .x$lat[3] %+% "," %+%
                .x$long[4] %+% "%20" %+% .x$lat[4] %+% "," %+%
                .x$long[5] %+% "%20" %+% .x$lat[5] %+% "," %+%
                .x$long[6] %+% "%20" %+% .x$lat[6] %+% "," %+%
                .x$long[7] %+% "%20" %+% .x$lat[7] %+% "))&"
    ) %>%
    flatten_chr() %>%
    stringr::str_replace_all(" ","")

  d = tibble(cell,geometry)

  grid = merge(grid,d,id="cell",all.x=TRUE) # merge with grid

  centerList = dgSEQNUM_to_GEO(dggs,as.numeric(grid$cell))
  cellCenters = tibble(cell=grid$cell,lonCenter = pluck(centerList,"lon_deg"), latCenter =  pluck(centerList,"lat_deg")) %>% select(lonCenter, latCenter)

  grid = cbind(grid,cellCenters)

  return(grid)
}
