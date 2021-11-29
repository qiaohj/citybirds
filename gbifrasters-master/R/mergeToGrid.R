mergeToGrid = function(grid,D) {

  es50Table = gbifrasters::getEs50Table(D,esNum=50)
  es20Table = gbifrasters::getEs50Table(D,esNum=20) %>% rename(es20 = es50)
  es100Table = gbifrasters::getEs50Table(D,esNum=100) %>% rename(es100 = es50)
  es200Table = gbifrasters::getEs50Table(D,esNum=200) %>% rename(es200 = es50)
  es1000Table = gbifrasters::getEs50Table(D,esNum=1000) %>% rename(es1000 = es50)
  cellSpCounts = gbifrasters::getSpCountsByCell(D)
  cellOccCounts = gbifrasters::getOccCountsByCell(D)

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

  return(grid)
}

