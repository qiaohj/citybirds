getOccCountLatTable = function(D,grid,round=2) {

  cellOccCounts = gbifrasters::getOccCountsByCell(D)

  grid = merge(grid,cellOccCounts,id="cell",all.x=TRUE) %>%
    arrange(cell, order)


  occCountLatTable = grid %>%
    select(latCenter,occCounts) %>%
    mutate(latCenter = plyr::round_any(latCenter,round)) %>%
    group_by(latCenter) %>%
    summarise(occCounts = mean(occCounts))

  return(occCountLatTable)
}
