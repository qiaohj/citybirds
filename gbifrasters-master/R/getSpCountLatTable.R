
getSpCountLatTable = function(grid,round=2) {

  spCountLatTable = grid %>%
    select(latCenter,spCount) %>%
    mutate(latCenter = plyr::round_any(latCenter,round)) %>%
    group_by(latCenter) %>%
    summarise(spCountMean = mean(spCount,na.rm=TRUE),spCountMax = max(spCount,na.rm=TRUE),spCountMin = min(spCount,na.rm=TRUE))

  return(spCountLatTable)
}


