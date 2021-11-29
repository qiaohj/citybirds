
getEs50LatTable = function(grid,round=2) {

  es50LatTable = grid %>%
    select(latCenter,es50) %>%
    mutate(latCenter = plyr::round_any(latCenter,round)) %>%
    group_by(latCenter) %>%
    summarise(es50Mean = mean(es50,na.rm=TRUE)) %>%
    na.omit()

  return(es50LatTable)
}

