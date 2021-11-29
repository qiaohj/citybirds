
getSpCountsByCell = function(D) {

  cellSpCounts = D %>%
    group_by(cell) %>%
    summarise(spCount=n()) %>%
    mutate(cell = as.character(cell))

  return(cellSpCounts)
}

