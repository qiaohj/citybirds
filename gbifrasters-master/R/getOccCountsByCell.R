getOccCountsByCell = function(D) { # a data.frame of taxonkey,count,cell

  cellOccCounts = D %>%
    group_by(cell) %>%
    summarise(occCounts = sum(count))

  return(cellOccCounts)
}
