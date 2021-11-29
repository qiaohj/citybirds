getGridData = function(stepSize) {

  latitudeBreaks = seq(-90,(90-stepSize),stepSize)
  latitudeBreaksOffset = seq(-1*(90-stepSize),90,stepSize)

  longitudeBreaks = seq(-180,(180-stepSize),stepSize)
  longitudeBreaksOffset = seq(-1*(180-stepSize),180,stepSize)

  gridGeometries = longitudeBreaks %>%
    map2(
      longitudeBreaksOffset, ~
        "geometry=POLYGON((" %+%
        .x %+% "%20" %+% latitudeBreaks %+%
        "," %+%
        .y %+% "%20" %+% latitudeBreaks %+%
        "," %+%
        .y %+% "%20" %+% latitudeBreaksOffset %+%
        "," %+%
        .x %+% "%20" %+% latitudeBreaksOffset %+%
        "," %+%
        .x %+% "%20" %+% latitudeBreaks%+%
        "))&"
    )

  longitudeMid = rep(longitudeBreaks, each = length(latitudeBreaks)) + (stepSize/2)
  latitudeMid = rep(latitudeBreaks,length(longitudeBreaks)) + (stepSize/2)

  # print(latitudeMid)

  gridData = gridGeometries %>%
    map(~
          enframe(.x) %>%
          select(geometry = value)
    ) %>%
    plyr::rbind.fill() %>%
    add_column(latitudeMid) %>%
    add_column(longitudeMid)

  return(gridData)
}
