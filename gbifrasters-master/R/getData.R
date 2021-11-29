getSpeciesOccCounts = function(url,Step=1000,maxPages=200) {

  L = gbifapi::page_api_facet(url,pluck="facets",Step=Step,maxPages=maxPages,verbose=TRUE) %>%
    map(~ .x$counts) %>%
    flatten()

  taxonkey = L %>%
    map(~ .x$name) %>%
    flatten_chr()

  count = L %>%
    map(~ .x$count) %>%
    flatten_dbl()

  spOcc = tibble(taxonkey,count)

  return(spOcc)
}

getData = function(grid,taxonKey,facet="speciesKey",datasetkey=NULL) {

  gridData = grid %>% select(cell,geometry) %>% unique()

  geometries = gridData %>% pull(geometry)

  # taxonKey = "212"
  noFossils = "basis_of_record=HUMAN_OBSERVATION&basis_of_record=MACHINE_OBSERVATION&basis_of_record=OBSERVATION&basis_of_record=MATERIAL_SAMPLE&basis_of_record=LITERATURE&basis_of_record=PRESERVED_SPECIMEN&"
  api = "http://api.gbif.org/v1/occurrence/search?limit=0&"

  queryParams = "has_coordinate=true" %+% "&" %+% "has_geospatial_issue=false&" %+% noFossils

  if(!is.null(datasetkey)) {
    datasetkey = "dataset_key=" %+% datasetkey %+% "&"
    queryParams = "has_coordinate=true" %+% "&" %+% "has_geospatial_issue=false&" %+% noFossils %+% datasetkey
  }

  facet = "facet=" %+% facet %+% "&taxonKey=" %+% taxonKey

  apiCalls = geometries %>% map(~ api %+% queryParams %+% .x %+% facet)

  # print(apiCalls)
  DL = apiCalls %>% map(~ getSpeciesOccCounts(.x,Step=100,maxPages=2000))

  cell = gridData %>% pull(cell)

  D = DL %>%
    map2(cell,~ .x %>% mutate(cell = .y)) %>% # add the cell bank into to the downloaded data
    plyr::rbind.fill() %>%
    as_tibble()

  grid = gbifrasters::mergeToGrid(grid,D)

  return(grid)
}
