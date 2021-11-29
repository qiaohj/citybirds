getEs50Table = function(D,esNum=50) { # expects a data.frame with cell, taxonkey, and count

  es50 = D %>%
    group_split(cell) %>%
    map(~ .x %>%
          group_by(taxonkey) %>%
          summarize(occCount = sum(count))
    ) %>%
    map(~ deframe(.x)) %>%
    modify_if(~ length(.x) <= esNum, ~ NA) %>% # run only if more than 50 species
    modify_if(~ !anyNA(.x), ~ entropart::Hurlbert(.x, esNum)) %>%
    map(~ unname(.x)) %>%
    flatten_dbl()

  cell = D %>%
    group_split(cell) %>%
    map(~ .x %>% mutate(cell = as.character(cell))) %>%
    map_chr(~ unique(.x$cell))

  es50Table = tibble(cell,es50)

  return(es50Table)
}
