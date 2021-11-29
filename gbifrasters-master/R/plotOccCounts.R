plotOccCounts = function(D,
                         grid,
                         breaks=c(0,2e03,4e03,4e03,6e03,8e03,1e04,1.2e04,1.4e04,1.6e04)
                         ,labels=NULL,
                         legend_title="",
                         zoom_x=c(-150,170),
                         zoom_y=c(-55,80),
                         legend.position = c(.50,-0.05),
                         margin = margin(0.7, 0.7, 0, 0, "cm"),
                         polygon_text_size = 1
                         ) { # a data.frame of taxonkey,count,cell

  library(ggplot2)
  library(ggthemes)

  cellOccCounts = gbifrasters::getOccCountsByCell(D)

  grid = merge(grid,cellOccCounts,id="cell",all.x=TRUE) %>%
    arrange(cell, order) %>%
    na.omit() %>%
    mutate(plotLabel = case_when(
      occCounts >= 1e6 ~ round(occCounts/1e6) %+% "M",
      occCounts >= 1e3 ~ round(occCounts/1e3) %+% "K",
      TRUE ~ as.character(occCounts)
    ))

  countries = map_data("world")

  p = ggplot() +
    coord_cartesian(zoom_x,zoom_y) +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="#D8DACF", color=NA,alpha=0.8) +
    geom_path(data=grid,aes(x=long, y=lat, group=group), alpha=1, color="#7d8085") +
    geom_polygon(data=grid, aes(x=long, y=lat, group=group, fill=cut(
      occCounts,
      breaks,
      labels
    ))) +
    scale_fill_brewer(palette = "Spectral",direction=-1) +
    guides(
      fill=guide_legend(
        title=legend_title,label.position="bottom",keywidth=0.01,keyheight=0.2,nrow=1
      )) +
    theme_bw() +
    theme(panel.background = element_rect(fill = '#F0F3F8')) +
    scale_y_continuous(name="",breaks = seq(-90,90,by = 20),labels=rep("",length(seq(-90,90,by = 20)))) +
    scale_x_continuous(name="",breaks = seq(-180,180,by = 20),labels=rep("",length(seq(-180,180,by = 20)))) +
    theme(legend.position = legend.position, legend.direction = "horizontal") +
    theme(plot.margin = margin) +
    geom_text(data=grid,aes(lonCenter,latCenter,label=plotLabel),size=polygon_text_size)

  return(p)
}
