plotEs50 = function(D,
                    grid,
                    breaks=c(1,27,38,42,44,46,48,50),
                    labels=NULL,
                    legend_title = "es50",
                    pretty_breaks=7,
                    zoom_x=c(-150,170),
                    zoom_y=c(-55,80),
                    legend.position = c(.50,-0.05),
                    margin=margin(0.7, 0.7, 0, 0, "cm"),
                    polygon_text_size=1,
                    polygon_alpha=1,
                    orthoProjection=FALSE,
                    orientation = c(-38.49831, -179.9223, 0),
                    coord_map_xlim=NULL,
                    coord_map_ylim=NULL) { # a data.frame of taxonkey,count,cell

  library(ggplot2)
  library(ggthemes)

  es50Table = gbifrasters::getEs50Table(D)

  grid = merge(grid,es50Table,id="cell",all.x=TRUE) %>%
    arrange(cell, order) %>%
    na.omit()

  grid %>% glimpse()

  countries = map_data("world")

  # breaks10 = c(0,10,20,30,40,50)
  # obisBreaks = c(1,27,38,42,44,46,48,50)
  #
  p = ggplot() +
    coord_cartesian(zoom_x,zoom_y) +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="#D8DACF", color=NA,alpha=0.8) +
    geom_path(data=grid,aes(x=long, y=lat, group=group), alpha=1, color="#7d8085") +
    geom_polygon(data=grid, aes(x=long, y=lat, group=group, fill=cut(
      es50,
      breaks,
      labels
      )),
      alpha=polygon_alpha
      ) +
    scale_fill_brewer(palette = "Spectral",direction=-1) +
    guides(
      fill=guide_legend(
        title=legend_title,label.position="bottom",keywidth=0.01,keyheight=0.2,nrow=1
      )) +
    theme_bw() +
    theme(panel.background = element_rect(fill = '#F0F3F8', colour = 'red')) +
    scale_y_continuous(name="",breaks = seq(-90,90,by = 20),labels=rep("",length(seq(-90,90,by = 20)))) +
    scale_x_continuous(name="",breaks = seq(-180,180,by = 20),labels=rep("",length(seq(-180,180,by = 20)))) +
    theme(legend.position = legend.position, legend.direction = "horizontal") +
    theme(plot.margin = margin) +
    geom_text(data=grid,aes(lonCenter,latCenter,label=round(es50,0)),size=polygon_text_size)


    if(orthoProjection) {
      p = p + coord_map("ortho",orientation=orientation,xlim=coord_map_xlim, ylim=coord_map_ylim)
    }

  return(p)

}


