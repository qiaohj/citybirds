
plotMiniMap = function(grid,
                       orientation = c(-0.238566,116.2990554,4),
                       xlim_ortho = c(100,150),
                       ylim_ortho = c(-18,14)) {

  countries = map_data("world")

  gridMiniMap =  grid %>%
    filter(lonCenter >= xlim_ortho[1] & lonCenter <= xlim_ortho[2]) %>%
    filter(latCenter >= ylim_ortho[1] & latCenter <= ylim_ortho[2])

  miniMap = ggplot() +
    coord_cartesian() +
    theme_void() +
    geom_polygon(data=countries,aes(x=long, y=lat, group=group), fill="#D8DACF", color=NA,alpha=1) +
    scale_y_continuous(name="",breaks = seq(-90,90,by = 20),labels=rep("",length(seq(-90,90,by = 20)))) +
    scale_x_continuous(name="",breaks = seq(-180,180,by = 20),labels=rep("",length(seq(-180,180,by = 20)))) +
    coord_map("ortho",orientation=orientation,xlim=NULL, ylim=NULL,clip="off") +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.ticks.y=element_blank()) +
    theme(plot.margin=margin(0, 0, 0, 0, "cm")) +
    theme(panel.grid.major = element_line(colour = "#D8DACF")) +
    geom_path(data=gridMiniMap,aes(x=long, y=lat, group=group),size=0.2, alpha=0.8, color="#939878")

  return(miniMap)
}



plotZoomMap = function(grid,
                       orientation = c(-0.238566,116.2990554,4),
                       xlim_ortho = c(100,150),
                       ylim_ortho = c(-18,14),
                       variable="occCounts",
                       legend_title="occ counts",
                       breaks = c(0,1e03,1e04,1e05,1e06,1e07,1e8),
                       labels = c("0-1K","1K-10K","10K-100K","100K-1M","1M-10M","10M-100M"),
                       polygon_text_size = 3,
                       labelType="") {


  # avoid this "promise already under evaluation: recursive default argument reference or earlier problems?"
  variable_ = variable
  legend_title_ = legend_title
  breaks_ = breaks
  labels_ = labels
  labelType_ = labelType
  polygon_text_size_ = polygon_text_size


  p = gbifrasters::plotPolyMap(grid,
                               variable=variable_,
                               legend_title=legend_title_,
                               legend.position = c(.50,-0.1),
                               breaks=breaks_,
                               labels=labels_,
                               polygon_text_size=polygon_text_size_,
                               labelType=labelType_
  ) +
    theme(plot.margin=margin(0, 0.7, 0, 0, "cm")) +
    coord_map("ortho",orientation=orientation,xlim=xlim_ortho, ylim=ylim_ortho) +
    theme(axis.ticks.x=element_blank())+
    theme(axis.ticks.y=element_blank())

  return(p)
}
