plotWireFrame = function(grid,textSize=1) {

  countries = map_data("world")

  p = ggplot() +
    coord_cartesian(c(-150,170),c(-55,80)) +
    theme_void() +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group,alpha=0.8), fill="gray", color="gray",alpha=0.8) +
    geom_path(data=grid,aes(x=long, y=lat, group=group), alpha=0.8, color="black") +
    geom_text(data=grid,aes(lonCenter,latCenter,label=cell),size=textSize)


  return(p)

}
