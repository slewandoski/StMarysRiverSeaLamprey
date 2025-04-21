library(sf)
library(sfnetworks)
library(ggplot2)

##spatial boundaries
regions=st_read('region_polys/Regions_2013.shp')
regions=st_transform(regions,crs=26916)

st_marys_shoreline=st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")

st_marys_shoreline<-st_transform(st_marys_shoreline,crs = 26916)

st_marys_shoreline_region=st_marys_shoreline[regions,]
st_marys_shoreline_lineString<-st_cast(st_marys_shoreline_region,"LINESTRING")
st_marys_nonsense_network = as_sfnetwork(st_marys_shoreline_lineString)
plot(st_marys_nonsense_network, main="nonsense St Marys River Network")

##overlay a spatial grid
spatial_corr_grid=st_as_sf(st_make_grid(st_marys_shoreline_region,cellsize = c(2000),what = "polygons", square = FALSE)) #dimension of grid in meters
spatial_corr_grid_in_river=spatial_corr_grid[st_marys_shoreline_region,]
spatial_corr_grid_in_river2shoreline=st_intersection(spatial_corr_grid_in_river,st_marys_shoreline_region)
  
adjecency_list=sf::st_intersects(spatial_corr_grid_in_river2shoreline,remove_self=T)


adjecency_matrix=matrix(0,length(adjecency_list),length(adjecency_list))

#make a list of lines between the centroids of adjacent grids
line_list=list()
for(j in 1:length(adjecency_list)){
  one_grid_line_list=list()
  for(q in 1:length(adjecency_list[[j]])){
    one_connection=st_as_sf(st_cast(st_combine(st_centroid(spatial_corr_grid_in_river2shoreline[c(j,adjecency_list[[j]][q]),])),"LINESTRING"))
    one_connection$line_id=paste(j,adjecency_list[[j]][q],sep="-") ##from grid id - to grid id
    one_grid_line_list[[q]]=st_as_sf(one_connection)
  }
  line_list[[j]]=do.call(rbind,one_grid_line_list)
  print(j)
}

stMarys_network_linestring=do.call(rbind,line_list)
ggplot() + geom_sf(data=st_marys_shoreline) + geom_sf(data=st_centroid(spatial_corr_grid_in_river2shoreline)) + 
  geom_sf(data=stMarys_network_linestring)

st_marys_network = as_sfnetwork(stMarys_network_linestring)
plot(st_marys_network, main="St Marys River Network")

