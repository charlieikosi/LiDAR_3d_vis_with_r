# Set crs
proj_crs <- "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000
+ellps=GRS80 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


##---------------------------------LIBRARIES------------------------------------
library(tidyverse)
library(sf)
library(tmap)
library(buffeRs)
library(raster)
library(lidR)
library(rlas)
library(rgl)
library(rgdal)
library(gdalUtils)
library(quadmesh)
library(alphashape3d)
library(tmaptools)
library(RStoolbox)
library(gstat)
library(akima)
library(Morpho)
##--------------------------BUILDING FOOTPRINTS---------------------------------


# Create vector object by using the st_read function.
nz_footprint <- st_read("nz-building-outlines.shp")
st_crs(nz_footprint)

nz_basemap_footprint <- tm_shape(nz_footprint) +
  tm_borders(alpha = 0.3)
nz_basemap_footprint


# Select one building footprint.
building_ID <- nz_footprint[, 1]
uc_library_footprint <- building_ID[building_ID$building_i == 2042309, ]
library_plot <- tm_shape(uc_library_footprint) +
  tm_fill("red") #creates plot object. Library footprint is in red
library_plot
nz_basemap_footprint + library_plot


# Create a 1m buffer 
buffer <- st_buffer(uc_library_footprint, dist = 1)
buffer_1m <- st_buffer(uc_library_footprint, dist = 1) %>%
  tm_shape() +
  tm_borders(alpha = 0.3) + 
  tm_fill("yellow")


# view buffer, overlay the library building footprint onto it
buffer_1m + library_plot
nz_basemap_footprint + buffer_1m + library_plot


#save the library and buffer building footprints as shapefile for use in gis
write_sf(uc_library_footprint, "main_library_footprint.shp")
write_sf(buffer, "main_library_footprint_buffer.shp")


##---------------------------LIDAR POINT CLOUD----------------------------------

# load the lidar point cloud
nz_pointcloud <- readLAS("points.las")
summary(nz_pointcloud)
las_check(nz_pointcloud)
plot(nz_pointcloud, color = "Z") 


# Clip point cloud and building footprint
uc_library_sp <- as_Spatial(buffer) # we convert this sf object to an spDataframe before we clip
uc_library_pointcloud <- clip_roi(nz_pointcloud, uc_library_sp) # lasclip is deprecated so I used clip_roi instead.
plot(uc_library_pointcloud)

# Save las file.
writeLAS(uc_library_pointcloud, "uc_library.las")

##----------------------------------DEM-----------------------------------------

#Load DEM
nzdem <- raster("dem/DEM_BX24_2018_1000_1305.tif")
plot(nzdem) # plot to check dem
crs(nzdem)

library_dem <- nzdem %>%
  crop(uc_library_footprint) %>%
  mask(uc_library_footprint) # this code clips the dem to the 1m buffer building footprint
plot(library_dem)
writeRaster(library_dem, "path/to/directory")
library_dem
getValues(library_dem)

#Check extents
extent(library_dem) 
#extent(only_uc_main_library_footprint_buffer)
class(uc_library_pointcloud)
class(library_dem)


##----------------------------POINTS EXTRACTION---------------------------------

#install.packages("skimr")
#library(skimr)
# Subtract lidar elevation from dem.

# We will do this in two parts. First we create a data frame from the lidar data and save it as an sf object
# 1. To do this, we extract the data information contained in the "datatable" of the las object

df <- uc_library_pointcloud@data # this extracts only the information contained in the data table of the point cloud object
# and creates a dataframe.
# lidar data is very noisy. We did see this in the early 3D plot. Lets plot this in a ggplot
# and identify which classes they fall in so we can filter them out.
df_plot <- df %>% ggplot(mapping = aes(x = X, y = Z, color = factor(Classification))) +
  geom_point()
df_plot

# We see that Class 7 needs to be filtered out. Lets filter them
# Before we filter them, lets select the columns we need and then we plot.
df_2 <- df %>% filter(Classification <= 1) %>% 
  dplyr::select(X, Y, Z, Classification) # filtered dataframe without class 2 and 7
# To see the plot use ggplot
df_plot2 <- ggplot(df_2, aes(x = X, y = Z, color = factor(Classification))) +
  geom_point()
df_plot2 

# 2. Then we convert it to an sf object by using the X and Y values. 
dfsf <- st_as_sf(df_2, coords = c("X","Y")) #converts the above df to an sf object by using the X and Y values

# The following code works with the dem object:
# We need the cell values of the dem to do the subtraction. Therefore we must use raster::extract.
dfsf.ground <- raster::extract(nzdem, dfsf) # Extracts the data values in raster.
# Note, the use of the double colon is to specify which extract function to use from a library package, in this case,
# it is from the raster package. If you look in the environment window pane, the code has created a list of numbers which are stored as values in R
# under the name dfsf.ground.

dfsf$ground <- dfsf.ground # this code will append the list of cell values we extracted from the dem as dfsf.ground to
# as a new column with the heading "ground" in the data frame dfsf. Note the use of $ is used to give this heading name.
# ground values in the dataframe


# Now the last step, we subtract the point cloud elevation and the dem elevation to get the corrected relative height.
dfsf <- dfsf %>% mutate(relheight = Z - ground) # performs the subtraction using mutate and 
class(dfsf)

# Also create a tibble with the corrected height
x <- as.numeric(df_2$X)
y <- as.numeric(df_2$Y) 
z <- as.numeric(dfsf$relheight) 
classification <- as.numeric(dfsf$Classification)

corrected <- tibble(x, y, z, classification)

##------------------------INTERPOLATE POINTCLOUD--------------------------------

# IDW Interpolation
library(gstat)

# create empty grid
grid = as(library_dem, "SpatialPixels")
crs(grid) <- proj_crs
grddf = as.data.frame(grid)

# create spatial points df
df_sp <- df %>% filter(Classification <= 1) %>% 
  dplyr::select(X, Y, Z, Classification) # filtered dataframe without class 2 and 7

coordinates(df_sp) = ~ X + Y
crs(df_sp) <-  proj_crs

# IDW
interp = idw(formula = Z~1, 
             locations = df_sp, 
             newdata = grid)
interp_df = as.data.frame(interp) %>% dplyr::select("x", "y", "var1.pred")

interp_raster <- rasterFromXYZ(interp_df, res = 1, crs=proj_crs, digits=3)
class(interp_raster)

plot(interp_raster, col = terrain.colors(20))

ggplot()+
  geom_tile(data = interp_df, aes(x = x, y = y, fill = var1.pred))+
  scale_fill_gradientn(colors = terrain.colors(10))+
  theme_bw()

# KRIGE


##---------------------------3D Plots------------------------------------------

# Points in 3D
plot3d(x, y, z, zlab = "height")

# 3D Mesh2 from IDW Interpolation
idw_3d <- quadmesh(interp_raster)
shade3d(idw_3d, col = "red")


# 3D Mesh3 from delaunay
dxyz <- deldir::deldir(x - mean(x), y - mean(y), z = z)
persp3d(dxyz, col = "red")

mesh_3d <- as.mesh3d(dxyz)
shade3d(mesh_3d, color = "red")
mesh2obj(mesh_3d, filename = dataname) #write 3d file
##------------------








