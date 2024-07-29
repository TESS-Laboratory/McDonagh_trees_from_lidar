# Load necessary libraries
library(raster)       # Geographic Data Analysis and Modelling
library(sf)  
library(lidR)
library(rgl)
library(terra)

# Set the working directory
setwd("F:/uni_msc/Dissertation_RStudio")
getwd()

# Read in Environment Agency point cloud
EAPointCloud2018 <- readLAS("TF7000_P_10787_20181113_20181117.laz")
if (is.null(EAPointCloud2018)) stop("Failed to read the point cloud file.")
print(EAPointCloud2018)
summary(EAPointCloud2018)

# Validating lidar data
las_check(EAPointCloud2018)

# Plotting - basic 3D rendering
plot(EAPointCloud2018)
plot(EAPointCloud2018, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Read the shapefile (AOI)
aoi <- st_read("F:/uni_msc/Dissertation_RStudio/oxburgh_aoi.shp")
if (is.null(aoi)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(EAPointCloud2018))
print(st_crs(aoi))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi) != st_crs(EAPointCloud2018)) {
  aoi <- st_transform(aoi, crs = st_crs(EAPointCloud2018))
}

# Clip the point cloud using the shapefile AOI
clipped_las_EA_2018 <- clip_roi(EAPointCloud2018, aoi)

# Check if the clipping was successful
if (is.null(clipped_las_EA_2018)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped data
las_check(clipped_las_EA_2018)
summary(clipped_las_EA_2018)

# Remove duplicate points
clipped_las_EA_2018 <- filter_duplicates(clipped_las_EA_2018)

# Validate the clipped data after removing duplicates
las_check(clipped_las_EA_2018)

# Plot the clipped point cloud
plot(clipped_las_EA_2018)
plot(clipped_las_EA_2018, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_las_EA_2018, "F:/uni_msc/Dissertation_RStudio/clipped_point_cloud_EA_2018.laz")

### Digital Terrain Model DTM ###
# DTM clipped 2018
#dtmclipped2018 <- grid_terrain(clipped_las_EA_2018, res = 1, algorithm = tin())
#plot(dtmclipped2018, main = "Digital Terrain Model (DTM)")
# Save the DTM clipped 2018 plot to a file
#png(filename = "F:/uni_msc/Dissertation_RStudio/dtmclipped2018_plot.png", width = 800, height = 600)
#plot(dtmclipped2018, main = "Digital Terrain Model (DTM)")
#dev.off()

## (should I use )
#dtm <- grid_terrain(EAPointCloud2018, res = 1, algorithm = tin())
#plot(dtm, main = "Digital Terrain Model (DTM)")

# Render shaded DTM using terra and tin method
#dtmclipped2018_terra <- rasterize_terrain(clipped_las_EA_2018, algorithm = tin(), pkg ="terra")
#dtm_prod_2018 <- terrain(dtmclipped2018_terra, v = c("slope", "aspect"), unit = "radians")
#dtm_hillshade_2018 <- shade(slope = dtm_prod_2018$slope, aspect = dtm_prod_2018$aspect)
#plot(dtm_hillshade_2018, col =gray(0:30/30), legend = FALSE, main = "Digital Terrain Model (DTM)")
# Save the DTM plot to a file
#png(filename = "F:/uni_msc/Dissertation_RStudio/dtmclipped2018_terra_plot.png", width = 800, height = 600)
#plot(dtm_hillshade_2018, col =gray(0:30/30), legend = FALSE, main = "Digital Terrain Model (DTM)")
#dev.off()

### Height normalisation ###
## Normalize the point cloud without DTM (5.2 in lidar lidRbook)
normclipped_las_EA_2018 <- normalize_height(clipped_las_EA_2018, knnidw())
# check ground points are exactly 0
hist(filter_ground(normclipped_las_EA_2018)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
las_check(normclipped_las_EA_2018)

# Handle negative outliers
negative_points_2018 <- filter_poi(normclipped_las_EA_2018, Z < 0)
plot(negative_points_2018, color = "Z")

# Filter out points with Z < 0
cleaned_las_2018 <- filter_poi(normclipped_las_EA_2018, Z >= 0)
plot(cleaned_las_2018)

# Validate the cleaned data
las_check(cleaned_las_2018)
summary(cleaned_las_2018)

### Individual tree segmentation (ITS) (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 1, dt2 = 2)
# segment trees
trees_2018 <- segment_trees(cleaned_las_2018, algo)
# Plot the segmented trees
png(filename = "F:/uni_msc/Dissertation_RStudio/trees2018.png", width = 800, height = 600)
plot(trees_2018, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(trees_2018, "F:/uni_msc/Dissertation_RStudio/segmented_trees_EA_2018.laz")

# Filters tree number 110
#tree110 <- filter_poi(trees_2018, treeID == 110)
#plot(tree110, size = 8, bg = "white")

# Generate crown metrics
crowns <- crown_metrics(trees_2018, func = .stdtreemetrics, geom = "convex")
# Plot the crown areas
png(filename = "F:/uni_msc/Dissertation_RStudio/crownarea_2018.png", width = 800, height = 600)
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the crowns as a shapefile
st_write(crowns, "F:/uni_msc/Dissertation_RStudio/crown_metrics.shp")

# Compute the maximum height for each tree
maxtree_heights_2018 <- tree_metrics(trees_2018, ~max(Z))
# Save tree heights to a CSV file
write.csv(maxtree_heights_2018, "F:/uni_msc/Dissertation_RStudio/maxtree_heights2018.csv")
# Print tree heights
print(maxtree_heights_2018)

##### Doing ITS using canopy height model (rasterising) ####
### Canopy Height Model (CHM) ###
# Generate a Canopy Height Model (CHM)
chm_2018 <- grid_canopy(cleaned_las_2018, res = 0.5, p2r())
# Plot the CHM
plot(chm_2018, main = "Canopy Height Model (CHM)")
# Plot the segmented trees
png(filename = "F:/uni_msc/Dissertation_RStudio/chm_2018.png", width = 800, height = 600)
plot(chm_2018, main = "Canopy Height Model (CHM) EA 2018")
dev.off()
### Individual tree segmentation (ITS) ###
# Set the segmentation algorithm (dalponte2016 in this case)
ttops <- locate_trees(cleaned_las_2018, lmf(5, 2)) # Detect treetops
algo2 <- dalponte2016(chm = chm_2018, treetops = ttops)
# Segment trees
trees_2018_2 <- segment_trees(cleaned_las_2018, algo2)
# Plot the segmented trees
plot(trees_2018_2, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
# Compute the maximum height for each tree
tree_heights_2 <- tree_metrics(trees_2018_2, ~max(Z))
# Print tree heights
print(tree_heights_2)


