## Oxborough drone LIDAR data

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
DronePointCloud <- readLAS("Oxburgh_ground_points.las")
if (is.null(DronePointCloud)) stop("Failed to read the point cloud file.")
print(DronePointCloud)
summary(DronePointCloud)
#oxbpointcloud2 <- readLAScatalog("Oxburgh_ground_points.las")

# Validating lidar data
las_check(DronePointCloud)

# Plotting - basic 3D rendering
plot(DronePointCloud)
plot(DronePointCloud, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Read the shapefile (AOI)
aoi <- st_read("F:/uni_msc/Dissertation_RStudio/oxburgh_aoi.shp")
if (is.null(aoi)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(DronePointCloud))
print(st_crs(aoi))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi) != st_crs(DronePointCloud)) {
  aoi <- st_transform(aoi, crs = st_crs(DronePointCloud))
}

# Clip the point cloud using the shapefile AOI
clipped_DronePointCloud <- clip_roi(DronePointCloud, aoi)

# Check if the clipping was successful
if (is.null(clipped_DronePointCloud)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped data
las_check(clipped_DronePointCloud)
summary(clipped_DronePointCloud)

# Remove duplicate points
clipped_DronePointCloud <- filter_duplicates(clipped_DronePointCloud)

# Plot the clipped point cloud
plot(clipped_DronePointCloud)
plot(clipped_DronePointCloud, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_DronePointCloud, "F:/uni_msc/Dissertation_RStudio/clipped_DronePointCloud.las")

#GROUND CLASSIFICATION
clipped_DronePointCloud <- classify_ground(clipped_DronePointCloud, algorithm = pmf(ws = 5, th = 3))

### Height normalisation ###
## Normalize the point cloud without DTM (5.2 in lidar lidRbook)
normclipped_DronePointCloud <- normalize_height(clipped_DronePointCloud, knnidw())
# check ground points are exactly 0
hist(filter_ground(normclipped_DronePointCloud)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
las_check(normclipped_DronePointCloud)

# Handle negative outliers
negative_points_drone <- filter_poi(normclipped_DronePointCloud, Z < 0)
plot(negative_points_drone, color = "Z")

# Filter out points with Z < 0
cleaned_dronepointcloud <- filter_poi(normclipped_DronePointCloud, Z >= 0)
plot(cleaned_dronepointcloud)

# Validate the cleaned data
las_check(cleaned_dronepointcloud)

### Individual tree segmentation (ITS) (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 1, dt2 = 2)
# segment trees
trees_drone <- segment_trees(cleaned_dronepointcloud, algo)
# Plot the segmented trees
png(filename = "F:/uni_msc/Dissertation_RStudio/treesdrone.png", width = 800, height = 600)
plot(trees_drone, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(trees_drone, "F:/uni_msc/Dissertation_RStudio/segmented_rees_drone.laz")
