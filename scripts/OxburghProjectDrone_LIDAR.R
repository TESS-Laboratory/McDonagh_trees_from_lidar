## Oxborough drone LIDAR data

# Load necessary libraries
library(raster)       # Geographic Data Analysis and Modelling
library(sf)  
library(lidR)
library(rgl)
library(terra)

# Read in Environment Agency point cloud
#DronePointCloud <- readLAS("F:/uni_msc/Dissertation/data/Oxburgh_ground_points.las")
DronePointCloud <- readLAScatalog("D:/oxburgh_pointcloudD/Oxburgh_full_lidar_cloud.las")
if (is.null(DronePointCloud)) stop("Failed to read the point cloud file.")
print(DronePointCloud)
summary(DronePointCloud)

# Validating lidar data
las_check(DronePointCloud)

# Read the shapefile (AOI)
aoi <- st_read("F:/uni_msc/Dissertation/oxburgh_aoi.shp")
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

clipped_DronePointCloud <-readLAScatalog("C:/workspace/McDonagh_trees_from_lidar/outputs/clipped_DronePointCloud.las")
clipped_DronePointCloud <-readLAS("C:/workspace/McDonagh_trees_from_lidar/outputs/clipped_DronePointCloud.las")
# Validate the clipped data
las_check(clipped_DronePointCloud)
summary(clipped_DronePointCloud)

# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_DronePointCloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/clipped_DronePointCloud.las")

##Density Reduction as point cloud is too large
# Reduce density by randomly subsampling 50% of the points
subsampled_DronePointCloud <- filter_poi(clipped_DronePointCloud, runif(npoints(clipped_DronePointCloud)) < 0.2)
subsampled_DronePointCloud_2 <- filter_poi(clipped_DronePointCloud, runif(npoints(clipped_DronePointCloud)) < 0.3)

# Optionally, save the output to a new LAS file
writeLAS(subsampled_DronePointCloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/subsampled_DronePointCloud.las")
writeLAS(subsampled_DronePointCloud_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/subsampled_DronePointCloud_2.las")

subsampled_DronePointCloud <-readLAS("C:/workspace/McDonagh_trees_from_lidar/outputs/subsampled_DronePointCloud.las")
subsampled_DronePointCloud_2 <-readLAS("C:/workspace/McDonagh_trees_from_lidar/outputs/subsampled_DronePointCloud_2.las")
# Check the size and density of the subsampled point cloud
summary(subsampled_DronePointCloud)
plot(subsampled_DronePointCloud)

#GROUND CLASSIFICATION
classified_subsampled_DronePointCloud <- classify_ground(subsampled_DronePointCloud, algorithm = pmf(ws = 5, th = 3))
writeLAS(classified_subsampled_DronePointCloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/classified_subsampled_DronePointCloud.las")

summary(classified_subsampled_DronePointCloud)
plot(classified_subsampled_DronePointCloud)
las_check(classified_subsampled_DronePointCloud)

### Height normalisation ###
## Normalize the point cloud without DTM (5.2 in lidar lidRbook)
normDronePointCloud <- normalize_height(classified_subsampled_DronePointCloud, knnidw())
las_check(normDronePointCloud)
writeLAS(normDronePointCloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/normDronePointCloud.las")
##Ground classification
# check ground points are exactly 0
ground_points <- filter_ground(normDronePointCloud)
range(ground_points$Z)
# Adjust the breaks to cover the full range of Z values
breaks <- seq(0, 2.2, by = 0.01)
# Plot the histogram with the adjusted breaks
hist(ground_points$Z, breaks = breaks, main = "", xlab = "Elevation")
las_check(normclipped_DronePointCloud)
plot(normclipped_DronePointCloud)

# Optionally, save the clipped point cloud to a new file
writeLAS(normclipped_DronePointCloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/normclipped_DronePointCloud.las")

# Handle negative outliers
negative_points_drone <- filter_poi(normDronePointCloud, Z < 0)
plot(negative_points_drone, color = "Z")

# Filter out points with Z < 0
cleaned_dronepointcloud <- filter_poi(normDronePointCloud, Z >= 0)
plot(cleaned_dronepointcloud)

# Optionally, save the clipped point cloud to a new file
writeLAS(cleaned_dronepointcloud, "C:/workspace/McDonagh_trees_from_lidar/outputs/cleaned_dronepointcloud.las")

# Extract the ground points
ground_points <- filter_ground(normDronePointCloud)

# View the degenerated points
print(degenerated_points)
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
writeLAS(trees_drone, "F:/uni_msc/Dissertation_RStudio/segmented_trees_drone.laz")
