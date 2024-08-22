library(raster)       # Geographic Data Analysis and Modelling
library(sf)  
library(lidR)
library(rgl)

# Set the working directory
setwd("D:/uni_msc/Dissertation_RStudio")
getwd()

### Environment Agency ###
## Read in Environment Agency point cloud
EAPointCloud <- readLAS("TF7000_P_10787_20181113_20181117.laz")
print(EAPointCloud)
#summary(EAPointCloud)

# Load XYZ and intensity only
EAPointCloudxyzi <- readLAS("TF7000_P_10787_20181113_20181117.laz", select = "xyzi")

# Read using LAS catalog
EAPointCloud2 <- readLAScatalog("TF7000_P_10787_20181113_20181117.laz")
print(EAPointCloud2)
#summary(EAPointCloud2)

## Validating lidar data
las_check(EAPointCloud)
las_check(EAPointCloud2)

## Plotting - basic 3D rendering
#plot(EAPointCloud)
#plot(EAPointCloud2)
#plot(EAPointCloudxyzi)
#plot(EAPointCloud, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Read the shapefile (AOI)
aoi <- st_read("D:/uni_msc/Dissertation_RStudio/weoodland2.shp")

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi) != st_crs(EAPointCloud)) {
  aoi <- st_transform(aoi, crs = st_crs(EAPointCloud))
}

# Clip the point cloud using the shapefile AOI
clipped_las <- clip_roi(EAPointCloud, aoi)

# Validate the clipped data
las_check(clipped_las)

# Plot the clipped point cloud
plot(clipped_las, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)
plot(clipped_las)
# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_las, "D:/uni_msc/Dissertation_RStudio/clipped_point_cloud.laz")

#### CODE I HAVENT USED ####
# Load XYZ and intensity only
#EAPointCloudxyzi <- readLAS("TF7000_P_10787_20181113_20181117.laz", select = "xyzi")
#if (is.null(EAPointCloudxyzi)) stop("Failed to read the point cloud file with selected attributes.")

# Read using LAS catalog
#EAPointCloud2 <- readLAScatalog("TF7000_P_10787_20181113_20181117.laz")
#print(EAPointCloud2)
#summary(EAPointCloud2)

# Handle degenerated ground points by averaging Z values
#coords <- as.data.frame(clipped_las_EA_2018@data)
#duplicate_indices <- which(duplicated(coords[, c("X", "Y")]))

#if (length(duplicate_indices) > 0) {
#message("Handling degenerated ground points...")

# Find unique (X, Y) pairs and average their Z values
#averaged_coords <- aggregate(Z ~ X + Y, data = coords, FUN = mean)

# Keep necessary columns
#common_columns <- intersect(names(coords), c("X", "Y", "Z", "Intensity", "ReturnNumber", "NumberOfReturns", "ScanDirectionFlag", "EdgeOfFlightLine", "Classification", "ScanAngleRank", "UserData", "PointSourceID", "GpsTime"))

# Merge the averaged Z values back into the original data frame
#coords <- merge(coords, averaged_coords, by = c("X", "Y"), suffixes = c("", "_avg"))
#coords$Z <- coords$Z_avg
#coords$Z_avg <- NULL

# Select only common columns for new LAS object
#coords <- coords[, common_columns]

# Create a new LAS object with the updated coordinates
#clipped_las_EA_2018 <- LAS(coords, header = EAPointCloud2018@header)
#projection(clipped_las_EA_2018) <- projection(EAPointCloud2018)
#}

## DTM
dtm <- grid_terrain(EAPointCloud2018, res = 1, algorithm = tin())
plot(dtm, main = "Digital Terrain Model (DTM)")
# tin method
dtm_tin <- rasterize_terrain(clipped_las_EA_2018, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin, bg = "white") 
# invert distance weighing for dtm
dtm_idw <- rasterize_terrain(clipped_las_EA_2018, algorithm = knnidw(k = 10L, p = 2))
plot_dtm3d(dtm_idw, bg = "white") 
# kriging for dtm
dtm_kriging <- rasterize_terrain(clipped_las_EA_2018, algorithm = kriging(k = 40))
plot_dtm3d(dtm_kriging, bg = "white")

## DTM (5.3 Hybrid method)
normclipped_las_EA_2018_2 <- normalize_height(clipped_las_EA_2018, tin(), dtm = dtm)
hist(filter_ground(normclipped_las_EA_2018_2)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")

##### Doing ITS using canopy height model (rasterising) ####
# Plot the segmented trees
#png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/chm_2018.png", width = 800, height = 600)
#plot(chm_2018, main = "Canopy Height Model (CHM) EA 2018")
#dev.off()
### Individual tree segmentation (ITS) ###
# Set the segmentation algorithm (dalponte2016 in this case)
#ttops <- locate_trees(cleaned_las_2018, lmf(5, 2)) # Detect treetops
#algo2 <- dalponte2016(chm = chm_2018, treetops = ttops)
# Segment trees
#trees_2018_2 <- segment_trees(cleaned_las_2018, algo2)
# Plot the segmented trees
#plot(trees_2018_2, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
# Compute the maximum height for each tree
#tree_heights_2 <- tree_metrics(trees_2018_2, ~max(Z))
# Print tree heights
#print(tree_heights_2)

### CODE I NEED TO MAKE

# Segment trees
#segment_trees

# Extract tree metrics
#crown_metrics 

# Use sf package for tree crown area and perimeter 

## Oxborough data
## Read in drone point cloud
#oxbpointcloud <- readLAS("Oxburgh_ground_points.las")
#print(oxbpointcloud)
#summary(oxbpointcloud)

#oxbpointcloud2 <- readLAScatalog("Oxburgh_ground_points.las")
#print(oxbpointcloud2)
#summary(oxbpointcloud2)

## Validating lidar data
#las_check(oxbpointcloud)
#las_check(oxbpointcloud2)

#plot(oxbpointcloud)
