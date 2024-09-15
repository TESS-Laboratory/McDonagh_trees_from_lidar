## Project area of interest: Home Covert woodland at Oxburgh Estate, Oxborough, Norfolk
## Using National Trust (NT) 2023 Drone LiDAR point cloud data for tree metrics extraction
## and estimation of Aboveground Biomass using allometric equations
## Author of code: Amber McDonagh am1287@exeter.ac.uk  

# Load necessary libraries
library(sf)  # Geographic Data Analysis and Modelling
library(lidR) # Manipulating and Visualising LiDAR data
library(rgl) # Interactive 3D plots and Modelling
library(terra) # Spatial Data Analysis
library(tidyverse) # Data manipulation and Visualisation
library(ggplot2) # Statistical Graphics for Data Analysis

## Plotting 'fancy' theme created by Andrew Cunliffe
theme_fancy <- function() {
  theme_bw() +
    theme(
      text = element_text(family = "Helvetica"),
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 6, color = "black"),
      axis.line.x = element_line(size = 0.3, color = "black"),
      axis.line.y = element_line(size = 0.3, color = "black"),
      axis.ticks = element_line(size = 0.3, color = "black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
      plot.title = element_text(
        size = 8,
        vjust = 1,
        hjust = 0.5,
        color = "black"
      ),
      legend.text = element_text(size = 6, color = "black"),
      legend.title = element_text(size = 6, color = "black"),
      legend.position = c(0.3, 1),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 2,
        linetype = "blank"
      )
    )
}
windowsFonts("Helvetica" = windowsFont("Helvetica")) # Ensure font is mapped correctly

## Read and Prepare NT Drone Data 
# Read in drone point cloud
DronePointCloud <- readLAScatalog("D:/oxburgh_pointcloudD/Oxburgh_full_lidar_cloud.las")
if (is.null(DronePointCloud)) stop("Failed to read the point cloud file.")
print(DronePointCloud)
summary(DronePointCloud)
las_check(DronePointCloud)

# Read the shapefile (Area of Interest - AOI, Home Covert woodland)
aoi <- st_read("F:/uni_msc/Dissertation/TopographicArea.shp")
if (is.null(aoi)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(DronePointCloud))
print(st_crs(aoi))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi) != st_crs(DronePointCloud)) {
  aoi <- st_transform(aoi, crs = st_crs(DronePointCloud))
}

## Clip and Validate point cloud
# Clip the point cloud using the shapefile AOI
clipped_DronePointCloud <- clip_roi(DronePointCloud, aoi)

# Check if the clipping was successful
if (is.null(clipped_DronePointCloud)) stop("Clipping the point cloud to the AOI failed.")

# Code to read clipped point cloud (use when want to avoid re-processing data)
clipped_DronePointCloud <-readLAScatalog("D:/datafilesfordissertation/interim_pointclouds/clipped_DronePointCloud.las")
clipped_DronePointCloud <-readLAS("D:/datafilesfordissertation/interim_pointclouds/clipped_DronePointCloud.las")

# Validate the clipped data
las_check(clipped_DronePointCloud)
summary(clipped_DronePointCloud)

# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_DronePointCloud, "D:/datafilesfordissertation/interim_pointclouds/clipped_DronePointCloud.las")

##Density Reduction as point cloud is too large
# Reduce density by randomly subsampling 80% of the points
subsampled_DronePointCloud <- filter_poi(clipped_DronePointCloud, runif(npoints(clipped_DronePointCloud)) < 0.2)
# Reduce density by randomly subsampling 70% of the points
subsampled_DronePointCloud_2 <- filter_poi(clipped_DronePointCloud, runif(npoints(clipped_DronePointCloud)) < 0.3)
# Reduce density by randomly subsampling 90% of the points
subsampled_DronePointCloud_3 <- filter_poi(clipped_DronePointCloud, runif(npoints(clipped_DronePointCloud)) < 0.1)

# Optionally, save the output to a new LAS file
writeLAS(subsampled_DronePointCloud, "D:/datafilesfordissertation/interim_pointclouds/subsampled_DronePointCloud.las")
writeLAS(subsampled_DronePointCloud_2, "D:/datafilesfordissertation/interim_pointclouds/outputs/subsampled_DronePointCloud_2.las")
writeLAS(subsampled_DronePointCloud_3, "C:/workspace/McDonagh_trees_from_lidar/outputs/subsampled_DronePointCloud_3.las")

# Code to read subsampled point cloud (use when want to avoid re-processing data)
subsampled_DronePointCloud <-readLAS("D:/datafilesfordissertation/interim_pointclouds/subsampled_DronePointCloud.las")
subsampled_DronePointCloud_2 <-readLAS("D:/datafilesfordissertation/interim_pointclouds/subsampled_DronePointCloud_2.las")

# Check the size and density of the subsampled point cloud
summary(subsampled_DronePointCloud)
plot(subsampled_DronePointCloud)

summary(subsampled_DronePointCloud_3)
plot(subsampled_DronePointCloud_3)

## Process Data
# Ground classification
classified_subsampled_DronePointCloud <- classify_ground(subsampled_DronePointCloud, algorithm = pmf(ws = 5, th = 3))
writeLAS(classified_subsampled_DronePointCloud, "D:/datafilesfordissertation/interim_pointclouds/classified_subsampled_DronePointCloud.las")

classified_subsampled_DronePointCloud_3<- classify_ground(subsampled_DronePointCloud_3, algorithm = pmf(ws = 5, th = 3))
writeLAS(classified_subsampled_DronePointCloud_3, "D:/datafilesfordissertation/interim_pointclouds/classified_subsampled_DronePointCloud_3.las")

# Validate classified subsampled drone point cloud
summary(classified_subsampled_DronePointCloud)
plot(classified_subsampled_DronePointCloud)
las_check(classified_subsampled_DronePointCloud)

summary(classified_subsampled_DronePointCloud_3)
plot(classified_subsampled_DronePointCloud_3)
las_check(classified_subsampled_DronePointCloud_3)

## Normalise the point cloud without DTM (5.2 in lidar lidRbook)
normDronePointCloud <- normalize_height(classified_subsampled_DronePointCloud, knnidw())
las_check(normDronePointCloud)
writeLAS(normDronePointCloud, "D:/datafilesfordissertation/interim_pointclouds/normDronePointCloud.las")

normDronePointCloud_3 <- normalize_height(classified_subsampled_DronePointCloud_3, knnidw())
las_check(normDronePointCloud_3)
writeLAS(normDronePointCloud_3, "D:/datafilesfordissertation/interim_pointclouds/normDronePointCloud_3.las")

# check ground points are exactly 0
ground_points <- filter_ground(normDronePointCloud)
range(ground_points$Z)

ground_points_3 <- filter_ground(normDronePointCloud_3)
range(ground_points_3$Z)

# Adjust the breaks to cover the full range of Z values
breaks <- seq(0, 2.2, by = 0.01)

# Plot the histogram with the adjusted breaks
hist(ground_points_3$Z, breaks = breaks, main = "", xlab = "Elevation")
las_check(normDronePointCloud_3)
plot(normDronePointCloud_3)

# Optionally, save the clipped point cloud to a new file
writeLAS(normDronePointCloud_3, "D:/datafilesfordissertation/interim_pointclouds/normDronePointCloud_3.las")

# Handle negative outliers
negative_points_drone_3 <- filter_poi(normDronePointCloud_3, Z < 0)
plot(negative_points_drone_3, color = "Z")

# Filter out points with Z < 0
cleaned_dronepointcloud_3 <- filter_poi(normDronePointCloud_3, Z >= 0)
plot(cleaned_dronepointcloud_3)
las_check(cleaned_dronepointcloud_3)

# Optionally, save the clipped point cloud to a new file
writeLAS(cleaned_dronepointcloud_3, "D:/datafilesfordissertation/interim_pointclouds/cleaned_dronepointcloud_3.las")
cleaned_dronepointcloud_3 <-readLAS("D:/datafilesfordissertation/interim_pointclouds/cleaned_dronepointcloud_3.las")

# Extract the ground points
ground_points <- filter_ground(cleaned_dronepointcloud)
range(ground_points$Z)

# View the degenerated points
print(degenerated_points)

# Validate the cleaned data
las_check(cleaned_dronepointcloud)

### Individual tree segmentation (ITS) (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 1, dt2 = 2)
# segment trees
trees_drone <- segment_trees(cleaned_dronepointcloud_3, algo)
# Plot the segmented trees
png(filename = "C:/workspace/McDonagh_trees_from_lidar/treesdrone.png", width = 800, height = 600)
plot(trees_drone, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(trees_drone, "D:/datafilesfordissertation/interim_pointclouds/segmented_trees_drone.las")
trees_drone <-readLAS("D:/datafilesfordissertation/interim_pointclouds/segmented_trees_drone.las")

plot(trees_drone)

##Counting drone segmented trees
# Extract the treeID attribute
drone_tree_ids <- trees_drone@data$treeID

# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
# Depending on the segmentation method, treeID might be negative for some non-tree points.
unique_tree_count <- length(unique(drone_tree_ids[drone_tree_ids > 0]))

# Output the number of detected trees
cat("Number of trees segmented and detected:", unique_tree_count, "\n")

# Filters tree number 110
#tree110 <- filter_poi(trees_2018, treeID == 110)
#plot(tree110, size = 8, bg = "white")

## Generate Crown Metrics
drone_crowns <- crown_metrics(trees_drone, func = .stdtreemetrics, geom = "convex")
# Plot the crown areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_drone.png", width = 600, height = 600)
plot(drone_crowns["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the crowns as a shapefile
st_write(drone_crowns, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics.shp", append=FALSE)

### Calculate additional tree metrics
# Compute the maximum height for each tree
maxtree_heights_drone <- tree_metrics(trees_drone, ~max(Z))
# Save tree heights to a CSV file
write.csv(maxtree_heights_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heightsdrone.csv")
# Print tree heights
print(maxtree_heights_drone)

# Compute maximum and average tree height
max_height_drone <- max(trees_drone$V1, na.rm = TRUE)
avg_height_drone <- mean(trees_dronw$V1, na.rm = TRUE)

print(max_height_drone)
print(avg_height_drone)

# Convert crowns to sf object
# crowns will be an sf object with geometries and associated metrics
drone_crowns_sf <- drone_crowns

## Calculate tree metrics using sf
# 1. Crown Area
drone_crowns_sf$area <- st_area(drone_crowns_sf)  # Calculate area for each crown polygon
# 2. Crown Perimeter
drone_crowns_sf$perimeter <- st_length(st_cast(drone_crowns_sf, "MULTILINESTRING"))
# 3. Tree Height
max_area_drone <- max(drone_crowns_sf$area, na.rm = TRUE)
avg_area_drone <- mean(drone_crowns_sf$area, na.rm = TRUE)
max_perimeter_drone <- max(drone_crowns_sf$perimeter, na.rm = TRUE)
avg_perimeter_drone <- mean(drone_crowns_sf$perimeter, na.rm = TRUE)
# Print summary tree metrics 
print(paste("Max crown area:", max_area_drone))
print(paste("Average crown area:", avg_area_drone))
print(paste("Max crown perimeter:", max_perimeter_drone))
print(paste("Average crown perimeter:", avg_perimeter_drone))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(drone_crowns_sf, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics_with_sf.shp", append=FALSE)

# Convert sf object to a data frame for CSV output
drone_crowns_df <- as.data.frame(drone_crowns_sf)
# Ensure the geometry column is not included in the CSV
drone_crowns_df$geometry <- NULL
# Add summary statistics to the data frame
drone_crowns_df$max_area_drone <- max_area_drone
drone_crowns_df$avg_area_drone <- avg_area_drone
drone_crowns_df$max_perimeter_drone <- max_perimeter_drone
drone_crowns_df$avg_perimeter_drone <- avg_perimeter_drone

## Calculating Crown diameter using max length of crown (convex hull)
drone_crowns_sf$crown_diameter_hull_drone <- st_length(st_cast(st_convex_hull(drone_crowns_sf), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_drone <- max(drone_crowns_sf$crown_diameter_hull_drone, na.rm = TRUE)
avg_diameter_hull_drone <- mean(drone_crowns_sf$crown_diameter_hull_drone, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_drone, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_drone, "\n")

# Add crown diameters to the data frame
drone_crowns_df$crown_diameter_drone <- drone_crowns_sf$crown_diameter_drone
drone_crowns_df$crown_diameter_hull_drone <- drone_crowns_sf$crown_diameter_hull_drone

# Save the metrics to a CSV file
write.csv(crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_tree_metrics.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics.csv'.\n")


### ITS using canopy height model (rasterising) ###
### Canopy Height Model (CHM) ###
# Generate a Canopy Height Model (CHM)
chm_drone <- grid_canopy(cleaned_dronepointcloud_3, res = 0.5, p2r())
# Plot the CHM
plot(chm_drone, main = "Canopy Height Model (CHM)")

### Aboveground Biomass (AGB) using mixed temperate woodland calculation - generalised calculation ###
## Implementing AGB calculations using height, crown area and perimeter
TreeHeightDrone <- maxtree_heights_drone$V1 # Use max tree height from data
print(TreeHeightDrone)

# Use `crowns_df' data frame that includes `avg_area` and `avg_perimeter`
CrownAreaDrone <- drone_crowns_sf$area # Use average crown area
CrownPerimeterDrone <- drone_crowns_sf$perimeter  # Use average crown perimeter

# Coefficients for mixed broadleaf and conifer forest
beta_0_mixed <- 2.5
beta_1_mixed <- 0.9
beta_2_mixed <- 0.5
beta_3_mixed <- 0.3

# Convert variables to numeric
TreeHeightDrone_numeric <- as.numeric(TreeHeightDrone)
CrownAreaDrone_numeric <- as.numeric(CrownAreaDrone)
CrownPerimeterDrone_numeric <- as.numeric(CrownPerimeterDrone)

# AGB calculation using the allometric model
AGB_mixed_drone <- beta_0_mixed + 
  beta_1_mixed * log(TreeHeightDrone_numeric) + 
  beta_2_mixed * log(CrownAreaDrone_numeric) + 
  beta_3_mixed * log(CrownPerimeterDrone_numeric)

# Exponentiate to get the actual AGB in the original units
AGB_mixed_drone <- exp(AGB_mixed_drone)

# Print the calculated AGB
print(AGB_mixed_drone)

# Assuming 'AGB' is a vector of AGB values for each tree
total_AGB_mixed_drone <- sum(AGB_mixed_drone)

# Print the total AGB for the forest
print(total_AGB_mixed_drone)

## Save tree metrics to a CSV file
# Create a data frame with the relevant metrics
tree_metrics_df_mixed_drone <- data.frame(
  TreeHeightDrone = TreeHeightDrone_numeric,
  CrownAreaDrone = CrownAreaDrone_numeric,
  CrownPerimeterDrone = CrownPerimeterDrone_numeric,
  AGB_mixed_drone = AGB_mixed_drone
)

# Print the calculated AGB for each tree
print(tree_metrics_df_mixed_drone)

# Write the tree metrics to a CSV file
write.csv(tree_metrics_df_mixed_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_tree_metrics_agb_mixed_forest.csv", row.names = FALSE)

drone_crowns_df$AGB_mixed_drone <- tree_metrics_df_mixed_drone$AGB_mixed_drone

### Aboveground Biomass (AGB) using deciduous-specific calculation from Jucker et al. 2017###
## Allometric calculation and their coefficients from Jucker et al. https://doi.org/10.1111/gcb.13388 ##
TreeHeightDrone_numeric <- as.numeric(TreeHeightDrone)
CrownAreaDrone_numeric <- as.numeric(CrownAreaDrone)
CrownPerimeterDrone_numeric <- as.numeric(CrownPerimeterDrone)

# Coefficients for deciduous AGB calculation (based on Angiosperm coefficients from Jucker et al. 2017)
alpha_G <- 0
beta_G <- 0

# Convert deciduous max tree height and crown diameter hull variables to numeric
drone_crowns_df$Z <- as.numeric(drone_crowns_df$Z)
drone_crowns_df$crown_diameter_hull_drone <- as.numeric(drone_crowns_df$crown_diameter_hull_drone)

# Allometric AGB Deciduous Calculation
drone_crowns_df$AGB_drone_mixed <- (0.016 + alpha_G) * 
  (drone_crowns_df$Z * drone_crowns_df$crown_diameter_hull_drone)^(2.013 + beta_G) * 
  exp(0.2042 / 2)

print(drone_crowns_df$AGB_drone_mixed)

# Save ABG for deciduous trees to a CSV file
write.csv(drone_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_drone_mixed_jucker.csv", row.names = FALSE)

# Assuming 'AGB_deciduous' is a vector of AGB values for each tree
total_AGB_drone_mixed <- sum(drone_crowns_df$AGB_drone_mixed)

# Print the total AGB for deciduous trees
print(total_AGB_drone_mixed)

# Calculate the product of tree height and crown diameter
drone_crowns_df['tree_height_crown_diameter'] = drone_crowns_df['Z'] * drone_crowns_df['crown_diameter_hull_drone']



### Splitting Home Covert into Deciduous and Coniferous ###
### DECIDUOUS ###
# Read the shapefile for deciduous area of interest (AOI)
aoi_deciduous <- st_read("F:/uni_msc/Dissertation/deciduous_aoi.shp")
if (is.null(aoi_deciduous)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(cleaned_dronepointcloud_3))
print(st_crs(aoi_deciduous))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi_deciduous) != st_crs(cleaned_dronepointcloud_3)) {
  aoi_deciduous <- st_transform(aoi_deciduous, crs = st_crs(cleaned_dronepointcloud_3))
}

## Clip and Validate point cloud
# Clip the segmented tree point cloud using the deciduous shapefile AOI
deciduous_trees_drone <- clip_roi(trees_drone, aoi_deciduous)

# Check if the clipping was successful
if (is.null(deciduous_trees_drone)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped deciduous data
las_check(deciduous_trees_drone)
summary(deciduous_trees_drone)
plot(deciduous_trees_drone)

# Extract the treeID attribute
deciduous_tree_ids_drone <- deciduous_trees_drone@data$treeID
# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
unique_tree_count_decid_drone <- length(unique(deciduous_tree_ids_drone[deciduous_tree_ids_drone > 0]))

# Output the number of detected deciduous trees
cat("Number of trees segmented and detected:", unique_tree_count_decid_drone, "\n")

## Generate deciduous tree crown metrics
deciduous_crowns_drone <- crown_metrics(deciduous_trees_drone, func = .stdtreemetrics, geom = "convex")
# Plot the deciduous crown areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_drone_deciduous.png", width = 600, height = 600)
plot(deciduous_crowns_drone["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the deciduous crowns as a shapefile
st_write(deciduous_crowns_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics_deciduous.shp", append=FALSE)

### Calculate additional tree metrics
# Compute the maximum height for each deciduous tree
maxtree_heights_drone_deciduous <- tree_metrics(deciduous_trees_drone, ~max(Z))
# Save deciduous tree heights to a CSV file
write.csv(maxtree_heights_drone_deciduous, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights_drone_deciduous.csv", append=FALSE)
# Print deciduous tree heights
print(maxtree_heights_drone_deciduous)

max_height_deciduous_drone <- max(maxtree_heights_drone_deciduous$V1, na.rm = TRUE)
avg_height_deciduous_drone <- mean(deciduous_trees_drone$Z, na.rm = TRUE)

print(max_height_deciduous_drone)
print(avg_height_deciduous_drone)

# Convert deciduous crowns to sf object
# Crowns will be an sf object with geometries and associated metrics
deciduous_crowns_sf_drone <- deciduous_crowns_drone

## Calculate deciduous tree metrics
# 1. Crown Area
deciduous_crowns_sf_drone$area <- st_area(deciduous_crowns_sf_drone)  # Calculate area for each crown polygon
# 2. Crown Perimeter
deciduous_crowns_sf_drone$perimeter <- st_length(st_cast(deciduous_crowns_sf_drone, "MULTILINESTRING"))
# 3. Tree Height
max_area_deciduous_drone <- max(deciduous_crowns_sf_drone$area, na.rm = TRUE)
avg_area_deciduous_drone <- mean(deciduous_crowns_sf_drone$area, na.rm = TRUE)
max_perimeter_deciduous_drone <- max(deciduous_crowns_sf_drone$perimeter, na.rm = TRUE)
avg_perimeter_deciduous_drone <- mean(deciduous_crowns_sf_drone$perimeter, na.rm = TRUE)
# Print summary deciduous tree metrics 
print(paste("Max crown area:", max_area_deciduous_drone))
print(paste("Average crown area:", avg_area_deciduous_drone))
print(paste("Max crown perimeter:", max_perimeter_deciduous_drone))
print(paste("Average crown perimeter:", avg_perimeter_deciduous_drone))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(deciduous_crowns_sf_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics_with_sf_deciduous.csv", append = FALSE)

# Convert sf object to a data frame for CSV output
deciduous_crowns_df_drone <- as.data.frame(deciduous_crowns_sf_drone)
# Ensure the geometry column is not included in the CSV
deciduous_crowns_df_drone$geometry <- NULL
# Add summary statistics to the data frame
deciduous_crowns_df_drone$max_area_deciduous_drone <- max_area_deciduous_drone
deciduous_crowns_df_drone$avg_area_deciduous_drone <- avg_area_deciduous_drone
deciduous_crowns_df_drone$max_perimeter_deciduous_drone <- max_perimeter_deciduous_drone
deciduous_crowns_df_drone$avg_perimeter_deciduous_drone <- avg_perimeter_deciduous_drone

#Calculating Crown diameter using max length of crown (convex hull)
# Assumes you want the longest straight line within the convex hull
deciduous_crowns_sf_drone$crown_diameter_hull_deciduous_drone <- st_length(st_cast(st_convex_hull(deciduous_crowns_sf_drone), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_deciduous_drone <- max(deciduous_crowns_sf_drone$crown_diameter_hull_deciduous_drone, na.rm = TRUE)
avg_diameter_hull_deciduous_drone <- mean(deciduous_crowns_sf_drone$crown_diameter_hull_deciduous_drone, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_deciduous_drone, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_deciduous_drone, "\n")

# Add crown diameters to the data frame
deciduous_crowns_df_drone$crown_diameter_deciduous_drone <- deciduous_crowns_sf_drone$crown_diameter_deciduous_drone
deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone <- deciduous_crowns_sf_drone$crown_diameter_hull_deciduous_drone

# Save the metrics to a CSV file
write.csv(deciduous_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_deciduous_drone.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics_deciduous_drone.csv'.\n")


### Aboveground Biomass (AGB) using deciduous-specific calculation from Jucker et al. 2017###
## Allometric calculation and their coefficients from Jucker et al. https://doi.org/10.1111/gcb.13388 ##
## Generate maximum deciduous tree height
TreeHeightDeciduousDrone <- maxtree_heights_drone_deciduous$V1 # Use max tree height from your data
print(TreeHeightDeciduousDrone)

# Use `deciduous_crowns_df' data frame that includes `avg_area` and `avg_perimeter`
CrownAreaDeciduousDrone <- deciduous_crowns_sf_drone$area # Use average crown area
CrownPerimeterDeciduousDrone <- deciduous_crowns_sf_drone$perimeter  # Use average crown perimeter

# Convert variables to numeric values for deciduous allometric calculation
TreeHeightDeciduousDrone_numeric <- as.numeric(TreeHeightDeciduousDrone)
CrownAreaDeciduousDrone_numeric <- as.numeric(CrownAreaDeciduousDrone)
CrownPerimeterDeciduousDrone_numeric <- as.numeric(CrownPerimeterDeciduousDrone)

# Coefficients for deciduous AGB calculation (based on Angiosperm coefficients from Jucker et al. 2017)
alpha_G <- 0
beta_G <- 0

# Convert deciduous max tree height and crown diameter hull variables to numeric
deciduous_crowns_df_drone$Z <- as.numeric(deciduous_crowns_df_drone$Z)
deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone <- as.numeric(deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone)

# Allometric AGB Deciduous Calculation
deciduous_crowns_df_drone$AGB_deciduous_drone <- (0.016 + alpha_G) * 
  (deciduous_crowns_df_drone$Z * deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone)^(2.013 + beta_G) * 
  exp(0.2042 / 2)

print(deciduous_crowns_df_drone$AGB_deciduous_drone)

# Save ABG for deciduous trees to a CSV file
write.csv(deciduous_crowns_df_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_deciduous_drone.csv", row.names = FALSE)

# Assuming 'AGB_deciduous' is a vector of AGB values for each tree
total_AGB_deciduous_drone <- sum(deciduous_crowns_df_drone$AGB_deciduous_drone)

# Print the total AGB for deciduous trees
print(total_AGB_deciduous_drone)
avg_AGB_deciduous <- mean(deciduous_crowns_df_drone$AGB_deciduous_drone, na.rm = TRUE)
print(paste("Average AGB for deciduous trees:", avg_AGB_deciduous))

# Calculate the product of tree height and crown diameter
deciduous_crowns_df_drone['tree_height_crown_diameter'] = deciduous_crowns_df_drone['Z'] * deciduous_crowns_df_drone['crown_diameter_hull_deciduous_drone']



### Splitting Home Covert into Deciduous and Coniferous ###
### CONIFEROUS###
# Read the shapefile for deciduous area of interest (AOI)
aoi_conifer <- st_read("F:/uni_msc/Dissertation/conifer_aoi.shp")
if (is.null(aoi_conifer)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(trees_drone))
print(st_crs(aoi_conifer))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi_conifer) != st_crs(cleaned_dronepointcloud_3)) {
  aoi_conifer <- st_transform(aoi_conifer, crs = st_crs(cleaned_dronepointcloud_3))
}

## Clip and Validate point cloud
# Clip the segmented tree point cloud using the conifer shapefile AOI
conifer_trees_drone <- clip_roi(trees_drone, aoi_conifer)

# Check if the clipping was successful
if (is.null(conifer_trees_drone)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped conifer data
las_check(conifer_trees_drone)
summary(conifer_trees_drone)
plot(conifer_trees_drone)

# Extract the treeID attribute
conifer_tree_ids_drone <- conifer_trees_drone@data$treeID
# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
unique_tree_count_conifer_drone <- length(unique(conifer_tree_ids_drone[conifer_tree_ids_drone > 0]))

# Output the number of detected conifer trees
cat("Number of trees segmented and detected:", unique_tree_count_conifer_drone, "\n")

## Generate conifer tree crown metrics
conifer_crowns_drone <- crown_metrics(conifer_trees_drone, func = .stdtreemetrics, geom = "convex")
# Plot the conifer crown areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_drone_conifer.png", width = 600, height = 600)
plot(conifer_crowns_drone["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the conifer crowns as a shapefile
st_write(conifer_crowns_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics_conifer.shp", append=FALSE)

### Calculate additional tree metrics
# Compute the maximum height for each conifer tree
maxtree_heights_drone_conifer <- tree_metrics(conifer_trees_drone, ~max(Z))
# Save conifer tree heights to a CSV file
write.csv(maxtree_heights_drone_conifer, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights_drone_conifer.csv", append=FALSE)
# Print conifer tree heights
print(maxtree_heights_drone_conifer)

max_height_conifer_drone <- max(maxtree_heights_drone_conifer$V1, na.rm = TRUE)
avg_height_conifer_drone <- mean(conifer_trees_drone$Z, na.rm = TRUE)

print(max_height_conifer_drone)
print(avg_height_conifer_drone)

# Convert conifer crowns to sf object
# Crowns will be an sf object with geometries and associated metrics
conifer_crowns_sf_drone <- conifer_crowns_drone

## Calculate conifer tree metrics
# 1. Crown Area
conifer_crowns_sf_drone$area <- st_area(conifer_crowns_sf_drone)  # Calculate area for each crown polygon
# 2. Crown Perimeter
conifer_crowns_sf_drone$perimeter <- st_length(st_cast(conifer_crowns_sf_drone, "MULTILINESTRING"))
# 3. Tree Height
max_area_conifer_drone <- max(conifer_crowns_sf_drone$area, na.rm = TRUE)
avg_area_conifer_drone <- mean(conifer_crowns_sf_drone$area, na.rm = TRUE)
max_perimeter_conifer_drone <- max(conifer_crowns_sf_drone$perimeter, na.rm = TRUE)
avg_perimeter_conifer_drone <- mean(conifer_crowns_sf_drone$perimeter, na.rm = TRUE)
# Print summary deciduous tree metrics 
print(paste("Max crown area:", max_area_conifer_drone))
print(paste("Average crown area:", avg_area_conifer_drone))
print(paste("Max crown perimeter:", max_perimeter_conifer_drone))
print(paste("Average crown perimeter:", avg_perimeter_conifer_drone))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(conifer_crowns_sf_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/drone_crown_metrics_with_sf_conifer.csv", append = FALSE)

# Convert sf object to a data frame for CSV output
conifer_crowns_df_drone <- as.data.frame(conifer_crowns_sf_drone)
# Ensure the geometry column is not included in the CSV
conifer_crowns_df_drone$geometry <- NULL
# Add summary statistics to the data frame
conifer_crowns_df_drone$max_area_conifer_drone <- max_area_conifer_drone
conifer_crowns_df_drone$avg_area_conifer_drone <- avg_area_conifer_drone
conifer_crowns_df_drone$max_perimeter_conifer_drone <- max_perimeter_conifer_drone
conifer_crowns_df_drone$avg_perimeter_conifer_drone <- avg_perimeter_conifer_drone

#Calculating Crown diameter using max length of crown (convex hull)
# Assumes you want the longest straight line within the convex hull
conifer_crowns_sf_drone$crown_diameter_hull_conifer_drone <- st_length(st_cast(st_convex_hull(conifer_crowns_sf_drone), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_conifer_drone <- max(conifer_crowns_sf_drone$crown_diameter_hull_conifer_drone, na.rm = TRUE)
avg_diameter_hull_conifer_drone <- mean(conifer_crowns_sf_drone$crown_diameter_hull_conifer_drone, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_conifer_drone, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_conifer_drone, "\n")

# Add crown diameters to the data frame
conifer_crowns_df_drone$crown_diameter_conifer_drone <- conifer_crowns_sf_drone$crown_diameter_conifer_drone
conifer_crowns_df_drone$crown_diameter_hull_conifer_drone <- conifer_crowns_sf_drone$crown_diameter_hull_conifer_drone

# Save the metrics to a CSV file
write.csv(conifer_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer_drone.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics_conifer_drone.csv'.\n")


### Aboveground Biomass (AGB) using conifer-specific calculation from Jucker et al. 2017###
## Allometric calculation and their coefficients from Jucker et al. https://doi.org/10.1111/gcb.13388 ##
## Generate maximum conifer tree height
TreeHeightConiferDrone <- maxtree_heights_drone_conifer$V1 # Use max tree height from your data
print(TreeHeightConiferDrone)

# Use `deciduous_crowns_df' data frame that includes `avg_area` and `avg_perimeter`
CrownAreaConiferDrone <- conifer_crowns_sf_drone$area # Use average crown area
CrownPerimeterConiferDrone <- conifer_crowns_sf_drone$perimeter  # Use average crown perimeter

# Convert variables to numeric values for conifer allometric calculation
TreeHeightConiferDrone_numeric <- as.numeric(TreeHeightConiferDrone)
CrownAreaConiferDrone_numeric <- as.numeric(CrownAreaConiferDrone)
CrownPerimeterConiferDrone_numeric <- as.numeric(CrownPerimeterConiferDrone)

# Coefficients for conifer AGB calculation (based on Angiosperm coefficients from Jucker et al. 2017)
alpha_G <- 0.093
beta_G <- -0.223

# Convert deciduous max tree height and crown diameter hull variables to numeric
conifer_crowns_df_drone$Z <- as.numeric(conifer_crowns_df_drone$Z)
conifer_crowns_df_drone$crown_diameter_hull_conifer_drone <- as.numeric(conifer_crowns_df_drone$crown_diameter_hull_conifer_drone)

# Allometric AGB coniferous Calculation
conifer_crowns_df_drone$AGB_conifer_drone <- (0.016 + alpha_G) * 
  (conifer_crowns_df_drone$Z * conifer_crowns_df_drone$crown_diameter_hull_conifer_drone)^(2.013 + beta_G) * 
  exp(0.2042 / 2)

print(conifer_crowns_df_drone$AGB_conifer_drone)

# Save ABG for deciduous trees to a CSV file
write.csv(conifer_crowns_df_drone, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer_drone.csv", row.names = FALSE)

# Assuming 'AGB_coniferous' is a vector of AGB values for each tree
total_AGB_conifer_drone <- sum(conifer_crowns_df_drone$AGB_conifer_drone)

# Print the total AGB for conifer trees
print(total_AGB_conifer_drone)
avg_AGB_conifer <- mean(conifer_crowns_df_drone$AGB_conifer_drone, na.rm = TRUE)
print(paste("Average AGB for coniferous trees:", avg_AGB_conifer))

# Calculate the product of tree height and crown diameter
conifer_crowns_df_drone['tree_height_crown_diameter'] = conifer_crowns_df_drone['Z'] * conifer_crowns_df_drone['crown_diameter_hull_conifer_drone']


### Creating plots for analysis
## Create the scatter plot with both coniferous and deciduous data
# Tree height x Crown Diameter vs AGB
AGB_scatterplot_drone <- ggplot() +
  # Add coniferous data points
  geom_point(data = conifer_crowns_df_drone, 
             aes(x = tree_height_crown_diameter, y = AGB_conifer_drone, color = "Coniferous"), 
             alpha = 0.7) +
  # Add deciduous data points
  geom_point(data = deciduous_crowns_df_drone, 
             aes(x = tree_height_crown_diameter, y = AGB_deciduous_drone, color = "Deciduous"), 
             alpha = 0.7) +
  # Manually specify the colors for the legend
  scale_color_manual(values = c("Coniferous" = "blue", "Deciduous" = "green")) +
  # Add labels and title
  labs(x = "Tree Height x Crown Diameter", 
       y = "Above Ground Biomass (AGB)",
       color = "Tree Type") +  # This label will appear as the legend title
  # Use fancy theme
  theme_fancy(
  # Make x and y axis titles bigger and bold
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold")
  
)
# Save the plot as a PNG file to the specified path
ggsave("plots/dronescatter_plot_tree_height_vs_AGB.png", 
       AGB_scatterplot_drone, units = "cm", width = 16, height = 10)

print(AGB_scatterplot_drone)
class(ggplot)


## Create violin plots for tree metrics
## Tree type vs Crown Diameter
# Ensure crown diameter columns are numeric
conifer_crowns_df_drone$crown_diameter_hull_conifer_drone <- as.numeric(conifer_crowns_df_drone$crown_diameter_hull_conifer_drone)
deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone <- as.numeric(deciduous_crowns_df_drone$crown_diameter_hull_deciduous_drone)
drone_crowns_df$crown_diameter_hull_drone <- as.numeric(drone_crowns_df$crown_diameter_hull_drone)

# Create violin plot for Woodland Type against Crown Diameter (CD)
violinCD_drone <- ggplot() +
  # Add violin plot for coniferous trees
  geom_violin(data = conifer_crowns_df_drone, 
              aes(x = crown_diameter_hull_conifer_drone, y = factor("Coniferous"), fill = "Coniferous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for deciduous trees
  geom_violin(data = deciduous_crowns_df_drone, 
              aes(x = crown_diameter_hull_deciduous_drone, y = factor("Deciduous"), fill = "Deciduous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for the whole forest
  geom_violin(data = drone_crowns_df, 
              aes(x = crown_diameter_hull_drone, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Add labels and title
  labs(x = "Crown Diameter", 
       y = "Tree Type") +
  # Customize the theme
  theme_minimal() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Coniferous" = "blue", "Deciduous" = "green", "Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_CD_drone.png", 
       violinCD_drone, width = 7, height = 5)

# Print the violin CD plot
print(violinCD_drone)


## Tree type vs AGB Violin Plot
# Ensure crown diameter columns are numeric
conifer_crowns_df_drone$AGB_conifer_drone <- as.numeric(conifer_crowns_df_drone$AGB_conifer_drone)
deciduous_crowns_df_drone$AGB_deciduous_drone <- as.numeric(deciduous_crowns_df_drone$AGB_deciduous_drone)
drone_crowns_df$AGB_mixed_drone <- as.numeric(drone_crowns_df$AGB_mixed_drone)

# Create violin plot for Woodland type against AGB
violinAGB_drone <- ggplot() +
  # Add violin plot for coniferous trees
  geom_violin(data = conifer_crowns_df_drone, 
              aes(x = AGB_conifer_drone, y = factor("Coniferous"), fill = "Coniferous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for deciduous trees
  geom_violin(data = deciduous_crowns_df_drone, 
              aes(x = AGB_deciduous_drone, y = factor("Deciduous"), fill = "Deciduous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for the whole forest
  geom_violin(data = drone_crowns_df, 
              aes(x = AGB_mixed_drone, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Cap x-axis to 5000
  xlim(0, 10000) +
  # Add labels and title
  labs(x = "Aboveground Biomass (AGB)", 
       y = "Tree Type") +
  # Customize the theme
  theme_fancy() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Coniferous" = "blue", "Deciduous" = "green", "Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_drone.png", 
       violinAGB_drone, width = 7, height = 5)

# Print the violin AGB plot
print(violinAGB_drone)



### Creating plots for analysis (using AGB with Jucker et al. (2017) deciduous equation)
## Create the scatter plot with whole woodland data
# Tree height x Crown Diameter vs AGB
AGB_scatterplot_drone_2 <- ggplot() +
  # Add coniferous data points
  geom_point(data = drone_crowns_df, 
             aes(x = tree_height_crown_diameter, y = AGB_drone_mixed, color = "Whole Forest"), 
             alpha = 0.7) +
  # Manually specify the colors for the legend
  scale_color_manual(values = c("Whole Forest" = "purple")) +
  # Add labels and title
  labs(x = "Tree Height x Crown Diameter", 
       y = "Above Ground Biomass (AGB)",
       color = "Tree Type") +  # This label will appear as the legend title
  # Use fancy theme
  theme_fancy()

# Save the plot as a PNG file to the specified path
ggsave("plots/dronescatter_plot_tree_height_vs_AGB_2.png", 
       AGB_scatterplot_drone_2, units = "cm", width = 16, height = 10)

print(AGB_scatterplot_drone_2)
class(ggplot)

## Create violin plots for tree metrics
## Tree type vs Crown Diameter
# Ensure crown diameter columns are numeric
drone_crowns_df$crown_diameter_hull_drone <- as.numeric(drone_crowns_df$crown_diameter_hull_drone)

# Create violin plot for Woodland Type against Crown Diameter (CD)
violinCD_drone_2 <- ggplot() +
  # Add violin plot for the whole forest
  geom_violin(data = drone_crowns_df, 
              aes(x = crown_diameter_hull_drone, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Add labels and title
  labs(x = "Crown Diameter", 
       y = "Tree Type") +
  # Customize the theme
  theme_fancy() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_CD_drone_2.png", 
       violinCD_drone_2, width = 7, height = 5)

# Print the violin CD plot
print(violinCD_drone_2)


# Create violin plot for Woodland type against AGB
violinAGB_drone_2 <- ggplot() +
  # Add violin plot for the whole forest
  geom_violin(data = drone_crowns_df, 
              aes(x = AGB_drone_mixed, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Cap x-axis to 5000
  xlim(0, 5000) +
  # Add labels and title
  labs(x = "Aboveground Biomass (AGB)", 
       y = "Tree Type") +
  # Customize the theme
  theme_fancy() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_drone_2.png", 
       violinAGB_drone_2, width = 7, height = 5)

# Print the violin AGB plot
print(violinAGB_drone_2)


