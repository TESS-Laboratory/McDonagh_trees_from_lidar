## Project area of interest: Home Covert woodland at Oxburgh Estate, Oxborough, Norfolk
## Using Environment Agency 2018 LiDAR point cloud data for tree metrics extraction
## and estimation of Aboveground Biomass using allometric equations
## Author of code: Amber McDonagh am1287@exeter.ac.uk  
## Environment Agency 2018 LiDAR point cloud data available at: 
## https://www.data.gov.uk/dataset/f0db0249-f17b-4036-9e65-309148c97ce4/national-lidar-programme 

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

## Read and Prepare Environment Agency (EA) Data 
# Read in EA point cloud
EAPointCloud2018 <- readLAS("F:/uni_msc/Dissertation/data/TF7000_P_10787_20181113_20181117.laz")
if (is.null(EAPointCloud2018)) stop("Failed to read the point cloud file.")
print(EAPointCloud2018)
summary(EAPointCloud2018)
las_check(EAPointCloud2018)

# Plotting - basic 3D rendering
plot(EAPointCloud2018)
plot(EAPointCloud2018, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Read the shapefile (Area of Interest - AOI, Home Covert woodland)
aoi <- st_read("F:/uni_msc/Dissertation/TopographicArea.shp")
if (is.null(aoi)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(EAPointCloud2018))
print(st_crs(aoi))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi) != st_crs(EAPointCloud2018)) {
  aoi <- st_transform(aoi, crs = st_crs(EAPointCloud2018))
}

## Clip and Validate point cloud
# Clip the point cloud using the shapefile AOI
clipped_las_EA_2018 <- clip_roi(EAPointCloud2018, aoi)


# Check if the clipping was successful
if (is.null(clipped_las_EA_2018)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped data
las_check(clipped_las_EA_2018)
summary(clipped_las_EA_2018)

## Process Data
# Remove duplicate points
clipped_las_EA_2018 <- filter_duplicates(clipped_las_EA_2018)

# Validate the clipped data after removing duplicates
las_check(clipped_las_EA_2018)

# Plot the clipped point cloud
plot(clipped_las_EA_2018)
plot(clipped_las_EA_2018, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)

# Optionally, save the clipped point cloud to a new file
writeLAS(clipped_las_EA_2018, "D:/datafilesfordissertation/interim_pointclouds/clipped_point_cloud_EA_2018.laz")

## Normalise the point cloud without DTM (5.2 in lidar lidRbook)
normclipped_las_EA_2018 <- normalize_height(clipped_las_EA_2018, knnidw())

# Check ground points are exactly 0
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

writeLAS(clipped_las_EA_2018, "D:/datafilesfordissertation/interim_pointclouds/cleaned_las_2018.las")
cleaned_las_2018 <- readLAS("D:/datafilesfordissertation/interim_pointclouds/cleaned_las_2018.las")
plot(cleaned_las_2018)

### Digital Terrain Model DTM  ###
# DTM clipped 2018
dtmclipped2018 <- grid_terrain(clipped_las_EA_2018, res = 1, algorithm = tin(), theme_fancy())
plot(dtmclipped2018)

# Save the DTM clipped 2018 plot to a file
plot(dtmclipped2018, main = "Digital Terrain Model (DTM)")
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/dtmclipped2018_plot_2.png", width = 16, height = 10)
dev.off()

# Render shaded DTM using terra and tin method
dtmclipped2018_terra <- rasterize_terrain(clipped_las_EA_2018, algorithm = tin(), pkg ="terra")
dtm_prod_2018 <- terrain(dtmclipped2018_terra, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade_2018 <- shade(slope = dtm_prod_2018$slope, aspect = dtm_prod_2018$aspect)
plot(dtm_hillshade_2018, col =gray(0:30/30), legend = FALSE, main = "Digital Terrain Model (DTM)")

# Save the DTM plot to a file
png(filename = "F:/uni_msc/Dissertation/dtmclipped2018_terra_plot.png", width = 800, height = 600)
plot(dtm_hillshade_2018, col =gray(0:30/30), legend = FALSE, main = "Digital Terrain Model (DTM)")
dev.off()

### Individual tree segmentation (ITS) (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 1, dt2 = 2)
# segment trees
trees_2018 <- segment_trees(cleaned_las_2018, algo)
# Plot the segmented trees
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/trees2018.png", width = 800, height = 600)
plot(trees_2018, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(trees_2018, "D:/datafilesfordissertation/interim_pointclouds/segmented_trees_EA_2018.laz")
trees_2018 <- readLAS("C:/workspace/McDonagh_trees_from_lidar/outputs/segmented_trees_EA_2018.laz")


##Counting 2018 segmented trees
# Extract the treeID attribute
tree_ids <- trees_2018@data$treeID

# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
# Depending on the segmentation method, treeID might be negative for some non-tree points.
unique_tree_count <- length(unique(tree_ids[tree_ids > 0]))

# Output the number of detected trees
cat("Number of trees segmented and detected:", unique_tree_count, "\n")

# Filters tree number 110
#tree110 <- filter_poi(trees_2018, treeID == 110)
#plot(tree110, size = 8, bg = "white")

## Generate Crown Metrics
crowns <- crown_metrics(trees_2018, func = .stdtreemetrics, geom = "convex")
# Plot the crown areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_2018.png", width = 600, height = 600)
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the crowns as a shapefile
st_write(crowns, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics.shp", append=FALSE)

### Calculate additional tree metrics
# Compute the maximum height for each tree
maxtree_heights_2018 <- tree_metrics(trees_2018, ~max(Z))
# Save tree heights to a CSV file
write.csv(maxtree_heights_2018, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights2018.csv")
# Print tree heights
print(maxtree_heights_2018)

# Compute maximum and average tree height
max_height <- max(trees_2018$height, na.rm = TRUE)
avg_height <- mean(trees_sf$height, na.rm = TRUE)

print(max_height)
print(avg_height)

# Convert crowns to sf object
# crowns will be an sf object with geometries and associated metrics
crowns_sf <- crowns

## Calculate tree metrics using sf
# 1. Crown Area
crowns_sf$area <- st_area(crowns_sf)  # Calculate area for each crown polygon
# 2. Crown Perimeter
crowns_sf$perimeter <- st_length(st_cast(crowns_sf, "MULTILINESTRING"))
# 3. Tree Height
max_area <- max(crowns_sf$area, na.rm = TRUE)
avg_area <- mean(crowns_sf$area, na.rm = TRUE)
max_perimeter <- max(crowns_sf$perimeter, na.rm = TRUE)
avg_perimeter <- mean(crowns_sf$perimeter, na.rm = TRUE)
# Print summary tree metrics 
print(paste("Max crown area:", max_area))
print(paste("Average crown area:", avg_area))
print(paste("Max crown perimeter:", max_perimeter))
print(paste("Average crown perimeter:", avg_perimeter))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(crowns_sf, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_with_sf.shp", append=FALSE)

# Convert sf object to a data frame for CSV output
crowns_df <- as.data.frame(crowns_sf)
# Ensure the geometry column is not included in the CSV
crowns_df$geometry <- NULL
# Add summary statistics to the data frame
crowns_df$max_area <- max_area
crowns_df$avg_area <- avg_area
crowns_df$max_perimeter <- max_perimeter
crowns_df$avg_perimeter <- avg_perimeter

## Calculating Crown diameter using max length of crown (convex hull)
crowns_sf$crown_diameter_hull <- st_length(st_cast(st_convex_hull(crowns_sf), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull <- max(crowns_sf$crown_diameter_hull, na.rm = TRUE)
avg_diameter_hull <- mean(crowns_sf$crown_diameter_hull, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull, "\n")

# Add crown diameters to the data frame
crowns_df$crown_diameter <- crowns_sf$crown_diameter
crowns_df$crown_diameter_hull <- crowns_sf$crown_diameter_hull

# Save the metrics to a CSV file
write.csv(crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics.csv'.\n")

### ITS using canopy height model (rasterising) ###
### Canopy Height Model (CHM) ###
# Generate a Canopy Height Model (CHM)
chm_2018 <- grid_canopy(cleaned_las_2018, res = 0.5, p2r())
# Plot the CHM
plot(chm_2018, main = "Canopy Height Model (CHM)")

### Aboveground Biomass (AGB) using mixed temperate woodland calculation ###
## Implementing AGB calculations using height, crown area and perimeter
TreeHeight <- maxtree_heights_2018$V1 # Use max tree height from data
print(TreeHeight)

# Use `crowns_df' data frame that includes `avg_area` and `avg_perimeter`
CrownArea <- crowns_df$area # Use average crown area
CrownPerimeter <- crowns_df$perimeter  # Use average crown perimeter

# Coefficients for mixed broadleaf and conifer forest
beta_0_mixed <- 2.5
beta_1_mixed <- 0.9
beta_2_mixed <- 0.5
beta_3_mixed <- 0.3

# Convert variables to numeric
TreeHeight_numeric <- as.numeric(TreeHeight)
CrownArea_numeric <- as.numeric(CrownArea)
CrownPerimeter_numeric <- as.numeric(CrownPerimeter)

# AGB calculation using the allometric model
AGB_mixed <- beta_0_mixed + 
  beta_1_mixed * log(TreeHeight_numeric) + 
  beta_2_mixed * log(CrownArea_numeric) + 
  beta_3_mixed * log(CrownPerimeter_numeric)

# Exponentiate to get the actual AGB in the original units
AGB_mixed <- exp(AGB_mixed)

# Print the calculated AGB
print(AGB_mixed)

# Assuming 'AGB' is a vector of AGB values for each tree
total_AGB_mixed <- sum(AGB_mixed)

# Print the total AGB for the forest
print(total_AGB_mixed)

## Save tree metrics to a CSV file
# Create a data frame with the relevant metrics
tree_metrics_df_mixed <- data.frame(
  TreeHeight = TreeHeight_numeric,
  CrownArea = CrownArea_numeric,
  CrownPerimeter = CrownPerimeter_numeric,
  AGB_mixed = AGB_mixed
)

# Print the calculated AGB for each tree
print(tree_metrics_df_mixed)

# Write the tree metrics to a CSV file
write.csv(tree_metrics_df_mixed, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_agb_mixed_forest.csv", row.names = FALSE)

crowns_df$AGB_mixed <- tree_metrics_df_mixed$AGB_mixed


### Splitting Home Covert into Deciduous and Coniferous ###
### DECIDUOUS ###
# Read the shapefile for deciduous area of interest (AOI)
aoi_deciduous <- st_read("F:/uni_msc/Dissertation/deciduous_aoi.shp")
if (is.null(aoi_deciduous)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(cleaned_las_2018))
print(st_crs(aoi_deciduous))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi_deciduous) != st_crs(cleaned_las_2018)) {
  aoi_deciduous <- st_transform(aoi_deciduous, crs = st_crs(cleaned_las_2018))
}

## Clip and Validate point cloud
# Clip the point cloud using the shapefile AOI
deciduous_forest <- clip_roi(cleaned_las_2018, aoi_deciduous)

# Check if the clipping was successful
if (is.null(deciduous_forest)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped deciduous data
las_check(deciduous_forest)
summary(deciduous_forest)
plot(deciduous_forest)

### Individual tree segmentation (ITS) (7.2.1) for Deciduous area ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 1, dt2 = 2)
# segment deciduous trees
deciduous_trees_2018 <- segment_trees(deciduous_forest, algo)
# Plot the segmented trees
plot(deciduous_trees_2018, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/deciduoustrees2018.png", width = 800, height = 600)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(deciduous_trees_2018, "C:/workspace/McDonagh_trees_from_lidar/outputs/segmented_trees_deciduous_2018.laz")
deciduous_trees_2018 <- readLAS("C:/workspace/McDonagh_trees_from_lidar/outputs/segmented_trees_deciduous_2018.laz")

##Counting deciduous segmented trees
# Extract the treeID attribute
deciduous_tree_ids <- deciduous_trees_2018@data$treeID
# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
unique_tree_count <- length(unique(deciduous_tree_ids[deciduous_tree_ids > 0]))

# Output the number of detected deciduous trees
cat("Number of trees segmented and detected:", unique_tree_count, "\n")

## Generate deciduous tree crown metrics
deciduous_crowns <- crown_metrics(deciduous_trees_2018, func = .stdtreemetrics, geom = "convex")
# Plot the deciduous crown areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_2018_deciduous.png", width = 600, height = 600)
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the deciduous crowns as a shapefile
st_write(deciduous_crowns, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_deciduous.shp", append=FALSE)

### Calculate additional tree metrics
# Compute the maximum height for each deciduous tree
maxtree_heights_2018_deciduous <- tree_metrics(deciduous_trees_2018, ~max(Z))
# Save deciduous tree heights to a CSV file
write.csv(maxtree_heights_2018_deciduous, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights2018_deciduous.csv", append=FALSE)
# Print deciduous tree heights
print(maxtree_heights_2018_deciduous)

max_height_deciduous <- max(maxtree_heights_2018_deciduous$V1, na.rm = TRUE)
avg_height_deciduous <- mean(deciduous_trees_2018$Z, na.rm = TRUE)

print(max_height_deciduous)
print(avg_height_deciduous)

# Convert deciduous crowns to sf object
# Crowns will be an sf object with geometries and associated metrics
deciduous_crowns_sf <- deciduous_crowns

## Calculate deciduous tree metrics
# 1. Crown Area
deciduous_crowns_sf$area <- st_area(deciduous_crowns_sf)  # Calculate area for each crown polygon
# 2. Crown Perimeter
deciduous_crowns_sf$perimeter <- st_length(st_cast(deciduous_crowns_sf, "MULTILINESTRING"))
# 3. Tree Height
max_area_deciduous <- max(deciduous_crowns_sf$area, na.rm = TRUE)
avg_area_deciduous <- mean(deciduous_crowns_sf$area, na.rm = TRUE)
max_perimeter_deciduous <- max(deciduous_crowns_sf$perimeter, na.rm = TRUE)
avg_perimeter_deciduous <- mean(deciduous_crowns_sf$perimeter, na.rm = TRUE)
# Print summary deciduous tree metrics 
print(paste("Max crown area:", max_area_deciduous))
print(paste("Average crown area:", avg_area_deciduous))
print(paste("Max crown perimeter:", max_perimeter_deciduous))
print(paste("Average crown perimeter:", avg_perimeter_deciduous))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(deciduous_crowns_sf, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_with_sf_deciduous.csv", append = FALSE)

# Convert sf object to a data frame for CSV output
deciduous_crowns_df <- as.data.frame(deciduous_crowns_sf)
# Ensure the geometry column is not included in the CSV
deciduous_crowns_df$geometry <- NULL
# Add summary statistics to the data frame
deciduous_crowns_df$max_area_deciduous <- max_area_deciduous
deciduous_crowns_df$avg_area_deciduous <- avg_area_deciduous
deciduous_crowns_df$max_perimeter_deciduous <- max_perimeter_deciduous
deciduous_crowns_df$avg_perimeter_deciduous <- avg_perimeter_deciduous

#Calculating Crown diameter using max length of crown (convex hull)
# Assumes you want the longest straight line within the convex hull
deciduous_crowns_sf$crown_diameter_hull_deciduous <- st_length(st_cast(st_convex_hull(deciduous_crowns_sf), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_deciduous <- max(deciduous_crowns_sf$crown_diameter_hull_deciduous, na.rm = TRUE)
avg_diameter_hull_deciduous <- mean(deciduous_crowns_sf$crown_diameter_hull_deciduous, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_deciduous, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_deciduous, "\n")

# Add crown diameters to the data frame
deciduous_crowns_df$crown_diameter_deciduous <- deciduous_crowns_sf$crown_diameter_deciduous
deciduous_crowns_df$crown_diameter_hull_deciduous <- deciduous_crowns_sf$crown_diameter_hull_deciduous

# Save the metrics to a CSV file
write.csv(deciduous_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_deciduous.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics_deciduous.csv'.\n")

### Aboveground Biomass (AGB) using deciduous-specific calculation from Jucker et al. 2017###
## Allometric calculation and their coefficients from Jucker et al. https://doi.org/10.1111/gcb.13388 ##
## Generate maximum deciduous tree height
TreeHeightDeciduous <- maxtree_heights_2018_deciduous$V1 # Use max tree height from your data
print(TreeHeightDeciduous)

# Use `deciduous_crowns_df' data frame that includes `avg_area` and `avg_perimeter`
CrownAreaDeciduous <- deciduous_crowns_sf$area # Use average crown area
CrownPerimeterDeciduous <- deciduous_crowns_sf$perimeter  # Use average crown perimeter

# Convert variables to numeric values for deciduous allometric calculation
TreeHeightDeciduous_numeric <- as.numeric(TreeHeightDeciduous)
CrownArea_numeric_Deciduous <- as.numeric(CrownAreaDeciduous)
CrownPerimeter_numeric_Deciduous <- as.numeric(CrownPerimeterDeciduous)

# Coefficients for deciduous AGB calculation (based on Angiosperm coefficients from Jucker et al. 2017)
alpha_G <- 0
beta_G <- 0

# Convert deciduous max tree height and crown diameter hull variables to numeric
drone_crowns_df$Z <- as.numeric(deciduous_crowns_df$Z)
drone_crowns_df$crown_diameter_hull_deciduous <- as.numeric(deciduous_crowns_df$crown_diameter_hull_deciduous)

# Allometric AGB Deciduous Calculation
deciduous_crowns_df$AGB_deciduous <- (0.016 + alpha_G) * 
  (deciduous_crowns_df$Z * deciduous_crowns_df$crown_diameter_hull_deciduous)^(2.013 + beta_G) * 
  exp(0.2042 / 2)

print(deciduous_crowns_df$AGB_deciduous)

# Save ABG for deciduous trees to a CSV file
write.csv(deciduous_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_deciduous.csv", row.names = FALSE)

# Assuming 'AGB_deciduous' is a vector of AGB values for each tree
total_AGB_deciduous <- sum(deciduous_crowns_df$AGB_deciduous)

# Print the total AGB for deciduous trees
print(total_AGB_deciduous)

# Calculate the product of tree height and crown diameter
deciduous_crowns_df['tree_height_crown_diameter'] = deciduous_crowns_df['Z'] * deciduous_crowns_df['crown_diameter_hull_deciduous']





###CONIFEROUS###
# Read the shapefile for coniferous area of interest (AOI)
aoi_conifer <- st_read("F:/uni_msc/Dissertation/conifer_aoi.shp")
if (is.null(aoi)) stop("Failed to read the shapefile.")

# Check CRS of the point cloud and shapefile
print(st_crs(cleaned_las_2018))
print(st_crs(aoi_conifer))

# Ensure the shapefile and point cloud have the same CRS
if (st_crs(aoi_conifer) != st_crs(cleaned_las_2018)) {
  aoi <- st_transform(aoi_conifer, crs = st_crs(cleaned_las_2018))
}

## Clip and Validate point cloud
# Clip the point cloud using the shapefile AOI
conifer_forest <- clip_roi(cleaned_las_2018, aoi_conifer)

# Check if the clipping was successful
if (is.null(conifer_forest)) stop("Clipping the point cloud to the AOI failed.")

# Validate the clipped coniferous data
las_check(conifer_forest)
summary(conifer_forest)
plot(conifer_forest)

### Individual tree segmentation (ITS) for Coniferous area (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 0.5, dt2 = 1) # Changed li2012 parameters
# segment trees
conifer_trees_2018 <- segment_trees(conifer_forest, algo)
# Plot the segmented trees
plot(conifer_trees_2018, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/conifertrees2018.png", width = 800, height = 600)
dev.off()
# Save the segmented trees point cloud to a new file
writeLAS(conifer_trees_2018, "C:/workspace/McDonagh_trees_from_lidar/outputs/segmented_trees_conifer_2018.laz")

##Counting coniferous segmented trees
# Extract the treeID attribute
conifer_tree_ids <- conifer_trees_2018@data$treeID
# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
unique_tree_count <- length(unique(conifer_tree_ids[conifer_tree_ids > 0]))

# Output the number of detected coniferous trees
cat("Number of trees segmented and detected:", unique_tree_count, "\n")

## Generate coniferous crown metrics
conifer_crowns <- crown_metrics(conifer_trees_2018, func = .stdtreemetrics, geom = "convex")
# Plot the crown coniferous areas
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_2018_conifer.png", width = 600, height = 600)
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
dev.off()
# Save the coniferous crowns as a shapefile
st_write(conifer_crowns, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_conifer.shp", append = FALSE)

## Calculate additional tree metrics
# Compute the maximum height for each coniferous tree
maxtree_heights_2018_conifer <- tree_metrics(conifer_trees_2018, ~max(Z))
# Save coniferous maximum tree heights to a CSV file
write.csv(maxtree_heights_2018_conifer, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights2018_conifer.csv")
# Print coniferous tree heights
print(maxtree_heights_2018_conifer)

max_height_conifer <- max(maxtree_heights_2018_conifer$V1, na.rm = TRUE)
avg_height_conifer <- mean(conifer_trees_2018$Z, na.rm = TRUE)

print(max_height_conifer)
print(avg_height_conifer)

# Convert coniferous crowns to sf object
# crowns will be an sf object with geometries and associated metrics
conifer_crowns_sf <- conifer_crowns

## Calculate coniferous tree metrics
# 1. Crown Area
conifer_crowns_sf$area <- st_area(conifer_crowns_sf)  # Calculate area for each crown polygon
# 2. Crown Perimeter
conifer_crowns_sf$perimeter <- st_length(st_cast(conifer_crowns_sf, "MULTILINESTRING"))
# 3. Tree Height
max_area_conifer <- max(conifer_crowns_sf$area, na.rm = TRUE)
avg_area_conifer <- mean(conifer_crowns_sf$area, na.rm = TRUE)
max_perimeter_conifer <- max(conifer_crowns_sf$perimeter, na.rm = TRUE)
avg_perimeter_conifer <- mean(conifer_crowns_sf$perimeter, na.rm = TRUE)
# Print summary coniferous tree metrics 
print(paste("Max crown area:", max_area_conifer))
print(paste("Average crown area:", avg_area_conifer))
print(paste("Max crown perimeter:", max_perimeter_conifer))
print(paste("Average crown perimeter:", avg_perimeter_conifer))
# Save the crowns with metrics as a shapefile for further use or visualisation
st_write(conifer_crowns_sf, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_with_sf_conifer.shp", append = FALSE)

# Convert sf object to a data frame for CSV output
conifer_crowns_df <- as.data.frame(conifer_crowns_sf)
# Ensure the geometry column is not included in the CSV
conifer_crowns_df$geometry <- NULL
# Add summary statistics to the data frame
conifer_crowns_df$max_area_conifer <- max_area_conifer
conifer_crowns_df$avg_area_conifer <- avg_area_conifer
conifer_crowns_df$max_perimeter_conifer <- max_perimeter_conifer
conifer_crowns_df$avg_perimeter_conifer <- avg_perimeter_conifer

#Calculating Crown diameter using max length of crown (convex hull)
# Assumes you want the longest straight line within the convex hull
conifer_crowns_sf$crown_diameter_hull_conifer <- st_length(st_cast(st_convex_hull(conifer_crowns_sf), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_conifer <- max(conifer_crowns_sf$crown_diameter_hull_conifer, na.rm = TRUE)
avg_diameter_hull_conifer <- mean(conifer_crowns_sf$crown_diameter_hull_conifer, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_conifer, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_conifer, "\n")

# Add crown diameters to the data frame
conifer_crowns_df$crown_diameter_conifer <- conifer_crowns_sf$crown_diameter_conifer
conifer_crowns_df$crown_diameter_hull_conifer <- conifer_crowns_sf$crown_diameter_hull_conifer

# Save the metrics to a CSV file
write.csv(conifer_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics_conifer.csv'.\n")

### Aboveground Biomass (AGB) using coniferous-specific calculation from Jucker et al. 2017 ###
#TreeHeightConifer <- maxtree_heights_2018_conifer$V1 # Use max tree height from data
#print(TreeHeightConifer)

# Assuming `conifer_crowns_df` is a data frame that includes `avg_area` and `avg_perimeter`
CrownAreaConifer <- conifer_crowns_sf$area # Use average crown area
CrownPerimeterConifer <- conifer_crowns_sf$perimeter  # Use average crown perimeter

# Coefficients for coniferous AGB calculation (based on Gymnosperm coefficients from Jucker et al. 2017)
alpha_G_2 <- 0.093
beta_G_2 <- -0.223

# Convert coniferous max tree height and crown diameter hull variables to numeric
conifer_crowns_df$Z <- as.numeric(conifer_crowns_df$Z)
conifer_crowns_df$crown_diameter_hull_conifer <- as.numeric(conifer_crowns_df$crown_diameter_hull_conifer)

# Allometric AGB Coniferous Calculation
conifer_crowns_df$AGB_conifer <- (0.016 + alpha_G_2) * 
  (conifer_crowns_df$Z * conifer_crowns_df$crown_diameter_hull_conifer)^(2.013 + beta_G_2) * 
  exp(0.2042 / 2)

print(conifer_crowns_df$AGB_conifer)

# Save ABG for coniferous trees to a CSV file
write.csv(conifer_crowns_df, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer.csv", row.names = FALSE)

# Assuming 'AGB_conifer' is a vector of AGB values for each tree
total_AGB_conifer <- sum(conifer_crowns_df$AGB_conifer)

# Print the total AGB for coniferous trees
print(total_AGB_conifer)

# Calculate the product of tree height and crown diameter
conifer_crowns_df['tree_height_crown_diameter'] = conifer_crowns_df['Z'] * conifer_crowns_df['crown_diameter_hull_conifer']

#coniferttops <- locate_trees(conifer_forest, lmf(ws = 5))

#x <- plot(conifer_forest, bg = "white", size = 4)
#add_treetops3d(x, coniferttops)



## CONIFEROUS  - with density reduction and changed LI2012 parameters
# Conifer Density reduction
# Decimate the point cloud using a random point sub sampling method
# Define the proportion of points you want to retain
desired_proportion <- 0.5  # Retain 50% of the points

# Randomly sample points to reduce density
conifer_forest_reduced <- filter_poi(conifer_forest, runif(npoints(conifer_forest)) < desired_proportion)

# Validate and summarize the reduced coniferous point cloud
las_check(conifer_forest_reduced)
summary(conifer_forest_reduced)
plot(conifer_forest_reduced)

### Height normalization ###
## Normalize the point cloud without DTM (5.2 in lidar lidRbook)
normconifer_forest_reduced <- normalize_height(conifer_forest_reduced, knnidw())
# check ground points are exactly 0
hist(filter_ground(normconifer_forest_reduced)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
las_check(normconifer_forest_reduced)

# Filter out points with Z < 0
conifer_forest_cleaned <- filter_poi(normconifer_forest_reduced, Z >= 0)
plot(conifer_forest_cleaned)

# Validate the cleaned data
las_check(conifer_forest_cleaned)
summary(conifer_forest_cleaned)

# Validate the cleaned data
las_check(conifer_forest_cleaned)

### Individual tree segmentation (ITS) for Coniferous (7.2.1) ###
# Set the segmentation algorithm (li2012 in this case)
algo <- li2012(dt1 = 0.5, dt2 = 1) # changed parameters
# segment coniferous trees
conifer_trees_cleaned <- segment_trees(conifer_forest_cleaned, algo)
# Plot the coniferous segmented trees
plot(conifer_trees_cleaned, color = "treeID", bg = "white", axis = TRUE, legend = FALSE)
png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/conifertrees2018.png", width = 800, height = 600)
dev.off()

##Counting coniferous segmented trees
# Extract the treeID attribute
conifer_tree_ids_cleaned <- conifer_trees_cleaned@data$treeID
# Count unique tree IDs, excluding non-tree points (assuming treeID = 0 or NA for non-tree points)
unique_tree_count <- length(unique(conifer_tree_ids_cleaned[conifer_tree_ids_cleaned > 0]))

# Output the number of detected coniferous trees
cat("Number of trees segmented and detected:", unique_tree_count, "\n")

## Generate coniferous tree crown metrics
conifer_crowns_2 <- crown_metrics(conifer_trees_cleaned, func = .stdtreemetrics, geom = "convex")
# Plot the crown areas
#png(filename = "C:/workspace/McDonagh_trees_from_lidar/plots/crownarea_2018_conifer.png", width = 600, height = 600)
#plot(crowns["convhull_area"], main = "Crown area (convex hull)")
#dev.off()
# Save the coniferous crowns as a shapefile
st_write(conifer_crowns_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_conifer_2.shp", append = FALSE)

## Calculate additional tree metrics
# Compute the maximum height for each tree
maxtree_heights_2018_conifer_2 <- tree_metrics(conifer_trees_cleaned, ~max(Z))
# Save tree heights to a CSV file
write.csv(maxtree_heights_2018_conifer_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/maxtree_heights2018_conifer_2.csv")
# Print tree heights
print(maxtree_heights_2018_conifer_2)

max_height_conifer_2 <- max(maxtree_heights_2018_conifer_2$V1, na.rm = TRUE)
avg_height_conifer_2 <- mean(conifer_trees_cleaned$Z, na.rm = TRUE)

print(max_height_conifer_2)
print(avg_height_conifer_2)

# Convert to an sf object
conifer_crowns_2 <- crown_metrics(conifer_trees_cleaned, func = .stdtreemetrics, geom = "convex")
# Convert coniferous crowns to sf object
# crowns will be an sf object with geometries and associated metrics
conifer_crowns_sf_2 <- conifer_crowns_2

## Calculate coniferous tree metrics
# 1. Crown Area
conifer_crowns_sf_2$area <- st_area(conifer_crowns_sf_2)  # Calculate area for each crown polygon
# 2. Crown Perimeter
conifer_crowns_sf_2$perimeter <- st_length(st_cast(conifer_crowns_sf_2, "MULTILINESTRING"))
# 3. Tree Height
max_area_conifer_2 <- max(conifer_crowns_sf_2$area, na.rm = TRUE)
avg_area_conifer_2 <- mean(conifer_crowns_sf_2$area, na.rm = TRUE)
max_perimeter_conifer_2 <- max(conifer_crowns_sf_2$perimeter, na.rm = TRUE)
avg_perimeter_conifer_2 <- mean(conifer_crowns_sf_2$perimeter, na.rm = TRUE)
# Print summary tree metrics 
print(paste("Max crown area:", max_area_conifer_2))
print(paste("Average crown area:", avg_area_conifer_2))
print(paste("Max crown perimeter:", max_perimeter_conifer_2))
print(paste("Average crown perimeter:", avg_perimeter_conifer_2))
# Save the crowns with metrics as a shapefile for further use or visualization
st_write(conifer_crowns_sf_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/crown_metrics_with_sf_conifer_2.shp", append = FALSE)

# Convert sf object to a data frame for CSV output
conifer_crowns_df_2 <- as.data.frame(conifer_crowns_sf_2)
# Ensure the geometry column is not included in the CSV
conifer_crowns_df_2$geometry <- NULL
# Add summary statistics to the data frame
conifer_crowns_df_2$max_area_conifer_2 <- max_area_conifer_2
conifer_crowns_df_2$avg_area_conifer_2 <- avg_area_conifer_2
conifer_crowns_df_2$max_perimeter_conifer_2 <- max_perimeter_conifer_2
conifer_crowns_df_2$avg_perimeter_conifer_2 <- avg_perimeter_conifer_2

#Calculating Crown diameter using max length of crown (convex hull)
# Assumes you want the longest straight line within the convex hull
conifer_crowns_sf_2$crown_diameter_hull_conifer_2 <- st_length(st_cast(st_convex_hull(conifer_crowns_sf_2), "LINESTRING"))

# Print the max and average crown diameters using convex hull
max_diameter_hull_conifer_2 <- max(conifer_crowns_sf_2$crown_diameter_hull_conifer_2, na.rm = TRUE)
avg_diameter_hull_conifer_2 <- mean(conifer_crowns_sf_2$crown_diameter_hull_conifer_2, na.rm = TRUE)

cat("Max crown diameter (hull):", max_diameter_hull_conifer_2, "\n")
cat("Average crown diameter (hull):", avg_diameter_hull_conifer_2, "\n")

# Add crown diameters to the data frame
conifer_crowns_df_2$crown_diameter_conifer_2 <- conifer_crowns_sf_2$crown_diameter_conifer_2
conifer_crowns_df_2$crown_diameter_hull_conifer_2 <- conifer_crowns_sf_2$crown_diameter_hull_conifer_2

# Save the metrics to a CSV file
write.csv(conifer_crowns_df_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer_2.csv", row.names = FALSE)
# Print a message to confirm the file was saved
cat("Tree metrics have been saved to 'tree_metrics_conifer_2.csv'.\n")


### Aboveground Biomass (AGB) using coniferous-specific calculation from Jucker et al. 2017###
# Assuming `coniferous_crowns_df_2` is a data frame that includes `avg_area` and `avg_perimeter`
CrownAreaConifer_2 <- conifer_crowns_sf_2$area # Use average crown area
CrownPerimeterConifer_2 <- conifer_crowns_sf_2$perimeter  # Use average crown perimeter

# Coefficients for coniferous AGB calculation (based on Gymnosperm coefficients from Jucker et al. 2017)
alpha_G_2 <- 0.093
beta_G_2 <- -0.223

# Convert coniferous max tree height and crown diameter hull variables to numeric
conifer_crowns_df_2$Z <- as.numeric(conifer_crowns_df_2$Z)
conifer_crowns_df_2$crown_diameter_hull_conifer_2 <- as.numeric(conifer_crowns_df_2$crown_diameter_hull_conifer_2)

# Allometric AGB Coniferous Calculation
conifer_crowns_df_2$AGB_conifer_2 <- (0.016 + alpha_G_2) * 
  (conifer_crowns_df_2$Z * conifer_crowns_df_2$crown_diameter_hull_conifer_2)^(2.013 + beta_G_2) * 
  exp(0.2042 / 2)

print(conifer_crowns_df_2$AGB_conifer_2)

# Save ABG for coniferous trees to a CSV file
write.csv(conifer_crowns_df_2, "C:/workspace/McDonagh_trees_from_lidar/outputs/tree_metrics_conifer_2.csv", row.names = FALSE)

# Assuming 'AGB_coniferous_2' is a vector of AGB values for each tree
total_AGB_conifer_2 <- sum(conifer_crowns_df_2$AGB_conifer_2)

# Print the total AGB for coniferous trees
print(total_AGB_conifer_2)

# Calculate the product of tree height and crown diameter
conifer_crowns_df_2['tree_height_crown_diameter'] = conifer_crowns_df_2['Z'] * conifer_crowns_df_2['crown_diameter_hull_conifer_2']



### Creating plots for analysis
## Create the scatter plot with both coniferous and deciduous data
# Tree height x Crown Diameter vs AGB
AGB_scatterplot <- ggplot() +
  # Add coniferous data points
  geom_point(data = conifer_crowns_df, 
             aes(x = tree_height_crown_diameter, y = AGB_conifer, color = "Coniferous"), 
             alpha = 0.7) +
  # Add deciduous data points
  geom_point(data = deciduous_crowns_df, 
             aes(x = tree_height_crown_diameter, y = AGB_deciduous, color = "Deciduous"), 
             alpha = 0.7) +
  # Manually specify the colors for the legend
  scale_color_manual(values = c("Coniferous" = "blue", "Deciduous" = "green")) +
  # Add labels and title
  labs(x = "Tree Height x Crown Diameter", 
       y = "Above Ground Biomass (AGB)",
       color = "Tree Type") +  # This label will appear as the legend title
  # Use fancy theme
  theme_fancy()
  theme(
  # Make x and y axis titles bigger and bold
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold")
    
  )
# Save the plot as a PNG file to the specified path
ggsave("plots/2018_scatter_plot_tree_height_vs_AGB.png", 
       AGB_scatterplot, units = "cm", width = 16, height = 10)

class(ggplot)


## Create violin plots for tree metrics
## Tree type vs Crown Diameter
# Ensure crown diameter columns are numeric
conifer_crowns_df_2$crown_diameter_hull_conifer_2 <- as.numeric(conifer_crowns_df_2$crown_diameter_hull_conifer_2)
deciduous_crowns_df$crown_diameter_hull_deciduous <- as.numeric(deciduous_crowns_df$crown_diameter_hull_deciduous)
crowns_df$crown_diameter <- as.numeric(crowns_df$crown_diameter)

# Create violin plot for Woodland Type against Crown Diameter (CD)
violinCD <- ggplot() +
  # Add violin plot for coniferous trees
  geom_violin(data = conifer_crowns_df_2, 
              aes(x = crown_diameter_hull_conifer_2, y = factor("Coniferous"), fill = "Coniferous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for deciduous trees
  geom_violin(data = deciduous_crowns_df, 
              aes(x = crown_diameter_hull_deciduous, y = factor("Deciduous"), fill = "Deciduous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for the whole forest
  geom_violin(data = crowns_df, 
              aes(x = crown_diameter, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Add labels and title
  labs(x = "Crown Diameter", 
       y = "Tree Type") +
  # Customize the theme
  theme_fancy() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Coniferous" = "blue", "Deciduous" = "green", "Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_CD_2018.png", 
       plot = plot, width = 7, height = 5)

# Print the violin CD plot
print(violinCD)


## Tree type vs AGB
# Ensure crown diameter columns are numeric
conifer_crowns_df_2$AGB_conifer_2 <- as.numeric(conifer_crowns_df_2$AGB_conifer_2)
deciduous_crowns_df$AGB_deciduous <- as.numeric(deciduous_crowns_df$AGB_deciduous)
crowns_df$AGB_mixed <- as.numeric(crowns_df$AGB_mixed)

# Create violin plot for Woodland type against AGB
violinAGB <- ggplot() +
  # Add violin plot for coniferous trees
  geom_violin(data = conifer_crowns_df_2, 
              aes(x = AGB_conifer_2, y = factor("Coniferous"), fill = "Coniferous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for deciduous trees
  geom_violin(data = deciduous_crowns_df, 
              aes(x = AGB_deciduous, y = factor("Deciduous"), fill = "Deciduous"), 
              trim = FALSE, alpha = 0.7) +
  # Add violin plot for the whole forest
  geom_violin(data = crowns_df, 
              aes(x = AGB_mixed, y = factor("Whole Forest"), fill = "Whole Forest"), 
              trim = FALSE, alpha = 0.7) +
  # Cap x-axis to 5000
  xlim(0, 2000) +
  # Add labels and title
  labs(x = "Aboveground Biomass (AGB)", 
       y = "Tree Type") +
  # Customize the theme
  theme_fancy() +
  # Optionally, you can adjust the color scheme
  scale_fill_manual(values = c("Coniferous" = "blue", "Deciduous" = "green", "Whole Forest" = "purple"))

ggsave("C:/workspace/McDonagh_trees_from_lidar/plots/violetplot_agb_2018.png", 
       plot = plot, width = 7, height = 5)

# Print the violin AGB plot
print(violinAGB)

