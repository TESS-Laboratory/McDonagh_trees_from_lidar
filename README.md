# McDonagh_trees_from_lidar-
This repository contains R code for processing, analysing and visualising the data presented in *"Quantifying woodland aboveground biomass in the Oxborough Estate in Norfolk using airborne LiDAR"*by Amber McDonagh under the supervision of Dr Andrew Cunliffe and Keith Challis for her MSc in Geographical Information Science (Environment and Sustainability) at the University of Exeter
Data Source: This study uses a combination of data:
1. Airborne LiDAR data in point cloud LAS format (2018) from Environment Agency.
2. Airbone drone LiDAR data in point cloud LAS format (2023), supplied by Keith Challis.
3. National Trust Toographic data supplied from Emapsite.

Scripts
This repository contains the following scripts within the "Script" folder:

OxburghProject_EA2018_LIDAR.R ~  This script contains the importing and analysis of EA 2018 point cloud data. In particular, normalisation, filtering, clipping of point cloud .Extracting tree metrics, applying allometric equations for abovergroound biomass and plotting the results.

OxburghProjectDrone_LIDAR.R ~ This script contains the importing and analysis of National Trust 2023 drone point cloud. In particular, normalisation, filtering, clipping of point cloud. Extracting tree metrics, applying allometric equations for abovergroound biomass and plotting the results. Includes density reduction methods as drone data is very lage file.

Outputs
This repository contains output files that were created using the scripts. In formats of CSV, shapefiles and LAS.

Plots
This repository contains plots and outputs used to visualise the data. In PNG format.
