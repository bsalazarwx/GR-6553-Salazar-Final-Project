# GR-6553-Salazar-Final-Project

This is the final project for Dr. Johna Rudzin's Computer Methods Course. 

My Final Project focuses on the March 31st through April 1st, 2023 tornado outbreak. I am honing in on two tornadoes from the bi-modal high risk event: the Little Rock, AR EF-3 and the Keota, IA EF-4.

To do this, I have six scripts that plot 300mb, 500mb, 850mb charts, SPC Convective Outlooks, SPC Storm Reports, Plan-Position Indicator (PPI) scans, WPC Surface Analysis, and Skew-Ts with calculated meteorological parameters. 

###
300mb_500mb_850mg_Charts.py

I used NOAA NCEI data archive to obtain Rapid Refresh (RAP) grib file for the event. The two data sets are for 18Z and 21Z. The link to grab the data is https://www.ncei.noaa.gov/data/rapid-refresh/access/rap-130-13km/analysis/202303/20230331/

The Python scripts are designed to analyze and visualize RAP (Rapid Refresh) grib data, focusing on atmospheric variables such as geopotential height, wind speed, and absolute vorticity at different pressure levels.

The script utilizes the PyGrib library to open RAP grib files and extract necessary atmospheric variables.

Various plots are generated for different pressure levels and timestamps, including 300mb heights and wind speeds, 500mb heights and absolute vorticity, and 850mb heights and wind speeds.

Each plot is projected onto a Lambert Conformal projection and features map elements such as coastlines, borders, states, and lakes.

Wind barbs and contour plots are used to visualize wind speed and geopotential height, respectively.

###
NOAASPCConvectiveOutlook.py
I used thge SPC Product and Report Archives site to obtain the necessary files, including shapefile, for the SPC Outlooks.

This python file contains Python scripts for plotting convective outlooks issued by the Storm Prediction Center (SPC) using Geopandas and MetPy libraries. Convective outlooks provide forecasts for severe weather potential across the United States.

The script reads shapefiles for SPC Day 1 convective outlooks valid at 1200Z and 1630Z. It uses Geopandas to read shapefiles and extract geometries, fills, strokes, and labels. MetPy is utilized for plotting the outlooks on map panels with specified projections.

###
PPI_Radar.py
I used AWS S3 Explorer archive to download radar data for KLZK and KDVN: https://s3.amazonaws.com/noaa-nexrad-level2/index.html.

This python file contains Python scripts for plotting Plan Position Indicators (PPIs) of radar reflectivity data using the Py-ART library. PPIs are commonly used in meteorology to display radar data in a two-dimensional format.

This script reads multiple radar data files from the KLZK and KDVN radar site and plots PPIs for each file.

###
SPCStormReports.py

This python file contains Python scripts for mapping severe weather reports (tornadoes, hail, and wind) from NOAA using geospatial data visualization tools.

The script retrieves tornado, hail, and wind reports CSV data from NOAA's Storm Prediction Center (SPC). It replaces unknown ratings or sizes with 0 and converts data to numeric format for consistent plotting. Using Cartopy and Matplotlib, it creates separate maps for each severe weather type (tornadoes, hail, and wind) with different colors for each type of report.

###
Soundings.py

This python file contains Python scripts for analyzing upper air soundings from Little Rock, Arkansas and Lincoln, IL at 1200Z and 0000Z. The scripts utilize data retrieval, processing, and plotting techniques to visualize atmospheric conditions, including thermodynamic and kinematic parameters.

The script imports necessary packages and modules for data analysis and plotting.
It retrieves upper air sounding data from the Wyoming Upper Air archive using the Siphon library.

The data is processed to handle missing values and converted into appropriate units for plotting.

The script generates Skew-T plots with hodographs for both 12Z and 00Z soundings from Little Rock, AR and Lincoln, IL. Various thermodynamic and kinematic parameters are calculated and overlaid on the plots for analysis.

###
WPC_Surface_Analysis.py

This python file contains Python scripts for plotting surface features on a map using data from the Weather Prediction Center (WPC) surface bulletin. The scripts utilize Cartopy for map projections and MetPy for parsing and plotting surface features.

The script imports necessary packages and modules for data analysis and plotting. It defines a function plot_bulletin to plot surface features on a map using Cartopy.

Surface analysis data is retrieved from WPC surface bulletins using MetPy's parse_wpc_surface_bulletin function.

The script generates maps with surface features plotted for different timestamps.
Features such as high and low-pressure systems, warm fronts, cold fronts, etc., are plotted with appropriate styling.
