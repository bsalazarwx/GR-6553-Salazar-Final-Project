# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:01:48 2024

@author: brees
"""

#Import necessary libraries
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import pandas as pd
import geopandas
import matplotlib.pyplot as plt
from metpy.plots import USCOUNTIES

#Geocode the location 'Little Rock, AR'
my_town = geopandas.tools.geocode('Little Rock, AR')

#Extract the first point (latitude and longitude) from geocoded location
my_town = my_town['geometry'][0]

#Print the latitude and longitude of the location
print(my_town)

#Read tornado reports CSV data from NOAA
tornado_reports = pd.read_csv('https://www.spc.noaa.gov/climo/reports/230331_rpts_torn.csv')

#Replace unknown F_Scale ratings with 0
tornado_reports['F_Scale'].replace(to_replace='UNK', value=0, inplace=True)
tornado_reports['F_Scale'] = pd.to_numeric(tornado_reports['F_Scale'], errors='coerce').fillna(0)

#Define Lambert Conformal projection
proj = ccrs.LambertConformal(central_longitude=-95, central_latitude=35, standard_parallels=[35])

#Create a figure and subplot for plotting tornado reports
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,1,1,projection=proj)

#Add map features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)

#Set the geographical extent
ax.set_extent([-98, -83, 30, 46])

#Plot tornado reports
ax.scatter(tornado_reports['Lon'], tornado_reports['Lat'],
           color='red',  # Set the color to red
           transform=ccrs.PlateCarree(),
           s=20)

plt.savefig('SPCTornadoReports.png')

#Read hail reports CSV data from NOAA
hail_reports = pd.read_csv('https://www.spc.noaa.gov/climo/reports/230331_rpts_hail.csv')

#Replace unknown Size ratings with 0
hail_reports['Size'].replace(to_replace='UNK', value=0, inplace=True)
hail_reports['Size'] = pd.to_numeric(hail_reports['Size'], errors='coerce').fillna(0)

#Create a figure and subplot for plotting hail reports
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,1,1,projection=proj)

#Add map features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)

#Set the geographical extent
ax.set_extent([-98, -83, 30, 46])

#Plot hail reports
ax.scatter(hail_reports['Lon'], hail_reports['Lat'],
           color='green',  # Set the color to green
           transform=ccrs.PlateCarree(),
           s=20)

plt.savefig('SPCHailReports.png')

#Read wind reports CSV data from NOAA
wind_reports = pd.read_csv('https://www.spc.noaa.gov/climo/reports/230331_rpts_filtered_wind.csv')

#Replace unknown Speed ratings with 0
wind_reports['Speed'].replace(to_replace='UNK', value=0, inplace=True)
wind_reports['Speed'] = pd.to_numeric(wind_reports['Speed'], errors='coerce').fillna(0)

#Create a figure and subplot for plotting wind reports
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,1,1,projection=proj)

#Add map features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)

#Set the geographical extent
ax.set_extent([-98, -83, 30, 46])

#Plot wind reports
ax.scatter(wind_reports['Lon'], wind_reports['Lat'],
           color='blue',  # Set the color to blue
           transform=ccrs.PlateCarree(),
           s=20)

plt.savefig('SPCWindReports.png')
