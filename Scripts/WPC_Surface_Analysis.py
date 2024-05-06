# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 14:32:45 2024

@author: brees
"""


#Import necessary packages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from metpy.io import parse_wpc_surface_bulletin
from metpy.plots import (ColdFront, OccludedFront, StationaryFront,
                         StationPlot, WarmFront)

#Function to plot surface features on a map
def plot_bulletin(ax, data):
    """Plot a dataframe of surface features on a map."""
    size = 4
    fontsize = 9
    complete_style = {'HIGH': {'color': 'blue', 'fontsize': fontsize},
                      'LOW': {'color': 'red', 'fontsize': fontsize},
                      'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
                      'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
                      'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
                      'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
                      'TROF': {'linewidth': 2, 'linestyle': 'dashed', 'edgecolor': 'darkorange'}}
    
    for field in ('HIGH', 'LOW'):
        rows = data[data.feature == field]
        x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
        sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
        sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
        sp.plot_parameter('S', rows.strength, **complete_style[field])

    for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
        rows = data[data.feature == field]
        ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field], facecolor='none')

#Create figure and axis for plotting
fig = plt.figure(figsize=(10, 6), dpi=150)  
ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=-100))

#Set the geographical extent
ax.set_extent([-125, -66.5, 24.5, 49.5], crs=ccrs.PlateCarree())  

#Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAKES)

#Update the file path to the file
file_path = 'CODSUS_202303311331.txt'  
df = parse_wpc_surface_bulletin(file_path)
plot_bulletin(ax, df)

#Set title
ax.set_title(f'WPC Surface Analysis Valid {df.valid.dt.strftime("%HZ %d %b %Y")[0]}')

#Save figure and plot
plt.savefig('WPCSurfaceAnalysis_12Z.png')
plt.show()

#Import necessary packages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from metpy.io import parse_wpc_surface_bulletin
from metpy.plots import (ColdFront, OccludedFront, StationaryFront,
                         StationPlot, WarmFront)
#Function to plot surface features on a map
def plot_bulletin(ax, data):
    """Plot a dataframe of surface features on a map."""
    size = 4
    fontsize = 9
    complete_style = {'HIGH': {'color': 'blue', 'fontsize': fontsize},
                      'LOW': {'color': 'red', 'fontsize': fontsize},
                      'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
                      'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
                      'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
                      'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
                      'TROF': {'linewidth': 2, 'linestyle': 'dashed', 'edgecolor': 'darkorange'}}
    
    for field in ('HIGH', 'LOW'):
        rows = data[data.feature == field]
        x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
        sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
        sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
        sp.plot_parameter('S', rows.strength, **complete_style[field])

    for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
        rows = data[data.feature == field]
        ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field], facecolor='none')

#Create figure and axis for plotting
fig = plt.figure(figsize=(10, 6), dpi=150)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=-100))

#Set geographical extent
ax.set_extent([-125, -66.5, 24.5, 49.5], crs=ccrs.PlateCarree())

#Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAKES)

#Update file path
file_path = 'CODSUS_202303311626.txt'
df = parse_wpc_surface_bulletin(file_path)
plot_bulletin(ax, df)

#Set title
ax.set_title(f'WPC Surface Analysis Valid {df.valid.dt.strftime("%HZ %d %b %Y")[0]}')

#Save figure and show
plt.savefig('WPCSurfaceAnalysis_15Z.png')
plt.show()

#Function to plor surface features on map
def plot_bulletin(ax, data):
    """Plot a dataframe of surface features on a map."""
    size = 4
    fontsize = 9
    complete_style = {'HIGH': {'color': 'blue', 'fontsize': fontsize},
                      'LOW': {'color': 'red', 'fontsize': fontsize},
                      'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
                      'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
                      'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
                      'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
                      'TROF': {'linewidth': 2, 'linestyle': 'dashed', 'edgecolor': 'darkorange'}}
    
    for field in ('HIGH', 'LOW'):
        rows = data[data.feature == field]
        x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
        sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
        sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
        sp.plot_parameter('S', rows.strength, **complete_style[field])

    for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
        rows = data[data.feature == field]
        ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field], facecolor='none')

#Create figure and axis for plotting
fig = plt.figure(figsize=(10, 6), dpi=150)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=-100))

#Set geographical extent
ax.set_extent([-125, -66.5, 24.5, 49.5], crs=ccrs.PlateCarree())

#Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAKES)

#Update file path
file_path = 'CODSUS_202303311925.txt'
df = parse_wpc_surface_bulletin(file_path)
plot_bulletin(ax, df)

#Set title
ax.set_title(f'WPC Surface Analysis Valid {df.valid.dt.strftime("%HZ %d %b %Y")[0]}')

#Save figure and show
plt.savefig('WPCSurfaceAnalysis_18Z.png')
plt.show()

#Function to plot surface features on a map
def plot_bulletin(ax, data):
    """Plot a dataframe of surface features on a map."""
    size = 4
    fontsize = 9
    complete_style = {'HIGH': {'color': 'blue', 'fontsize': fontsize},
                      'LOW': {'color': 'red', 'fontsize': fontsize},
                      'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
                      'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
                      'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
                      'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
                      'TROF': {'linewidth': 2, 'linestyle': 'dashed', 'edgecolor': 'darkorange'}}
    
    for field in ('HIGH', 'LOW'):
        rows = data[data.feature == field]
        x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
        sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
        sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
        sp.plot_parameter('S', rows.strength, **complete_style[field])

    for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
        rows = data[data.feature == field]
        ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field], facecolor='none')

#Create figure and axis for plotting
fig = plt.figure(figsize=(10, 6), dpi=150)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=-100))

#Set geographical extent
ax.set_extent([-125, -66.5, 24.5, 49.5], crs=ccrs.PlateCarree())

#Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAKES)

#Update file path
file_path = 'CODSUS_202303312224.txt'
df = parse_wpc_surface_bulletin(file_path)
plot_bulletin(ax, df)

#Set title
ax.set_title(f'WPC Surface Analysis Valid {df.valid.dt.strftime("%HZ %d %b %Y")[0]}')

#Save figure and plot
plt.savefig('WPCSurfaceAnalysis_21Z.png')
plt.show()
