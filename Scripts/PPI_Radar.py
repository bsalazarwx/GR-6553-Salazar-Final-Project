# -*- coding: utf-8 -*-
"""
Created on Sun May  5 23:09:50 2024

@author: brees
"""

#Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyart

#Looping over multiple radar files to plot PPIs
for time in ['192215', '192732', '193250', '193806', '194310', '194827', '195326', '195852']:
    #File path for each radar data file
    filename = f'KLZK20230331_{time}_V06'
    
    #Reading radar data from the file
    radar3 = pyart.io.read(filename)
    
    #Creating radar display object
    displayc = pyart.graph.RadarDisplay(radar3)
    
    #Creating a new figure
    fig = plt.figure()
    
    #Setting axis limits
    ax = plt.axis([-20, 20, -20, 20])
    
    #Plotting PPI for reflectivity
    displayc.plot_ppi('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef',
                      title=f'KLZK 0.5\u00B0 2023-03-31 {time[:-2]}:{time[-2:]}Z\nEquivalent Reflectivity Factor',
                      axislabels=('N-S Distance from Radar (km)', 'E-W Distance from Radar (km)'),
                      colorbar_label='Equivalent Reflectivity Factor (dBZ)')
    #Save figure
    plt.savefig(f'KLZK_PPI_{time}.png')


#Importing necessary packages
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyart

#Plotting PPIs for multiple radar scans
#Looping over each radar file and plotting PPI for each
for time in ['210603', '211218', '211833', '212500', '213236', '213905', '214533', '215201', '215829']:
    # File path for each radar data file
    filename2 = f'KDVN20230331_{time}_V06'
    
    #Reading radar data from the file
    radar3 = pyart.io.read(filename2)
    
    #Creating radar display object
    displayc = pyart.graph.RadarDisplay(radar3)
    
    #Creating a new figure
    fig = plt.figure()
    
    #Setting axis limits
    ax = plt.axis([-150, 5, -40, 60])
    
    #Plotting PPI for reflectivity
    displayc.plot_ppi('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef',
                      title=f'KDVN 0.5\u00B0 2023-03-31 {time[:2]}:{time[2:]}Z\nEquivalent Reflectivity Factor',
                      axislabels=('N-S Distance from Radar (km)', 'E-W Distance from Radar (km)'),
                      colorbar_label='Equivalent Reflectivity Factor (dBZ)')
    
    #Save figure
    plt.savefig(f'KDVN_PPI_{time}.png')
