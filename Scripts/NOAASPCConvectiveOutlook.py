# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 12:37:00 2024

@author: brees
"""
#Import necessary packages
import geopandas
from metpy.cbook import get_test_data
from metpy.plots import MapPanel, PanelContainer, PlotGeometry

#SPC Day 1 Outlook Valid 1200Z 
#Read shapefile
day1_outlook = geopandas.read_file('C:/Users/brees/Desktop/CompMethods/Final_Project/Final_Codes/Data/day1otlk_20230331_1200_cat.shp')

#Define plot geometry object
geo = PlotGeometry()

#Set geometry, fill, stroke, and labels for PlotGeometry
geo.geometry = day1_outlook['geometry']
geo.fill = day1_outlook['fill']
geo.stroke = day1_outlook['stroke']
geo.labels = day1_outlook['LABEL']
geo.label_fontsize = 'large'

#Define MapPanel object
panel = MapPanel()
panel.title = 'SPC Day 1 Convective Outlook (Valid 1200z Mar 31 2023)'
panel.plots = [geo]
panel.area = [-120, -75, 25, 50]
panel.projection = 'lcc'
panel.layers = ['lakes', 'land', 'ocean', 'states', 'coastline', 'borders']

#Define PanelContainer object
pc = PanelContainer()

#Set size and panels for the PanelContainer
pc.size = (12, 8)
pc.panels = [panel]

#Save and show plot
pc.save('SPCDay1ConvectiveOutlook1200z.png')
pc.show()


#SPC Day 1 Outlook Valild 1630Z
#Read shapefile
day1_outlook = geopandas.read_file('C:/Users/brees/Desktop/CompMethods/Final_Project/Final_Codes/Data/day1otlk_20230331_1630_cat.shp')

#Define plot geometry object
geo = PlotGeometry()

#Set geometry, fill, stroke, and labels for PlotGeometry
geo.geometry = day1_outlook['geometry']
geo.fill = day1_outlook['fill']
geo.stroke = day1_outlook['stroke']
geo.labels = day1_outlook['LABEL']
geo.label_fontsize = 'large'

#Define MapPanel object
panel = MapPanel()

#Set title and properties
panel.title = 'SPC Day 1 Convective Outlook (Valid 1630z Mar 31 2023)'
panel.plots = [geo]
panel.area = [-120, -75, 25, 50]
panel.projection = 'lcc'
panel.layers = ['lakes', 'land', 'ocean', 'states', 'coastline', 'borders']

#Define PanelContainer object
pc = PanelContainer()

#Set size and panels
pc.size = (12, 8)
pc.panels = [panel]

#Save and show plot
pc.save('SPCDay1ConvectiveOutlook1630z.png')
pc.show()
