# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:54:49 2024

@author: brees
"""


import numpy as np
import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import pygrib

#Open RAP grib file 
grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_1800_000.grb2')

#PLOT 1

############300mb Heights and Wind Speeds

#Extract 300mb Heights
gp300_msg=grbRAP.message(56)
gp_300=gp300_msg.values
[lats,lons]=gp300_msg.latlons()

#Extract 300mb wind magnitude (knots)
u = grbRAP.select(name='U component of wind')
v = grbRAP.select(name='V component of wind')
uwnd300 = grbRAP.message(60)
uwnd300_arr=uwnd300.values
vwnd300=grbRAP.message(61)
vwnd300_arr=vwnd300.values

#Calculate wind speed in knots
wndspd=np.sqrt(uwnd300_arr**2+vwnd300_arr**2)*1.94

#Set projection
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='wheat')
ax.add_feature(cf.OCEAN,color='lightsteelblue')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='gray')
ax.add_feature(cf.BORDERS,edgecolor='gray',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue', alpha=0.5)

#Define wind speed contour levels
wndspd_bounds=[30,40,50,60,80,100,125,150,200]

#Plot filled contours of wind speed
plt.contourf(lons,lats,wndspd, wndspd_bounds, cmap=plt.cm.hot_r,transform=ccrs.PlateCarree())
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label ('knots')

#Plot 300mb Height Contours
h=plt.contour(lons,lats,gp_300/10,np.arange(np.min(gp_300/10),np.max(gp_300/10),6),linestyles='-',linewidths=2,colors='black',transform=ccrs.PlateCarree())

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='white', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Add wind barbs
plt.barbs(lons[::48,::48],lats[::48,::48],uwnd300_arr[::48,::48],vwnd300_arr[::48,::48],transform=ccrs.PlateCarree())

#Set title for plot
plt.title ('300 mb Heights (dm) / Isotachs (knots)')
plt.savefig('300mb_Heights_Winds_18Z.png')
plt.show()

#Open RAP grib file 
grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_2100_000.grb2')

#PLOT 2

############300mb Heights and Wind Speeds

#Extract 300mb Heights
gp300_msg=grbRAP.message(56)
gp_300=gp300_msg.values
[lats,lons]=gp300_msg.latlons()

#Extract 300mb wind magnitude (knots)
u = grbRAP.select(name='U component of wind')
v = grbRAP.select(name='V component of wind')
uwnd300 = grbRAP.message(60)
uwnd300_arr=uwnd300.values
vwnd300=grbRAP.message(61)
vwnd300_arr=vwnd300.values

#Calculate wind speed in knots
wndspd=np.sqrt(uwnd300_arr**2+vwnd300_arr**2)*1.94

#Set projection
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='wheat')
ax.add_feature(cf.OCEAN,color='lightsteelblue')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='gray')
ax.add_feature(cf.BORDERS,edgecolor='gray',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue', alpha=0.5)

#Define wind speed contour levels
wndspd_bounds=[30,40,50,60,80,100,125,150,200]

#Plot filled contours of wind speed
plt.contourf(lons,lats,wndspd, wndspd_bounds, cmap=plt.cm.hot_r,transform=ccrs.PlateCarree())
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label ('knots')

#Plot 300mb Height Contours
h=plt.contour(lons,lats,gp_300/10,np.arange(np.min(gp_300/10),np.max(gp_300/10),6),linestyles='-',linewidths=2,colors='black',transform=ccrs.PlateCarree())

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='white', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Add wind barbs
plt.barbs(lons[::48,::48],lats[::48,::48],uwnd300_arr[::48,::48],vwnd300_arr[::48,::48],transform=ccrs.PlateCarree())

#Set title for plot
plt.title ('300 mb Heights (dm) / Isotachs (knots)')
plt.savefig('300mb_Heights_Winds_21Z.png')
plt.show()

#PLOT 3 -- 500mb Plots

grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_1800_000.grb2')

#Extract 500mb geopotential height
gpRAP=grbRAP.select(name='Geopotential height')
gp_msg=grbRAP.message(104)
[lats,lons]=gp_msg.latlons()
gp500=gp_msg.values

#Extract 500mb absolute vorticity 
vortHRRR=grbRAP.select(name='Absolute vorticity')
vort_msg=grbRAP.message(110)
vort500=vort_msg.values

#Set up projection and plot
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='wheat')
ax.add_feature(cf.OCEAN,color='lightsteelblue')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='black')
ax.add_feature(cf.BORDERS,edgecolor='black',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue', alpha=0.5)

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='white', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Plot 500mb geopotential height contours
h=plt.contour(lons,lats,gp500/10,np.arange(np.min(gp500/10),np.max(gp500/10),6),linestyles='-',linewidths=2,colors='black',transform=ccrs.PlateCarree())

#Define and plot contour levels for absolute vorticity
bounds1 = [-6, 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60]
plt.contourf(lons,lats,vort500*100000,bounds1,cmap=plt.cm.YlOrRd,transform=ccrs.PlateCarree())

#Add colorbar for absolute vorticity
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label ('10^-5 s^-1')

#Set title for plot
plt.title ('500 mb Heights (dm) / Absolute Vorticity (10^-5)')

plt.savefig('500mb_Heights_Vorticity_18Z.png')
plt.show()

#PLOT 4

#Open RAP grib file 
grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_2100_000.grb2')

#Extract 500mb geopotential height
gpRAP=grbRAP.select(name='Geopotential height')
gp_msg=grbRAP.message(104)
[lats,lons]=gp_msg.latlons()
gp500=gp_msg.values

#Extract 500mb absolute vorticity
vortHRRR=grbRAP.select(name='Absolute vorticity')
vort_msg=grbRAP.message(110)
vort500=vort_msg.values

#500mb absolute vorticity 
vortHRRR=grbRAP.select(name='Absolute vorticity')
vort_msg=grbRAP.message(110)
vort500=vort_msg.values

#Set up projection and plot
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='wheat')
ax.add_feature(cf.OCEAN,color='lightsteelblue')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='black')
ax.add_feature(cf.BORDERS,edgecolor='black',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue', alpha=0.5)

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='white', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Plot 500mb geopotential height contours
h=plt.contour(lons,lats,gp500/10,np.arange(np.min(gp500/10),np.max(gp500/10),6),linestyles='-',linewidths=2,colors='black',transform=ccrs.PlateCarree())

#Define and plot contour levels for absolute vorticity
bounds1 = [-6, 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60]
plt.contourf(lons,lats,vort500*100000,bounds1,cmap=plt.cm.YlOrRd,transform=ccrs.PlateCarree())

#Add colorbar for absolute vorticity
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label ('10^-5 s^-1')

#Set title for plot
plt.title ('500 mb Heights (dm) / Absolute Vorticity (10^-5)')
plt.savefig('500mb_Heights_Vorticity_21Z.png')
plt.show()

###########850mb Heights and Wind Speeds

#Plot 5 - 850mb

grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_1800_000.grb2')

#Extract geopotential/wind magnitude at 850mb 
gp_850=grbRAP.select(name='Geopotential height', level=850)[0]
gp_msg_850=grbRAP.message(189)
gp_data=gp_850.values
gp_units=gp_msg_850.units

uwnd_850=grbRAP.select(name='U component of wind', level=850)[0]
uwnd_data=uwnd_850.values
vwnd_850=grbRAP.select(name='V component of wind', level=850)[0]
vwnd_data=vwnd_850.values

#Extract latitude and longitude
[lats,lons]=gp_msg_850.latlons()

#Set up plotting 
wndspd_850 = np.sqrt(uwnd_data**2 + vwnd_data**2) * 1.94384
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='lightgrey')
ax.add_feature(cf.OCEAN,color='white')
ax.add_feature(cf.COASTLINE,edgecolor='grey')
ax.add_feature(cf.STATES,edgecolor='grey')
ax.add_feature(cf.BORDERS,edgecolor='grey',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue')

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              	linewidth=2, color='white', linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Set wind speed contour bounds
bounds2 = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

#Plot geopotential height contour
gp_contour = ax.contour(lons, lats, gp_data, colors='black', transform=ccrs.PlateCarree())
plt.clabel(gp_contour, inline=1, fontsize=12, fmt='%1.0f')

#Plot wind speed contours
wndspd_contour_850 = ax.contourf(lons, lats, wndspd_850, bounds2, cmap=plt.cm.BuPu, transform=ccrs.PlateCarree())
cbar = plt.colorbar(wndspd_contour_850, orientation='horizontal')
cbar.set_label('(knots)')

#Add wind barbs
ax.barbs(lons[::25, ::25], lats[::25, ::25], uwnd_data[::25, ::25], vwnd_data[::25, ::25], transform=ccrs.PlateCarree())

#Set title for plot
plt.title('850 mb Heights (dm) / Wind Speed (knots)')
plt.savefig('850mb_Height_Winds_18Z.png')
plt.show()

#PLOT 6

grbRAP=pygrib.open('C:/Users/brees/Desktop/CompMethods/Final_Project/SFC_UPPERAIR/rap_130_20230331_2100_000.grb2')

#Extract geopotential/wind magnitude at 850mb 
gp_850=grbRAP.select(name='Geopotential height', level=850)[0]
gp_msg_850=grbRAP.message(189)
gp_data=gp_850.values
gp_units=gp_msg_850.units

uwnd_850=grbRAP.select(name='U component of wind', level=850)[0]
uwnd_data=uwnd_850.values
vwnd_850=grbRAP.select(name='V component of wind', level=850)[0]
vwnd_data=vwnd_850.values

#Extract latitude and longitude
[lats,lons]=gp_msg_850.latlons()

#Set up plotting 
wndspd_850 = np.sqrt(uwnd_data**2 + vwnd_data**2) * 1.94384
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
fig=plt.figure(figsize=(8,8))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add map features
ax.add_feature(cf.LAND,color='lightgrey')
ax.add_feature(cf.OCEAN,color='white')
ax.add_feature(cf.COASTLINE,edgecolor='grey')
ax.add_feature(cf.STATES,edgecolor='grey')
ax.add_feature(cf.BORDERS,edgecolor='grey',linestyle='-')
ax.add_feature(cf.LAKES,edgecolor='lightsteelblue')

#Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              	linewidth=2, color='white', linestyle='--')
gl.top_labels = False
gl.left_labels = False

#Set wind speed contour bounds
bounds2 = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

#Plot geopotential height contour
gp_contour = ax.contour(lons, lats, gp_data, colors='black', transform=ccrs.PlateCarree())
plt.clabel(gp_contour, inline=1, fontsize=12, fmt='%1.0f')

#Plot wind speed contours
wndspd_contour_850 = ax.contourf(lons, lats, wndspd_850, bounds2, cmap=plt.cm.BuPu, transform=ccrs.PlateCarree())
cbar = plt.colorbar(wndspd_contour_850, orientation='horizontal')
cbar.set_label('(knots)')

#Add wind barbs
ax.barbs(lons[::25, ::25], lats[::25, ::25], uwnd_data[::25, ::25], vwnd_data[::25, ::25], transform=ccrs.PlateCarree())

#Set title for plot
plt.title('850 mb Heights (dm) / Wind Speed (knots)')
plt.savefig('850mb_Height_Winds_21Z.png')
plt.show()


