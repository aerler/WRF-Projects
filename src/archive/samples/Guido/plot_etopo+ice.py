from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm as bcm
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import rcParams
import numpy.ma as ma
import matplotlib
from matplotlib.colors import colorConverter
import sys

ftopo = Dataset('../topo10min.nc','r')
#fitopo = Dataset('../pct_glacier_yr1850_0.23x0.31.nc','r')
fitopo = Dataset('../pct_glacier_yr2000_0.23x0.31.nc','r')

#topoin = ftopo.variables['htopo'] 
#lons = ftopo.variables['lon']
#lats = ftopo.variables['lat']
#print topoin.shape, lons.shape, lats.shape
#print lons[:]

etopotmp = ftopo.variables['htopo'][:,:]
lonm1 = ftopo.variables['lon'][:]
lat = ftopo.variables['lat'][:]
etopo=np.zeros((etopotmp.shape[0],etopotmp.shape[1]+1),np.float64)
etopo[:,0:-1] = etopotmp; etopo[:,-1] = etopotmp[:,0]
lon = np.zeros(lonm1.shape[0]+1,np.float64)
lon[0:-1] = lonm1[:]; lon[-1] = lonm1[0]+360.


itopotmp = fitopo.variables['PCT_GLACIER'][:,:]
ilonm1 = fitopo.variables['LONGXY'][0,:]
ilat = fitopo.variables['LATIXY'][:,0]
itopo=np.zeros((itopotmp.shape[0],itopotmp.shape[1]+1),np.float64)
itopo[:,0:-1] = itopotmp; itopo[:,-1] = itopotmp[:,0]
ilon = np.zeros(ilonm1.shape[0]+1,np.float64)
ilon[0:-1] = ilonm1[:]; ilon[-1] = ilonm1[0]+360.

#print ilat
#print ilon

clon=-30.0
clat=90.0
lon_0 = -90.
lat_0 = 90.
bounding_lat = 60.
projection="npaeqd"
ncont=10

#label sublot character size and position (e.g. a,b,c,d)
lbfs=20
lbpx=0.04
lbpy=0.96

#fig = plt.figure(facecolor='white', figsize = (8.5, 11))
fig = plt.figure(facecolor='white', figsize = (11., 8.5))

baroffset=0.023


# --------- plots
plt.subplot(111)

ma = Basemap(boundinglat=bounding_lat,lon_0=lon_0,\
            area_thresh=10000.,projection=projection)
xa, ya = ma(*np.meshgrid(lon, lat))

ma.drawcoastlines()
ma.drawparallels(np.arange(-90.,90.,30.))    
ma.drawmeridians(np.arange(-180.,180.,30.))
#revhaxby=cmap_map(lambda x: x[::-1], cm.GMT_haxby)

#tlvs=np.arange(-20.,22.,1.)
#tlvs=(-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
#tlvs=(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
#tlvs=np.arange(0,5000,100)
tlvs=np.arange(0,2000,100)

#cs1 = m1.contourf(x1, y1, obsnaug, lvs, cmap=bcm.GMT_haxby_r)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_no_green)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_globe)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_haxby_r)
csa = ma.contourf(xa, ya, etopo, tlvs, cmap=cm.RdYlGn_r)


ma.drawlsmask(land_color='lightgrey',ocean_color='steelblue',lakes=False)

# draw colorbar.
axa = plt.gca()
#pos = axa.get_position()
#l, b, w, h = pos.bounds
#caxa = plt.axes([l+w+0.01+baroffset, b, 0.02, h]) # setup colorbar axes
#plt.colorbar(csa, caxa,ticks=csa.levels[::2], format='%g', drawedges=False) 

#plt.axes(axa)  # make the original axes current again
#plt.axes(rect1)
plt.title("ETOPO 10 min (meters)", horizontalalignment='left', position=(0.0,1.0))


#make subfigure label for publication (a,b,c,d etc)
axa.text(lbpx, lbpy, "a", ha='left', va='top', size=lbfs,
         bbox=dict(boxstyle="round", fc="w", ec="k"), transform=axa.transAxes)


#plt.subplot(212)
#plt.hold()


mai = Basemap(boundinglat=bounding_lat,lon_0=lon_0,\
            area_thresh=10000.,projection=projection)
xai, yai = mai(*np.meshgrid(ilon, ilat))

#mai.drawcoastlines()
#mai.drawparallels(np.arange(-90.,90.,30.))    
#mai.drawmeridians(np.arange(-180.,180.,30.))
#revhaxby=cmap_map(lambda x: x[::-1], cm.GMT_haxby)

#tlvs=np.arange(-20.,22.,1.)
#tlvs=(-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
#tlvs=(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
#tlvs=np.arange(0,5000,100)
tlvsi=np.arange(0,105,5)

#cs1 = m1.contourf(x1, y1, obsnaug, lvs, cmap=bcm.GMT_haxby_r)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_no_green)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_globe)
#csa = ma.contourf(xa, ya, etopo, tlvs, cmap=bcm.GMT_haxby_r)
#csai = mai.contourf(xai, yai, itopo, tlvsi, cmap=cm.PuBu_r)
#csai = mai.contourf(xai, yai, itopo, tlvsi, cmap=cm.jet)
csai = mai.contourf(xai, yai, itopo, tlvsi, cmap=cm.bone, alpha = 1.0)


#mai.drawlsmask(land_color='lightgrey',ocean_color='steelblue',lakes=False)

# draw colorbar.
axai = plt.gca()
posi = axai.get_position()
li, bi, wi, hi = posi.bounds
caxai = plt.axes([li+wi-0.05+baroffset, bi, 0.02, hi]) # setup colorbar axes
plt.colorbar(csai, caxai,ticks=csai.levels[::2], format='%g', drawedges=False) 

plt.axes(axai)  # make the original axes current again
plt.title("Percent Glacier (%)", horizontalalignment='left', position=(0.0,1.0))


#make subfigure label for publication (a,b,c,d etc)
#axai.text(lbpx, lbpy, "a", ha='left', va='top', size=lbfs,
#         bbox=dict(boxstyle="round", fc="w", ec="k"), transform=axa.transAxes)



fig.subplots_adjust(hspace=0.2)
fig.subplots_adjust(wspace=-0.05)
#plt.figtext(0.5, 0.925,  'August', 
#               ha='center', color='black', weight='bold', size='large')        

ftopo.close()

plt.show()
