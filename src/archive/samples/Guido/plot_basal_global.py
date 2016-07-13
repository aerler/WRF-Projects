from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import rcParams
import numpy.ma as ma
import matplotlib
from matplotlib.colors import colorConverter
import sys

#nchf = open("./hfmap.asc", 'r')
#ncyd = Dataset("./T85ydts_dec.nc", 'r')
nchf = Dataset("./hfmap_global_shift.nc", 'r')

# partition data into (lon, lat, v1 ,v2) array
#Format: lon  lat  HF (mW/m**2)  sigmaHF(mW/m**2)
#data = np.array(nchf.read().split())
lat = nchf.variables["lat"][:]
lon = nchf.variables["lon"][:]
hf = nchf.variables["bheatflx"][:,:]
shf = nchf.variables["sdevbheatflx"][:,:]
 
#print data[0][0],data[0][1],data[0][2],data[0][3]
#print data.shape
#lon = np.array(range(0,360,1),dtype="float")
#lat = np.array(range(-90,91,1),dtype="float")
#print lon,lat

#flon, flat, thf, tshf = data[0::4].reshape(360,181),data[1::4].reshape(360,181),data[2::4].reshape(360,181),data[3::4].reshape(360,181)

#flon, flat, thf, tshf = data[0::4].reshape(181,360),data[1::4].reshape(181,360),data[2::4].reshape(181,360),data[3::4].reshape(181,360)

#hf=np.transpose(thf)
#shf=np.transpose(tshf)
#hf=data
print hf.shape
#print shf.shape
#print lon.shape
#print lat.shape

#print np.min(hf),np.max(hf)
#print min(shf),max(shf)



nchf.close()

#add cyclic
#hflon = np.zeros(lon.shape[0]+1,np.float64)
#hflon[0:-1] = lon[:]; hflon[-1] = lon[0]+360.

# YD (20th Decade)
#ydsttmp = ncyd.variables['TS'][20,:,:]-zkelvin
#ydstlonm1 = ncyd.variables['lon'][:]
#ydstlat = ncyd.variables['lat'][:]

#ydstavg=np.zeros((ydsttmp.shape[0],ydsttmp.shape[1]+1),np.float64)
#ydstavg[:,0:-1] = ydsttmp; ydstavg[:,-1] = ydsttmp[:,0]
#ydstlon = np.zeros(ydstlonm1.shape[0]+1,np.float64)
#ydstlon[0:-1] = ydstlonm1[:]; ydstlon[-1] = ydstlonm1[0]+360.


clon=0.0
clat=0.0
lon_0 = -90.
lat_0 = 90.
bounding_lat = 0.
projection="merc"
ncont=10

#label sublot character size and position (e.g. a,b,c,d)
lbfs=20
lbpx=0.04
lbpy=0.96

#fig = plt.figure(facecolor='white', figsize = (8.5, 11))
fig = plt.figure(facecolor='white', figsize = (11., 8.5))

#baroffset=0.023
baroffset=0.0

# --------- plots
plt.subplot(211)

#ma = Basemap(boundinglat=bounding_lat,lon_0=lon_0,\
#            resolution='l',area_thresh=10000.,projection=projection)

ma = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')

xa, ya = ma(*np.meshgrid(lon, lat))

ma.drawcoastlines()
ma.drawparallels(np.arange(-90.,90.,30.))    
ma.drawmeridians(np.arange(-180.,180.,30.))
#revhaxby=cmap_map(lambda x: x[::-1], cm.GMT_haxby)

#tlvs=np.arange(0.,160.,10.)
tlvs=np.arange(-160.,0.,10.)

#tlvs=(-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30)

#cs1 = m1.contourf(x1, y1, obsnaug, lvs, cmap=cm.GMT_haxby_r)
csa = ma.contourf(xa, ya, hf, tlvs, cmap=cm.GMT_no_green_r)
#csa = ma.contourf(xa, ya, obstaugavg, tlvs, cmap=cm.GMT_globe)
ma.drawlsmask(land_color='lightgrey',ocean_color='steelblue',lakes=False)

# draw colorbar.
axa = plt.gca()
pos = axa.get_position()
l, b, w, h = pos.bounds
#caxa = plt.axes([l+w+0.01, b+0.18, 0.02, h-0.35]) # setup colorbar axes
caxa = plt.axes([l+w+0.01, b, 0.02, h]) # setup colorbar axes

plt.colorbar(csa, caxa,ticks=csa.levels[::2], format='%g', drawedges=False) 

plt.axes(axa)  # make the original axes current again
#plt.axes(rect1)
plt.title("Heat Flux (mW/m**2)", horizontalalignment='left', position=(0.0,1.0))


#make subfigure label for publication (a,b,c,d etc)
axa.text(lbpx, lbpy, "a", ha='left', va='top', size=lbfs,
         bbox=dict(boxstyle="round", fc="w", ec="k"), transform=axa.transAxes)

# ---- model ground temperature
plt.subplot(212)

#mb = Basemap(boundinglat=bounding_lat,lon_0=lon_0,\
#            resolution='l',area_thresh=10000.,projection=projection)

mb = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')


xb, yb = mb(*np.meshgrid(lon, lat))
mb.drawcoastlines()
mb.drawparallels(np.arange(-90.,90.,30.))    
mb.drawmeridians(np.arange(-180.,180.,30.))

dtlvs=np.arange(0,200,20.)
#csb = mb.contourf(xb,yb,ydstavg-obstavg,dtlvs,cmap=cm.GMT_no_green)
csb = mb.contourf(xb,yb,shf,dtlvs,cmap=cm.sstanom)
mb.drawlsmask(land_color='lightgrey',ocean_color='steelblue',lakes=False)
#czero = mb.contour(xb,yb, shf, (0,),
#                colors = 'black',
#                linewidths = 2,
#                hold='on')

# draw colorbar.
axb =plt.gca()
pos = axb.get_position()
l, b, w, h = pos.bounds
#caxb = plt.axes([l+w+0.01, b+0.18, 0.02, h-0.35]) # setup colorbar axes
caxb = plt.axes([l+w+0.01, b, 0.02, h]) # setup colorbar axes
plt.colorbar(csb, caxb,ticks=csb.levels[::2], format='%g', drawedges=False) 

plt.axes(axb)  # make the original axes current again
#plt.axes(rect2)
plt.title('$\sigma$(Heat Flux) (mW/m**2)', horizontalalignment='left', position=(0.0,1.0))

#make subfigure label for publication (a,b,c,d etc)
axb.text(lbpx, lbpy, "b", ha='left', va='top', size=lbfs,
         bbox=dict(boxstyle="round", fc="w", ec="k"), transform=axb.transAxes)



#fig.subplots_adjust(hspace=0.2)
#fig.subplots_adjust(wspace=-0.05)
#plt.figtext(0.5, 0.925,  'August', 
#               ha='center', color='black', weight='bold', size='large')        


plt.show()
