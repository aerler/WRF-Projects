'''
Created on 2012-09-29

A simple script that reads a WRF netcdf-4 file and displays a 2D field in a proper geographic projection;
primary application is plotting the domain with topography.

@author: Andre R. Erler
'''

## includes
# matplotlib config: size etc.
import numpy as np
import matplotlib.pylab as pyl
import matplotlib as mpl
mpl.rc('lines', linewidth=1.)
mpl.rc('font', size=10)
# pygeode stuff
from myDatasets.loadWRF import openWRF
from myPlots.plots import surfacePlot
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm, maskoceans

#from pygeode.plot import plot_v1 as pl
#from pygeode.plot import basemap as bm

## settings
sf = dict(dpi=150) # print properties
folder = '/home/me/Research/Dynamical Downscaling/figures/' # figure directory


if __name__ == '__main__':
  
  ## read data
  data1, data2 = openWRF('ctrl-1',[1981],[1,2,12])
  print(data1)
  ## compute data
#  precip1 = data1.rain(time=) 
#  precip2 = 
  
  ## setup projection
  f = pyl.figure(facecolor='white', figsize = (6.25,4.25))
  ax = f.add_subplot(1,1,1)
  f.subplots_adjust(bottom=0.05, left=0.05, right=.93, top=.93) # hspace, wspace
  # setup lambert conformal basemap.
  # lat_1 is first standard parallel.
  # lat_2 is second standard parallel (defaults to lat_1).
  # lon_0,lat_0 is central point.
  # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
  # area_thresh=1000 means don't plot coastline features less
  # than 1000 km^2 in area.
  lcc = dict(projection='lcc', lat_0=57, lon_0=-137, lat_1=53, rsphere=(6378137.00,6356752.3142),#
              width=280*30e3, height=210*30e3, area_thresh = 3000., resolution='l')
  map = Basemap(ax=ax,**lcc)
  
  ## Plot data
  grid = 1.25; res = 'l'; 
  clev = np.linspace(0,3500,36); cmap = mpl.cm.gist_earth # s3pcpn(_l)
  # coordinates
  lat1 = data1.lat.get(); lat2 = data2.lat.get() 
  lon1 = data1.lon.get(); lon2 = data2.lon.get()
  x1, y1 = list(map(lon1,lat1)) # convert to map-native coordinates
  x2, y2 = list(map(lon2,lat2)) # convert second domain to first domain map coordinates
  # draw boundaries of inner and outer domains
  bdy1 = np.ones_like(lat1); bdy1[0,:]=0; bdy1[-1,:]=0; bdy1[:,0]=0; bdy1[:,-1]=0
  bdy2 = np.ones_like(lat2); bdy2[0,:]=0; bdy2[-1,:]=0; bdy2[:,0]=0; bdy2[:,-1]=0
  map.contour(x1,y1,bdy1,[0],ax=ax, colors='k') # draw boundary of outer domain
  map.contour(x2,y2,bdy2,[0],ax=ax, colors='k') # draw boundary of inner domain
  # data pre-processing: mask out ocean
  zs1 = maskoceans(lon1,lat1,data1.zs.get(),resolution=res,grid=grid) # coarse terrain
  zs2 = maskoceans(lon2,lat2,data2.zs.get(),resolution=res,grid=grid) # hires terrain
  # draw data
  cmap.set_under('blue')
  cmap.set_over('white')
  c1 = map.contourf(x1,y1,zs1,clev,ax=ax,cmap=cmap,extend='both')
  c2 = map.contourf(x2,y2,zs2,clev,ax=ax,cmap=cmap,extend='both')
  
  # add colorbar
  c1.set_clim(vmin=-1500,vmax=3500); c2.set_clim(vmin=-1500,vmax=3500)
  cbar = map.colorbar(mappable=c1,location='right',size='3%',pad='2%')
  cbl = np.linspace(0,3500,8)
  cbar.set_ticks(cbl); cbar.set_ticklabels(['%02.1f km'%(lev/1e3) for lev in cbl])

  ## Annotation
  # land/sea mask
  map.drawlsmask(ocean_color='blue', land_color='green',resolution=res,grid=grid)
  # add map stuff
  map.drawmapscale(-162, 31.5, -137, 57, 2000, barstyle='fancy', yoffset=0.01*(map.ymax-map.ymin))
  map.drawcoastlines(linewidth=0.5)
  map.drawcountries(linewidth=0.5)
#  map.drawrivers(linewidth=0.5)
#  map.fillcontinents(color = 'coral')
  map.drawmapboundary(fill_color='k',linewidth=2)
  # labels = [left,right,top,bottom]
  map.drawparallels([30,60],linewidth=1, labels=[True,False,False,False])
  map.drawparallels([45,75],linewidth=0.5, labels=[True,False,False,False])
  map.drawmeridians([-180,-150,-120,-90],linewidth=1, labels=[False,False,False,True])
  map.drawmeridians([150],linewidth=1, labels=[True,False,False,False])
  map.drawmeridians([165,-165,-135,-105,-75],linewidth=0.5, labels=[False,False,False,True])
  # add labels
  ax.set_title('Topography of Outer and Inner WRF Domain')
#  ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    
  # save figure to disk
  f.savefig(folder+'Topography.pdf', **sf) # save figure to pdf
  print(('\nSaved figure in '+folder+'Topography.pdf'))
  # show plots
  pyl.show()

  ## more projections  
  # setup lambert azimuthal equal area basemap.
  # lat_ts is latitude of true scale.
  # lon_0,lat_0 is central point.
#  laea = dict(projection='laea', lat_0=57, lon_0=-137, lat_ts=53, resolution='l', #
#              width=259*30e3, height=179*30e3, rsphere=(6378137.00,6356752.3142), area_thresh = 1000.)
  # lon_0, lat_0 are the center point of the projection.
  # resolution = 'l' means use low resolution coastlines.
#  ortho    = dict(projection='ortho', lat_0 = 57, lon_0 = -137, resolution = 'l', area_thresh = 1000.)
  #              'parallels':[30,50,70], 'meridians':[-180,-150,-120,-90], 'labels':[1,0,0,1]} 
