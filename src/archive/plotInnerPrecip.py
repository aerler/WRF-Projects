'''
Created on 2012-09-29

A simple script that reads a WRF netcdf-4 file and displays a 2D field in a proper geographic projection;
application here is plotting precipitation in the inner WRF domain.

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
nax = 2 # number of panels
ndom = 2
sf = dict(dpi=150) # print properties
folder = '/home/me/Research/Dynamical Downscaling/figures/' # figure directory

if __name__ == '__main__':
  
  ## read data
  data = openWRF('ctrl-1',[1982],list(range(11,12)))
  print(data[ndom-1])
  ## compute data
  precip = []; ndays = []
  for n in range(ndom):
    nrec = data[n].time.values[-1]+1
    ndays = data[n].xtime(time=nrec-1).get() /24/60 # xtime is in minutes, need days    
    dailyrain = data[n].rain(time=nrec-1).get() / ndays
#    ndays = ( data[n].xtime(time=nrec-1).get() - data[n].xtime(time=0).get() )/24/60 # xtime is in minutes, need days    
#    dailyrain = ( data[n].rain(time=nrec-1).get() - data[n].rain(time=0).get() ) / ndays 
    precip.append(dailyrain.squeeze()) 
  
  ## setup projection
  f = pyl.figure(facecolor='white', figsize = (6.25,4.25))
  ax = []
  for n in range(nax):
    ax.append(f.add_subplot(1,2,n+1))
  f.subplots_adjust(bottom=0.12, left=0.06, right=.97, top=.94, hspace=0.05, wspace=0.05) # hspace, wspace
  # setup lambert conformal basemap.
  # lat_1 is first standard parallel.
  # lat_2 is second standard parallel (defaults to lat_1).
  # lon_0,lat_0 is central point.
  # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
  # area_thresh=1000 means don't plot coastline features less
  # than 1000 km^2 in area.
  lcc = dict(projection='lcc', lat_0=59, lon_0=-123, lat_1=53, rsphere=(6378137.00,6356752.3142),#
              width=310*10e3, height=315*10e3, area_thresh = 1000., resolution='l')
  # map projection boundaries for inner WRF domain
  map = [] 
  for n in range(nax):
    map.append(Basemap(ax=ax[n],**lcc)) # one map for each panel!!
  
  ## Plot data
  grid = 10; res = 'l'
  clevs = np.linspace(0,25,51)
  norm = mpl.colors.Normalize(vmin=min(clevs),vmax=max(clevs),clip=True) 
  cmap = mpl.cm.gist_ncar #s3pcpn
  cmap.set_over('purple'); cmap.set_under('blue')
  # coordinates
  lat = []; lon = []; x = []; y = []
  for n in range(ndom):
    lat.append(data[n].lat.get()) 
    lon.append(data[n].lon.get())
    xx, yy = map[0](lon[n],lat[n]) # convert to map-native coordinates
    x.append(xx); y.append(yy)
  # draw boundaries of inner and outer domains
  bdy2 = np.ones_like(lat[1]); bdy2[0,:]=0; bdy2[-1,:]=0; bdy2[:,0]=0; bdy2[:,-1]=0
  for n in range(nax):
    # N.B.: bdy2 depends on inner domain coordinates x[1],y[1]
    map[n].contour(x[1],y[1],bdy2,[0],ax=ax[n], colors='k') # draw boundary of inner domain
#  # terrain data: mask out ocean
#  zs = []
#  for n in xrange(ndom):
#    zs.append(maskoceans(lon[n],lat[n],data[n].zs.get(),resolution=res,grid=grid))
  # draw data
  cd = []  
  for n in range(nax): # only plot first domain in first panel
    for m in range(n+1): # but also plot first domain in second panel (as background) 
      print('panel %i / domain %i'%(n,m))
      print('precip: min %f / max %f / mean %f'%(precip[m].min(),precip[m].max(),precip[m].mean()))
      cd.append(map[n].contourf(x[m],y[m],precip[m],clevs,ax=ax[n],cmap=cmap, norm=norm,extend='both'))
  
  # add colorbar
  cax = f.add_axes([0.1, 0.06, 0.8, 0.03])
  for cn in cd: # [c1d1, c1d2, c2d2]:
    cn.set_clim(vmin=min(clevs),vmax=max(clevs))
  cbar = f.colorbar(cax=cax,mappable=cd[0],orientation='h',extend='both') # ,size='3%',pad='2%' 
  cbl = np.linspace(min(clevs),max(clevs),6)
  cbar.set_ticks(cbl); cbar.set_ticklabels(['%02.1f mm'%(lev) for lev in cbl])

  ## Annotation
  # add labels
  f.suptitle('Average Daily Precipitation',fontsize=12)
  ax[0].set_title('Outer Domain (30 km)',fontsize=11)
  ax[1].set_title('Inner Domain (10 km)',fontsize=11)
#  ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
  map[0].drawmapscale(-135, 49, -137, 57, 800, barstyle='fancy', yoffset=0.01*(map[n].ymax-map[n].ymin))
  for n in range(nax):
    if n == 0 or n == 1: Bottom = True
    else: Bottom = False 
    if n == 0: Left = True
    else: Left = False
    # land/sea mask
    map[n].drawlsmask(ocean_color='blue', land_color='green',resolution=res,grid=grid)
    # add map stuff
    map[n].drawcoastlines(linewidth=0.5)
    map[n].drawcountries(linewidth=0.5)
  #  map[n].drawrivers(linewidth=0.5)
  #  map[n].fillcontinents(color = 'coral')
    map[n].drawmapboundary(fill_color='k',linewidth=2)
    # labels = [left,right,top,bottom]
    map[n].drawparallels([45,65],linewidth=1, labels=[Left,False,False,False])
    map[n].drawparallels([55,75],linewidth=0.5, labels=[Left,False,False,False])
    map[n].drawmeridians([-140,-120,-100],linewidth=1, labels=[False,False,False,Bottom])
    map[n].drawmeridians([-150,-130,-110],linewidth=0.5, labels=[False,False,False,Bottom])
    
  # save figure to disk
  f.savefig(folder+'AnnualPrecip.pdf', **sf) # save figure to pdf
  print(('\nSaved figure in '+folder+'AnnualPrecip.pdf'))
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
