'''
Created on 2012-11-05

A simple script that reads a WRF and NARR data and displays it in a proper geographic projection;
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
from myDatasets.loadNARR import openNARR
from myPlots.plots import surfacePlot
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm, maskoceans

#from pygeode.plot import plot_v1 as pl
#from pygeode.plot import basemap as bm

## figure settings
#subplot = (2,2)
#figsize = (6.25,6.25)
figtitle = 'Annual Average Precipitation [mm/day]'
axtitles = ['WRF Ctrl-1 (GPC)', 'WRF Ctrl-1 (TCS)']
sf = dict(dpi=150) # print properties
folder = '/home/me/Research/Dynamical Downscaling/figures/' # figure directory

if __name__ == '__main__':
  
  # load WRF data
  domains = [2]; ndom = len(domains) # WRF domain
#  way2, = openWRF(exp='2way-ctrl',years=[1982],months=[9], domains=domains)
#  grell3, = openWRF(exp='grell3-ctrl',years=[1980],months=[8], domains=domains) 
#  ctrl1, = openWRF(exp='ctrl-1',years=[1984],months=[7], domains=domains)
  ctrl1, = openWRF(exp='ctrl-1',years=[1980],months=[12], domains=domains)
  ctrl2, = openWRF(exp='ctrl-2',years=[1980],months=[12], domains=domains)
  WRF = (ctrl1, ctrl2)
  # load NARR data 
  NARR = openNARR() # just load whole climatology
  print(NARR)
  
  ## compute data
  data = []; lon = []; lat=[]  # list of data and coordinate fields to be plotted 
  # compute average WRF precip
  for exp in WRF:
    nrec = exp.time.values[-1]+1
    ndays = exp.xtime(time=nrec-1).get() /24/60 # xtime is in minutes, need days    
    dailyrain = exp.rain(time=nrec-1).get().squeeze() / ndays
    data.append(dailyrain)
    lon.append(exp.lon.get()); lat.append(exp.lat.get())
  # compute annual mean for NARR
  if len(data) < len(axtitles):
    months = [31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]; days = 0.
  #  rng = [0,1,2,9,10,11]
    rng = list(range(12))
    months = [months[n] for n in rng] 
    for mon in months: days += mon  
    weight = [ tmp/days for tmp in months ] # weights
    dailyrain = np.zeros((len(NARR.y.values),len(NARR.x.values))) # allocate array
    # sum up contributions from each month
    for n in range(len(rng)):
      time = NARR.time.values[rng[n]]
      dailyrain += NARR.precip(time=time).get().squeeze() * weight[n]
    data.append(dailyrain)
    lon.append(NARR.lon.get()); lat.append(NARR.lat.get())
  
  ## setup projection
  figsize = (6.25,4.25)
  (ie,je) = (1,2)
  nax = ie*je # number of panels
  # make figure and axes
  f = pyl.figure(facecolor='white', figsize=figsize)
  ax = []
  for n in range(nax):
    ax.append(f.add_subplot(ie,je,n+1))
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
  grid = 10; res = lcc['resolution']
  # map projection boundaries for inner WRF domain
  map = [] 
  for n in range(nax):
    map.append(Basemap(ax=ax[n],**lcc)) # one map for each panel!!  
  # transform coordinates (on per-map basis)
  x = []; y = []
  for n in range(nax):
    xx, yy = map[0](lon[n],lat[n]) # convert to map-native coordinates
    x.append(xx); y.append(yy)
    
  ## Plot data
  # color map
  clevs = np.linspace(0,25,51)
  norm = mpl.colors.Normalize(vmin=min(clevs),vmax=max(clevs),clip=True) 
  cmap = mpl.cm.gist_ncar #s3pcpn
  cmap.set_over('purple'); cmap.set_under('blue')
  # draw boundaries of inner domain
  bdy2 = np.ones_like(x[1]); bdy2[0,:]=0; bdy2[-1,:]=0; bdy2[:,0]=0; bdy2[:,-1]=0
  for n in range(nax):
    # N.B.: bdy2 depends on inner domain coordinates x[1],y[1]
    map[n].contour(x[1],y[1],bdy2,[0],ax=ax[n], colors='k') # draw boundary of inner domain
  # draw data
  cd = []  
  for n in range(nax): # only plot first domain in first panel 
    print('panel %i: min %f / max %f / mean %f'%(n,data[n].min(),data[n].max(),data[n].mean()))
    cd.append(map[n].contourf(x[n],y[n],data[n],clevs,ax=ax[n],cmap=cmap, norm=norm,extend='both'))  
  # add colorbar
  cax = f.add_axes([0.1, 0.06, 0.8, 0.03])
  for cn in cd: # [c1d1, c1d2, c2d2]:
    cn.set_clim(vmin=min(clevs),vmax=max(clevs))
  cbar = f.colorbar(cax=cax,mappable=cd[0],orientation='h',extend='both') # ,size='3%',pad='2%' 
  cbl = np.linspace(min(clevs),max(clevs),6)
  cbar.set_ticks(cbl); cbar.set_ticklabels(['%02.1f'%(lev) for lev in cbl])

  ## Annotation
  # add labels
  f.suptitle(figtitle,fontsize=12)
#  ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
  map[0].drawmapscale(-135, 49, -137, 57, 800, barstyle='fancy', yoffset=0.01*(map[n].ymax-map[n].ymin))
  n = -1 # axes counter
  for i in range(ie):
    for j in range(je):
      n += 1 # count up
      ax[n].set_title(axtitles[n],fontsize=11) # axes title
      if j == 0 : Left = True
      else: Left = False 
      if i == ie-1: Bottom = True
      else: Bottom = False
      # land/sea mask
      map[n].drawlsmask(ocean_color='blue', land_color='green',resolution=res,grid=grid)
      # add map stuff
      map[n].drawcoastlines(linewidth=0.5)
      map[n].drawcountries(linewidth=0.5)
      map[n].drawmapboundary(fill_color='k',linewidth=2)
      # labels = [left,right,top,bottom]
      map[n].drawparallels([45,65],linewidth=1, labels=[Left,False,False,False])
      map[n].drawparallels([55,75],linewidth=0.5, labels=[Left,False,False,False])
      map[n].drawmeridians([-140,-120,-100],linewidth=1, labels=[False,False,False,Bottom])
      map[n].drawmeridians([-150,-130,-110],linewidth=0.5, labels=[False,False,False,Bottom])
    
  # save figure to disk
#  f.savefig(folder+'AnnualPrecip.pdf', **sf) # save figure to pdf
#  print('\nSaved figure in '+folder+'AnnualPrecip.pdf')
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
