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
from myDatasets.loadWRF import openWRFclim
from myDatasets.loadCESM import openCESMclim
from myDatasets.loadCFSR import openCFSRclim
from myDatasets.loadNARR import openNARRclim
from myPlots.plots import surfacePlot
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm, maskoceans

#from pygeode.plot import plot_v1 as pl
#from pygeode.plot import basemap as bm

## figure settings
## 1 panel
#subplot = (1,1)
#figsize = (7,5.5)
#figsize = (6.25,6.25)
## 2 panel
subplot = (1,2)
figsize = (6.25,5.5)
# 4 panel
#subplot = (2,2)
#figsize = (6.25,6.25)
#margins = dict(bottom=0.025, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
#caxpos = [0.91, 0.05, 0.03, 0.9]
sf = dict(dpi=150) # print properties
folder = '/home/me/Research/Dynamical Downscaling/Figures/' # figure directory
lprint = True

if __name__ == '__main__':
  
  
  ## 1979-1981 comparison
#  dom = (1,2)
#  filename = 'wrfsrfc_d%02i_clim_1979-1985.nc' # domain substitution
#  ctrl1 = openWRFclim(exp='ctrl-1', filetypes=[filename], domains=dom)
#  ensA = openWRFclim(exp='ens-A-ctrl', filetypes=[filename], domains=dom)
#  ensB = openWRFclim(exp='ens-B-ctrl', filetypes=[filename], domains=dom)
#  ensC = openWRFclim(exp='ens-C-ctrl', filetypes=[filename], domains=dom)
##  ctrl2 = openWRFclim(exp='ctrl-2', filetypes=[filename], domains=dom)
##  way2 = openWRFclim(exp='2way-ctrl', filetypes=[filename], domains=dom)
##  grell3 = openWRFclim(exp='grell3-ctrl', filetypes=[filename], domains=dom)
##  exps = [ctrl1, way2, ctrl2, grell3]; case = 'config_d%02i'%dom[-1]
##  axtitles = ['WRF Control-1 (GPC)','WRF 2-way Nesting','WRF Control-2 (TCS)', 'WRF Grell-3']
#  exps = [ctrl1, ensA, ensB, ensC]; case = 'ensemble_d%02i'%dom[-1] 
#  axtitles = ['WRF Control-1', 'Ensemble-A', 'Ensemble-C', 'Ensemble-C']
#  print ctrl1[-1]

  ## Projection RCP 8.5
  dom = (1,2)
  ctrl1 = openWRFclim(exp='ctrl-1',filetypes=['wrfsrfc_d%02i_clim_1979-1985.nc'], domains=dom)
  rcp85 = openWRFclim(exp='ctrl-1-rcp85', domains=dom); print rcp85[-1]
#  CFSR = openCFSRclim(filename='CFSRclimFineRes1979-1986_6.nc'); print CFSR
#  NARR = openNARRclim(); print NARR
  exps = [ctrl1, rcp85]; case = '2045_d%02i'%dom[-1]
  axtitles = ['WRF Control','2045 - 2055, RCP 8.5']     
 
  ## 1979 only
#  dom = (1,)
#  fdda = openWRFclim(exp='ctrl-1',filetypes=['wrfsrfc_d%02i_monthly.nc'], domains=dom)
#  hitop = openWRFclim(exp='hitop-ctrl', domains=dom)
#  nofdda = openWRFclim(exp='nofdda-ctrl', domains=dom)
#  nodahi = openWRFclim(exp='nofdda-hitop', domains=dom)
#  exps = [fdda, hitop, nofdda, nodahi]; case = 'hitop_d%02i'%dom[-1]
#  axtitles = ['WRF Control (1979)','High-top (with Nudging)','No Spectral Nudging', 'High-top, no Nudging']

  ## Hi-Top test
#  dom = (1,2)
#  ctrl1 = openWRFclim(exp='ctrl-1', domains=dom)
#  hitop = openWRFclim(exp='hitop-ctrl', domains=dom); print hitop
#  CFSR = openCFSRclim(filename='CFSRclimFineRes1979-1986_6.nc'); print CFSR
#  NARR = openNARRclim(); print NARR
#  exps = [ctrl1, hitop, CFSR, NARR]; case = 'hitop_d%02i'%dom[-1]
#  axtitles = ['WRF Control','High-top','CFSR 1979 - 1986','NARR Climatology']     
  
  ## Reanalysis and CESM comparison
#  ctrl1 = openWRFclim(exp='ctrl-1', domains=(1,2)); print ctrl1
#  CESM = openCESMclim(); print CESM
#  CFSR = openCFSRclim(filename='CFSRclimFineRes1979-1986_6.nc'); print CFSR
#  NARR = openNARRclim(); print NARR
#  exps = [CESM]; case = 'ortho' # CESM only
#  axtitles = ['CESM 1979 - 1986']
#  exps = [[CESM,]]*4; case = 'CESM_small' #  ctrl1[0], ctrl1[1]
#  axtitles = ['CESM']*4
#  exps = [ctrl1[0],CESM,CFSR,ctrl1]; case = 'OGCI_intermed' #  ctrl1[0], ctrl1[1]
#  axtitles = ['WRF Outer Domain','CESM 1979 - 1986','CFSR Climatology', 'WRF Inner Domain']
#  exps = [[CESM,],]; case = 'CESM_large' #  ctrl1[0], ctrl1[1]
#  axtitles = ['']
#  exps = [CESM,ctrl1]; case = 'noWRF' #  ctrl1[0], ctrl1[1]
#  axtitles = ['CESM ~80km', '10km Resolution']
#  exps = [[CESM,ctrl1[0], ctrl1[1]]]; case = 'nested_large' #   
#  axtitles = ['']
#  exps = [ctrl1[0], CESM, CFSR, NARR]; case = 'OGCN4' # only outer domain
#  axtitles = ['WRF Outer Domain','CESM 1979 - 1986','CFSR 1979 - 1986', 'NARR Climatology']  
#  exps = [ctrl1, CESM, CFSR, NARR]; case = 'IGCN4' # with inner domain inset
#  axtitles = ['WRF Inner Domain','CESM 1979 - 1986','CFSR 1979 - 1986', 'NARR Climatology']

  
  ## make tuples (or lists actually)
  nexps = []
  for n in xrange(len(exps)):
    if not isinstance(exps[n],(tuple,list)):
      exps[n] = (exps[n],)
    nexps.append(len(exps[n]))    
  
  ## setup projection: lambert conformal
  # lon_0,lat_0 is central point.
  ## Lambert Conic Conformal - Fine Domain
#  projection = dict(projection='lcc', lat_0=58, lon_0=-132, lat_1=53, rsphere=(6378137.00,6356752.3142),#
#              width=200*10e3, height=300*10e3, area_thresh = 1000., resolution='l')
#  projtype = 'lcc-fine'
  ## Lambert Conic Conformal - Small Domain
  projection = dict(projection='lcc', lat_0=56, lon_0=-130, lat_1=53, rsphere=(6378137.00,6356752.3142),#
              width=2500e3, height=2650e3, area_thresh = 1000., resolution='l')
  projtype = 'lcc-small'
  ## Lambert Conic Conformal - Intermed Domain
#  projection = dict(projection='lcc', lat_0=57, lon_0=-140, lat_1=53, rsphere=(6378137.00,6356752.3142),#
#              width=4000e3, height=3400e3, area_thresh = 1000., resolution='l')
#  projtype = 'lcc-intermed'
  ## Lambert Conic Conformal - Large Domain
#  projection = dict(projection='lcc', lat_0=54.5, lon_0=-140, lat_1=53, #rsphere=(6378137.00,6356752.3142),#
#              width=11000e3, height=7500e3, area_thresh = 10e3, resolution='l')
#  projtype = 'lcc-large'
  ## Lambert Azimuthal Equal Area
  # lat_ts is latitude of true scale.
#  laea = dict(projection='laea', lat_0=57, lon_0=-137, lat_ts=53, resolution='l', #
#              width=259*30e3, height=179*30e3, rsphere=(6378137.00,6356752.3142), area_thresh = 1000.)
  ## Orthographic Projection
#  projection = dict(projection='ortho', lat_0 = 57, lon_0 = -137, resolution = 'l', area_thresh = 1000.)
#  projtype = 'ortho-NA'
  grid = 10; res = projection['resolution']
  
  ## loop over vars and seasons
  vars = ['snow']; seasons = ['spring']
#  vars = ['SST','T2','rain']
#  vars = ['T2','rain']
  seasons = ['winter', 'summer', 'annual']
#  vars = ['rain']; seasons = ['annual']
#  vars = ['zs']; seasons = ['hidef']
  # start loop
  for var in vars:
    oldvar = var
    for season in seasons:
      
      ## settings
      # plot variable and averaging 
      cbl = None; clim = None; cbo = 'horizontal' # vertical 
      lframe = True # draw domain boundary
      lmskocn = False; lmsklnd = False # mask ocean or land?
      # color maps and   scale (contour levels)
      cmap = mpl.cm.gist_ncar; cmap.set_over('white'); cmap.set_under('black')
      if var == 'snow': # snow (liquid water equivalent) 
        clevs = np.linspace(0,5,51); clbl = '%3.2f' # m
      elif var == 'rain': # total precipitation 
        clevs = np.linspace(0,25,51); clbl = '%02.1f' # mm/day
      elif oldvar=='SST' or var=='SST': # skin temperature (SST)
        clevs = np.linspace(220,310,36); clbl = '%03.0f' # K
        var = 'Ts'; lmsklnd = True # mask land
      elif var=='T2' or var=='Ts': # 2m or skin temperature (SST)
        clevs = np.linspace(240,310,36); clbl = '%03.0f' # K
      elif var == 'zs': # surface elevation / topography
        if season == 'topo':
          lmskocn = True; clim = (-1.,2.5); # nice geographic map feel
          clevs = np.hstack((np.array((-1.5,)), np.linspace(0,2.5,26))); clbl = '%02.1f' # km
          cmap = mpl.cm.gist_earth; cmap.set_over('white'); cmap.set_under('blue') # topography
        elif season == 'hidef': 
          lmskocn = True; clim = (-0.5,2.5); # good contrast for high elevation
          clevs = np.hstack((np.array((-.5,)), np.linspace(0,2.5,26))); clbl = '%02.1f' # km
          cmap = mpl.cm.gist_ncar; cmap.set_over('white'); cmap.set_under('blue')
        cbl = np.linspace(0,clim[-1],6)
      # time frame / season
      if season == 'annual':  # all month
        month = range(1,13); plottype = 'Annual Average'
      elif season == 'winter':# DJF
        month = [12, 1, 2]; plottype = 'Winter Average'
      elif season == 'spring': # MAM
        month = [3, 4, 5]; plottype = 'Spring Average'
      elif season == 'summer': # JJA
        month = [6, 7, 8]; plottype = 'Summer Average'
      else:
        plottype = '' # for static fields      
      # assemble plot title
      filename = '%s_%s_%s.pdf'%(var,season,case)
      plat = exps[0][0].vardict[var].plotatts 
      figtitle = '%s %s [%s]'%(plottype,plat['plottitle'],plat['plotunits']) # Annual Average
      
      # feedback
      print('\n\n   ***  %s %s (%s)   ***   \n'%(plottype,plat['plottitle'],var))
      
      ## compute data
      data = []; lons = []; lats=[]  # list of data and coordinate fields to be plotted 
      # compute average WRF precip
      wrfmon = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
      stdmon = np.array([31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
      for exptpl in exps:
        lontpl = []; lattpl = []; datatpl = []                
        for exp in exptpl:
          # handle dimensions
          assert (exp.lon.naxes == 2) and (exp.lat.naxes == 2), '\nWARNING: no coordinate fields found!'
          lon = exp.lon.get(); lat = exp.lat.get()
          lontpl.append(lon); lattpl.append(lat) # append to data list
          # figure out calendar
          if 'WRF' in exp.atts.get('description',''): mon = wrfmon
          else: mon = stdmon
          # extract data field
          vardata = np.zeros((exp.y.size,exp.x.size)) # allocate array
          # compute average over seasonal range
          days = 0
          expvar = exp.vardict[var]
          if expvar.hasaxis('time'):
            for m in month:
              n = m-1 
              vardata += expvar(time=exp.time.values[n]).get().squeeze() * mon[n]
              days += mon[n]
            vardata /=  days # normalize
          else:
            vardata = expvar.get().squeeze()
          vardata = vardata * expvar.plotatts.get('scalefactor',1) # apply unit conversion          
          if lmskocn: 
            if exp.vardict.has_key('lnd'): # CESM and CFSR 
              vardata[exp.lnd.get()<0.5] = -2. # use land fraction
            elif exp.vardict.has_key('lndidx'): 
              mask = exp.lndidx.get()
              vardata[mask==16] = -2. # use land use index (ocean)  
              vardata[mask==24] = -2. # use land use index (lake)
            else : vardata = maskoceans(lon,lat,vardata,resolution=res,grid=grid)
          if lmsklnd: 
            if exp.vardict.has_key('lnd'): # CESM and CFSR 
              vardata[exp.lnd.get()>0.5] = 0 # use land fraction
            elif exp.vardict.has_key('lndidx'): # use land use index (ocean and lake)
              mask = exp.lndidx.get(); tmp = vardata.copy(); vardata[:] = 0.
              vardata[mask==16] = tmp[mask==16]; vardata[mask==24] = tmp[mask==24]
          datatpl.append(vardata) # append to data list
        # add tuples to master list
        lons.append(lontpl); lats.append(lattpl); data.append(datatpl)
        
      ## setup projection
      nax = subplot[0]*subplot[1] # number of panels
      # figure out colorbar placement
      if cbo == 'vertical':
        margins = dict(bottom=0.02, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
        caxpos = [0.91, 0.05, 0.03, 0.9]
      else: 
        margins = dict(bottom=0.1, left=0.065, right=.9725, top=.925, hspace=0.05, wspace=0.05)
        caxpos = [0.05, 0.05, 0.9, 0.03]
      # make figure and axes
      f = pyl.figure(facecolor='white', figsize=figsize)
      ax = []
      for n in xrange(nax):
        ax.append(f.add_subplot(subplot[0],subplot[1],n+1))
      f.subplots_adjust(**margins) # hspace, wspace
        # lat_1 is first standard parallel.
      # lat_2 is second standard parallel (defaults to lat_1).
      # lon_0,lat_0 is central point.
      # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
      # area_thresh=1000 means don't plot coastline features less
      # than 1000 km^2 in area.
      # map projection boundaries for inner WRF domain
      maps = [] 
      for n in xrange(nax):
        maps.append(Basemap(ax=ax[n],**projection)) # one map for each panel!!  
      # transform coordinates (on per-map basis)
      x = []; y = []
      for n in xrange(nax):
        xtpl = []; ytpl = []
        for m in xrange(nexps[n]):
          xx, yy = maps[n](lons[n][m],lats[n][m]) # convert to map-native coordinates
          xtpl.append(xx); ytpl.append(yy)
        x.append(xtpl); y.append(ytpl) 
        
      ## Plot data
      # draw boundaries of inner domain
      if lframe:
        for n in xrange(nax):
          for m in xrange(nexps[n]):   
            bdy = np.ones_like(x[n][m]); bdy[0,:]=0; bdy[-1,:]=0; bdy[:,0]=0; bdy[:,-1]=0
            maps[n].contour(x[n][m],y[n][m],bdy,[0],ax=ax[n], colors='k') # draw boundary of inner domain
      # draw data
      norm = mpl.colors.Normalize(vmin=min(clevs),vmax=max(clevs),clip=True) # for colormap
      cd = []  
      for n in xrange(nax): 
        for m in xrange(nexps[n]):
          print 'panel %i: min %f / max %f / mean %f'%(n,data[n][m].min(),data[n][m].max(),data[n][m].mean())
          cd.append(maps[n].contourf(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],cmap=cmap, norm=norm,extend='both'))  
      # add colorbar
      cax = f.add_axes(caxpos)
      for cn in cd: # [c1d1, c1d2, c2d2]:
        if clim: cn.set_clim(vmin=clim[0],vmax=clim[1])
        else: cn.set_clim(vmin=min(clevs),vmax=max(clevs))
      cbar = f.colorbar(cax=cax,mappable=cd[0],orientation=cbo,extend='both') # ,size='3%',pad='2%'       
      if cbl is None: cbl = np.linspace(min(clevs),max(clevs),6)
      cbar.set_ticks(cbl); cbar.set_ticklabels([clbl%(lev) for lev in cbl])
    
      ## Annotation
      # add labels
      f.suptitle(figtitle,fontsize=12)
    #  ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
      if projtype == 'lcc-small':
        maps[0].drawmapscale(-136, 49, -137, 57, 800, barstyle='fancy', yoffset=0.01*(maps[n].ymax-maps[n].ymin))
      elif projtype == 'lcc-large':
        maps[0].drawmapscale(-171, 21, -137, 57, 2000, barstyle='fancy', yoffset=0.01*(maps[n].ymax-maps[n].ymin))
      n = -1 # axes counter
      for i in xrange(subplot[0]):
        for j in xrange(subplot[1]):
          n += 1 # count up
          ax[n].set_title(axtitles[n],fontsize=11) # axes title
          if j == 0 : Left = True
          else: Left = False 
          if i == subplot[0]-1: Bottom = True
          else: Bottom = False
          # land/sea mask
          maps[n].drawlsmask(ocean_color='blue', land_color='green',resolution=res,grid=grid)
          # black-out continents, if we have no proper land mask 
          if lmsklnd and not (exps[n][0].vardict.has_key('lnd') or exps[n][0].vardict.has_key('lndidx')): 
            maps[n].fillcontinents(color='black',lake_color='black') 
          # add maps stuff
          maps[n].drawcoastlines(linewidth=0.5)
          maps[n].drawcountries(linewidth=0.5)
          maps[n].drawmapboundary(fill_color='k',linewidth=2)
          # labels = [left,right,top,bottom]
          if projtype=='lcc-fine' or projtype=='lcc-small' or projtype=='lcc-intermed':
            maps[n].drawparallels([45,65],linewidth=1, labels=[Left,False,False,False])
            maps[n].drawparallels([55,75],linewidth=0.5, labels=[Left,False,False,False])
            maps[n].drawmeridians([-180,-160,-140,-120,-100],linewidth=1, labels=[False,False,False,Bottom])
            maps[n].drawmeridians([-170,-150,-130,-110],linewidth=0.5, labels=[False,False,False,Bottom])
          elif projtype == 'lcc-large':
            maps[n].drawparallels(range(0,90,30),linewidth=1, labels=[Left,False,False,False])
            maps[n].drawparallels(range(15,90,30),linewidth=0.5, labels=[Left,False,False,False])
            maps[n].drawmeridians(range(-180,180,30),linewidth=1, labels=[False,False,False,Bottom])
            maps[n].drawmeridians(range(-165,180,30),linewidth=0.5, labels=[False,False,False,Bottom])
          elif projtype == 'ortho':
            maps[n].drawparallels(range(-90,90,30),linewidth=1)
            maps[n].drawmeridians(range(-180,180,30),linewidth=1)
        
      # save figure to disk
      if lprint:
        f.savefig(folder+filename, **sf) # save figure to pdf
        print('\nSaved figure in '+filename)
        print(folder)
  
  ## show plots after all iterations
  pyl.show()

