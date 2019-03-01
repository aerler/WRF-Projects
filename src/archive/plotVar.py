'''
Created on 2012-11-05

A simple script that reads a WRF and NARR data and displays it in a proper geographic projection;
application here is plotting precipitation in the inner WRF domain.

@author: Andre R. Erler
'''

## includes
from copy import copy # to copy map projection objects
# matplotlib config: size etc.
import numpy as np
import matplotlib.pylab as pyl
import matplotlib as mpl
mpl.rc('lines', linewidth=1.)
mpl.rc('font', size=10)
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
# pygeode stuff
from myDatasets.loadWRF import openWRFclim, WRFtitle
from myDatasets.loadCESM import openCESMclim, CESMtitle
from myDatasets.loadCFSR import openCFSRclim
from myDatasets.loadNARR import openNARRclim
from myDatasets.loadGPCC import openGPCCclim
from myDatasets.loadCRU import openCRUclim
from myDatasets.loadPRISM import openPRISMclim

#from pygeode.plot import plot_v1 as pl
#from pygeode.plot import basemap as bm

## figure settings
def getFigureSettings(nexp, cbo):
  sf = dict(dpi=150) # print properties
  figformat = 'png'
  folder = '/home/me/Research/Dynamical Downscaling/Figures/' # figure directory
  # figure out colorbar placement
  if cbo == 'vertical':
    margins = dict(bottom=0.02, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
    caxpos = [0.91, 0.05, 0.03, 0.9]
  else:# 'horizontal'
    margins = dict(bottom=0.1, left=0.065, right=.9725, top=.925, hspace=0.05, wspace=0.05)
    caxpos = [0.05, 0.05, 0.9, 0.03]        
  # pane settings
  if nexp == 1:
    ## 1 panel
    subplot = (1,1)
    figsize = (3.75,3.75) #figsize = (6.25,6.25)  #figsize = (7,5.5)
    margins = dict(bottom=0.025, left=0.075, right=0.875, top=0.875, hspace=0.0, wspace=0.0)
#     margins = dict(bottom=0.12, left=0.075, right=.9725, top=.95, hspace=0.05, wspace=0.05)
#    margins = dict(bottom=0.025, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
  elif nexp == 2:
    ## 2 panel
    subplot = (1,2)
    figsize = (6.25,5.5)
  elif nexp == 4:
    # 4 panel
    subplot = (2,2)
    figsize = (6.25,6.25)
    margins = dict(bottom=0.025, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
  elif nexp == 4:
    # 4 panel
    subplot = (2,2)
    figsize = (6.25,6.25)
    margins = dict(bottom=0.025, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
  elif nexp == 6:
    # 6 panel
    subplot = (2,3) # rows, columns
    figsize = (9.25,6.5) # width, height (inches)
    cbo = 'horizontal'
    margins = dict(bottom=0.09, left=0.05, right=.97, top=.92, hspace=0.1, wspace=0.05)
    caxpos = [0.05, 0.025, 0.9, 0.03]
  #    margins = dict(bottom=0.025, left=0.065, right=.885, top=.925, hspace=0.05, wspace=0.05)
  # return values
  return sf, figformat, folder, margins, caxpos, subplot, figsize, cbo


## setup projection: lambert conformal
def getProjectionSettings(projtype):
  # lon_0,lat_0 is central point. lat_ts is latitude of true scale.
  if projtype == 'lcc-new':
    ## Lambert Conic Conformal - New Fine Domain
    projection = dict(projection='lcc', lat_0=55, lon_0=-120, lat_1=52, rsphere=(6378137.00,6356752.3142),#
                width=180*10e3, height=180*10e3, area_thresh = 1000., resolution='l')
  elif projtype == 'lcc-fine':
    ## Lambert Conic Conformal - Fine Domain
    projection = dict(projection='lcc', lat_0=58, lon_0=-132, lat_1=53, rsphere=(6378137.00,6356752.3142),#
                width=200*10e3, height=300*10e3, area_thresh = 1000., resolution='l')
  elif projtype == 'lcc-small':
    ## Lambert Conic Conformal - Small Domain
    projection = dict(projection='lcc', lat_0=56, lon_0=-130, lat_1=53, rsphere=(6378137.00,6356752.3142),#
                width=2500e3, height=2650e3, area_thresh = 1000., resolution='l')
  elif projtype == 'lcc-intermed':
    ## Lambert Conic Conformal - Intermed Domain
    projection = dict(projection='lcc', lat_0=57, lon_0=-140, lat_1=53, rsphere=(6378137.00,6356752.3142),#
                width=4000e3, height=3400e3, area_thresh = 1000., resolution='l')
  elif projtype == 'lcc-large':
    ## Lambert Conic Conformal - Large Domain
    projection = dict(projection='lcc', lat_0=54.5, lon_0=-140, lat_1=53, #rsphere=(6378137.00,6356752.3142),#
                width=11000e3, height=7500e3, area_thresh = 10e3, resolution='l')
  elif projtype == 'laea': 
    ## Lambert Azimuthal Equal Area
    projection = dict(projection='laea', lat_0=57, lon_0=-137, lat_ts=53, resolution='l', #
                width=259*30e3, height=179*30e3, rsphere=(6378137.00,6356752.3142), area_thresh = 1000.)  
  elif projtype == 'ortho-NA':
    ## Orthographic Projection
    projection = dict(projection='ortho', lat_0 = 75, lon_0 = -137, resolution = 'l', area_thresh = 1000.)
  # resolution of coast lines
  grid = 10; res = projection['resolution']
  # return values
  return projection, grid, res
 
# used for annotation
monthnames = ['January  ', 'February ', 'March    ', 'April    ', 'May      ', 'June     ', #
              'July     ', 'August   ', 'September', 'October  ', 'November ', 'December ']

if __name__ == '__main__':
  
#  filename = 'wrfsrfc_d%02i_clim.nc' # domain substitution
#  CFSR = openCFSRclim(filename='CFSRclimFineRes1979-2009.nc')
#  RRTMG = openWRFclim(exp='rrtmg-arb1', filetypes=['wrfsrfc_d%02i_clim.nc'], domains=dom)
#  axtitles = ['CRU Climatology', 'WRF Control', 'WRF RRTMG', 'Polar WRF']


  ## general settings and shortcuts
  H01 = '1979'; H02 = '1979-1980'; H03 = '1979-1981'; H30 = '1979-2009' # for tests 
  H05 = '1979-1983'; H10 = '1979-1988'; H15 = '1979-1993' # historical validation periods
  G10 = '1969-1978'; I10 = '1989-1998'; J10 = '1999-2008' # additional historical periods
  A03 = '2045-2047'; A05 = '2045-2049'; A10 = '2045-2054'; A15 = '2045-2059' # mid-21st century
  B03 = '2095-2097'; B05 = '2095-2099'; B10 = '2095-2104'; B15 = '2095-2109' # late 21st century
  lprint = True # write plots to disk
  ltitle = True # plot/figure title
  lcontour = False # contour or pcolor plot
  lframe = True # draw domain boundary
  cbo = 'vertical' # vertical horizontal
  resolution=None # only for GPCC (None = default/highest)
 
  ## case settings
  
  # observations
  case = 'new' # name tag
  projtype = 'lcc-new' # 'lcc-new'  
  period = H01; dom = (1,2,)
#   explist = ['CRU']*3 + ['NARR', 'CRU', 'CRU']; period = [H30, G10, H10, None, I10, J10]
#   explist = ['PRISM-10km','ctrl-1','NARR','PRISM','max','CRU']
  explist = ['PRISM-10km','new','noah','nogulf','max','CRU']; period = [H01]*5 + [H10]  
#  explist = ['GPCC']; vars = ['stns']; seasons = ['annual']
#   explist = ['cfsr', 'ctrl-1', 'max', 'NARR', 'PRISM', 'CRU']; # period = [H10]*5 + [None]
#   explist = ['cfsr', 'ens-Z', 'max', 'ctrl-1']
#   explist = ['tom', 'ens-Z', 'max', 'ctrl-1']
#  explist = ['CFSR', 'ctrl-1', 'CRU', 'NARR', 'CESM', 'GPCC']
#   explist = ['ctrl-1', 'PRISM', 'Ctrl-1', 'CRU']
#   explist = ['ctrl-1', 'grell', 'tiedt', 'PRISM', 'CRU', 'GPCC']
#   explist = ['milb', 'PRISM', 'wdm6', 'tom', 'ctrl-1', 'max']
#   explist = ['wdm6','tom', 'ctrl-1', 'max']
#  explist = ['ctrl-1', 'cam3', 'noahmp', 'pwrf']
#   explist = ['PRISM', 'max', 'ctrl-1', 'GPCC']
#   explist = ['PRISM', 'max', 'ctrl-1', 'CFSR']
#   explist = ['nmpbar', 'nmpsnw', 'nmpdef', 'ctrl-1']
#   explist = ['nmpbar', 'clm4', 'max', 'ctrl-1']
#   explist = ['PRISM', 'tom', 'tiedt', 'ctrl-1']
#   explist = ['PRISM', 'tom', 'nmpbar', 'ctrl-1']
#   explist = ['grell','PRISM', 'tiedt', 'nmpbar', 'ctrl-1', 'max']
#   explist = ['ctrl-1', 'grell', 'tiedt', 'pbl4']
#   explist = ['PRISM', 'ctrl-1-2000', 'cam3', 'Ctrl-1']
#   explist = ['PRISM', 'nmpbar', 'tom', 'ctrl-1']
#   explist = ['ctrl-1', 'tom', 'tiedt', 'nmpbar']
#   explist = ['ens-Z', 'ens-B', 'ens-A', 'ens-C'];
#   explist = ['Ens-Z', 'Ens-B', 'Ens-A', 'Ens-C']; # CESM historical
#   explist = ['SeaIce', 'Ens-B-rcp85', 'Ens-A-rcp85', 'Ens-Z-rcp85']; # CESM RCP 8.5 projections
#   explist = ['SeaIce']; # CESM RCP 8.5 projections
#   explist = ['Ctrl-1', 'Ctrl-1-rcp85', 'Ctrl-1-rcp85', 'SeaIce']; period = (H10, A10, B10, A10) 
#  explist = ['ens-Z', 'CRU', 'ens-B', 'ens-A', 'GPCC', 'ens-C']; dom = (1,)
#   explist = ['ctrl-1', 'ctrl-1-2000','ctrl-1-2050','ctrl-2-2050']; period = (H10, H10, A10, A10)
#  explist = ['ctrl-1', 'CRU', 'NARR', 'CESM']
#   explist = ['max', 'PRISM', 'grell', 'ctrl-1']
#   explist = ['ctrl-1', 'PRISM', 'GPCC', 'NARR']
#   explist = ['ctrl-1']
#   explist = ['modis']
#   explist = ['PRISM']; period = None
  
  ## select variables and seasons
#   vars = ['rainnc', 'rainc', 'T2']
#  vars = ['snowh'];  seasons = [8]
#   vars = ['rain']
#   vars = ['evap']
#   vars = ['snow']
#   vars = ['rain', 'T2', 'p-et','evap']
#   vars = ['p-et','rain','snow']
#   vars = ['GLW','OLR','qtfx']
#   vars = ['SWDOWN','GLW','OLR']
#   vars = ['hfx','lhfx']
#   vars = ['qtfx','lhfr']
  vars = ['rain','T2']
#   vars = ['T2']
#   vars = ['seaice']; seasons = [8] # September seaice
#  vars = ['rain','T2','snow']
#   vars = ['snow', 'snowh']
#  vars = ['SST','T2','rain','snow','snowh']
#   seasons = [ [i] for i in xrange(12) ] # monthly
#   seasons = ['annual']
#   seasons = ['summer']
#   seasons = ['winter']    
  seasons = ['winter', 'summer', 'annual']
#   vars = ['snow']; seasons = ['fall','winter','spring']
#  vars = ['rain']; seasons = ['annual']
#  vars = ['zs']; seasons = ['hidef']
#  vars = ['stns']; seasons = ['annual']
#   vars = ['lndcls']; seasons = [''] # static

  

  ## load data 
  if not isinstance(period,(tuple,list)): period = (period,)*len(explist)
  exps = []; axtitles = []
  for exp,prd in zip(explist,period): 
    ext = exp; axt = ''
    if isinstance(exp,str):
      if exp[0].isupper():
        if exp == 'GPCC': ext = (openGPCCclim(resolution=resolution,period=prd),); axt = 'GPCC Observations' # ,period=prd
        elif exp == 'CRU': ext = (openCRUclim(period=prd),); axt = 'CRU Observations' 
        elif exp[0:5] == 'PRISM': # all PRISM derivatives
          if exp == 'PRISM': prismfile = 'prism_clim.nc'
          elif exp == 'PRISM-10km': prismfile = 'prism_10km.nc'
          if len(vars) ==1 and vars[0] == 'rain': 
            ext = (openGPCCclim(resolution='0.25'), openPRISMclim(filename=prismfile)); axt = 'PRISM (and GPCC)'
          else: ext = (openCRUclim(period='1979-2009'), openPRISMclim(filename=prismfile)); axt = 'PRISM (and CRU)'
          # ext = (openPRISMclim(),)          
        elif exp == 'CFSR': ext = (openCFSRclim(period=prd),); axt = 'CFSR Reanalysis' 
        elif exp == 'NARR': ext = (openNARRclim(),); axt = 'NARR Reanalysis'
        else: # all other uppercase names are CESM runs  
          ext = (openCESMclim(exp=exp, period=prd),)
          axt = CESMtitle.get(exp,exp)
      else: # WRF runs are all in lower case
        ext = openWRFclim(exp=exp, period=prd, domains=dom)
        axt = WRFtitle.get(exp,exp)
    exps.append(ext); axtitles.append(axt)  
  print(exps[-1][-1])
  # count experiment tuples (layers per panel)
  nexps = []; nlen = len(exps)
  for n in range(nlen):
    if not isinstance(exps[n],(tuple,list)): # should not be necessary
      exps[n] = (exps[n],)
    nexps.append(len(exps[n])) # layer counter for each panel
  
  # get figure settings
  sf, figformat, folder, margins, caxpos, subplot, figsize, cbo = getFigureSettings(nexp=nlen, cbo=cbo)
  
  # get projections settings
  projection, grid, res = getProjectionSettings(projtype=projtype)
  
  ## loop over vars and seasons
  maps = []; x = []; y = [] # projection objects and coordinate fields (only computed once)
  # start loop
  for var in vars:
    oldvar = var
    for season in seasons:
      
      ## settings
      # plot variable and averaging 
      cbl = None; clim = None       
      lmskocn = False; lmsklnd = False # mask ocean or land?
      # color maps and   scale (contour levels)
      cmap = mpl.cm.gist_ncar; cmap.set_over('white'); cmap.set_under('black')
      if var == 'snow': # snow (liquid water equivalent) 
        lmskocn = True; clbl = '%2.0f' # kg/m^2
        clevs = np.linspace(0,200,41)
      elif var == 'snowh': # snow (depth/height) 
        lmskocn = True; clbl = '%2.1f' # m
        clevs = np.linspace(0,2,41)
      elif var=='hfx' or var=='lhfx' or var=='qtfx': # heat fluxes (W / m^2)
        clevs = np.linspace(-20,100,41); clbl = '%03.0f'
        if var == 'qtfx': clevs = clevs * 2
        if season == 'winter': clevs = clevs - 30
        elif season == 'summer': clevs = clevs + 30
      elif var=='GLW': # heat fluxes (W / m^2)
        clevs = np.linspace(200,320,41); clbl = '%03.0f'
        if season == 'winter': clevs = clevs - 40
        elif season == 'summer': clevs = clevs + 40
      elif var=='OLR': # heat fluxes (W / m^2)
        clevs = np.linspace(190,240,31); clbl = '%03.0f'
        if season == 'winter': clevs = clevs - 20
        elif season == 'summer': clevs = clevs + 30
      elif var=='rfx': # heat fluxes (W / m^2)
        clevs = np.linspace(320,470,51); clbl = '%03.0f'
        if season == 'winter': clevs = clevs - 100
        elif season == 'summer': clevs = clevs + 80
      elif var=='SWDOWN' or var=='SWNORM': # heat fluxes (W / m^2)
        clevs = np.linspace(80,220,51); clbl = '%03.0f'
        if season == 'winter': clevs = clevs - 80
        elif season == 'summer': clevs = clevs + 120
      elif var == 'lhfr': # relative latent heat flux (fraction)        
        clevs = np.linspace(0,1,26); clbl = '%2.1f' # fraction
      elif var == 'evap': # moisture fluxes (kg /(m^2 s))
        clevs = np.linspace(-4,4,25); clbl = '%02.1f'
        cmap = mpl.cm.PuOr
      elif var == 'p-et': # moisture fluxes (kg /(m^2 s))
        # clevs = np.linspace(-3,22,51); clbl = '%02.1f'
        clevs = np.linspace(-2,2,25); cmap = mpl.cm.PuOr; clbl = '%02.1f'
      elif var == 'rain' or var == 'rainnc': # total precipitation 
        clevs = np.linspace(0,20,41); clbl = '%02.1f' # mm/day
      elif var == 'rainc': # convective precipitation 
        clevs = np.linspace(0,5,26); clbl = '%02.1f' # mm/day
      elif oldvar=='SST' or var=='SST': # skin temperature (SST)
        clevs = np.linspace(240,300,61); clbl = '%03.0f' # K
        var = 'Ts'; lmsklnd = True # mask land
      elif var=='T2' or var=='Ts': # 2m or skin temperature (SST)
        clevs = np.linspace(255,290,36); clbl = '%03.0f' # K
        if season == 'winter': clevs = clevs - 10
        elif season == 'summer': clevs = clevs + 10
      elif var == 'seaice': # sea ice fraction
        lmsklnd = True # mask land        
        clevs = np.linspace(0.04,1,25); clbl = '%2.1f' # fraction
        cmap.set_under('white')
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
      elif var=='stns': # station density
        clevs = np.linspace(0,5,6); clbl = '%2i' # stations per grid points  
        cmap.set_over('purple'); cmap.set_under('white')      
      elif var=='lndcls': # land use classes (works best with contour plot)
        clevs = np.linspace(0.5,24.5,25); cbl = np.linspace(4,24,6)  
        clbl = '%2i'; cmap.set_over('purple'); cmap.set_under('white')
      # time frame / season
      if isinstance(season,str):
        if season == 'annual':  # all month
          month = list(range(1,13)); plottype = 'Annual Average'
        elif season == 'winter':# DJF
          month = [12, 1, 2]; plottype = 'Winter Average'
        elif season == 'spring': # MAM
          month = [3, 4, 5]; plottype = 'Spring Average'
        elif season == 'summer': # JJA
          month = [6, 7, 8]; plottype = 'Summer Average'
        elif season == 'fall': # SON
          month = [9, 10, 11]; plottype = 'Fall Average'
        else:
          plottype = '' # for static fields
          month = [1]
      else:                
        month = season      
        if len(season) == 1 and isinstance(season[0],int):
          plottype =  '%s Average'%monthnames[season[0]].strip()
          season = '%02i'%(season[0]+1) # number of month, used for file name
        else: plottype = 'Average'
      # assemble plot title
      filename = '%s_%s_%s.%s'%(var,season,case,figformat)
      plat = exps[0][0].vardict[var].plotatts 
      if plat['plotunits']: figtitle = '%s %s [%s]'%(plottype,plat['plottitle'],plat['plotunits'])
      else: figtitle = '%s %s'%(plottype,plat['plottitle'])
      
      # feedback
      print(('\n\n   ***  %s %s (%s)   ***   \n'%(plottype,plat['plottitle'],var)))
      
      ## compute data
      data = []; lons = []; lats=[]  # list of data and coordinate fields to be plotted 
      # compute average WRF precip
      wrfmon = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
      stdmon = np.array([31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
      print(' - loading data\n')
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
            if 'lnd' in exp.vardict: # CESM and CFSR 
              vardata[exp.lnd.get()<0.5] = -2. # use land fraction
            elif 'lndidx' in exp.vardict: 
              mask = exp.lndidx.get()
              vardata[mask==16] = -2. # use land use index (ocean)  
              vardata[mask==24] = -2. # use land use index (lake)
            else : vardata = maskoceans(lon,lat,vardata,resolution=res,grid=grid)
          if lmsklnd: 
            if 'lnd' in exp.vardict: # CESM and CFSR 
              vardata[exp.lnd.get()>0.5] = 0 # use land fraction
            elif 'lndidx' in exp.vardict: # use land use index (ocean and lake)
              mask = exp.lndidx.get(); tmp = vardata.copy(); vardata[:] = 0.
              vardata[mask==16] = tmp[mask==16]; vardata[mask==24] = tmp[mask==24]
          datatpl.append(vardata) # append to data list
        # add tuples to master list
        lons.append(lontpl); lats.append(lattpl); data.append(datatpl)
              
      ## setup projection
      #print(' - setting up figure\n') 
      nax = subplot[0]*subplot[1] # number of panels
      # make figure and axes
      f = pyl.figure(facecolor='white', figsize=figsize)
      ax = []
      for n in range(nax):
        ax.append(f.add_subplot(subplot[0],subplot[1],n+1))
      f.subplots_adjust(**margins) # hspace, wspace
      # lat_1 is first standard parallel.
      # lat_2 is second standard parallel (defaults to lat_1).
      # lon_0,lat_0 is central point.
      # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
      # area_thresh=1000 means don't plot coastline features less
      # than 1000 km^2 in area.
      if not maps:
        print(' - setting up map projection\n') 
        mastermap = Basemap(ax=ax[n],**projection)
        for axi in ax:          
          tmp = copy(mastermap)
          tmp.ax = axi  
          maps.append(tmp) # one map for each panel!!  
      else:
        print(' - resetting map projection\n') 
        for n in range(nax):
          maps[n].ax=ax[n] # assign new axes to old projection
      # transform coordinates (on per-map basis)
      if not (x and y):
        print(' - transforming coordinate fields\n')
        for n in range(nax):
          xtpl = []; ytpl = []
          for m in range(nexps[n]):
            xx, yy = maps[n](lons[n][m],lats[n][m]) # convert to map-native coordinates
            xtpl.append(xx); ytpl.append(yy)
          x.append(xtpl); y.append(ytpl) 
        
      ## Plot data
      # draw boundaries of inner domain
      if lframe:
        print(' - drawing data frames\n')
        for n in range(nax):
          for m in range(nexps[n]):   
            bdy = np.ones_like(x[n][m]); bdy[0,:]=0; bdy[-1,:]=0; bdy[:,0]=0; bdy[:,-1]=0
            maps[n].contour(x[n][m],y[n][m],bdy,[0],ax=ax[n], colors='k') # draw boundary of inner domain
      # draw data
      norm = mpl.colors.Normalize(vmin=min(clevs),vmax=max(clevs),clip=True) # for colormap
      cd = []
#       print(' - creating plots\n')  
      for n in range(nax): 
        for m in range(nexps[n]):
          print(('panel %i: min %f / max %f / mean %f'%(n,data[n][m].min(),data[n][m].max(),data[n][m].mean())))
          if lcontour: cd.append(maps[n].contourf(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],cmap=cmap, norm=norm,extend='both'))  
          else: cd.append(maps[n].pcolormesh(x[n][m],y[n][m],data[n][m],cmap=cmap,shading='gouraud'))
      # add colorbar
      cax = f.add_axes(caxpos)
      for cn in cd: # [c1d1, c1d2, c2d2]:
        if clim: cn.set_clim(vmin=clim[0],vmax=clim[1])
        else: cn.set_clim(vmin=min(clevs),vmax=max(clevs))
      cbar = f.colorbar(cax=cax,mappable=cd[0],orientation=cbo,extend='both') # ,size='3%',pad='2%'       
      if cbl is None: cbl = np.linspace(min(clevs),max(clevs),6)
      cbar.set_ticks(cbl); cbar.set_ticklabels([clbl%(lev) for lev in cbl])
    
      ## Annotation
      #print('\n - annotating plots\n')      
      # add labels
      if ltitle: f.suptitle(figtitle,fontsize=12)
    #  ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
      msn = len(maps)/2 # place scale 
      if projtype == 'lcc-new':
        maps[msn].drawmapscale(-128, 48, -120, 55, 400, barstyle='fancy', 
                             fontsize=8, yoffset=0.01*(maps[n].ymax-maps[n].ymin))
      elif projtype == 'lcc-small':
        maps[msn].drawmapscale(-136, 49, -137, 57, 800, barstyle='fancy', yoffset=0.01*(maps[n].ymax-maps[n].ymin))
      elif projtype == 'lcc-large':
        maps[msn].drawmapscale(-171, 21, -137, 57, 2000, barstyle='fancy', yoffset=0.01*(maps[n].ymax-maps[n].ymin))
      n = -1 # axes counter
      for i in range(subplot[0]):
        for j in range(subplot[1]):
          n += 1 # count up
          ax[n].set_title(axtitles[n],fontsize=11) # axes title
          if j == 0 : Left = True
          else: Left = False 
          if i == subplot[0]-1: Bottom = True
          else: Bottom = False
          # land/sea mask
          maps[n].drawlsmask(ocean_color='blue', land_color='green',resolution=res,grid=grid)
          # black-out continents, if we have no proper land mask 
          if lmsklnd and not ('lnd' in exps[n][0].vardict or 'lndidx' in exps[n][0].vardict): 
            maps[n].fillcontinents(color='black',lake_color='black') 
          # add maps stuff
          maps[n].drawcoastlines(linewidth=0.5)
          maps[n].drawcountries(linewidth=0.5)
          maps[n].drawmapboundary(fill_color='k',linewidth=2)
          # labels = [left,right,top,bottom]
          if projtype=='lcc-new':
            maps[n].drawparallels([40,50,60,70],linewidth=1, labels=[Left,False,False,False])
            maps[n].drawparallels([45,55,65],linewidth=0.5, labels=[Left,False,False,False])
            maps[n].drawmeridians([-180,-160,-140,-120,-100],linewidth=1, labels=[False,False,False,Bottom])
            maps[n].drawmeridians([-170,-150,-130,-110],linewidth=0.5, labels=[False,False,False,Bottom])          
          elif projtype=='lcc-fine' or projtype=='lcc-small' or projtype=='lcc-intermed':
            maps[n].drawparallels([45,65],linewidth=1, labels=[Left,False,False,False])
            maps[n].drawparallels([55,75],linewidth=0.5, labels=[Left,False,False,False])
            maps[n].drawmeridians([-180,-160,-140,-120,-100],linewidth=1, labels=[False,False,False,Bottom])
            maps[n].drawmeridians([-170,-150,-130,-110],linewidth=0.5, labels=[False,False,False,Bottom])
          elif projtype == 'lcc-large':
            maps[n].drawparallels(list(range(0,90,30)),linewidth=1, labels=[Left,False,False,False])
            maps[n].drawparallels(list(range(15,90,30)),linewidth=0.5, labels=[Left,False,False,False])
            maps[n].drawmeridians(list(range(-180,180,30)),linewidth=1, labels=[False,False,False,Bottom])
            maps[n].drawmeridians(list(range(-165,180,30)),linewidth=0.5, labels=[False,False,False,Bottom])
          elif projtype == 'ortho':
            maps[n].drawparallels(list(range(-90,90,30)),linewidth=1)
            maps[n].drawmeridians(list(range(-180,180,30)),linewidth=1)
        
      # mark stations
      # ST_LINA, WESTLOCK_LITKE, JASPER, FORT_MCMURRAY, SHINING_BANK
      sn = ['SL', 'WL', 'J', 'FM', 'SB']
      slon = [-111.45,-113.85,-118.07,-111.22,-115.97]
      slat = [54.3,54.15,52.88,56.65,53.85]
      for (axn,mapt) in zip(ax,maps):
        for (name,lon,lat) in zip(sn,slon,slat):
          xx,yy = mapt(lon, lat)
          mapt.plot(xx,yy,'ko',markersize=3)
          axn.text(xx+1.5e4,yy-1.5e4,name,ha='left',va='top',fontsize=8)
      
      # save figure to disk
      if lprint:
        print(('\nSaving figure in '+filename))
        f.savefig(folder+filename, **sf) # save figure to pdf
        print(folder)
  
  ## show plots after all iterations
  pyl.show()

