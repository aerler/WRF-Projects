'''
Created on 2012-11-05, adapted for PyGeoDat on 2013-10-10

A simple script that mushroomed into a complex module... reads a Datasets and displays them in a proper 
geographic projection.

@author: Andre R. Erler, GPL v3
'''

## includes
from copy import copy # to copy map projection objects
# matplotlib config: size etc.
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
       
print "Importing 'pyplot' from 'matplotlib'\n"   
  
import matplotlib.pyplot as plt
# N.B.: importing pyplot actually takes quite long! 
# prevent figures from closing: don't run in interactive mode, or plt.show() will not block
plt.ioff()
# mpl.use('Agg') # enforce QT4
mpl.rc('lines', linewidth=1.)
mpl.rc('font', size=10)
# PyGeoDat stuff
from mpl_toolkits.basemap import maskoceans # used for masking data
  
from geodata.base import DatasetError
from datasets.WSC import basins_info
from datasets.EC import province_info
from datasets.common import stn_params
from legacy_plotting.legacy import loadDatasets, checkItemList
# project related stuff
# Western Canada
# from projects.WesternCanada import getSetup, figure_folder, map_folder, getFigureSettings, getVariableSettings
# Great Lakes
from projects.GreatLakes import getSetup, getFigureSettings, getVariableSettings
from projects.GreatLakes import figure_folder, map_folder, WRF_exps, CESM_exps

station_constraints = dict()
station_constraints['min_len'] = 15 # for valid climatology
station_constraints['lat'] = (40,55)
station_constraints['end_after'] = 1980
#station_constraints['max_zerr'] = 100 # can't use this, because we are loading EC data separately from WRF
station_constraints['prov'] = ('ON')

if __name__ == '__main__':

  ## general settings and shortcuts
  # period shortcuts
  H01 = '1979-1980'; H02 = '1979-1981'; H03 = '1979-1982'; H30 = '1979-2009' # for tests 
  H05 = '1979-1984'; H10 = '1979-1989'; H15 = '1979-1994'; H60 = '1949-2009' # historical validation periods
  G10 = '1969-1979'; I10 = '1989-1999'; J10 = '1999-2009' # additional historical periods
  A03 = '2045-2048'; A05 = '2045-2050'; A09 = '2045-2054'; A10 = '2045-2055'; A15 = '2045-2060' # mid-21st century
  B03 = '2085-2088'; B05 = '2085-2900'; B10 = '2085-2095'; B15 = '2085-2100' # late 21st century  
  ltitle = True # plot/figure title
  figtitles = None; comments = None
  subplot = None # subplot layout (or defaults based on number of plots)
  lbackground = True
  lcontour = False # contour or pcolor plot
  shading = 'gouraud' # shading for pixel plot: 'flat' | 'gouraud'
  laddContour = False # add black contour lines
  lframe = True # draw domain boundary
  loutline = True # draw boundaries around valid (non-masked) data
  framewidths = 2.5
  cbn = None # colorbar levels
  figuretype = None
  lsamesize = True
  l3pan = False # make a single column of a three-column figure
  lminor = True # draw minor tick mark labels
  locean = False # mask continent in white and omit country borders
  lstations = False; stations = 'EC'; cluster_symbols = {2:'o',5:'^',8:'s'}; cluster_name = 'cluster_projection'
  cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.iteritems()}
  lbasins = False; basinlist = ('ARB','FRB','GLB'); primary_basins = []; subbasins = {} #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
  lprovinces = True; provlist = ('BC','AB','ON')
  cbo = None # default based on figure type
  resolution = None # only for GPCC (None = default/highest)
  exptitles = None
  grid = None
  lWRFnative = False
  reflist = None # an additional list of experiments, that can be used to compute differences
  refprd = None; refdom = None; refvars = None
  ldiff = False # compute differences
  lfrac = False # compute fraction
  domain = (2,)
  variable_settings = None
  season_settings = None
  aggregation = 'mean'
  level_agg = dict(p=2,s='mean',soil_layers_stag='mean')
    
  # WRF file types
  WRFfiletypes = [] # WRF data source
  WRFfiletypes += ['aux']
  WRFfiletypes += ['hydro']
#   WRFfiletypes += ['lsm']
  WRFfiletypes += ['srfc']
#  WRFfiletypes += ['xtrm']
#   WRFfiletypes += ['plev3d']
  ## select variables and seasons
  variables = [] # variables
#   variables += ['Ts']
#  variables += ['T2']
#  variables += ['Tmean']
#  variables += ['Tmin', 'Tmax']
#   variables += ['MaxPrecip_1d']; aggregation = 'mean'
#   variables += ['MaxPrecip_1d']; aggregation = 'max'
#   variables += ['MaxPreccu_1d']; aggregation = 'max'
#   variables += ['MaxPrecnc_1d']; aggregation = 'max'
  variables += ['dryprec_010']
  variables += ['precip'] 
#  variables += ['precip','pet'] 
#   variables += ['preccu']
#   variables += ['precnc']
#   variables += ['wetfrq']
#   variables += ['cqwu']
#   variables += ['cqwv']
#   variables += ['cqw']
#   variables += ['RH']; level_agg['p'] = 1
#   variables += ['T']; level_agg['p'] = 2
#   variables += ['u']; level_agg['p'] = 1
#   variables += ['Z']; level_agg['p'] = 2; laddContour = True
#   variables += ['waterflx']
#   variables += ['p-et']
#   variables += ['OIPX']
#   variables += ['OrographicIndex']
#  variables += ['Q2']
#   variables += ['evap']
#  variables += ['pet']
#  variables += ['runoff']
#  variables += ['sfroff']
#  variables += ['ugroff']
#   variables += ['snwmlt']
#   variables += ['aSM']
#   variables += ['snow']
#   variables += ['snowh']
#   variables += ['GLW','OLR','qtfx']
#   variables += ['SWDOWN','GLW','OLR']
#   variables += ['hfx','lhfx']
#   variables += ['qtfx','lhfr']
#   variables += ['SST']
#   variables += ['lat2D','lon2D']
  seasons = [] # seasons
#   seasons += ['cold']
#   seasons += ['warm']
#   seasons += ['melt']
#   seasons += ['annual']
#   seasons += ['summer']
#   seasons += ['winter']
#   seasons += ['spring']    
#   seasons += ['fall']
  # special variable/season combinations
#   variables = ['seaice']; seasons = [8] # September seaice
#  variables = ['snowh'];  seasons = [8] # September snow height
#   variables = ['stns']; seasons = ['annual']
#   variables = ['lndcls']; seasons = [''] # static
#   variables = ['zs']; seasons = ['topo']; lcontour = True; lframe = False 
#   WRFfiletypes = ['const'] if grid is None else ['const','srfc'] # static
#   variables = ['zs']; seasons = ['hidef']; WRFfiletypes=['const']; lcontour = True # static
  
  ## case settings
    
  folder = figure_folder
  lpickle = True # load projection from file or recompute
  lprint = True # write plots to disk using case as a name tag
  maptype = 'lcc-grw'; lstations = False; lbasins = True; domain = 2
#   maptype = 'lcc-can'; lstations = False; domain = 1
#   lbasins = True; basinlist = ('ARB','FRB','CRB','NRB','PSB'); lprovinces = False; provlist = ['BC','AB','ON']
#   lbasins = False; basinlist = ('ARB','FRB','GLB'); lprovinces = False; provlist = ['BC','AB','ON']
  lbasins = True; basinlist = ('GLB','GRW'); lprovinces = False; provlist = ['ON']

## PET validation for GRW
  lprint = True; lpickle = False
  explist = ['g-ens']*3+['GPCC']*3
  seasons = [['annual','summer','winter']*2]
  period  = [H15]*3+[H15]*3; grid = 'grw2'; case = 'grw'
  variables = [['dryprec_010']*3+['precip']*3, 'precip']
# # 2-panel map with CESM and Obs: precip and topo
#   lprint = True; lpickle = True; lbackground = True
# #   explist = ['Ens','PRISM']
# #   exptitles = ['CESM ~80 km','Observations']; case = 'resobs'
# #   seasons = ['winter']; variables = ['precip']
#   explist = ['Ens','max-ens']; domain = (1,2); case = 'reswrf'
#   exptitles = ['CESM ~80 km','Topography ~10 km']
#   seasons = ['hidef']; variables = ['zs']; lcontour = True # static
# #   exptitles = ['CESM ~80 km','WRF 10km']
# #   seasons = ['annual']; variables = ['precip']; lcontour = True # static
#   maptype = 'lcc-fine'; period = H15
#   lframe = False; loutline = False; lbasins = False; lprovinces = False

# Precip Extremes in Ensemble Members (PanAm, progrssion)
#   seasons = [('summer',)*3+('winter',)*3]; period = [H15,A15,B15]*2
#   period = [H15,A15,B15]*2; domain = 1
#   explist = ['max-ctrl','max-ctrl-2050','max-ctrl-2100','max-ens','max-ens-2050','max-ens-2100']
#   expstrs = ('WRF-1','Ensemble Mean'); case = 'panam_ens_10'
# #   explist = ['max-ctrl','max-ctrl-2050','max-ctrl-2100','max-ens-A','max-ens-A-2050','max-ens-A-2100']
# #   expstrs = ('WRF-1','WRF-2'); case = 'panam_ens_12'
# #   explist = ['max-ens-B','max-ens-B-2050','max-ens-B-2100','max-ens-C','max-ens-C-2050','max-ens-C-2100']
# #   expstrs = ('WRF-3','WRF-4'); case = 'panam_ens_34'
#   periodstrs = ('Historical','Mid-Century','End-Century'); 
#   exptitles = ['{:s}, {:s}'.format(expstr,prdstr) for expstr in expstrs for prdstr in periodstrs]
#   maptype = 'lcc-can'; lstations = False; lsamesize = False
#   lprovinces = False; provlist = ['BC','AB']
#   lbasins = False; basinlist = ['FRB','ARB']
# #   variables = ['precip']
#   variables = ['MaxPrecip_1d']; aggregation = 'max'; seasons = ['summer']
#   figtitles = ['Maxima of Daily Precipitation Totals in Summer [mm/day]']
  

# map with river basins
#   variables = ['zs']; seasons = ['topo']; lcontour = True; lframe = False 
#   maptype = 'lcc-prairies'; lstations = False; stations = 'EC'
#   period = H15; lWRFnative = True; loutline = True
#   explist = ['max-ctrl']; exptitles = ' '; domain = (1,2)
#   case = 'SSR'; figtitles = 'Basin Outlines and Topography [km]'
#   lbasins = True; basinlist = ('ARB','FRB','SSR')[:]; primary_basins = basinlist; subbasins = {} #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
#   lprovinces = True; provlist = ('BC','AB','SK')

# high resolution map
#   maptype = 'lcc-new'; lstations = True; stations = 'EC'; lbasins = True
#   period = H01; lWRFnative = True; lframe = False; loutline = False
#   explist = ['col1-ctrl']; exptitles = ' '; domain = (2,3)
#   case = 'BCAB_stns'; figtitles = 'EC Stations (BC & Alberta) and Terrain Height [km]'
# map with most of Canada
#   maptype = 'lcc-can'; lstations = False; lbasins = True
#   period = H05; lWRFnative = True; lframe = False; loutline = False
#   explist = ['max-ctrl']; exptitles = ' '; domain = (1,2)
#   case = 'topo_can'

# observations
#   explist = ['PRISM']; maptype = 'lcc-new'; period = H15
#   ldiff = True; reflist = ['Unity']; maptype = 'lcc-small'
#   exptitles = ['Merged Observations (10 km)']
#   case = 'prism'; lsamesize = True; grid = 'arb2_d02'

# # CESM + WRF domain, global
# #   explist = ['Ens']; case = 'cesm'
# #   explist = ['max-ens']; domain = (0,1,); case = 'wrf1'; lbackground = False
#   explist = ['max-ens']; domain = (0,1,2); case = 'wrf2'; lbackground = False
#   exptitles = ''; maptype = 'ortho-NA'; period = H15; lbasins = False; lprovinces = False   
#   lsamesize = True; lcontour = True; lframe = True; loutline = False
# #   exptitles = ['Merged Observations (10 km)']
#   variables = ['Ts']; seasons = ['annual']; WRFfiletypes = ['srfc'] 

# # observations
#   variables = ['precip']; seasons = ['annual']
#   explist = ['Unity']; maptype = 'lcc-bcab'; period = H15
# #   ldiff = True; reflist = ['Unity']; maptype = 'lcc-small'
#   exptitles = 'Annual Total Precipitation [mm/day]'; figtitles = '' # ['Merged Observations (10 km)']
#   case = 'unity'; lsamesize = False; grid = 'arb2_d02'
#   cluster_name = 'cluster_historical'; cluster_symbols = {i:'o' for i in xrange(10)} # '^','s'
#   cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.iteritems()}
#   lbasins = False; lstations = True; stations = 'EC'; case += '_stations'
  
# single panel plot
#   explist = ['max-ens-2100']; maptype = 'lcc-new'; period = B15
#   lfrac = True; reflist = ['max-ens']; refprd = H15
#   case = 'sum'; lsamesize = False; figtitles = ''; primary_basins = ['FRB','ARB']
#   variables = ['aSM']; seasons = ['jas']; exptitles = 'Soil Moisture Change [%]'
#   variable_settings = ['asm_red']

# GPCC stations
#   variables = ['stations']; seasons = ['annual']
#   explist = ['GPCC']; maptype = 'lcc-new'; period = H30
#   exptitles = ['GPCC Station Density']; figtitles = [''] 
#   case = 'gpcc'; lsamesize = False; #grid = 'arb2_d02'
#   variables = ['stations']; seasons = ['annual']; shading = 'flat'

# # ERA-Interim validation
# #   explist = ['max-ens']; seasons = ['annual']; period = H15; domain = 2
#   explist = ['Unity','erai-max','max-ens']*2
#   seasons = [('summer',)*3+('winter',)*3]
#   period = H15; domain = 2
# #   maptype = 'lcc-can'; lstations = False; lbasins = True; domain = 1
# #   lfrac = True; reflist = ['max-ens',]; refprd = H15;
#   case = 'erai'; lsamesize = True

# # differences to historical period
# #   reflist = ['max-ens-2050','max-ens-2100']*2; case = 'max-seaice'; # l3pan = True (left column of 6-panel figure)
# #   explist = [case+'-2050',case+'-2100']*2
# #   seasons = [('summer',)*2+('winter',)*2]; period = [A15,B15]*2
#   seasons = [('summer',)*3+('winter',)*3]; period = B15; refprd = H15
#   exptitles = ['{:s}, CESM', '{:s}, WRF (IC)', '{:s}, WRF (AE)']*2  
#   exptitles = [title.format(season.title())  for season,title in zip(seasons[0],exptitles)]
#   explist = ['Ens-2100','max-ens-2100','ctrl-ens-2100']*2; case = 'xtrm'
#   reflist = ['Ens','max-ens','ctrl-ens']*2; #grid = 'arb2_d02'
# #  exptitles = ['Historical, {:s}', 'Mid-century, {:s}', 'End-century, {:s}']*2
# #  exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons[0])]
#   domain = 2; maptype = 'lcc-new'; lstations = True; lbasins = False
# #   lsamesize = False; cbo = 'horizontal'  
# #   variables = ['precip']; ldiff = True
# #   variables = ['MaxPrecip_1d']; aggregation = 'max'; lfrac = True
# #   variables = ['precip']; lfrac = True 
#   variables = ['T2']; variable_settings = ['T2_prj']; ldiff = True 

# # differences to Obs
# #   reflist = ['max-ens-2050','max-ens-2100']*2; case = 'max-seaice'; # l3pan = True (left column of 6-panel figure)
# #   explist = [case+'-2050',case+'-2100']*2
# #   seasons = [('summer',)*2+('winter',)*2]; period = [A15,B15]*2
# #   exptitles = ['Mid-century, {:s}', 'End-century, {:s}']*2  
#   case = 'mc_val'; explist = ['ctrl-ens','max-ens']*2; reflist = 'Unity'
#   seasons = [('summer',)*2+('winter',)*2]; period = H15
#   domain = 1; maptype = 'lcc-can'; lstations = False; lbasins = True; grid = 'arb2_d01'
# #   lsamesize = False; cbo = 'horizontal'  
#   variables = ['precip']; ldiff = True
# #   variables = ['MaxPrecip_1d']; domain = 2; cbo = 'vertical'; lfrac = True

# # historical state (continental; vertical orientation)
#   refexp = 'max-ctrl'; case = refexp; exptitles = ['Historical, {:s}',]*2; l3pan = True
#   explist = [refexp]*2; seasons = [('summer','winter',)]; period = H15
#   domain = 1; maptype = 'lcc-can'; lstations = False; lbasins = True; primary_basins = ('ARB','FRB','GLB')
#   lsamesize = False; cbo = 'horizontal'; subplot = (2,1) # vertical
# #   lfrac = True; reflist = ['Unity']*2; grid = 'arb2_d01'
# #   variables = ['MaxPrecip_1d']; domain = 2; cbo = 'vertical'
# #   variables = ['Z']; level_agg['p'] = 2; laddContour = True
# #   variables = ['RH']; level_agg['p'] = 1
# #   variables = ['aSM']
# #   variables = ['p-et']; seasons = [('summer','annual',)]  
# #   variables = ['precip']; variable_settings = ['precip_hist']
# #   exptitles = [exptitle.format(season.title()) for exptitle,season in zip(exptitles,seasons[0])]
#   variables = [('cqwu','cqwv')]; vartitles = ['zonal flux','meridional flux']; seasons = ['summer','winter']  
#   exptitles = [exptitle.format(vartitle.title()) for exptitle,vartitle in zip(exptitles,vartitles)]
#   figtitles = ['Water Flux in {:s} [$kg m^{{-1}} s^{{-1}}$]'.format(season.title()) for season in seasons]

# # differences to historical (continental; "panam")
#   domain = 1; maptype = 'lcc-can'; lstations = False; lbasins = True; primary_basins = ('ARB','FRB','GLB')
#   lsamesize = False; cbo = 'horizontal'  
#   refexp = 'max-ens'; case = refexp; l3pan = True # (left column of 6-panel figure)
#   explist = [refexp+'-2050',refexp+'-2100']*2; reflist = [refexp,]
#   seasons = [('summer',)*2+('winter',)*2]; period = [A15,B15]*2; refprd = H15
#   exptitles = ['Mid-century, {:s}', 'End-century, {:s}']*2
# #   variables = ['MaxPrecip_1d']; domain = 2; cbo = 'vertical'; lfrac = True 
# #   variables = ['Z']; level_agg['p'] = 2; laddContour = True; ldiff = True 
# #   variables = ['precip']; variable_settings = ['precip_prj']; lfrac = True #ldiff = True 
# #   variables = ['RH']; level_agg['p'] = 1; ldiff = True
# #   variables = ['evap']; ldiff = True; variable_settings = ['precip_prj']
# #   variables = ['aSM']; lfrac = True
# #   variables = ['p-et']; variable_settings = ['precip_prj']; seasons = [('summer',)*2+('annual',)*2]; ldiff = True 
# #   exptitles = [exptitle.format(season.title()) for exptitle,season in zip(exptitles,seasons[0])]
#   variables = [('cqwu','cqwu','cqwv','cqwv')]; seasons = ['summer','winter']; ldiff = True
#   vartitles = ['zonal flux','zonal flux','meridional flux','meridional flux']
#   exptitles = [exptitle.format(vartitle.title()) for exptitle,vartitle in zip(exptitles,vartitles)]
#   figtitles = ['Water Flux Differences in {:s} [$kg m^{{-1}} s^{{-1}}$]'.format(season.title()) for season in seasons]


# summer and winter progression (continental)
#   refexp = 'ctrl-ens'; case = refexp; cbo = 'horizontal'
#   explist = [refexp,refexp+'-2050',refexp+'-2100']*2
#   seasons = [('summer',)*3+('winter',)*3]; period = [H15,A15,B15]*2
#   exptitles = ['Historical, {:s}', 'Mid-century, {:s}', 'End-century, {:s}']*2
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons[0])]
#   maptype = 'lcc-can'; lstations = False; lbasins = True 
#   domain = 1; lsamesize = True; cbo = 'horizontal' 
#   variables = ['Z']; level_agg['p'] = 2; laddContour = True 
#   variables = ['precip']; variable_settings = ['precip_prj'] 
#   variables = ['RH']; level_agg['p'] = 1 
#   variables = ['cqwu']
#   variables = ['p-et']; variable_settings = ['precip_prj']; seasons = [('summer',)*2+('annual',)*2]; ldiff = True 

# water transport
#   explist = ['max-ens-2050','max-ens-2100']*2
#   seasons = [('summer',)*2+('winter',)*2]; period = [A15,B15]*2
#   lfrac = True; reflist = ['max-ens',]; refprd = H15
#   explist = ['max-ens','max-ens-2050','max-ens-2100']*2
#   explist = ['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100']*2
#   seasons = [('summer',)*3+('winter',)*3]; period = [H15,A15,B15]*2
#   explist = ['gg-ctrl','gg-ctrl-2050','gg-ctrl-2100',
#              'g-ctrl','g-ctrl-2050','g-ctrl-2100',]
#   seasons = [('summer',)*6]; period = [H15,A15,B15]*2
#   explist = ['gg-ctrl-2050','gg-ctrl-2100',
#              'g-ctrl-2050','g-ctrl-2100',]
#   reflist = ['gg-ctrl','gg-ctrl', 'g-ctrl','g-ctrl',]; ldiff = True
#   seasons = [('summer',)*4]; period = [A15,B15]*2; refprd = H15
#   maptype = 'lcc-can'; lstations = False; lbasins = True; domain = 1
#   case = 'water'; lsamesize = True

# wet-day validation
# #   explist = ['max-ens']; seasons = ['annual']; period = H15; domain = 2
# #   explist = ['Unity','max-ctrl-dry','max-ctrl']*2; seasons = [('summer',)*3+('winter',)*3]
# #   variables = [('precip','dryprec')*2]; refvars = [('precip',)*4]
#   explist = ['max-ens']*4; seasons = [('summer',)*2+('winter',)*2]; period = H15
# #   explist = ['old-ctrl','ctrl-1','new-ctrl']*2; seasons = [('summer',)*3+('winter',)*3]
# #   explist = ['old-ctrl','ctrl-1','new-ctrl','new-v361-ctrl']; seasons = ['winter','summer']
#   ldiff = True; reflist = ['Unity']; refprd = None
# #   explist = ['max-ctrl','max-ens']*2; seasons = [('summer',)*2+('winter',)*2]  
# #   period = H15; domain = 2; grid = 'arb2_d02'; 
#   period = H15; domain = [1,2,]*2; grid = ['arb2_d01','arb2_d02']*2
# #   period = H15; domain = 1; grid = 'arb2_d01'
# #   maptype = 'lcc-can'; lstations = False; lbasins = True; basinlist = ['FRB','ARB','GLB']
# #   lfrac = True; reflist = ['max-ens',]; refprd = H15;
#   case = 'wetdays'; lsamesize = True

# ensemble projection
#   seasons = [('summer',)*2+('winter',)*2]; period = [A15,B15]*2
# #   explist = ['phys-ens-2050','phys-ens-2100']*2; reflist = ['phys-ens']; case = 'phys-prj'; period = [A15,B10]*2
# #   explist = ['ctrl-2050','ctrl-2100']*2; reflist = ['ctrl-1']; case = 'ctrl-prj'
#   explist = ['max-ens-2050','max-ens-2100']*2; reflist = ['max-ens']; case = 'ens-prj'
# #   explist = ['max-ctrl-2050','max-ctrl-2100']*2; reflist = ['max-ctrl']; case = 'max-prj'
# #   explist = ['max-ens-A-2050','max-ens-A-2100']*2; reflist = ['max-ens-A']; case = 'ens-A-prj';
# #   explist = ['max-ens-B-2050','max-ens-B-2100']*2; reflist = ['max-ens-B']; case = 'ens-B-prj';
# #   explist = ['max-ens-C-2050','max-ens-C-2100']*2; reflist = ['max-ens-C']; case = 'ens-C-prj';   
#   periodstrs = ('Mid-Century','End-Century')
#   exptitles = ['{:s}, {:s}'.format(season.title(),prdstr) for season in seasons[0][::2] for prdstr in periodstrs]
#   maptype = 'lcc-bcab'; lstations = True; stations = 'EC'; 
#   lprovinces = True; provlist = ['BC','AB'][1:]
#   lbasins = True; lsamesize = False; basinlist = ['FRB','ARB']
#   lfrac = True; refprd = H15
# #   variables = ['precip']
#   variables = ['MaxPrecip_1d']
# #   ldiff = True; refprd = H15
# #   variables = ['SST']; variable_settings = ['T2_prj'] # parallel execution

# surface sensitivity test
#   maptype = 'lcc-intermed'; lstations = False; lbasins = True
#   explist = ['max-grass','max-ens','max-ctrl','max-1deg']; period = H05; domain = 1
#   ldiff = True; reflist = ['Unity']; refprd = None; 
#   grid = ['arb2_d01']*4 
#   case = 'srfc'; lsamesize = True

# # resolution sensitivity test
#   explist = ['max-ctrl','max-ctrl','max-1deg','max-1deg']; period = H15; domain = [1,2]*2
#   ldiff = True; reflist = ['Unity']; refprd = None; 
#   grid = ['arb2_d01','arb2_d02']*2 
#   case = '1deg'; lsamesize = True

# # max sensitivity experiments
#   explist = ['max-nmp','max-ctrl','max-nosub','max-kf','max-hilev','new-ctrl']; period = H15
#   ldiff = True; reflist = ['Unity']; refprd = None; grid = ['arb2_d02']*5+['arb3_d02'] 
#   case = 'max'; lsamesize = True

# # v361 validation
#   explist = ['new-v361','new-ctrl','max-ctrl','erai-v361','erai-v361-noah','erai-max']; period = H10 
#   domain = 1; grid = ['arb3_d01','arb3_d01','arb2_d01']*2
# #   domain = 2; grid = ['arb3_d02','arb3_d02','arb2_d02']*2
#   ldiff = True; reflist = ['Unity']; refprd = None 
#   case = 'v361'; lsamesize = True

# # v361 projection
#   explist = ['new-v361-2050','new-ctrl-2050','max-ctrl-2050','new-v361-2100','new-ctrl-2100','max-ctrl-2100']  
#   ldiff = True; reflist = ['new-v361','new-ctrl','max-ctrl']*2; refprd = H10 
#   domain = 1; grid = ['arb3_d01','arb3_d01','arb2_d01']*2; period = [A10]*3+[B10]*3
# #   domain = 2; grid = ['arb3_d02','arb3_d02','arb2_d02']*2
# #   ldiff = True; reflist = ['Unity']; refprd = None 
#   case = 'p361'; lsamesize = True

# physics ensemble validation
#   explist = ['new-ctrl','old-ctrl','ctrl-1','max-ctrl']; period = H15; domain = 2
#   ldiff = True; reflist = ['Unity']; refprd = None; grid = ['arb3_d02', 'arb2_d02', 'arb2_d02', 'arb2_d02'] 
#   case = 'phys'; lsamesize = True
#   explist = ['erai-v361','erai-v361-noah','erai-3km','erai-max']; period = H05; domain = [2,2,3,2]
#   ldiff = True; reflist = ['Unity']; refprd = None; grid = ['arb3_d02']*2 + ['arb2_d02']*2 
#   case = 'erai'; lsamesize = True

# hires validation
#   maptype = 'lcc-col'; lstations = False; lbasins = True
# #   explist = ['erai-wc2-rocks',]; period = [1,]; domain = [2,]
#   ldiff = True; reflist = ['PCIC']; refprd = None; grid = 'wc2_d02' 
# #   grid = ['wc2_d02','arb2_d02','col2_d02','wc2_d02']
#   explist = ['erai-wc2-2013','erai-max','erai-3km','erai-wc2-2010']
#   period = [1,15,5,1]; domain = [2, 2, 3, 2]
# #   explist = ['erai-wc2-bugaboo','PCIC','erai-3km','erai-wc2-rocks']
# #   period = [1,None,3,1]; domain = [(1,2), None, (2,3), (1,2)]
# #   exptitles = ['CESM (80 km)','Merged Observations (10 km)', 'Outer WRF Domain (30 km)', 'Inner WRF Domain (10 km)']
#   case = 'wc2'; lsamesize = True

# water transport
#   explist = ['max-1deg']; domain = [1,]; period = H15
#   explist = ['max-1deg', 'max-1deg', 'max-ctrl', 'max-ctrl']; period = H15; domain = [1, 2, 1, 2]
#   exptitles = ['CESM (80 km)','WRF Max-1deg (10 km)', 'WRF Max-Ctrl (30 km)', 'WRF Max-Ctrl (10 km)']
#   explist = ['max-1deg', 'max-1deg', 'max-ctrl', 'max-ctrl']; period = H15; domain = [1, 2, 1, 2]
#   exptitles = ['WRF Max-1deg (30 km)','WRF Max-1deg (10 km)', 'WRF Max-Ctrl (30 km)', 'WRF Max-Ctrl (10 km)']
#   case = 'val1deg'; lsamesize = True; # grid = 'arb2_d02'

# Validation: differences to obs (T2, precip, annual, summer, winter)
#   explist = ['max-ens','Ens',]*2; grid = ['arb2_d02','cesm1x1',]*2
#   seasons = ['summer']*2+['winter']*2; period = H15
#   exptitles = ['WRF, 10 km ({:s} Average)','CESM ({:s} Average)']*2
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons)]
#   case = 'val'; reflist = 'Unity'; refprd = H15; lsamesize = True
#   ldiff = True;  variables = ['T2','precip']; seasons = [seasons] # only make one plot with all seasons!
# #   lfrac = True; variables = ['precip']; seasons = [seasons] # only make one plot with all seasons!

# Projection: T2 and pecip diffs
#   explist = ['max-ens-2100','Ens-2100',]*2; period = B15
#   seasons = ['summer']*2+['winter']*2
#   exptitles = ['WRF, 10 km ({:s} Average)','CESM ({:s} Average)']*2
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons)]
#   case = 'prj'; lbasins = True; lsamesize = True
#   reflist = ['max-ens','Ens',]*2; refprd = H15
#   seasons = [seasons] # only make one plot with all seasons!
#   ldiff = True; variables = ['T2']; variable_settings = ['T2_prj'] # parallel execution
#   lfrac = True; variables = ['precip']; variable_settings = ['precip_prj']
#   ldiff = True; variables = ['precip']; variable_settings = ['precip_prj']

# Fig. 2 Annual Precip: as it is
#   explist = ['Ens']; period = H15; lcontour = True
#   explist = ['Ens', 'Unity', 'max-ens', 'max-ens']; period = H15; domain = [None, None, 1, 2]
#   exptitles = ['CESM (80 km)','Merged Observations (10 km)', 'Outer WRF Domain (30 km)', 'Inner WRF Domain (10 km)']
#   case = 'valobs'; lsamesize = False # grid = [None,'arb2_d02',None,None]

# Fig. 3/4 Validation: differences to obs (T2, precip, annual, summer, winter)
#   explist = ['Ens']; period = H15; grid = ['cesm1x1']
#   explist = ['max-ens']*3+['Ens']*3; grid = ['arb2_d02']*3+['cesm1x1']*3
#   seasons = ['annual', 'summer', 'winter']*2; period = H15
#   exptitles = ['WRF, 10 km ({:s} Average)']*3+['CESM ({:s} Average)']*3
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons)]
#   case = 'val'; lsamesize = True; cbo = 'horizontal'
#   ldiff = True; reflist = ['Unity']*6; refprd = H15
#   variables = ['T2','precip']; seasons = [seasons] # only make one plot with all seasons!
#   lfrac = True; reflist = ['Unity']*6; refprd = H15
#   variables = ['precip']; seasons = [seasons] # only make one plot with all seasons!

# # Fig. 3/4 Validation: differences to obs (T2, precip, annual, summer, winter)
#   explist = ['Ens']; period = H15; grid = ['cesm1x1']
#   explist = ['max-ens-A']*3+['Ens']*3; grid = ['arb2_d02']*3+['cesm1x1']*3
#   seasons = ['annual', 'summer', 'winter']*2; period = H15
#   exptitles = ['WRF, 10 km ({:s} Average)']*3+['CESM ({:s} Average)']*3
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons)]
#   case = 'val'; lsamesize = True; cbo = 'horizontal'
# #   ldiff = True; reflist = ['Unity']*6; refprd = H15
# #   variables = ['T2','precip']; seasons = [seasons] # only make one plot with all seasons!
# #   lfrac = True; reflist = ['Unity']*6; refprd = H15
# #   variables = ['precip']; seasons = [seasons] # only make one plot with all seasons!
#   variables = ['MaxPrecip_1d']; seasons = [seasons] # only make one plot with all seasons!

# Fig. 5 Summer Ensemble
# #   explist = ['Unity']; period = H15
#   explist = ['max', 'max-A', 'Unity', 'max-B', 'max-C', 'NARR']; period = H15
#   exptitles = ['WRF-1', 'WRF-2','Merged Observations', 'WRF-3', 'WRF-4', 'NARR (Reanalysis)']
#   case = 'val-ens'; lsamesize = False; #grid = 'arb2_d02'
#   variables = ['precip']; seasons = ['annual','summer','winter']

## ensemble differences
#   explist = ['Unity']; period = H15
#   explist = ['max-2100', 'max-A-2100', 'max-B-2100', 'max-C-2100']
#   reflist = ['max', 'max-A', 'max-B', 'max-C']; case = 'val-ens'
#   exptitles = ['WRF-1', 'WRF-2', 'WRF-3', 'WRF-4']
#   explist = ['Ctrl-1-2100', 'Ctrl-A-2100', 'Ctrl-B-2100', 'Ctrl-C-2100']
#   reflist = ['Ctrl-1', 'Ctrl-A', 'Ctrl-B', 'Ctrl-C']; case = 'val-Ens'
#   period = B15; refprd = H15; lfrac = True; lsamesize = False; grid = 'arb2_d02'
#   variables = ['precip']; seasons = ['annual','summer','winter'][:] 
 

# Fig. 6/7 Climate Change: T2 and pecip diffs
#   explist = ['max-ens-2100']*3+['Ens-2100']*3; period = B15
# #   explist = ['max-ens-2050']*3+['Ens-2050']*3; period = A15
# #   explist = ['max-1deg-2100']*3+['Ens-2100']*3; period = B15
#   seasons = ['annual', 'summer', 'winter']*2; #grid = ['arb2_d02']*3+['cesm1x1']*3
# # #   seasons = ['annual', 'spring', 'fall']*2; period = B15
#   exptitles = ['WRF, 10 km ({:s} Average)']*3+['CESM ({:s} Average)']*3
#   exptitles = [model.format(season.title()) for model,season in zip(exptitles,seasons)]
#   case = 'prj'; lbasins = True; lsamesize = False; cbo = 'horizontal'
#   ldiff = True; reflist = ['max-ens']*3+['Ens']*3; refprd = H15
# #   ldiff = True; reflist = ['max-1deg']*3+['Ens']*3; refprd = H15  
#   seasons = [seasons] # only make one plot with all seasons!
# #   grid = 'arb2_d02' # necessary for properly area averaged stats
# #   variables = ['SST']; locean = True; variable_settings = ['T2_prj'] # parallel execution
# # #   variables = ['T2']; variable_settings = ['T2_prj'] # parallel execution
#   variables = ['precip']; variable_settings = ['precip_prj']
#   ldiff = False; lfrac = True 

# Fig. 8 Net Precip: the hydro plot
#   case = 'hydro'; lsamesize = False; cbo = 'vertical'; ltitle = True
#   variables = ['p-et']; seasons = [['annual', 'summer']]
#   exptitles = [r'Annual Average', r'Summer Average']
# top row
#   figtitles = r'WRF Ensemble Mean Net Precipitation $(P - ET)$' 
#   explist = ['max-ens']*2; period = H15  
# bottom row
#   figtitles = r'Change in Net Precipitation $\Delta(P - ET)$' 
#   explist = ['max-ens-2100']*2; period = B15
#   ldiff = True; reflist = ['max-ens']; refprd = H15

# Fig. 13 (PDO, and now also AMO, because it is wrong)
#   maptype = 'robinson'; lstations = False; lbasins = False; lprovinces=False; lminor = False; locean = True  
#   case = 'cvdp'; lsamesize = False; cbo = 'horizontal'; ltitle = True; seasons = [None]  
#   exptitles = [r'HadISST', r'CESM Ensemble']; subplot = (2,1)
#   variables = ['PDO_eof','AMO_eof',]; explist = ['SST_CVDP','Ctrl-1_CVDP']; period = H15
#   figtitles = ['Pacific Decadal Oscillation SST Pattern', 'Atlantic Multi-Decadal Oscillation Pattern']   
#   exptitles = [r'20th Cent. Reanalysis V2', r'CESM Ensemble']
#   explist = ['PSL_CVDP','Ctrl-1_CVDP']; period = H15
#   variables = ['NAM_eof','NAO_eof','NPO_eof','PNA_eof',]
#   figtitles = ['Northern Annular Mode Pattern', 'North Atlantic Oscillation Pattern',
#                'North Pacific Oscillation Pattern', 'Pacific North America Pattern']

#   case = '3km'; stations = 'cities'
#   maptype = 'lcc-col'; lstations = True; lbasins = True # 'lcc-new'  
#   period = [H01]*4; period[1] = H15 
#   domain = [3,2,1,2]; lbackground = False
#   ldiff = True; reflist = ['Unity']; refprd = H30
#   grid = ['col2_d03','arb2_d02','col2_d01','col2_d02'] 
#   explist = ['max-3km','max-ctrl','max-3km','max-3km'] 
#   exptitles = ['WRF 3km','WRF 10km (15 yrs)','WRF 30km','WRF 10km']

## large map for all domains
#   maptype = 'lcc-large'; figuretype = 'largemap'; loutline = False; lframe = True
#   lstations = False; lbasins = True; lprovinces = False
#   explist = ['max']; exptitles = ' '; domain = (0,1,2); lWRFnative = True; period = H15 
#   case = 'arb2_basins'; basinlist = ('FRB','ARB','CRB','NRB'); primary_basins = ('FRB','ARB')
# #   lprovinces = True; provlist = ('BC','AB') 
#   variables = ['zs']; seasons = ['topo']; WRFfiletypes += ['const']; lcontour = True
#   figtitles = ['Topographic Height [km]' + ' and Domain Outlines' if lframe else '']
## smaller map for western Canada
#   maptype = 'lcc-new'; lstations = False; lbasins = True
#   case = 'arb'; period = None; lWRFnative = True; lframe = False; loutline = False
#   explist = ['columbia']; exptitles = ' '; domain = (2,3)
#   case = 'frb'; basins = ('FRB',)
#   case = 'arb'; basins = ('ARB',)

    
  if not case: raise ValueError, 'Need to define a \'case\' name!'

  # setup projection and map
  mapSetup = getSetup(maptype, lpickle=lpickle, folder=map_folder)
  
  ## load data
  if not lfrac and not ldiff: reflist = None

  if reflist is not None:
    if isinstance(reflist,basestring): reflist = [reflist]
    elif not isinstance(reflist,(list,tuple)): raise TypeError
    if len(explist) > len(reflist):
      if len(reflist) == 1: reflist *= len(explist)  
      else: raise DatasetError 
    lref = True    
  else: lref = False
  lbackground = not lref and lbackground
  exps, axtitles, nexps = loadDatasets(explist, n=None, varlist=variables, titles=exptitles, periods=period, 
                                       domains=domain, grids=grid, resolutions=resolution, 
                                       filetypes=WRFfiletypes, lWRFnative=lWRFnative, ltuple=True, 
                                       lbackground=lbackground, lautoregrid=True,
                                       WRF_exps=WRF_exps, CESM_exps=CESM_exps)
  nlen = len(exps)
#   print exps[-1][-1]
  # load reference list
  if refvars is None: refvars = variables
  if lref:
    if refprd is None: refprd = period    
    if refdom is None: refdom = domain
    refs, a, b = loadDatasets(reflist, n=None, varlist=refvars, titles=None, periods=refprd, 
                              domains=refdom, grids=grid, resolutions=resolution, filetypes=WRFfiletypes, 
                              lWRFnative=lWRFnative, ltuple=True, lbackground=lbackground,
                              WRF_exps=WRF_exps, CESM_exps=CESM_exps)
    # merge lists
    if len(exps) != len(refs): 
      raise DatasetError, 'Experiments and reference list need to have the same length!'
    for i in xrange(len(exps)):
      if not isinstance(exps[i],tuple): raise TypeError 
      if not isinstance(refs[i],tuple): raise TypeError
      if len(exps[i]) != len(refs[i]): 
        DatasetError, 'Experiments and reference tuples need to have the same length!'
      exps[i] = exps[i] + refs[i] # merge lists/tuples
  
  # get figure settings
  subplot = subplot or nlen
  tmp = getFigureSettings(subplot, cbar=True, cbo=cbo, figuretype=figuretype, 
                          sameSize=lsamesize, l3pan=l3pan)
  sf, figformat, margins, caxpos, subplot, figsize, cbo = tmp # distribute output to variables  
  if not ltitle: margins['top'] += 0.05
  
  # get projections settings
  projection, grid, res = mapSetup.getProjectionSettings()
  
  ## loop over variables and seasons
  maps = []; x = []; y = [] # projection objects and coordinate fields (only computed once)
  fn = -1 # figure counter
  N = len(variables)*len(seasons) # number of figures
  comments = checkItemList(comments, N, basestring, default=None)
  figtitles = checkItemList(figtitles, N, basestring, default=None)
  variable_settings = checkItemList(variable_settings, N, basestring, default=None)
  season_settings = checkItemList(season_settings, N, basestring, default=None)
  
#   if figtitles is not None:
#     if not isinstance(figtitles,(tuple,list)): figtitles = (figtitles,)*N
#     elif len(figtitles) != N: raise ValueError

  # start loop
  for varlist,ravlist in zip(variables,refvars):
    
    for sealist in seasons:
      
      # increment counter
      fn += 1
      M = len(exps) # number of panels
      
      # expand variables
      if isinstance(varlist,basestring): varstr = varlist
      elif isinstance(varlist,(list,tuple)):
        if all([var==varlist[0] for var in varlist]): varstr = varlist[0]
        else: varstr = ''.join([s[0] for s in varlist])
      else: varstr = ''
      varlist = checkItemList(varlist, M, basestring)
      ravlist = checkItemList(ravlist, M, basestring)
#       if ldiff or lfrac:
#       else: 
#         varlist = checkItemList(varlist, M, basestring)
      # expand seasons
      if isinstance(sealist,basestring): seastr = '_'+sealist
      elif isinstance(sealist,(list,tuple)): seastr = '_'+''.join([s[0] for s in sealist])
      else: seastr = ''
      sealist = checkItemList(sealist, M, basestring)
      
      # get smart defaults for variables and seasons
      varlist_settings = variable_settings[fn] or varlist[0] # default: first variable
      if season_settings[fn]: sealist_settings = season_settings[fn]
      elif all([sea==sealist[0] for sea in sealist]): sealist_settings = sealist[0]
      else: sealist_settings = ''      
      # get variable properties and additional settings
      tmp = getVariableSettings(varlist_settings, sealist_settings, ldiff=ldiff, lfrac=lfrac)
      clevs, clim, cbl, clbl, cmap, lmskocn, lmsklnd, plottype = tmp

      # assemble filename      
      filename = varstr
      if ldiff: filename += '_diff'
      elif lfrac: filename += '_frac'
      filename += '{:s}_{:s}.{:s}'.format(seastr,case,figformat)
      # assemble plot title
      plat = exps[0][0].variables[varlist[0]].plot
      figtitle = figtitles[fn]
      if figtitle is None:
        figtitle = plottype + ' ' + plat.title
        if lfrac: figtitle += ' Fractions'
        if ldiff: figtitle += ' Differences'
        if plat.units: 
          if lfrac: figtitle += ' [%]'
          else: figtitle += ' [{:s}]'.format(plat.units) 
      if comments[fn]: figtitle += comments[fn]
      
      # feedback
      print('\n\n   ***  %s %s (%s)   ***   \n'%(plottype,plat.title,varstr))
      
      ## compute data
      data = []; lons = []; lats=[]  # list of data and coordinate fields to be plotted 
      # compute average WRF precip            
      print(' - loading data ({0:s})'.format(varstr))
      for var,rav,season,exptpl in zip(varlist,ravlist,sealist,exps):
        lontpl = []; lattpl = []; datatpl = []
        for i,exp in enumerate(exptpl):
          varname = rav if ( lfrac or ldiff ) and i >= len(exptpl)//2 else var 
          if varname not in exp: 
            raise DatasetError, "Variable '{:s}' not found in Dataset '{:s}!".format(varname,exp.name)
          expvar = exp.variables[varname]
          expvar.load()
          #print expvar.name, exp.name, expvar.masked
          #print(exp.name)
          assert expvar.gdal
          # handle dimensions
          if expvar.isProjected: 
            assert (exp.lon2D.ndim == 2) and (exp.lat2D.ndim == 2), 'No coordinate fields found!'
            if not exp.lon2D.data: exp.lon2D.load()
            if not exp.lat2D.data: exp.lat2D.load()
            lon = exp.lon2D.getArray(); lat = exp.lat2D.getArray()          
          else: 
            assert expvar.hasAxis('lon') and expvar.hasAxis('lat'), 'No geographic axes found!'
            lon, lat = np.meshgrid(expvar.lon.getArray(),expvar.lat.getArray())
          lontpl.append(lon); lattpl.append(lat) # append to data list
          # reduce certain dimensions
          for ax,la in level_agg.iteritems():
            if expvar.hasAxis(ax):
              if isinstance(la, basestring): # aggregate over axis 
                if ax == 'p' and la[:4] == 'low_':
                  la = la[4:]
                  expvar = expvar(asVar=True, p=(0,3), lidx=True) # slice axis (default rules)
                expvar = getattr(expvar,la)(axis=ax, asVar=True)
              else: 
                expvar = expvar(asVar=True, **{ax:la}) # slice axis (default rules)
                if ax == 'p' and not figtitle.endswith(' hPa'): figtitle += " at {:3.0f} hPa".format((850,700,500,250,100)[la])
                if ax == 's' and not figtitle.endswith(' soil layer)'): figtitle += " ({:s} soil layer)".format(('1st','2nd','3rd','4th')[la])
          # extract data field
          if expvar.hasAxis('time'):
            method = aggregation if aggregation.isupper() else aggregation.title() 
            vardata = getattr(expvar,'seasonal'+method)(season, asVar=False, lclim=True)            
          else:
            vardata = expvar[:].squeeze()
#           if expvar.masked: vardata.set_fill_value(np.NaN) # fill with NaN
          vardata = vardata.squeeze() # make sure it is 2D
          if expvar.units != expvar.plot.units:
            vardata = vardata * expvar.plot.scalefactor # apply plot unit conversion
          # figure out ocean mask          
          if lmskocn:
            if exp.variables.has_key('landmask') and False:
              vardata[exp.landmask.getArray()] = -2.
            elif exp.variables.has_key('landfrac'): # CESM mostly 
              vardata[exp.landfrac.getArray(unmask=True,fillValue=0)<0.75] = -2. # use land fraction
            elif exp.variables.has_key('lndidx'): 
              mask = exp.lndidx.getArray()
              vardata[mask==16] = -2. # use land use index (ocean)  
              vardata[mask==24] = -2. # use land use index (lake)
            else :
              vardata = maskoceans(lon,lat,vardata,resolution=res,grid=grid)
          # figure out land mask
          if lmsklnd: 
            if exp.variables.has_key('landfrac'): # CESM and CFSR 
              vardata[exp.lnd.getArray(unmask=True,fillValue=0)>0.75] = 0 # use land fraction
            elif exp.variables.has_key('lndidx'): # use land use index (ocean and lake)
              mask = exp.lndidx.getArray(); tmp = vardata.copy(); vardata[:] = 0.
              vardata[mask==16] = tmp[mask==16]; vardata[mask==24] = tmp[mask==24]
          datatpl.append(vardata) # append to data list
        ## compute differences, if desired
        if ldiff or lfrac:
          assert len(datatpl)%2 == 0, 'needs to be divisible by 2'
          ntpl = len(datatpl)/2 # assuming (exp1, exp2, ..., ref1, ref2, ...)
          for i in xrange(ntpl):
            if ldiff: datatpl[i] = datatpl[i] - datatpl[i+ntpl] # compute differences in place
            elif lfrac: datatpl[i] = (datatpl[i]/datatpl[i+ntpl]-1)*100 # compute fractions in place
          del datatpl[ntpl+1:] # delete the rest 
        # add tuples to master list
        lons.append(lontpl); lats.append(lattpl); data.append(datatpl)
      print("")
              
      ## setup projection
      #print(' - setting up figure\n') 
      nax = subplot[0]*subplot[1] # number of panels
      # make figure and axes
      f = plt.figure(facecolor='white', figsize=figsize)
      ax = []
      for n in xrange(nax):
        ax.append(f.add_subplot(subplot[0],subplot[1],n+1, axisbg='blue'))
      f.subplots_adjust(**margins) # hspace, wspace
      if not maps:
        print(' - setting up map projection\n') 
        #mastermap = Basemap(ax=ax[n],**projection) # make just one basemap with dummy axes handle
        mastermap = mapSetup.basemap
        for axi in ax: # replace dummy axes handle with correct axes handle
          tmp = copy(mastermap)
          tmp.ax = axi  
          maps.append(tmp) # one map for each panel!!  
      else:
        print(' - resetting map projection\n') 
        for n in xrange(nax):
          maps[n].ax=ax[n] # assign new axes to old projection
      # transform coordinates (on per-map basis)
      if not (x and y):
        print(' - transforming coordinate fields\n')
        for n in xrange(nax):
          xtpl = []; ytpl = []
          for m in xrange(nexps[n]):
            xx, yy = maps[n](lons[n][m],lats[n][m]) # convert to map-native coordinates
            xtpl.append(xx); ytpl.append(yy)
          x.append(xtpl); y.append(ytpl) 
        
      ## Plot data
      # draw boundaries of inner domain
      if loutline or lframe:
        print(' - drawing domain outlines\n')
        for n in xrange(nax):
          for m in xrange(nexps[n]):   
            if loutline:
              bdy = ma.ones(data[n][m].shape); bdy[ma.getmaskarray(data[n][m])] = 0
              # N.B.: for some reason, using np.ones_like() causes a masked data array to fill with zeros  
              #print bdy.mean(), data[n][m].__class__.__name__, data[n][m].fill_value 
              bdy[0,:]=0; bdy[-1,:]=0; bdy[:,0]=0; bdy[:,-1]=0 # demarcate domain boundaries        
              maps[n].contour(x[n][m],y[n][m],bdy,[-1,0,1],ax=ax[n], colors='k', 
                              linewidths=framewidths, fill=False) # draw boundary of data
            if lframe:
              if isinstance(domain,(tuple,list)) and not ( domain[0] == 0 and m == 0):
                bdy = ma.ones(x[n][m].shape)   
                bdy[0,:]=0; bdy[-1,:]=0; bdy[:,0]=0; bdy[:,-1]=0 # demarcate domain boundaries        
                maps[n].contour(x[n][m],y[n][m],bdy,[-1,0,1],ax=ax[n], colors='k', 
                                linewidths=framewidths, fill=False) # draw boundary of domain
      # draw data
      norm = mpl.colors.Normalize(vmin=min(clevs),vmax=max(clevs),clip=True) # for colormap
      cd = []
#       print(' - creating plots\n')  
      toString = lambda v: v if isinstance(v,basestring) else clbl%v
      for n in xrange(nax): 
        for m in xrange(nexps[n]):
          vmean = toString(np.nanmean(data[n][m])) 
          vmin = toString(np.nanmin(data[n][m]))
          vmax = toString(np.nanmax(data[n][m]))
          if ldiff or lfrac: 
            vrms = toString(np.sqrt(np.nanmean(data[n][m]**2)))
            print('panel {:d}: bias {:s} / rms {:s} / min {:s} / max {:s}'.format(
                            n, vmean, vrms, vmin, vmax))
          else: 
            vstd = toString(np.nanstd(data[n][m]))
            print('panel {:d}: mean {:s} / std {:s} / min {:s} / max {:s}'.format(
                            n, vmean, vstd, vmin, vmax))
          if lcontour: 
            cd.append(maps[n].contourf(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],cmap=cmap, 
                                       norm=norm,extend='both'))  
          else: cd.append(maps[n].pcolormesh(x[n][m],y[n][m],data[n][m], cmap=cmap, 
                                             shading=shading))
          # add black contour lines (outlines) 
          if laddContour:
            cd.append(maps[n].contour(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],
                                      colors='k', linewidths=0.5))
      # add colorbar
      #TODO: use utils.sharedColorbar
      cax = f.add_axes(caxpos)
      for cn in cd: # [c1d1, c1d2, c2d2]:
        if clim: cn.set_clim(vmin=clim[0],vmax=clim[1])
        else: cn.set_clim(vmin=min(clevs),vmax=max(clevs))
      cbar = f.colorbar(cax=cax,mappable=cd[0],orientation=cbo,extend='both') # ,size='3%',pad='2%'       
      if cbl is None:
        if cbn is None:
          if ( cbo == 'horizontal' and subplot[1] == 1 ): cbn = 5
          elif ( cbo == 'vertical' and subplot[0] == 1 ): cbn = 7
          elif ( cbo == 'horizontal' and subplot[1] == 2 ): cbn = 7
          elif ( cbo == 'vertical' and subplot[0] == 2 ): cbn = 9
          else: cbn = 9
        cbl = np.linspace(min(clevs),max(clevs),cbn)
      cbar.set_ticks(cbl); cbar.set_ticklabels([clbl%lev for lev in cbl])
    
      ## Annotation
      #print('\n - annotating plots\n')      
      # add labels
      if ltitle: f.suptitle(figtitle,fontsize=16 if len(maps) == 1 else 12)
      # add a map scale to lower left axes
      msn = len(maps)/2 # place scale 
      mapSetup.drawScale(maps[msn])
      n = -1 # axes counter
      for i in xrange(subplot[0]):
        for j in xrange(subplot[1]):
          n += 1 # count up
          axn = ax[n] 
          axn.set_title(axtitles[n],fontsize=11) # axes title
          if j == 0 : left = True
          else: left = False 
          if i == subplot[0]-1: bottom = True
          else: bottom = False
          # begin annotation
          bmap = maps[n]
          kwargs = dict()
          # white-out continents, if we have no proper land mask 
          if locean or ( lmsklnd and not (exps[n][0].variables.has_key('lndmsk') ) 
                         or exps[n][0].variables.has_key('lndidx')): 
            kwargs['maskland'] = True          
          if ldiff or lfrac or locean: 
            kwargs['ocean_color'] = 'white' ; kwargs['land_color'] = 'white'
          # misc annotations
          mapSetup.miscAnnotation(bmap, **kwargs)
          # add parallels and meridians
          mapSetup.drawGrid(bmap, left=left, bottom=bottom, minor=lminor)
          # mark stations
          if lstations: 
            if stations == 'EC':
              # import station data
              from datasets.EC import loadEC_StnTS, selectStations
              from datasets.WRF import loadWRF_StnTS
              varlist = stn_params + [cluster_name]; station_type = 'ecprecip'
              ecstns = loadEC_StnTS(station=station_type, varlist=varlist)
              wrfstns = loadWRF_StnTS(experiment='max-ctrl', varlist=varlist, station=station_type, 
                                      filetypes='hydro', domains=2)

              ecstns,wrfstns = selectStations([ecstns, wrfstns] , stnaxis='station', linplace=False, lall=True, 
                                              **station_constraints)
              # loop over points
              for lon,lat,zerr,cln in zip(ecstns.stn_lon, ecstns.stn_lat, wrfstns.zs_err, ecstns[cluster_name]):
                if cln in cluster_symbols: # and zerr <= 300
                  tmp = cluster_symbols[cln].copy()
                  xx,yy = bmap(lon, lat)
                  marker = tmp.pop('marker')
                  bmap.plot(xx,yy,marker,**tmp)
            else: mapSetup.markPoints(ax[n], bmap, pointset=stations)     
          # add basin outlines          
          if lbasins:
            for basin in basinlist:      
              basininfo = basins_info[basin]
              try:
                if basin in subbasins:
                  for subbasin in subbasins[basin]:		  
                    bmap.readshapefile(basininfo.shapefiles[subbasin][:-4], subbasin, ax=axn, 
                                       drawbounds=True, linewidth = 0.5, color='k')          
                elif basin in primary_basins:
                  bmap.readshapefile(basininfo.shapefiles['Whole'+basin][:-4], basin, ax=axn, 
                                     drawbounds=True, linewidth = 1., color='k')            
                else:
                  bmap.readshapefile(basininfo.shapefiles['Whole'+basin][:-4], basin, ax=axn, 
                                     drawbounds=True, linewidth = 0.5, color='k')
              except:
                print(basin)
                raise # raise previous exception           
          # add certain provinces
          if lprovinces: 
            for province in provlist:    
              provinfo = province_info[province]
              bmap.readshapefile(provinfo.shapefiles[provinfo.long_name][:-4], province, 
                                 drawbounds=True, linewidth = 0.5, color='k')            

      # save figure to disk
      if lprint:
        print('\nSaving figure in '+filename)
        f.savefig(folder+filename, **sf) # save figure to pdf
        print(folder)
  
  ## show plots after all iterations
  plt.show()

