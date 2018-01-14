'''
Created on 2012-11-05, adapted for GeoPy on 2013-10-10

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
from datasets.common import stn_params
from legacy_plotting.legacy import loadDatasets, checkItemList
# project related stuff
# Western Canada
# from projects.WesternCanada import getSetup, figure_folder, map_folder, getFigureSettings, getVariableSettings
# Great Lakes
from projects.GreatLakes import getSetup, getFigureSettings, getVariableSettings
from projects.GreatLakes import figure_folder, map_folder, WRF_exps, CESM_exps
from projects.GreatLakes import basins, provinces

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
  G10 = '1969-1979'; I10 = '1989-1999'; J10 = '1999-2009'; NRC70 = '1970-2000'; NRC80 = '1980-2010' # additional historical periods
  A03 = '2045-2048'; A05 = '2045-2050'; A09 = '2045-2054'; A10 = '2045-2055'; A15 = '2045-2060' # mid-21st century
  B03 = '2085-2088'; B05 = '2085-2900'; B10 = '2085-2095'; B15 = '2085-2100' # late 21st century  
  ltitle = True # plot/figure title
  figtitles = None; comments = None
  subplot = None # subplot layout (or defaults based on number of plots)
  figsize = None
  lbackground = True
  lcontour = True # contour or pcolor plot
  shading = 'gouraud' # shading for pixel plot: 'flat' | 'gouraud'
  laddContour = False # add black contour lines
  lframe = True # draw domain boundary
  framewidths = 1.5; framecolor = 'k' # line width for domain outline
  loutline = True # draw boundaries around valid (non-masked) data
  outlinewidth = 1.; outlinecolor = 'k'
  cbn = None # colorbar levels
  figuretype = None
  lsamesize = False # same figure size, regardless of panels (square)
  l3pan = False # make a single column of a three-column figure
  lminor = True # draw minor tick mark labels
  locean = False # mask continent in white and omit country borders
  lstations = False; stations = 'EC'; #cluster_symbols = {2:'o',5:'^',8:'s'}; cluster_name = 'cluster_projection'
  #cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.iteritems()}
  lbasins = False; primary_basins = ('GLB',); subbasins = ('GRW','SNW')
  basin_args = dict(linewidth = 1., color='k'); subbasin_args = basin_args.copy()
  lprovinces = True; provlist = ('ON',)
  prov_args = dict(linewidth = 0.5, color='k')
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
  level_agg = dict(i_p=2,i_s='mean')
    
  # WRF file types
  WRFfiletypes = [] # WRF data source
#   WRFfiletypes += ['aux']
#   WRFfiletypes += ['rad']
#   WRFfiletypes += ['hydro']
#   WRFfiletypes += ['lsm']
  WRFfiletypes += ['srfc']
#   WRFfiletypes += ['xtrm']
#   WRFfiletypes += ['plev3d']
  ## select variables and seasons
  variables = [] # variables
#   variables += ['Ts']
#   variables += ['T2']
#   variables += ['Tmean']
#  variables += ['Tmin', 'Tmax']
#   variables += ['MaxPrecip_1d']; aggregation = 'mean'
#   variables += ['MaxPrecip_1d']; aggregation = 'max'
#   variables += ['MaxPreccu_1d']; aggregation = 'max'
#   variables += ['MaxPrecnc_1d']; aggregation = 'max'
#   variables += ['dryprec_010']
#   variables += ['precip'] 
#   variables += ['precip','pet'] 
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
    
  folder = figure_folder; lsamesize = False
  lpickle = True # load projection from file or recompute
  lprint = True # write plots to disk using case as a name tag
#   maptype = 'lcc-grw'; lstations = False; lbasins = True; domain = 2
#   maptype = 'lcc-glb'; lstations = False; lbasins = True; domain = 2
#   maptype = 'lcc-glb'; lstations = False; lbasins = True; domain = None
  maptype = 'lcc-NA'; lstations = False; domain = 1
#   lbasins = True; basinlist = ('ARB','FRB','CRB','NRB','PSB'); lprovinces = False; provlist = ['BC','AB','ON']
#   lbasins = False; basinlist = ('ARB','FRB','GLB'); lprovinces = False; provlist = ['BC','AB','ON']
  lbasins = True; basinlist = ('GLB','GRW'); lprovinces = False; provlist = ['ON']
#   lbasins = True; basinlist = ('GLB',); lprovinces = False; provlist = ['ON']
  isoline = None


## PET validation for GRW
#   lprint = True; lpickle = False
#   explist = ['g-ens']*3+['GPCC']*3
#   seasons = [['annual','summer','winter']*2]
#   period  = [H15]*3+[H15]*3; grid = 'grw2'; case = 'grw'
#   variables = [['dryprec_010']*3+['precip']*3, 'precip']


# # validation and projection for the Great Lakes region
# #   res = None; case = '{RES}'; domain = 1; wrftypes = None
#   maptype = 'lcc-NA'; lstations = False; lbasins = True; basinlist = ('GLB',); case = 'na'
#   domain = 1; grid = 'glb1_d{DOM:02d}'.format(DOM=domain); res = ['10km' if domain == 2 else '30km']*4
# #   explist = ['NRCan','GPCC',]*2; seasons = [['summer']*2+['winter']*2]; domain = None; period = [NRC70,None]*2; case = 'gpcc025'
# #   explist = ['NRCan','GPCC',]*2; seasons = [['summer']*2+['winter']*2]; domain = None; period = NRC70; case = 'gpcc05'
# #   explist = ['NRCan','CRU',]*2; seasons = [['summer']*2+['winter']*2]; domain = None; period = NRC70; case = 'cru'
# #   exptitles = [exptitle+', {S:s}' for exptitle in explist]; #grid = ['glb1_d02','glb1_d02']*2
# #   explist = ['CFSR','erai-g','erai-t',]*2; seasons = [['summer']*3+['winter']*3]; tag = 'd01'; domain = 2
# #   exptitles = ['CFSR','WRF G (30km, ERA-I)','WRF T (30km, ERA-I)',]*2; grid = 'glb1_d{:02d}'.format(domain)
# #   explist = ['g-ens','Ens','g-ens','Ens']; seasons = [['summer']*2+['winter']*2]; domain = 1; case = '{EXP0}_{RES}'
# #   exptitles = ['WRF Ensemble ({RES:s}), {S:s}','CESM Ensemble, {S:s}']*2; grid = ['glb1_d01','glb1_d01']*2
# #   explist = ['Ens','g-ens','t-ens',]*2; seasons = [['summer']*3+['winter']*3]; wrftypes = ['']*6
# #   domain = 1; grid = 'glb1_d{DOM:02d}'.format(DOM=domain); res = ['10km' if domain == 2 else '30km']*6
# #   explist = ['Ens','g3-ens','t3-ens',]*2; seasons = [['summer']*3+['winter']*3]
# #   domain = 1; grid = 'glb1-90km_d01'; res = '90km'
# #   exptitles = ['CESM Ensemble, {S:s}','WRF G Ens. ({RES:s}), {S:s}','WRF T Ens. ({RES:s}), {S:s}',]*2
# #   explist = ['g-ctrl','gg-ctrl','gg2-ctrl']; seasons = [['summer']*len(explist)]; wrftypes = ['']*len(explist)
#   explist = ['g-ens','t-ens']*2; seasons = [['summer']*2+['annual']*2]; wrftypes = ['']*len(explist); case = 'gt-ens'
#   exptitles = ['WRF G Ensemble, {S:s}','WRF T Ensemble, {S:s}',]*2
# #   wrftypes = ['G']*3+['T']*3; seasons = [['summer']*6]; case = '{:s}-ens_res'.format('gt')
# #   wrftypes = ['G']*6; seasons = [['summer']*3+['winter']*3]; case = '{:s}-ens_res'.format(wrftypes[0].lower())
# #   explist = [ wt.lower()+exp for wt,exp in zip(wrftypes,['3-ens','-ens','-ens',]*2) ]
# #   reflist = [ 'g'+exp for exp in ['3-ens','-ens','-ens',]*2 ]; wrftype = 'T/G'; case = 't_over_g'
# #   domain = [1,1,2]*2; grid = ['glb1-90km_d01','glb1_d01','glb1_d02']*2 
# #   res = ['90km','30km','10km']*2; exptitles = ['WRF {TYPE:s} Ens. ({RES:s}), {S:s}',]*6
# #   seasons = [['summer']*3+['annual']*3]
#   exptitles = [ t.format(RES=r,S=s.title(),TYPE=wt.upper()) for t,s,r,wt in zip(exptitles,seasons[0],res,wrftypes) ]
# #   explist = ['g-ens','t-ens','g-ens','t-ens']; seasons = [['summer']*2+['winter']*2]
# #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag; 
# #   explist = ['g-ens','t-ens','g-ens','t-ens']; seasons = [['summer']*2+['annual']*2]
# #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag; 
# #   exptitles = ['G Ensemble','T Ensemble']*2; res = '30km' if domain == 1 else '10km'
# #   explist = ['g-ens','g3-ens','g-ens','g3-ens']; seasons = [['summer']*2+['winter']*2]; tag = 'g3'
# #   exptitles = ['WRF Ensemble (30km)','WRF Ensemble (90km)']*2; grid = ['glb1_d01','glb1-90km_d01']*2
# #   if res is None: res = '30km' if domain == 1 else '10km' if domain == 2 else ''
# #   case = case.format(EXP0=explist[0],EXP1=explist[1],RES=res)
# #   exptitles = [ t.format(RES=res,S=s.title(),TYPE=wrftype) for t,s in zip(exptitles,seasons[0]) ]
# #   variables = ['SWDNB']; cbn = 5; lfrac = True; WRFfiletypes = ['rad']
# #   variables = ['LWDNB']; cbn = 5; lfrac = False; WRFfiletypes = ['rad']
# #   variables = ['ps']; cbn = 5; ldiff = True; WRFfiletypes = ['srfc']; laddContour = False
# #   variables = ['PBLH']; cbn = 5; lfrac = True; WRFfiletypes = ['plev3d']
# #   variables = ['Z']; cbn = 5; ldiff = True; WRFfiletypes = ['plev3d']; level_agg['i_p'] = 2; laddContour = True
# #   variables = ['RH']; cbn = 5; ldiff = True; WRFfiletypes = ['plev3d']; level_agg['p'] = 85000
# #   variables = ['T']; cbn = 5; ldiff = True; WRFfiletypes = ['plev3d']; variable_settings = ['T_prj']; level_agg['p'] = 85000
# #   variables = ['T2']; cbn = 5 # T2 without diffs
# #   variables = ['precip']; cbn = 6 # precip without diffs
#   variables = ['runoff']; cbn = 7; WRFfiletypes = ['lsm','hydro'] 
#   variable_settings = 'runoff_fraction'; refvars = ['precip']
#   lfrac = True; reflist = explist; period = H15; refprd = H15; case += '_frac'
# #   variables = ['preccu']; cbn = 6
# #   variables = ['precnc']; cbn = 6
# #   variables = ['preccu',]; cbn = 6; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['T2']; cbn = 6; ldiff = True; variable_settings = ['T2_prj'] # T2
# #   variables = ['precip']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['evap']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['preccu']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['precnc']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['zs']; seasons = [['hidef']*6]; WRFfiletypes += ['const']; lcontour = True
# #   variables = ['MaxPrecip_1d']; aggregation = 'max'; cbn = 7; lfrac = True; variable_settings = ['MaxPrecip_prj']
# #   variables = ['aSM']; aggregation = 'mean'; cbn = 7; lfrac = True
# #   period = B15; refprd = H15; reflist = reflist or explist; case += '_prj' # projection 
# #   period = H15; refprd = H15; case += '_val'; variable_settings = None; reflist = reflist or 'Unity' # validation
# #   period = H15; refprd = NRC70; variable_settings = None; reflist = 'NRCan' # validation  
# #   case = tag+'val_narr'; reflist = 'NARR'
# #   period = H15; ldiff = lfrac = False; variable_settings = None

# ## comparison of NRCan and GPCC/CRU
#   case = 'nrcan_val'; maptype = 'lcc-can'; grid = 'glb1_d01'; 
#   basinlist = []; lprovinces = True; provlist = ['AB','SK','MB','ON']
# #   seasons = [['summer','winter','spring','fall']]; exptitles = [s.title() for s in seasons[0]]
# #   explist = ['NRCan']*len(exptitles); reflist = ['Unity']; period = H30
# #   explist = ['NRCan']*2; period = H30; seasons = ['annual']
# #   variables = [('precip','pet')]; variable_settings = 'precip'; cbn = 9
# #   figtitles = ['Annual Averages from {} [mm/day]'.format(explist[0])]
# #   exptitles = ['Precipitation','PET']
#   explist = ['CRU']; figtitles = ['{} Precipitation - PET [mm/day]'.format(explist[0])]; case = explist[0].lower() 
#   exptitles = ['']; lsamesize = True; variable_settings = 'precip'; cbn = 9; lfrac = True
#   period = H30; seasons = ['annual']; reflist = explist; variables = ['precip',]; refvars = ['pet'] 
# #   reflist = ['CRU','GPCC']; exptitles = reflist
# #   variables = ['precip']; variable_settings = 'precip_obs'; cbn = 9; lfrac = True
# #   variables = ['precip']; variable_settings = 'precip_obs'; cbn = 7; ldiff = True
# #   variables = ['pet']; variable_settings = 'pet_obs'; cbn = 7; lfrac = True
# #   variables = ['pet']; variable_settings = 'pet_obs'; cbn = 7; ldiff = True
# #   variables = ['T2']; variable_settings = 'T2_obs'; cbn = 11; ldiff = True
# #   reflist = ['CRU']; case += '_cru'; figtitles = ['Relative Differences of Total Perecipitation w.r.t. CRU [%]']
# #   reflist = None; ldiff = False; lfrac = False; variable_settings = None; cbn = 6

# ## GRW maps
#   case = 'GRW'; maptype = 'lcc-grw'; cbo = 'vertical'; lcontour = True
#   seasons = [['November','December','January','February','March','April']]; exptitles = [s.title() for s in seasons[0]]
# #   explist = ['NRCan']*len(exptitles); period = NRC70
#   reflist = ['NRCan']*len(exptitles); refprd = NRC70; grid = 'glb1_d01'; domain = 1
# #   explist = ['erai-t']*len(exptitles); period = H30
#   explist = ['t3-ensemble']*len(exptitles); period = H15
#   variables = ['T2']; isoline = 273.5; cbn = 11; ldiff = True
# #   variables = ['Tslb']; variable_settings = 'T_freeze'; level_agg = dict(i_s=0); isoline = 273.5; cbn = 11
# #   variables = ['snwmlt',]; refvars = ['liqwatflx',]; reflist = explist; lfrac = True; variable_settings = 'negative_fraction'
# #   variables = ['snwmlt']; isoline = 1.; cbn = 11
# #   variables = ['ratio']; isoline = 1.; cbn = 11 
# #   case = 'ephemeral'
# #   case = 'maritime'
# #   case = 'prairies'
 
#   ## PET check
# #   explist = ['NRCan','CRU']; variables = ['pet']; period = [NRC70,H30]
#   explist = ['erai-g','erai-t','g-ens','t-ens',]; variables = ['pet_wrf']; period = H15
#   seasons = ['annual']; grid = 'can1'
#   maptype = 'lcc-can'; lstations = False; domain = 1; case = 'pet'


# # ERA-Interim validation
#   explist = ['erai-v36','erai-g','erai-t',]*2; seasons = [['summer']*3+['winter']*3]
# #   exptitles = ['CFSR','WRF G (30km, ERA-I)','WRF T (30km, ERA-I)',]*2; grid = 'glb1_d{:02d}'.format(domain)
#   domain = 1; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag
# #   explist = ['g-ens','t-ens','g-ens','t-ens']; seasons = [['summer']*2+['annual']*2]
# #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag; 
#   exptitles = ['ERA-I V3.6','ERA-I G', 'ERA-I T']*2; res = '30km' if domain == 1 else '10km'
# #   explist = ['g-ens','g3-ens','g-ens','g3-ens']; seasons = [['summer']*2+['winter']*2]; tag = 'g3'
# #   exptitles = ['WRF Ensemble (30km)','WRF Ensemble (90km)']*2; grid = ['glb1_d01','glb1-90km_d01']*2
#   exptitles = [ '{} ({}), {}'.format(e,res,s.title()) for e,s in zip(exptitles,seasons[0]) ]
# #   variables = ['T2']; cbn = 5; ldiff = True; variable_settings = ['T2_prj'] # T2
# #   variables = ['precip']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
# #   variables = ['MaxPrecip_1d']; aggregation = 'max'; cbn = 7; lfrac = True; variable_settings = ['MaxPrecip_prj']
# #   variables = ['aSM']; aggregation = 'mean'; cbn = 7; lfrac = True
#   period = H10; refprd = H15; case = tag+'val_new'; variable_settings = None; reflist = 'Unity' # validation 
# #   case = tag+'val_narr'; reflist = 'NARR'

# # snow observations
#   explist = ['NRCan']*2; maptype = 'lcc-can'; period = '1980-2010'
#   #ldiff = True; reflist = ['Unity']; 
#   exptitles = ['NRCan Snow 1960-1990','CMC SWE 1998-2015',]
#   case = 'CMC'; lbasins = False; lprovinces = False
#   grid = 'glb1_d01'; case += '_30km'
#   variables = [('snow','snow_CMC',)]; seasons = ['annual']; 
#   aggregation = 'max'; variable_settings = ['snow_max',]*2; case += '_max'
# #   variables = [('liqwatflx','liqwatflx_CMC')]; seasons = ['spring']
# #   seasons = ['summer','fall','winter','spring']

# # snow observations (single-panel)
#   explist = ['NRCan']; maptype = 'lcc-can'; period = '1980-2010'
#   #ldiff = True; reflist = ['Unity']; 
#   exptitles = ['NRCan Snow 1960-1990',]
#   case = 'CMC'; lbasins = False; lprovinces = False
# #   grid = 'glb1_d01'; case += '_30km'
#   variables = ['ratio',]; seasons = ['annual']; 
# #   aggregation = 'max'; variable_settings = ['snow_max',]*2; case += '_max'
# #   variables = [('liqwatflx','liqwatflx_CMC')]; seasons = ['spring']
# #   seasons = ['summer','fall','winter','spring']

# # single-panel validation with larger map
#   lsamesize = True
# #   case = 'ongl'; maptype = 'lcc-ongl'; lstations = False; lbasins = False; lprovinces = False
# #   variables = ['aSM']; WRFfiletypes = ['lsm']; seasons = ['jas']; figtitles = ['Summer Soil Moisture Change [%]']
# #   variables = ['T2',]; WRFfiletypes = ['srfc']; seasons = ['winter']; season_settings = 'annual'
# #   figtitle = '{:s} Surface Air Temperature [K]'
# #   variables = ['precip',]; WRFfiletypes = ['hydro']; seasons = ['winter']; season_settings = 'None'
# #   figtitle = '{:s} Precipitation Rate [mm/day]'
#   variables = ['snow',]; WRFfiletypes = ['hydro']; seasons = ['winter']; #aggregation = 'max'
#   figtitle = '{:s} Snow Water Equivalent [$kg/m^2$]'
# #   variables = ['snwmlt',]; WRFfiletypes = ['hydro']; seasons = ['winter']; season_settings = 'None'
# #   figtitle = '{:s} Snowmelt Rate [mm/day]'
# #   variables = ['Tlake',]; WRFfiletypes = ['srfc']; seasons = ['summer']; figtitles = ['Lake Surface Temperature']; case = 'lake'
# #   variables = ['MaxPrecip_1d']; WRFfiletypes = ['srfc']; seasons = ['annual']; aggregation = 'max'
# #   figtitles = ['Annual Max. Daily Precip. Extremes [%]']; cbn = 7; variable_settings = ['MaxPrecip_prj']
# #   variables = ['MaxPrecip_5d']; WRFfiletypes = ['hydro']; seasons = ['annual']; aggregation = 'max'
# #   figtitles = ['Annual Max. Pendat Precip. Extremes [%]']; cbn = 7; variable_settings = ['MaxPrecip_prj']
#   lWRFnative = True; loutline = False; lframe = True; lcontour = True; period = H15
# #   explist = ['erai-t']; case = 'tval'; figtitles = None; domain = (1,2)
# #   explist = ['g-ens',]; domain = 1; case = 'g-ens'
# #   explist = ['Ens',]; case = 'cesm'
#   explist = ['NRCan',]; period = NRC70; case = 'nrcan'
# #   explist = [('g-ens','g-ens',)]; exptitles = 'WRF Ensemble (10km)'; case = 'prj12'; domain = (1,2)
# #   explist = ['g-ens']; exptitles = 'WRF Ensemble (30km)'; case = 'prj'; domain = 1
# #   explist = ['g3-ens']; exptitles = 'WRF Ensemble (90km)'; case = 'g3prj'; domain = 1
# #   lfrac = True; refprd = H15; period = B15; reflist = explist
# #   seasons = ['summer','fall','winter','spring','annual']; figtitles = [figtitle.format(season.title()) for season in seasons]
#   seasons = ['annual']; figtitles = [figtitle.format(season.title()) for season in seasons]
  
# # validation over Great Lakes region
#   explist = ['g-ctrl','g-ens','g-ens-A','g-ens-B','NARR','g-ens-C']
#   seasons = ['annual','summer','winter']; variables = ['precip']
#   period  = H15; domain = 1; case = 'ens_d0{:d}'.format(domain); grid = 'glb1_d01'
# #   variables = ['precip'];  reflist = 'Unity'; lfrac = True
#   variables = ['T2']; reflist = 'Unity'; refprd = H15; ldiff = True

# # validation over Great Lakes region
#   explist = ['Ens','GPCC','g-ens','NARR']; seasons = ['annual','summer','winter']
#   exptitles = [None, None, 'WRF Ensemble Mean (30km)', None]
#   period  = H15; grid = None; case = 'val'; lframe = False
#   variables = ['precip']; WRFfiletypes = ['srfc']; domain = 1
# #   variables = ['T2']; explist[1] = 'CRU'

# # (topographic) map with river basins (single panel)
# #   lpickle = False; lprint = False
# #   maptype = 'lcc-glb'; case = 'glb'; lcontour = True; loutline = False
#   maptype = 'lcc-grw'; case = 'grw'; lcontour = True; loutline = False
#   lstations = False; lprovinces = True; provlist = ['ON'] 
#   lbasins = True; basinlist = ['GLB','GRW','SNW']; subbasin_args = dict(linewidth = 2., color='k')  
# #   variables = ['precip']; seasons = ['annual']; figtitles = 'Precipitation [mm/day]'
# #   variables = ['pet']; seasons = ['annual']; figtitles = 'Potential Evapotranspiration [mm/day]'
#   variables = ['snwmlt']; seasons = ['annual']; figtitles = 'Snowmelt [mm/day]'
# #   variables = ['stations']; seasons = ['annual']; figtitles = 'Station Density'
# #   explist = ['g-ens']; period = H15; domain = (1,2); lframe = True; lWRFnative = True 
#   explist = ['NRCan']; period = NRC70
# #   explist = ['GPCC']; period = None
#   case = explist[0].lower(); exptitles = ' '; lsamesize = True
# #   variables = ['zs']; seasons = ['hidef']; figtitles = 'Topography [km]'
# #   lsamesize = True; seasons = ['topo']
# #   basinlist = ['GLB','GRW',]; case = 'grw'
# #   basinlist = ['GLB','GRW',]; subbasin_args = dict(linewidth = 1.5, color='w'); case = 'grw_white'
# #   maptype = 'lcc-grw'; basinlist = ['GRW',]; subbasin_args = dict(linewidth = 1.5, color='w'); case = 'grw_local_white'
  
# # larger map with river basins
# #   lpickle = False; lprint = False
#   case = 'ongl'; lbasins = False; figtitles = 'Basin Outline and Topography [km]'
#   variables = ['zs']; seasons = ['topo']; lcontour = True
#   maptype = 'lcc-ongl'; lstations = False; stations = 'EC'
#   period = H15; lWRFnative = True; loutline = False; lframe = True 
#   explist = ['g-ctrl']; exptitles = ' '; domain = (1,2)
#   basinlist = ('GLB',); primary_basins = basinlist; subbasins = {} #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
# #   lbasins = True; basin_args = dict(linewidth = 1.5, color='k'); case = 'ongl_glb'
# #   lbasins = True; basin_args = dict(linewidth = 1.5, color='m'); case = 'ongl_glb_m'
#   lbasins = True; basin_args = dict(linewidth = 1.5, color='w'); case = 'ongl_glb_w'
#   lprovinces = False; provlist = ('ON',); lsamesize = True
  
# continental-scale map with domains
#   lpickle = False; lprint = False
#   variables = ['zs']; seasons = ['topo']; framecolor = 'w' 
  variables = ['T2']; seasons = ['annual']; framecolor = 'k'; variable_settings = 'T2_global'
  lframe = True; framewidths = 2.; lcontour = True; lsamesize = True
  loutline = False; outlinewidth = 1.; outlinecolor = 'k'
  maptype = 'ortho-can'; lstations = False; stations = 'EC'
  period = H15; lWRFnative = True; lbackground = False
  case = 'can_cesm'; #figtitles = 'Topography and Domain Outline [km]'
  explist = [('Ctrl-1',)]; exptitles = ' '; domain = (0,)
#   case = 'can_wrf1'; #figtitles = 'Topography and Domain Outline [km]'
#   explist = [('Ctrl-1','erai-g3','erai-g',)]; exptitles = ' '; domain = (0,1,1,)
#   case = 'can_wrf3'; #figtitles = 'Topography and Domain Outline [km]'
#   explist = [('Ctrl-1','erai-g3','erai-g','erai-g','erai-max')]; exptitles = ' '; domain = (0,1,1,2,2)
  lbasins = False; basinlist = ('GRW',); primary_basins = basinlist; subbasins = {} #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
  lprovinces = False; provlist = ('ON',)

# # larger map with river basins
# #   lpickle = False; lprint = False
#   case = 'obs'; #lbasins = False; figtitles = 'Basin Outline and Topography [km]'
#   variables = ['precip']; seasons = ['annual','summer','winter']; lcontour = True
#   maptype = 'lcc-ongl'; lstations = False; lprovinces = False
#   period = H15; lWRFnative = True; loutline = False; lframe = False 
#   explist = ['GPCC']; exptitles = ' '; grid = None
#   lbasins = True; basinlist = ('GLB',); primary_basins = basinlist; subbasins = {} 

# # CESM + WRF domain, global
# #   lpickle = False; lprint = False
#   explist = ['Ens']; case = 'cesm'
# #   explist = ['g-ens']; domain = (0,1,); case = 'wrf1'; lbackground = False
# #   explist = ['g-ens']; domain = (0,1,2); case = 'wrf2'; lbackground = False
#   exptitles = ''; maptype = 'ortho-can'; period = H15; lbasins = False; lprovinces = False   
#   lsamesize = True; lcontour = True; lframe = True; loutline = False
# #   exptitles = ['Merged Observations (10 km)']
#   variables = ['Ts']; seasons = ['annual']; WRFfiletypes = ['srfc'] 

# # observations and models side-by-side
#   seasons = [['summer']*3+['winter']*3]; season_settings = ['annual']
#   variables = ['precip',]; figtitles = ['Average Total Precipitation [mm/day]']
# #   variables = ['T2',]; figtitles = ['Average 2m Air Temperature [K]']
#   explist = ['NRCan','Ens','g-ens']*2; domain = 1
#   exptitles = ['NRCan Observations', 'CESM Ensemble', 'G Ensemble (30 km)']*2
#   period = [NRC70,H15,H15]*2; case = 'nrcan'
#   exptitles = ['{}, {}'.format(exp,season.title()) for exp,season in zip(exptitles,seasons[0])]
#   lsamesize = False; lbasins = True; lstations = False; lcontour = False

# GPCC stations
#   variables = ['stations']; seasons = ['annual']
#   explist = ['GPCC']; maptype = 'lcc-new'; period = H30
#   exptitles = ['GPCC Station Density']; figtitles = [''] 
#   case = 'gpcc'; lsamesize = False; #grid = 'arb2_d02'
#   variables = ['stations']; seasons = ['annual']; shading = 'flat'

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

# # large map for all domains
#   maptype = 'lcc-can'; figuretype = 'largemap'; loutline = False; lframe = True
#   lstations = False; lbasins = True; lprovinces = False
#   explist = ['g-ctrl']; exptitles = ' '; domain = (0,1,2); lWRFnative = True; period = H15 
#   case = 'arb2_basins'; basinlist = ('FRB','ARB','CRB','NRB'); primary_basins = ('FRB','ARB')
# # lprovinces = True; provlist = ('BC','AB') 
#   variables = ['zs']; seasons = ['topo']; WRFfiletypes += ['const']; lcontour = True
#   figtitles = ['Topographic Height [km]' + ' and Domain Outlines' if lframe else '']
    
    
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
  sf, figformat, margins, caxpos, subplot, fs, cbo = tmp # distribute output to variables  
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
        # add aggregation
        if aggregation.lower() not in figtitle.lower():
          if aggregation.upper() == 'SEM': agg_str = "SEM of "  
          elif aggregation.lower() == 'mean': agg_str = "Average "
          else: agg_str = aggregation.title()+'. '
          figtitle = agg_str + figtitle
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
                if ax == 'i_p' and not figtitle.endswith(' hPa'): figtitle += " at {:3.0f} hPa".format((850,700,500,250,100)[la])
                if ax == 'i_s' and not figtitle.endswith(' soil layer)'): figtitle += " ({:s} soil layer)".format(('1st','2nd','3rd','4th')[la])
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
      f = plt.figure(facecolor='white', figsize=fs if figsize is None else figsize)
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
              maps[n].contour(x[n][m],y[n][m],bdy,[-1,0,1],ax=ax[n], colors=outlinecolor,
                              linewidths=outlinewidth, fill=False) # draw boundary of data
            if lframe:
              if isinstance(domain,(tuple,list)) and not ( domain[0] == 0 and m == 0):
                bdy = ma.ones(x[n][m].shape)   
                bdy[0,:]=0; bdy[-1,:]=0; bdy[:,0]=0; bdy[:,-1]=0 # demarcate domain boundaries        
                maps[n].contour(x[n][m],y[n][m],bdy,[-1,0,1],ax=ax[n], colors=framecolor, 
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
          if isoline is not None:
            cd.append(maps[n].contour(x[n][m],y[n][m],data[n][m],[isoline],ax=ax[n],
                                      colors='w', linewidths=1.))
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
          elif ( cbo == 'horizontal' and subplot[1] == 2 ): cbn = 5
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
              basininfo = basins[basin]
              try:
                if basin in subbasins:
                  bmap.readshapefile(basininfo.shapefiles[basininfo.outline][:-4], basin, ax=axn, 
                                     drawbounds=True, **subbasin_args)
#                   for subbasin in subbasins[basin]:		  
#                     bmap.readshapefile(basininfo.shapefiles[subbasin][:-4], subbasin, ax=axn, 
#                                        drawbounds=True, **subbasin_args)          
                elif basin in primary_basins:
                  bmap.readshapefile(basininfo.shapefiles[basininfo.outline][:-4], basin, ax=axn, 
                                     drawbounds=True, **basin_args)            
                else:
                  bmap.readshapefile(basininfo.shapefiles[basininfo.outline][:-4], basin, ax=axn, 
                                     drawbounds=True, linewidth = 1., color='k')
              except:
                print(basin)
                raise # raise previous exception           
          # add certain provinces
          if lprovinces: 
            for province in provlist:    
              provinfo = provinces[province]
              bmap.readshapefile(provinfo.shapefile[:-4], province, drawbounds=True, **prov_args)            

      # save figure to disk
      if lprint:
        print('\nSaving figure in '+filename)
        f.savefig(folder+filename, **sf) # save figure to pdf
        print(folder)
  
  ## show plots after all iterations
  plt.show()

