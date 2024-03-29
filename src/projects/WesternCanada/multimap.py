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

print("Importing 'pyplot' from 'matplotlib'\n")

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
from projects.WesternCanada import getSetup, getFigureSettings, getVariableSettings
from projects.WesternCanada import figure_folder, map_folder, WRF_exps, CESM_exps
from projects.WesternCanada import basins, provinces

station_constraints = dict()
station_constraints['min_len'] = 15 # for valid climatology
station_constraints['lat'] = (40,55)
station_constraints['end_after'] = 1980
#station_constraints['max_zerr'] = 100 # can't use this, because we are loading EC data separately from WRF
station_constraints['prov'] = ('BC','AB')

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
    lbackground = True
    lcontour = False # contour or pcolor plot
    shading = 'gouraud' # shading for pixel plot: 'flat' | 'gouraud'
    laddContour = False # add black contour lines
    lframe = True # draw domain boundary
    framewidths = 1.; framecolor = 'k' # line width for domain outline
    loutline = True # draw boundaries around valid (non-masked) data
    outlinewidth = 1.; outlinecolor = 'k'
    cbn = None # colorbar levels
    figuretype = None
    lsamesize = True
    l3pan = False # make a single column of a three-column figure
    lminor = True # draw minor tick mark labels
    locean = False # mask continent in white and omit country borders
    lstations = False; stations = 'EC'; cluster_symbols = {2:'d',5:'^',8:'s',-1:'o'}; cluster_name = 'cluster_projection'
    cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.items()}
    lbasins = False; primary_basins = ('ARB','FRB','GLB'); subbasins = ('WholeARB','UpperARB','LowerCentralARB',) #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
    basin_args = dict(linewidth = 1., color='k'); subbasin_args = dict(linewidth = 0.75, color='k')
    lprovinces = True; provlist = ('BC','AB','ON')
    prov_args = dict(linewidth = 1., color='k')
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
    isoline = None
    variable_settings = None
    season_settings = None
    aggregation = 'mean'
    level_agg = dict(p=2,i_p=2,s='mean',i_s='mean',soil_layers_stag='mean')

    # WRF file types
    WRFfiletypes = [] # WRF data source
  #   WRFfiletypes += ['aux']
    WRFfiletypes += ['hydro']
    WRFfiletypes += ['lsm']
    WRFfiletypes += ['srfc']
    WRFfiletypes += ['xtrm']
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

    folder = figure_folder
    lpickle = True # load projection from file or recompute
    lprint = True # write plots to disk using case as a name tag
    maptype = 'lcc-new'; lstations = False; lbasins = True; domain = None
  #   maptype = 'lcc-can'; lstations = False; domain = 1
  #   lbasins = True; basinlist = ('ARB','FRB','CRB','NRB','PSB'); lprovinces = False; provlist = ['BC','AB','ON']
  #   lbasins = False; basinlist = ('ARB','FRB','GLB'); lprovinces = False; provlist = ['BC','AB','ON']
    lbasins = True; basinlist = ('ARB',); lprovinces = True; provlist = ['AB']; l60 = True
    unity_grid = 'arb2_d02'

    # station markers
    lstations = False; stations = 'EC';
    cluster_symbols = {2:'d',5:'^',8:'s',-1:'o'}; cluster_name = 'cluster_projection'
    cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.items()}

    ## Columbia Ice Field analysis
    # maptype = 'lcc-col_out'; lbasins = False; lprovinces = False
    # lcontour = False
    # variables = ['zs']; seasons = ['topo']; lstations = True; case = 'stations'
    # explist = ['erai-wc2']; period = '2010-2016'
    # domain = (1,2); lframe = True


  # ## ARB maps
  #   # map setup
  #   maptype = 'lcc-arb'; basinlist = ('ARB',); lbasins = True
  #   # figure settings
  #   cbo = 'vertical'; lcontour = False; shading = 'flat'
  # #   seasons = [['November','December','January','February','March','April']]; exptitles = [s.title() for s in seasons[0]]
  #   explist = ['NRCan']; period = NRC70
  # #   reflist = ['NRCan']*len(exptitles); refprd = NRC70; domain = 1
  # #   explist = ['max-ens']*len(exptitles); period = H15; grid = 'arb2_d01'
  # #   explist = ['new-ctrl']*len(exptitles); period = H15; grid = 'arb3_d01'
  #   case = 'forcing'
  # #   variables = ['snow']; aggregation = 'max'; seasons = ['annual']; cbn = 5
  #   variables = ['precip']; seasons = ['annual']; cbn = 5; variable_settings = 'pet'
  # #   variables = ['pet']; seasons = ['annual']; cbn = 5
  # #   variables = ['T2']; isoline = 273.5; cbn = 11; ldiff = True
  # #   variables = ['Tslb']; variable_settings = 'T_freeze'; level_agg = dict(i_s=0); isoline = 273.5; cbn = 11
  # #   variables = ['snwmlt',]; refvars = ['liqwatflx',]; reflist = explist; lfrac = True; variable_settings = 'negative_fraction'
  # #   variables = ['snwmlt']; isoline = 1.; cbn = 11
  # #   variables = ['ratio']; isoline = 1.; cbn = 11
  # #   if variables[0] == 'Tslb': case += '_{:d}'.format(level_agg['i_s'])

  # ## comparison of different observations over ARB
  #   maptype = 'lcc-arb'; lstations = False; lbasins = True; lprovinces = True;
  #   case = 'obsval'; figtitles = [None]; cbo = 'horizontal'; lcontour = False; loutline = False
  # #  explist = ['PCIC']; figtitles = ['Precipitation Climatology (PRISM & GPCC)']; exptitles = ['']; period = None
  #   explist = ['CRU','NRCan','max-ens','GPCC','Unity','ctrl-ens']; domain = 2
  #   period = [H30,NRC70,H15,H30,H30,H15]; unity_grid = 'arb2_d02'
  # #   exptitles = ['WRF 3km (15yr, ERA-I)','PCIC PRISM','WRF 10km (15yr, ERA-I)','WRF 30km (15yr, ERA-I)']; case = 'valera'
  #   variables = ['precip']; seasons = ['annual',]; variable_settings = 'preccu'
  # #   explist = ['CRU','NRCan',]; period = [H30,NRC70,]; unity_grid = 'arb2_d02'
  # #   variables = ['pet']; seasons = ['annual']; #variable_settings = 'preccu'
  # #   variables = ['T2']; WRFfiletypes = ['srfc']; seasons = ['annual']; lcontour = False; loutline = False
  # #  #explist = ['max-3km']; domain = 3; figtitles = ['WRF 3km (1979-1992, CESM)']; exptitles = ['']; period = H12; case = 'valmax'
  # #   seasons = ['annual','winter','summer','fall','spring']


  ## validation and projection for the Paper
  #  maptype = 'lcc-arb'; basinlist = ('ARB',); lbasins = True; lstations = False; lprovinces = False
    maptype = 'lcc-new'; basinlist = ('ARB',); lbasins = True; lstations = False; lprovinces = True; provlist = ['AB']
    explist = ['Ens','max-ens','ctrl-ens']*2; seasons = [['summer']*3+['winter']*3]
    domain = 2; tag = 'd{:02d}'.format(domain); grid = 2*( ['cesm1x1']+2*['arb2_'+tag] ) # 'arb2_'+tag;
    grid = 'arb2_'+tag # for bias estimates across the inner domain
  #   seasons = [['summer']*3+['annual']*3]
  #   seasons = [['spring']*3+['fall']*3]
  #   seasons = ['annual']; # grid = 'arb2_d02'; # to get stats
  #   explist = ['g-ens','t-ens','g-ens','t-ens']; seasons = [['summer']*2+['annual']*2]
  #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag;
  #   exptitles = ['CESM','1st WRF Ens.','Alt. WRF Ens.']*2; res = '30km' if domain == 1 else '10km'
    exptitles = ['CESM','WRF G','WRF T']*2; res = '30km' if domain == 1 else '10km'
  #   explist = ['g-ens','g3-ens','g-ens','g3-ens']; seasons = [['summer']*2+['winter']*2]; tag = 'g3'
  #   exptitles = ['WRF Ensemble (30km)','WRF Ensemble (90km)']*2; grid = ['glb1_d01','glb1-90km_d01']*2
    exptitles = [ 'CESM Ensemble, {}'.format(s.title()) if e == 'CESM' else '{} ({}), {}'.format(e,res,s.title()) for e,s in zip(exptitles,seasons[0]) ]
  #   variables = ['SWDNB']; cbn = 5; lfrac = True
  #   variables = ['LWDNB']; cbn = 5; lfrac = False
    # variables = ['T2']; cbn = 5; ldiff = True; variable_settings = ['T2_prj'] # T2
    variables = ['precip']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
  #   variables = ['MaxPrecip_1d']; aggregation = 'max'; cbn = 7; lfrac = True; variable_settings = ['MaxPrecip_prj']
  #   variables = ['aSM']; aggregation = 'mean'; cbn = 7; lfrac = True
    period = B15; refprd = H15; reflist = explist; case = tag+'prj' # projection
  #   period = A15; refprd = H15; reflist = explist; case = tag+'prjA' # projection
    # period = H15; case = tag+'val'; variable_settings = None; refprd = H15; reflist = 'Unity' # validation
  #   refprd = NRC70; reflist = 'NRCan'; case += '_nrcan' # validation using NRCan instead of Unity
  ##   case = tag+'val_narr'; reflist = 'NARR'

  # ## 2x2-panel projection for the Paper
  #   maptype = 'lcc-arb'; basinlist = ('ARB',); lbasins = True
  #   domain = 2; tag = 'd{:02d}'.format(domain); res = '30km' if domain == 1 else '10km'
  # #   grid = 2*( ['cesm1x1']+['arb2_'+tag] ); explist = ['Ens','max-ens']*2
  # #   exptitles = ['CESM Ens.','WRF G Ens. ({})'.format(res)]*2
  #   grid = 'arb2_'+tag; explist = ['max-ens','ctrl-ens']*2
  #   exptitles = ['WRF G Ens. ({})'.format(res),'WRF T Ens. ({})'.format(res)]*2
  #   seasons = [['summer']*2+['winter']*2]
  #   exptitles = [ '{}, {}'.format(e,s.title()) for e,s in zip(exptitles,seasons[0]) ]
  # #   variables = [('preccu','preccu','precip','precip')]; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
  # #   seasons = ['annual']; subtitles = ['Cu Precip']*2+['Total Precip']*2
  # #   exptitles = [ '{}, {}'.format(e,s.title()) for e,s in zip(exptitles,subtitles) ]
  # #  variables = ['OI']; lfrac = True
  # #  variables = ['U']; lfrac = True
  # #  variables = ['COIP']; lfrac = True
  # #  WRFfiletypes = ['plev3d']; level_agg = dict(p=0)
  # #   variables = [['aSM']*2+['p-et']*2]; aggregation = 'mean'
  # #   WRFfiletypes = ['lsm','hydro']; cbn = 7; lfrac = True
  # #   variables = ['snow']; cbn = 5; lfrac = True; #variable_settings = ['T2_prj'] # T2
  # #   variables = ['T2']; cbn = 5; ldiff = True; variable_settings = ['T2_prj'] # T2
  # #   variables = ['precip']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
  #   period = B15; refprd = H15; reflist = explist; case = tag+'prj' # projection
  # #   period = H15; refprd = NRC70; reflist = 'NRCan'; case = tag+'_nrcan'; variable_settings = None # validation using NRCan instead of Unity
  # #   period = H15; refprd = H15; case = tag+'val'; variable_settings = None; reflist = 'Unity' # validation

  # ## 2x1-panel projection for the Paper
  #   explist = ['max-ens','ctrl-ens']; cbo = 'vertical'; lsamesize = False
  #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'arb2_'+tag;
  # #   seasons = ['annual']; figtitles = ['Annual Water Supply Change [mm/day]']; ldiff = True
  #   variables = ['p-et']; WRFfiletypes = ['hydro']
  # #   seasons = ['JAS']; figtitles = ['Late Summer Net Precipitation Change [%]']; lfrac = True
  #   seasons = ['JAS']; figtitles = ['Late Summer Net Precipitation Change [mm/day]']; ldiff = True
  # #   seasons = ['JAS']; variables = ['aSM']; WRFfiletypes = ['lsm']; cbn = 7; lfrac = True
  # #  figtitles = ['Late Summer Soil Moisture Change [%]']
  #   period = B15; refprd = H15; reflist = explist; case = tag+'prj' # projection
  #   exptitles = ['1st WRF Ens.','Alt. WRF Ens.']; res = '30km' if domain == 1 else '10km'
  #   exptitles = [ '{} ({})'.format(e,res) for e in exptitles ]
  #   lbasins = True; basinlist = ['FRB','ARB']; lprovinces = True; provlist = ['AB']

  # ## single-panel Observations
  #   maptype = 'lcc-arb3_d03'; lstations = False; lbasins = False; lprovinces = True
  #   explist = ['PCIC']; figtitles = ['Precipitation Climatology (PRISM & GPCC)']; exptitles = ['']; period = None
  # #   explist = ['3km-ens']; domain = 3; figtitles = ['WRF 3km (30 years)']; exptitles = ['']; period = H15
  # #   explist = ['erai-3km']; domain = 3; figtitles = ['WRF 3km (1979-1994, ERA-I)']; exptitles = ['']; period = H15
  # #   explist = ['Ens']; figtitles = ['CESM (100km, 60yr)']; exptitles = ['']; period = H15
  #   #explist = ['max-3km']; domain = 3; figtitles = ['WRF 3km (1979-1992, CESM)']; exptitles = ['']; period = H12
  #   variables = ['precip']; WRFfiletypes = ['hydro']; seasons = ['annual']
  # #   variables = ['MaxPrecip_1d']; aggregation = 'max'; WRFfiletypes = ['hydro']; lcontour = True
  #   grid = None; loutline = False; lcontour = False; case = explist[0].lower() + ( '_d0{:d}'.format(domain) if domain else '' )
  #   seasons = ['annual','winter','summer','fall','spring'][:1]
  #   #loutline = True; case += '_outline'; lcontour = False

  # ## four-panel validation
  #   maptype = 'lcc-arb3_d03'; lstations = False; lbasins = False; lprovinces = True; figtitles = [None]; cbo = 'horizontal'
  # #  explist = ['PCIC']; figtitles = ['Precipitation Climatology (PRISM & GPCC)']; exptitles = ['']; period = None
  #   explist = ['erai-3km','PCIC','erai-3km','erai-3km']; domain = [3,None,2,1]; period = H15
  #   exptitles = ['WRF 3km (15yr, ERA-I)','PCIC PRISM','WRF 10km (15yr, ERA-I)','WRF 30km (15yr, ERA-I)']; case = 'valera'
  # #   explist = ['3km-ens','NRCan','3km-ens','3km-ens']; domain = [3,None,2,1]; period = [H15,NRC70,H15,H15]
  # #   exptitles = ['WRF 3km (30yr)','NRCan','WRF 10km (30yr)','WRF 30km (30yr)']; case = 'valens'
  # #   variables = ['snowh']; WRFfiletypes = ['srfc']; seasons = ['winter']; lcontour = False; loutline = False
  # #   explist = ['3km-ens','PCIC','3km-ens','3km-ens']; domain = [3,None,2,1]; period = H15; #grid = [None,'arb3_d03',None,None]
  # #   exptitles = ['WRF 3km (30yr)','PCIC PRISM','WRF 10km (30yr)','WRF 30km (30yr)']; case = 'valens'
  #   variables = ['precip']; WRFfiletypes = ['hydro']; seasons = ['annual']; lcontour = False; loutline = False
  # #   variables = ['T2']; WRFfiletypes = ['srfc']; seasons = ['annual']; lcontour = False; loutline = False
  # #  #explist = ['max-3km']; domain = 3; figtitles = ['WRF 3km (1979-1992, CESM)']; exptitles = ['']; period = H12; case = 'valmax'
  # #   seasons = ['annual','winter','summer','fall','spring']

  # # four-panel validation
  #   maptype = 'lcc-arb3_d03'; lstations = False; lbasins = False; lprovinces = True;
  #   lcontour = False; loutline = False; cbo = 'horizontal'
  # #   explist = ['3km-ens','Ens','3km-ens','3km-ens']; domain = [3,None,2,1]; period = H15; #grid = [None,'arb3_d03',None,None]
  # #   exptitles = ['WRF 3km (30yr)','CESM (100km, 60yr)','WRF 10km (30yr)','WRF 30km (30yr)']; case = 'valens'
  #   explist = ['erai-3km','Ens','erai-3km','erai-3km']; domain = [3,None,2,1]; period = H15; #grid = [None,'arb3_d03',None,None]
  #   exptitles = ['WRF 3km (15yr, ERA-I)','CESM (100km, 60yr)','WRF 10km (15yr, ERA-I)','WRF 30km (15yr, ERA-I)']; case = 'valerai'
  # # explist = ['max-3km']; domain = 3; figtitles = ['WRF 3km (1979-1992, CESM)']; exptitles = ['']; period = H12
  # #   variables = ['precip']; WRFfiletypes = ['hydro']; seasons = ['annual','winter','summer','fall','spring'][:1]
  #   variables = ['MaxPrecip_1d']; aggregation = 'max'; WRFfiletypes = ['hydro']; seasons = ['annual',]; lcontour = False; loutline = False
  # #   lfrac = True; reflist = 'Unity'; grid = 'arb3_d03'; # grid = ['arb3_d03','cesm1x1','arb3_d02','arb3_d01',]
  # #   variable_settings = 'precip_red'; cbo = 'horizontal'

  # # ERA-Interim validation
  #   explist = ['erai-v361','erai-max','erai-v361-noah',]*2; seasons = [['summer']*3+['winter']*3]
  # #   exptitles = ['CFSR','WRF G (30km, ERA-I)','WRF T (30km, ERA-I)',]*2; grid = 'glb1_d{:02d}'.format(domain)
  #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'arb2_'+tag
  # #   explist = ['g-ens','t-ens','g-ens','t-ens']; seasons = [['summer']*2+['annual']*2]
  # #   domain = 2; tag = 'd{:02d}'.format(domain); grid = 'glb1_'+tag;
  #   exptitles = ['ERA-I V3.6','ERA-I G', 'ERA-I T']*2; res = '30km' if domain == 1 else '10km'
  # #   explist = ['g-ens','g3-ens','g-ens','g3-ens']; seasons = [['summer']*2+['winter']*2]; tag = 'g3'
  # #   exptitles = ['WRF Ensemble (30km)','WRF Ensemble (90km)']*2; grid = ['glb1_d01','glb1-90km_d01']*2
  #   exptitles = [ '{} ({}), {}'.format(e,res,s.title()) for e,s in zip(exptitles,seasons[0]) ]
  #   variables = ['T2']; cbn = 5; ldiff = True; variable_settings = ['T2_prj'] # T2
  # #   variables = ['precip']; cbn = 7; lfrac = True; variable_settings = ['precip_prj'] # precip
  # #   variables = ['MaxPrecip_1d']; aggregation = 'max'; cbn = 7; lfrac = True; variable_settings = ['MaxPrecip_prj']
  # #   variables = ['aSM']; aggregation = 'mean'; cbn = 7; lfrac = True
  #   period = H10; refprd = H10; case = tag+'val_new'; variable_settings = None; reflist = 'Unity' # validation
  # #   case = tag+'val_narr'; reflist = 'NARR'


  # single-panel validation with larger map
  # #   case = 'ongl'; maptype = 'lcc-ongl'; lstations = False; lbasins = False; lprovinces = False
  # #   variables = ['aSM']; WRFfiletypes = ['lsm']; seasons = ['jas']; figtitles = ['Summer Soil Moisture Change [%]']
  #   variables = ['Tlake',]; WRFfiletypes = ['srfc']; seasons = ['summer']; figtitles = ['Lake Surface']
  # #   variables = ['MaxPrecip_1d']; WRFfiletypes = ['srfc']; seasons = ['annual']; aggregation = 'max'
  # #   figtitles = ['Annual Max. Daily Precip. Extremes [%]']; cbn = 7; variable_settings = ['MaxPrecip_prj']
  # #   variables = ['MaxPrecip_5d']; WRFfiletypes = ['hydro']; seasons = ['annual']; aggregation = 'max'
  # #   figtitles = ['Annual Max. Pendat Precip. Extremes [%]']; cbn = 7; variable_settings = ['MaxPrecip_prj']
  #   lWRFnative = True; loutline = False; lframe = True; lcontour = True; period = H15
  # #   explist = ['erai-t']; case = 'tval'; figtitles = None; domain = (1,2)
  # #   explist = ['erai-t']; case = 'tval'; figtitles = None; domain = (1,2)
  #   explist = ['g-ens',]; case = 'lake'; domain = 1
  # #   explist = [('g-ens','g-ens',)]; exptitles = 'WRF Ensemble (10km)'; case = 'prj12'; domain = (1,2)
  # #   explist = ['g-ens']; exptitles = 'WRF Ensemble (30km)'; case = 'prj'; domain = 1
  # #   explist = ['g3-ens']; exptitles = 'WRF Ensemble (90km)'; case = 'g3prj'; domain = 1
  # #   lfrac = True; refprd = H15; period = B15; reflist = explist

  # # continental-scale map with domains
  # #   lpickle = False; lprint = False
  #   case = 'can_wrf2'; figtitles = 'Topography and Domain Outline [km]'
  #   variables = ['zs']; seasons = ['topo']; lcontour = True;
  #   lframe = True; framewidths = 2.; framecolor = 'w'
  #   loutline = True; outlinewidth = 1.; outlinecolor = 'k'
  #   maptype = 'ortho-can'; lstations = False; stations = 'EC'
  #   period = H15; lWRFnative = True; lbackground = False
  #   explist = [('Ctrl-1','erai-g3','erai-g','erai-g')]; exptitles = ' '; domain = (0,1,1,2)
  #   lbasins = False; basinlist = ('GRW',); primary_basins = basinlist; subbasins = {} #dict(ARB=('WholeARB','UpperARB','LowerCentralARB'))
  #   lprovinces = False; provlist = ('ON',)

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

  # ## observations
  #   variables = ['precip']; seasons = ['annual']
  #   explist = ['Unity']; period = H15
  # #   ldiff = True; reflist = ['Unity']; maptype = 'lcc-small'
  #   exptitles = 'Annual Total Precipitation [mm/day]'; figtitles = '' # ['Merged Observations (10 km)']
  #   case = 'unity'; grid = 'arb2_d02'
  #   lsamesize = False; lcontour = True
  #   lprovinces = True; provlist = ('AB',)
  #   lbasins = True; basinlist = ('FRB','SSR'); maptype = 'lcc-prairies'
  # #  cluster_name = 'cluster_historical'; cluster_symbols = {i:'o' for i in xrange(10)} # '^','s'
  # #  cluster_symbols = {clu:dict(marker=sym, markersize=4, mfc='w', mec='k') for clu,sym in cluster_symbols.iteritems()}
  # #  lbasins = False; lstations = True; stations = 'EC'; case += '_stations'

  ## single panel plot
  #  explist = ['max-ens-2100']; maptype = 'lcc-new'; period = B15
  #  lfrac = True; reflist = ['max-ens']; refprd = H15
  #  case = 'sum'; lsamesize = False; figtitles = ''; # primary_basins = ['FRB','ARB']
  #  WRFfiletypes = ['lsm']; variables = ['aSM']; seasons = ['jas']; exptitles = 'Soil Moisture Change [%]'
  #  variable_settings = ['asm_red']
  #  lsamesize = False; lcontour = True
  #  lprovinces = True; provlist = ('AB',)
  #  lbasins = True; basinlist = ('FRB','SSR'); maptype = 'lcc-prairies'

  ## simple six-panel plot
  #  explist = ['Unity','max-ens','ctrl-ens','Ens','max-ens','ctrl-ens']; maptype = 'lcc-new'; period = H15
  #  grid = ['arb2_d02']*3+['cesm1x1']+['arb2_d01']*2; domain = [None, 2, 2, None, 1, 1]
  #  exptitles = ['Merged Observations','1st WRF Ens., 10km','Alt. WRF Ens., 10km','CESM','1st WRF Ens., 30km','Alt. WRF Ens., 30km']
  #  case = 'res'; lbasins = True; lprovinces = True
  #  variables = ['precip']; seasons = ['annual']

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

  # # large map for all domains
  #   variables = ['zs']; seasons = ['topo']; lcontour = True
  # #   lbasins = True; lprovinces = True; provlist = ['AB']; basinlist = ['FRB','ARB']; case = 'basins'
  # #   framewidths = 2; basin_args = dict(linewidth = 1.5, color='k'); case += '_fat'
  #   #lstations = False; lbasins = False; lprovinces = True; # provlist = ('AB',)
  #   figtitles = ['Topography [km]' + ' and Domain Outlines' if lframe else '']; exptitles = ' '
  # #   maptype = 'lcc-large'; figuretype = 'largemap'; loutline = False; lframe = True
  # #   explist = ['max']; exptitles = ' '; domain = (0,1,2); lWRFnative = True; period = H15
  #   explist = ['erai-wc2-2010']; exptitles = ' '; domain = (1,2); lWRFnative = True; period = '2010-2011'
  #   lWRFnative = True; lframe = True; loutline = False; lbasins = False
  # #  case = 'arb2_basins'; #basinlist = ('FRB','ARB','CRB','NRB'); primary_basins = ('FRB','ARB')
  # # smaller map for western Canada
  # #   explist = ['max-3km']; domain = (1,2,3); exptitles = ['Terrain Height [km]']; figtitles = ' '
  # #   period = H10; lWRFnative = True; lframe = True; loutline = False
  # #   maptype = 'lcc-prairies'
  #   maptype = 'lcc-col_out'; case = 'col'
  # #   maptype = 'lcc-arb3_d02'; case = 'large'
  # #   maptype = 'lcc-arb3_d03'; case = 'small'
  #   #maptype = 'lcc-arb3_d03'; case = 'hidef'; seasons = ['hidef']
  # #   maptype = 'ortho-NA'; case = 'global'; basinlist = []
  #   #domain = (1,); case += '_d01'
  #   domain = (1,2); case += '_d02'
  # #   domain = (0,1,); case += '_d01'
  # #   domain = (0,1,2); case += '_d02'
  # #   domain = (0,1,2,3); case += '_d03'
  # #  maptype = 'lcc-arb3_d03'; case = 'hidef'; seasons = ['hidef']; domain = (3,) # ; explist = ['col1-const']
  # #  case = 'frb'; basins = ('FRB',)
  # #  case = 'arb'; basins = ('ARB',)
  # #   case = 'ssr'; basinlist = ('SSR',)


    if not case: raise ValueError('Need to define a \'case\' name!')

    if not variables: raise ValueError('Need to define a variable list!')

    # setup projection and map
    mapSetup = getSetup(maptype, lpickle=lpickle, folder=map_folder)

    ## load data
    if not lfrac and not ldiff: reflist = None

    if reflist is not None:
      if isinstance(reflist,str): reflist = [reflist]
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
                                         WRF_exps=WRF_exps, CESM_exps=CESM_exps, unity_grid=unity_grid)
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
        raise DatasetError('Experiments and reference list need to have the same length!')
      for i in range(len(exps)):
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
    comments = checkItemList(comments, N, str, default=None)
    figtitles = checkItemList(figtitles, N, str, default=None)
    variable_settings = checkItemList(variable_settings, N, str, default=None)
    season_settings = checkItemList(season_settings, N, str, default=None)

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
        if isinstance(varlist,str): varstr = varlist
        elif isinstance(varlist,(list,tuple)):
          if all([var==varlist[0] for var in varlist]): varstr = varlist[0]
          else: varstr = ''.join([s[0] for s in varlist])
        else: varstr = ''
        varlist = checkItemList(varlist, M, str)
        ravlist = checkItemList(ravlist, M, str)
  #       if ldiff or lfrac:
  #       else:
  #         varlist = checkItemList(varlist, M, basestring)
        # expand seasons
        if isinstance(sealist,str): seastr = '_'+sealist
        elif isinstance(sealist,(list,tuple)): seastr = '_'+''.join([s[0] for s in sealist])
        else: seastr = ''
        sealist = checkItemList(sealist, M, str)

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
        print(('\n\n   ***  %s %s (%s)   ***   \n'%(plottype,plat.title,varstr)))

        ## compute data
        data = []; lons = []; lats=[]  # list of data and coordinate fields to be plotted
        # compute average WRF precip
        print((' - loading data ({0:s})'.format(varstr)))
        for var,rav,season,exptpl in zip(varlist,ravlist,sealist,exps):
          lontpl = []; lattpl = []; datatpl = []
          for i,exp in enumerate(exptpl):
            varname = rav if ( lfrac or ldiff ) and i >= len(exptpl)//2 else var
            if varname not in exp:
              raise DatasetError("Variable '{:s}' not found in Dataset '{:s}!".format(varname,exp.name))
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
            for ax,la in level_agg.items():
              if expvar.hasAxis(ax):
                if isinstance(la, str): # aggregate over axis
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
              if 'landmask' in exp.variables and False:
                vardata[exp.landmask.getArray()] = -2.
              elif 'landfrac' in exp.variables: # CESM mostly
                vardata[exp.landfrac.getArray(unmask=True,fillValue=0)<0.75] = -2. # use land fraction
              elif 'lndidx' in exp.variables:
                mask = exp.lndidx.getArray()
                vardata[mask==16] = -2. # use land use index (ocean)
                vardata[mask==24] = -2. # use land use index (lake)
              else :
                vardata = maskoceans(lon,lat,vardata,resolution=res,grid=grid)
            # figure out land mask
            if lmsklnd:
              if 'landfrac' in exp.variables: # CESM and CFSR
                vardata[exp.lnd.getArray(unmask=True,fillValue=0)>0.75] = 0 # use land fraction
              elif 'lndidx' in exp.variables: # use land use index (ocean and lake)
                mask = exp.lndidx.getArray(); tmp = vardata.copy(); vardata[:] = 0.
                vardata[mask==16] = tmp[mask==16]; vardata[mask==24] = tmp[mask==24]
            datatpl.append(vardata) # append to data list
          ## compute differences, if desired
          if ldiff or lfrac:
            assert len(datatpl)%2 == 0, 'needs to be divisible by 2'
            ntpl = len(datatpl)//2 # assuming (exp1, exp2, ..., ref1, ref2, ...)
            for i in range(ntpl):
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
        for n in range(nax):
          ax.append(f.add_subplot(subplot[0],subplot[1],n+1, facecolor='blue'))
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
        if loutline or lframe:
          print(' - drawing domain outlines\n')
          for n in range(nax):
            for m in range(nexps[n]):
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
        toString = lambda v: v if isinstance(v,str) else clbl%v
        for n in range(nax):
          for m in range(nexps[n]):
            vmean = toString(np.nanmean(data[n][m]))
            vmin = toString(np.nanmin(data[n][m]))
            vmax = toString(np.nanmax(data[n][m]))
            if ldiff or lfrac:
              vrms = toString(np.sqrt(np.nanmean(data[n][m]**2)))
              print(('panel {:d}: bias {:s} / rms {:s} / min {:s} / max {:s}'.format(
                              n, vmean, vrms, vmin, vmax)))
            else:
              vstd = toString(np.nanstd(data[n][m]))
              print(('panel {:d}: mean {:s} / std {:s} / min {:s} / max {:s}'.format(
                              n, vmean, vstd, vmin, vmax)))
            if lcontour:
              cd.append(maps[n].contourf(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],cmap=cmap,
                                         norm=norm,extend='both'))
            else:
                cd.append(maps[n].pcolormesh(x[n][m],y[n][m],data[n][m], cmap=cmap, shading=shading))
            # add black contour lines (outlines)
            if laddContour:
              cd.append(maps[n].contour(x[n][m],y[n][m],data[n][m],clevs,ax=ax[n],
                                        colors='k', linewidths=0.5))
            if isoline is not None:
              cd.append(maps[n].contour(x[n][m],y[n][m],data[n][m],[isoline],ax=ax[n],
                                        colors='w', linewidths=1.))

        ## highligh 60th parallel to delineate provinces and territories
        if l60:
            for m in maps: m.drawparallels([60],labels=[False,]*4, dashes=[], **prov_args)

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
        msn = len(maps)//2 # place scale
        mapSetup.drawScale(maps[msn])
        n = -1 # axes counter
        for i in range(subplot[0]):
          for j in range(subplot[1]):
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
            if locean or ( lmsklnd and not ('lndmsk' in exps[n][0].variables )
                           or 'lndidx' in exps[n][0].variables):
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
                from projects.WesternCanada import loadWRF_StnTS
                varlist = stn_params + [cluster_name]; station_type = 'ecprecip'
                ecstns = loadEC_StnTS(station=station_type, varlist=varlist)
                wrfstns = loadWRF_StnTS(experiment='max-ctrl', varlist=varlist, station=station_type,
                                        filetypes='hydro', domains=2)

                ecstns,wrfstns = selectStations([ecstns, wrfstns] , stnaxis='station', linplace=False, lall=True,
                                                **station_constraints)
                # loop over points
                if cluster_name in ecstns: cluster_axis = ecstns[cluster_name]
                else: cluster_axis = [-1]*len(ecstns.stn_lat)
                for lon,lat,zerr,cln in zip(ecstns.stn_lon, ecstns.stn_lat, wrfstns.zs_err, cluster_axis):
                    tmp = cluster_symbols.get(cln,cluster_symbols[-1]).copy()
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
                    for subbasin in subbasins[basin]:
                      bmap.readshapefile(basininfo.shapefiles[subbasin][:-4], subbasin, ax=axn,
                                         drawbounds=True, **subbasin_args)
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
          print(('\nSaving figure in '+filename))
          f.savefig(folder+filename, **sf) # save figure to pdf
          print(folder)

    ## show plots after all iterations
    #plt.show()

