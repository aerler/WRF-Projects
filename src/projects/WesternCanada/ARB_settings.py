'''
Created on Sep 4, 2016

This module contains a meta data for HGS simulations for the GRW and a wrapper to load them. 

@author: Andre R. Erler, GPL v3
'''
import numpy as np
from collections import OrderedDict
from collections import namedtuple
# internal imports
import hgs.HGS as hgs # need to prevent name collisions here
from hgs.PGMN import getWellName
import projects.WSC_basins as wsc
from projects.WesternCanada.WRF_experiments import WRF_exps
from geodata.misc import ArgumentError, DatasetError
# imports from HGS_settings
import projects.HGS_settings as default


project_name = 'ARB'
project_prefix = 'ARB_SSM_Jan2Dec_23Lyr'
conservation_authority = 'COSIA'
# some parameters
gage_scalefactors = 1. # default scalefactor for discharge plots
main_gage   = 'McMurray' # gage station(see station_list for proper HGS/WSC names)
main_basin  = 'ARB' # main basin name
main_grid   = 'arb2' # climate data grid
binary_grid = 'arb2' # grid used to interpolated binary output
# project folders
project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project_name) # the dataset root folder
project_folder_pattern = '{PROJECT_FOLDER:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_DIR:s}/{TASK:s}/'
station_file = '{PREFIX:s}o.hydrograph.Station_{STATION:s}.dat' # general HGS naming convention
station_file_default = station_file # alias
# mapping of WSC station names to HGS hydrograph names
Station = namedtuple('Station', ('HGS','WSC','ylim'),)
station_list = OrderedDict() # this is the main gage of the GRW
station_list['Jasper'] = Station(HGS='01_JASPER', WSC='Athabasca_Jasper', ylim=None)
station_list['Hinton'] = Station(HGS='02_HINTON', WSC='Athabasca_Hinton', ylim=None)
station_list['Windfall']  = Station(HGS='03_WINDFALL',  WSC='Athabasca_Windfall', ylim=None)
station_list['Athabasca'] = Station(HGS='04_ATHABASCA', WSC='Athabasca_Athabasca', ylim=None)
station_list['Fort McMurray']    = Station(HGS='05_MCMURRAY', WSC='Athabasca_FortMcMurray', ylim=None)
station_list['Embarras Airport'] = Station(HGS='06_EMBARRAS', WSC='Athabasca_EmbarrasAirport', ylim=None)
station_list_default = station_list # just an alias...
# look-up tables for WSC/HGS station name conversion                           
WSC_station_list = {stn.WSC:stn.HGS for stn in list(station_list.values())}
HGS_station_list = {stn.HGS:stn.WSC for stn in list(station_list.values())}


# plotting parameters for HGS simulations
hgs_plotargs_rc = dict() 
# stream observations
hgs_plotargs_rc['Observations'] = dict(color='#959595') # gray
hgs_plotargs_rc['Obs.']         = dict(color='#959595') # gray
hgs_plotargs_rc['WSC Obs.']     = dict(color='#959595') # gray
hgs_plotargs_rc['WSC']          = dict(color='#959595') # gray
# Observational forcing
hgs_plotargs_rc['CRU']          = dict(color='green') #, linewidth=3)
hgs_plotargs_rc['NRCan']        = dict(color='coral') # orange
hgs_plotargs_rc['SnoDAS']       = dict(color='#62A1C6') # blue
# Temporal Aggregation
hgs_plotargs_rc['Steady-State (V1)'] = dict(color='#AAA2D8') # purple
hgs_plotargs_rc['Periodic (V1)']     = dict(color='green') # 
hgs_plotargs_rc['Steady-State (V2)'] = dict(color='#E24B34') # red
hgs_plotargs_rc['Periodic (V2)']     = dict(color='coral')
hgs_plotargs_rc['Transient']    = dict(color='#62A1C6') # blue
hgs_plotargs_rc['Steady-State'] = dict(color='#AAA2D8') # purple
hgs_plotargs_rc['Periodic']     = dict(color='#E24B34') # red
hgs_plotargs_rc['Normals']      = dict(color='green') # 
hgs_plotargs_rc['Monthly']      = dict(color='#62A1C6') # blue
hgs_plotargs_rc['Daily']        = dict(color='#E24B34') # red
# WRF forcing
hgs_plotargs_rc['WRF 10km']     = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['WRF 30km']     = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF (MyBC)']   = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF (AABC)']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['WRF (Delta)']  = dict(color='#AAA2D8') # purple
# hgs_plotargs_rc['WRF G 10km']   = dict(color='#62A1C6') # blue                                     
# hgs_plotargs_rc['WRF G 30km']   = dict(color='#E24B34') # red                                      
# hgs_plotargs_rc['WRF G 90km']   = dict(color='#AAA2D8') # purple                                   
# hgs_plotargs_rc['WRF T 10km']   = dict(color='#62A1C6') # blue    
# hgs_plotargs_rc['WRF T 30km']   = dict(color='#E24B34') # red     
# hgs_plotargs_rc['WRF T 90km']   = dict(color='#AAA2D8') # purple  
hgs_plotargs_rc['WRF T 10km']   = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF G 10km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['WRF T 30km']   = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF G 30km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['WRF T 90km']   = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF G 90km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['WRF T']        = dict(color='#E24B34') # red
hgs_plotargs_rc['WRF G']        = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['T Ensemble']   = dict(color='#E24B34') # red
hgs_plotargs_rc['G Ensemble']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['ERAI-T']       = dict(color='#E24B34') # red
hgs_plotargs_rc['ERAI-G']       = dict(color='#62A1C6') # blue                                     
hgs_plotargs_rc['T Mean']       = dict(color='red') # red
hgs_plotargs_rc['G Mean']       = dict(color='black') # blue                                     
hgs_plotargs_rc['Ensemble']     = dict(color='red')
hgs_plotargs_rc['Mean']         = dict(color='black')
hgs_plotargs_rc['1980']         = dict(color='#62A1C6') # blue
hgs_plotargs_rc['2050']         = dict(color='#AAA2D8') # purple
hgs_plotargs_rc['2100']         = dict(color='#E24B34') # red
hgs_plotargs_rc['1979-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs_rc['1984-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs_rc['2045-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs_rc['2050-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs_rc['2085-2100']    = dict(color='#E24B34') # red
hgs_plotargs_rc['2090-2100']    = dict(color='#E24B34') # red
# adjust line thickness
for plotargs in list(hgs_plotargs_rc.values()): 
  if 'linewidth' not in plotargs: plotargs['linewidth'] = 1.
# extended color scheme for ensemble
color_args = {'Ctrl':'blue', 'Ens-A':'purple', 'Ens-B':'green','Ens-C':'coral','90km':'purple','30km':'red','10km':'blue'}
for key,value in list(color_args.items()): hgs_plotargs_rc[key] = dict(color=value, linewidth=.75)
# add "old" versions
for expname,plotarg in list(hgs_plotargs_rc.items()):
    if expname[-1] != ')' and expname+' (old)' not in hgs_plotargs_rc:
        oldarg = plotarg.copy(); 
        oldarg['linestyle'] = '--' 
        hgs_plotargs_rc[expname+' (old)'] = oldarg


# experiment aliases (for more systematic access)
exp_aliases = dict()
for config,alias in dict(g='max', t='ctrl').items():
    for ens in ('-ctrl','-ens-A','-ens-B','-ens-C','-ensemble'):
        for prd in ('','-2050','-2100'):
            exp_aliases[config+ens+prd] = alias+ens+prd
exp_aliases['t-ctrl'] = 'ctrl-1'
exp_aliases['t-ctrl-2050'] = 'ctrl-2050'
exp_aliases['t-ctrl-2100'] = 'ctrl-2100'
exp_aliases_default = exp_aliases # just an alias...
obs_datasets = ('NRCan','CRU','PRISM')
gage_datasets = ('wsc','obs','observations')
# ensemble definitions for GRW project 
ensemble_list = {'g-mean':('g-ctrl','g-ens-A','g-ens-B','g-ens-C'),
                 't-mean':('t-ctrl','t-ens-A','t-ens-B','t-ens-C'),
                 'daily-mean':('g-ens-A','g-ens-B')}
for name,members in list(ensemble_list.items()):
    for prd in ('-2050','-2100'):
        ensemble_list[name+prd] = tuple(member+prd for member in members)
  

## parameter definition function

# helper function to determine experiment parameters
def experimentParameters(experiment=None, domain=None, clim_mode=None, clim_dir=None, 
                         lold=True, task=None, bias_correction=None, exp_aliases=None,
                         clim_period=None, run_period=None, period=None, **kwargs):
    # figure out experiment name and domain
    lWRF = False # likely not a WRF experiment
    old_name = experiment # save for later
    if domain is not None:
        lWRF = True 
        exp_name = experiment
        experiment = '{:s}_d{:02d}'.format(experiment,domain) # used to determine aliases
        # extract resolution for later use (before we convert aliases)
        if domain == 0: resolution = '90km'
        elif domain == 1: resolution = '30km'
        elif domain == 2: resolution = '10km'
        elif domain == 3: resolution = '3km'
        else: raise NotImplementedError("Unsupported domain number '{:d}'.".format(domain))
        if 'resolution' not in kwargs: kwargs['resolution'] = resolution # for name expansion (will be capitalized)
    if exp_aliases is None: exp_aliases = exp_aliases_default
    if experiment in exp_aliases:
        lWRF = True
        # resolve aliases (always return full string format)
        experiment = exp_aliases[experiment] # resolve alias
        if domain is None and experiment[-4:-2] == '_d':
            exp_name = experiment[:-4]; domain = int(experiment[-2:]) # always full string format
    elif experiment[-4:-2] == '_d' and experiment[:-4] in exp_aliases:
        lWRF = True
        # resolve aliases (always return full string format)
        experiment = exp_aliases[experiment[:-4]]+experiment[-4:] # resolve alias
        exp_name = experiment[:-4]; domain = int(experiment[-2:]) # always full string format
    if experiment in WRF_exps: # i.e. if domain is None
        lWRF = True
        exp_name = experiment
        domain = WRF_exps[experiment].domains # innermost domain
        experiment = '{:s}_d{:02d}'.format(experiment,domain) # used to determine aliases
    # resolve climate input mode (time aggregation; assuming period is just length, i.e. int;)
    # translate climate mode into path convention, while retaining clim_mode
    if isinstance(clim_dir,str): lts = True
    elif clim_mode.lower() in ('daily_transient'):
        clim_dir = 'transient_daily'; lts = True
    elif clim_mode.lower() in ('ts','trans','timeseries','transient'):
        clim_dir = 'timeseries'; lts = True
    elif clim_mode.lower().startswith(('ts_','trans_','timeseries_','transient_')): 
        clim_dir = clim_mode; lts = True; clim_mode = 'timeseries'
    elif clim_mode.lower() in ('mn','norm','clim','peri','normals','climatology','periodic'): 
        clim_dir = 'clim'; lts = False
    elif clim_mode.lower() in ('ss','const','mean','annual','steady-state'): 
        clim_dir = 'annual'; lts = False
    else: raise ArgumentError(clim_mode)
    # update climate input mode for name expansion
    kwargs['clim_mode'] = clim_mode
    # set some defaults based on version
    if lold:
        if run_period is None and not lts: run_period = 15 # in years
        if task is None: task = 'hgs_run_wrfpet' if bias_correction else 'hgs_run' 
        if clim_period is None and not lts: clim_period = 15 # in years
    else:
        if run_period is None and not lts: run_period = 10 if lWRF else 5 # in years
        if task is None: task = 'hgs_run_v3_wrfpet' if bias_correction else 'hgs_run_v3'
        if clim_period is None and not lts: clim_period = 15 if lWRF else 30 # in years
    if not lts: clim_dir += '_{:02d}'.format(clim_period)
    # add period extensions of necessary/possible
    # N.B.: period is the period we want to show, run_period is the actual begin/end/length of the run
    if run_period is None: run_period = period
    elif period is None: period = run_period
    if isinstance(run_period, (tuple,list)):
        start_year, end_year = run_period
        if isinstance(period, (int,np.integer)): period = (end_year-period,end_year)
    elif isinstance(run_period, (int,np.integer)):
        if isinstance(period, (tuple,list)):
            end_year = period[1]
            start_year = end_year - run_period
        elif isinstance(period, (int,np.integer)):
            # try to guess from naming convention
            if '-2100' in old_name: start_year = 2085
            elif '-2050' in old_name: start_year = 2045
            else: start_year = 1979
            end_year = start_year + run_period
            period = (end_year-period,end_year)
        else: raise ArgumentError(period,run_period)
    else: start_year = end_year = None # not used
    # append conventional period name to name, if it appears to be missing
    if start_year is not None:
        if start_year > 2080:
          if exp_name[-5:] != '-2100': 
              if exp_name[-2:] == '-1': exp_name = exp_name[:-2] # handle inconsistent naming... 
              exp_name += '-2100'
          if old_name[-5:] != '-2100': 
              if old_name[-2:] == '-1': old_name = old_name[:-2] # handle inconsistent naming...
              old_name += '-2100'
        elif start_year > 2040:
          if exp_name[-5:] != '-2050': 
              if exp_name[-2:] == '-1': exp_name = exp_name[:-2] # handle inconsistent naming...
              exp_name += '-2050'
          if old_name[-5:] != '-2050': 
              if old_name[-2:] == '-1': old_name = old_name[:-2] # handle inconsistent naming...
              old_name += '-2050'
    if 'exp_name' not in kwargs: kwargs['exp_name'] = old_name # save for later use in name
    if 'exp_title' not in kwargs: kwargs['exp_title'] = old_name.title() # save for later use in name
    if 'exp_upper' not in kwargs: kwargs['exp_upper'] = old_name.upper() # save for later use in name
    # validate WRF experiment
    if lWRF: 
        # reassemble experiment name with domain extension
        experiment = '{:s}_d{:02d}'.format(exp_name,domain)
        exp = WRF_exps[exp_name]
        if domain > exp.domains or domain < 1: raise DatasetError("Invalid domain number: {}".format(domain))
        # get start date from experiment (if not defined in period)
        if start_year is None: 
          start_date = tuple(int(d) for d in exp.begindate.split('-')) # year,month,day tuple
          start_year = start_date[0]
        else:
          start_date = start_year  
    else:
        # assume observationa data, with time axis origin in Jan 1979
        if start_year is None: start_year = 1979
        start_date = start_year
        exp = None
    if end_year is None: end_year = start_year + run_period # period length
    else: run_period = end_year - start_year
    assert isinstance(run_period,(np.integer,int)), run_period
    # construct period information  
    # append conventional period name to name, if it appears to be missing
    if start_year < 1994 and end_year < 1995: prdext = '1980'; prdstr = '1979-1994'
    elif start_year < 2015  and end_year < 2016: prdext = '2000'; prdstr = '1979-2015'
    elif start_year < 2060 and end_year < 2061: prdext = '2050'; prdstr = '2045-2060'
    elif start_year < 2100 and end_year < 2101: prdext = '2100'; prdstr = '2085-2100'
    else: raise NotImplementedError("Unable to determine period extension for period '{:d}-{:d}'.".format(start_year,end_year))
    if 'prdext' not in kwargs: kwargs['prdext'] = prdext # for name expansion (will be capitalized)
    if 'prdstr' not in kwargs: kwargs['prdstr'] = prdstr # for name expansion (will be capitalized)
    # select bias-correction method, if applicable
    if bias_correction is not None:
        clim_dir = '{:s}_{:s}'.format(bias_correction,clim_dir)
    
    # return parameters
    return ( lWRF, task,experiment,exp, clim_dir, start_date,start_year,end_year, 
             clim_period,run_period,period, kwargs )

  
## wrapper for load functions  

def loadHGS_StnTS(experiment=None, domain=None, period=None, varlist=None, varatts=None, name=None, 
                  title=None, run_period=None, clim_mode=None, clim_period=None, lold=False,
                  station=None, well=None, bias_correction=None, Obs_well=None,  
                  task=None, WSC_station=None, lpad=True, ENSEMBLE=None, lWSCID=False,
                  project_folder=project_folder, station_file=station_file_default,
                  basin_list=None, project=project_name, folder=project_folder_pattern, 
                  conservation_authority=conservation_authority, basin=main_basin, main_gage=main_gage,
                  station_list=None, experimentParameters=experimentParameters, 
                  grid=main_grid, scalefactors=gage_scalefactors, **kwargs):
    ''' a wrapper for the HGS_settings functions, which sets some default values '''
    
    # resolve station and well names
    if station and station.lower() in ('water_balance','newton_info'): pass 
    elif well and well.lower() in ('water_balance','newton_info'): pass
    elif well and station: raise ArgumentError(station,well)
    elif station:
        if station_list is None: station_list = station_list_default
        if station in station_list: 
            if WSC_station is None: WSC_station = station_list[station].WSC
            station = station_list[station].HGS if lold else '{WSC_ID0:s}'
        elif station in WSC_station_list:
            WSC_station = station; station = WSC_station_list[station]
        elif station in HGS_station_list: 
            WSC_station =  HGS_station_list[station]
    elif well:
        # decompose well name...
        well_id,well_no = getWellName(well)
        well = 'W{WELL_ID:07d}_{WELL_NO:1d}'.format(WELL_ID=well_id, WELL_NO=well_no)
        Obs_well = 'W{WELL_ID:07d}-{WELL_NO:1d}'.format(WELL_ID=well_id, WELL_NO=well_no)
    if basin_list is None: basin_list = wsc.basin_list # default basin list
        
    # now call primary function
    return default.loadHGS_StnTS(
                  experiment=experiment, domain=domain, period=period, varlist=varlist, varatts=varatts, 
                  name=name, title=title, run_period=run_period, clim_mode=clim_mode, clim_period=clim_period, 
                  station=station, well=well, bias_correction=bias_correction, project_folder=project_folder, 
                  task=task, WSC_station=WSC_station, Obs_well=Obs_well, lpad=lpad, ENSEMBLE=ENSEMBLE, 
                  lWSCID=lWSCID, project=project, folder=folder, station_file=station_file, 
                  grid=grid, conservation_authority=conservation_authority, basin_list=basin_list, 
                  basin=basin, scalefactors=scalefactors, main_gage=main_gage, station_list=station_list, 
                  experimentParameters=experimentParameters, lold=lold, **kwargs)


# custom gage station loader
def loadWSC_StnTS(station=main_gage, name=None, title=None, basin=main_basin, basin_list=None, varlist=None, 
                  varatts=None, period=None, filetype='monthly', scalefactors=gage_scalefactors, lkgs=False, 
                  **kwargs):
  ''' a wrapper to load gage stations for the GRW with appropriate default values'''
  # load gage station
  if station in station_list: station = station_list[station].WSC # get WSC name for gage station (different from HGS name...)
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  return wsc.loadWSC_StnTS(basin=basin, station=station, varlist=varlist, varatts=varatts, filetype=filetype,
                           name=name, title=title, basin_list=basin_list, period=period, scalefactors=scalefactors, 
                           lkgs=lkgs, **kwargs)
  

# wrapper to load HGS ensembles, otherwise the same
def loadHGS_StnEns(ensemble=None, station=None, varlist=None, varatts=None, name=None, title=None, 
                   period=None, domain=None, exp_aliases=exp_aliases, run_period=None, clim_mode=None,
                   folder=project_folder_pattern, project_folder=None, obs_period=None, clim_period=None, 
                   ensemble_list=ensemble_list, ensemble_args=None, observation_list=gage_datasets, # ensemble and obs lists for project
                   project=project_name, grid=main_grid, task=None, 
                   bias_correction=None, conservation_authority=conservation_authority,
                   WSC_station=None, basin=main_basin, basin_list=None, scalefactors=gage_scalefactors, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations and assemble ensembles '''  
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  return hgs.loadHGS_StnEns(ensemble=ensemble, station=station, varlist=varlist, varatts=varatts, name=name, 
                            title=title, period=period, run_period=run_period, folder=folder, domain=domain,  
                            obs_period=obs_period, clim_period=clim_period, ensemble_list=ensemble_list, 
                            ensemble_args=ensemble_args, observation_list=observation_list, 
                            loadHGS_StnTS=loadHGS_StnTS, loadWSC_StnTS=loadWSC_StnTS, # use local versions of loaders
                            WSC_station=WSC_station, basin=basin, basin_list=basin_list, 
                            bias_correction=bias_correction, conservation_authority=conservation_authority,
                            project_folder=project_folder, project=project, grid=grid, clim_mode=clim_mode, 
                            exp_aliases=exp_aliases, task=task, scalefactors=scalefactors, **kwargs)  


## function to interpolate HGS binary output to regular grid
def gridDataset(dataset, griddef=binary_grid, basin=main_basin, subbasin=None, shape_file=None,  
                basin_list=None, grid_folder=None, laddMask=True, **kwargs):
    ''' interpolate nodal/elemental datasets to a regular grid, add GDAL, and mask to basin outlines,
        using values for the GRW and the grw1 grid'''
    if basin_list is None: basin_list = wsc.basin_list # default basin list
    return hgs.gridDataset(dataset, griddef=griddef, basin=basin, subbasin=subbasin, shape_file=shape_file,  
                           basin_list=basin_list, grid_folder=grid_folder, laddMask=laddMask, **kwargs)
    
  
## function to load HGS binary data
def loadHGS(experiment=None, varlist=None, name=None, title=None, lstrip=True, lgrid=False, sheet=None,
            lkgs=False, lflipdgw=False, season=None, period=None,
            griddef=binary_grid, basin=main_basin, subbasin=None, grid_folder=None, shape_file=None, 
            domain=None, clim_mode=None, clim_period=None, bias_correction=None, task=None, grid=main_grid,
            mode='climatology', file_mode='last_12', file_pattern='{PREFIX}o.head_olf.????', t_list=None, 
            varatts=None, constatts=None, project_folder=project_folder, project=project_name, lxyt=True,
            folder=project_folder_pattern, metadata=None, conservation_authority=conservation_authority, 
            basin_list=None, WRF_exps=WRF_exps, experimentParameters=experimentParameters, **kwargs):
    ''' a wrapper for the regular HGS binary load function '''
    if basin_list is None: basin_list = wsc.basin_list # default basin list
    # load dataset with some default values
    return default.loadHGS(experiment=experiment, varlist=varlist, name=name, title=title, lstrip=lstrip,
                           lgrid=lgrid, griddef=griddef, subbasin=subbasin, season=season, period=period,
                           grid_folder=grid_folder, shape_file=shape_file, sheet=sheet, lxyt=lxyt,
                           basin=basin, lkgs=lkgs, domain=domain, clim_mode=clim_mode, grid=grid,
                           clim_period=clim_period, bias_correction=bias_correction, task=task,    
                           mode=mode, file_mode=file_mode, file_pattern=file_pattern, t_list=t_list, 
                           varatts=varatts, constatts=constatts, project_folder=project_folder, 
                           project=project, folder=folder, basin_list=basin_list, metadata=metadata, 
                           conservation_authority=conservation_authority, WRF_exps=WRF_exps, 
                           lflipdgw=lflipdgw, experimentParameters=experimentParameters, **kwargs)


# abuse for testing
if __name__ == '__main__':
    
#   test_mode = 'gage_station'
#   test_mode = 'create_grid'
#   test_mode = 'dataset_regrid'
#   test_mode = 'binary_dataset'
  test_mode = 'dataset'
#   test_mode = 'ensemble'

  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadWSC_StnTS(period=(1979,2009),) # station='Nith River at Canning')
    print(ds)
    
  ## create a new grid
  elif test_mode == 'create_grid':
    
    import os
    from geodata.gdal import GridDefinition, pickleGridDef, loadPickledGridDef, grid_folder
    
    convention='Proj4'
    ## parameters for Athabasca River Basin project (COSIA)
    projection = '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs'
    #projection = 'PROJCS["US_National_Atlas_Equal_Area",GEOGCS["GCS_Unspecified datum based upon the Clarke 1866 Authalic Sphere",DATUM["D_Sphere_Clarke_1866_Authalic",SPHEROID["Clarke_1866_Authalic_Sphere",6370997,0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",45],PARAMETER["central_meridian",-100],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
#     name = 'arb1' # 1km
#     geotransform = [-1460500,1e3,0,810500,0,1e3]; size = (1420,1290)
#     name = 'arb2' # 5km
#     geotransform = [-1460500,5e3,0,810500,0,5e3]; size = (284,258)
    name = 'arb3' # a cropped version of arb2
    geotransform = [-1280e3,5e3,0,900e3,0,5e3]; size = (172,144)
    ## parameters for UTM 17 Great Lakes grids
#     name = 'glb1' # 5km resolution
#     geotransform = [ -709489.58091004, 5.e3, 0, 4148523.7226861, 0, 5.e3]; size = (405,371)
#     projection = "+proj=utm +zone=17 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
    ## grid for Elisha
#     name = 'uph1' # 5km resolution
#     geotransform = [ 443000, 5.e3, 0, 4774000, 0, 5.e3]; size = (3,6)
#     projection = "+proj=utm +zone=17 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
    ## GRW grids (also in Great Lakes, UTM 17)
#     name = 'grw1' # 1km resolution
#     geotransform = [500.e3,1.e3,0,4740.e3,0,1.e3]; size = (132,162)
#     name = 'grw2' # 5km resolution
#     geotransform = [500.e3,5.e3,0,4740.e3,0,5.e3]; size = (27,33)
#     name = 'grw3' # 500m resolution
#     d = 500
#     geotransform = [500.e3-d/2,d,0,4740.e3-d/2,0,d]; size = (250+1,320+1)
#     projection = "+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ## Southern Ontario grids (also in Great Lakes, UTM 17)
#     name = 'son1' # 5km resolution
#     # X 320919.7943000002 Y 4624073.9199, C 5890 R 4062 x 100m
#     geotransform = [320920.,5.e3,0,4624073.,0,5.e3]; size = (118,82)
#     name = 'son2' # 1km resolution
#     geotransform = [320920.,1.e3,0,4624073.,0,1.e3]; size = (590,410)
#     projection = "+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ## Hugo's grid for Quebec
#     name = 'hd1' # 5km resolution
#     geotransform = [-479184.769227,5.e3,0,68508.4877898,0,5.e3]; size = (70,49)
#     projection = "+proj=lcc +ellps=NAD83 +datum=NAD83 +lat_0=44.0 +lat_1=46.0 +lat_2=60.0 +lon_0=-68.5  +x_0=0 +y_0=0 +units=m +no_defs"
    ## UTM 18 parameters for South Nation grids
#     name = 'snw1' # 9km resolution
#     geotransform = [401826.125365249,9.e3,0,4851533.71730136,0,9.e3]; size = (22,29)
#     projection = "+proj=utm +zone=18 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
#     convention='Wkt'; projection = 'PROJCS["NAD_1983_UTM_Zone_14N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["false_easting",500000.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-99.0],PARAMETER["scale_factor",0.9996],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
#     name = 'snw2' # 2km resolution
#     geotransform = [438.e3,2.e3,0,4940.e3,0,2.e3]; size = (44,55)
#     projection = "+proj=utm +zone=18 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ## parameters for UTM 14 Assiniboine River Basin grids
    ## Bird River, a subbasin of the Assiniboine River Basin
#     name = 'brd1' # 5km resolution 
#     geotransform = (245.e3, 5.e3, 0., 5524.e3, 0., 5.e3); size = (39,32)
#     projection = "+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    ## Assiniboine River Basin
#     name = 'asb1' # 5km resolution 
#     geotransform = (-159.e3, 5.e3, 0., 5202.e3, 0., 5.e3); size = (191,135)
#     projection = "+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
#     name = 'asb2' # 1km resolution 
#     geotransform = (-159.e3, 1.e3, 0., 5202.e3, 0., 1.e3); size = (955,675)
#     projection = "+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    ## parameters for Canada-wide Lambert Azimuthal Equal-area
#     name = 'can1' # 5km resolution
#     llx = -3500000; lly = -425000; urx = 3000000; ury = 4000000; dx = dy = 5.e3
#     geotransform = [llx, dx, 0., lly, 0., dy]; size = ((urx-llx)/dx,(ury-lly)/dy)
#     size = tuple(int(i) for i in size)
#     projection = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs"
    # N.B.: (x_0, dx, 0, y_0, 0, dy); (xl,yl)
    #       GT(0),GT(3) are the coordinates of the bottom left corner
    #       GT(1) & GT(5) are pixel width and height
    #       GT(2) & GT(4) are usually zero for North-up, non-rotated maps
    # create grid
    griddef = GridDefinition(name=name, projection=projection, geotransform=geotransform, size=size, 
                             xlon=None, ylat=None, lwrap360=False, geolocator=True, convention=convention)

    #print(griddef.projection.exportToWkt())
    # save pickle to standard location
    filepath = pickleGridDef(griddef, folder=grid_folder, filename=None, lfeedback=True)
    assert os.path.exists(filepath)
    print('')
    
    # load pickle to make sure it is right
    del griddef
    griddef = loadPickledGridDef(grid=name, res=None, folder=grid_folder)
    print(griddef)

  elif test_mode == 'dataset_regrid':

    # load single dataset
#     ds = loadHGS(varlist=[], experiment='erai-g', domain=2,  
#                  clim_mode='periodic', bias_correction='AABC')
    ds = loadHGS(experiment='NRCan', task='hgs_run_v3', lgrid=True, 
                 clim_mode='periodic', varlist=['zs'])
    print('\n')
    print(ds)
    print('\n')
    
  elif test_mode == 'binary_dataset':

    # load single dataset
    ds = loadHGS(varlist=[], experiment='erai-g', domain=2,  
                 clim_mode='periodic', bias_correction='AABC')
#     ds = loadHGS(experiment='NRCan', task='hgs_run_v3', 
#                  clim_mode='periodic', varlist=[])
    print('\n')
    print(ds)
    if 'model_time' in ds:
        print('\n')
        print((ds.model_time))
        print((ds.model_time[:]))
    
  elif test_mode == 'dataset':

    # load single dataset
    ds = loadHGS_StnTS(experiment='max-ctrl', domain=2, basin='ARB', lkgs=False, lold=True,
                       station='Fort McMurray', period=(1979,1994), run_period=1, lskipNaN=False, 
                       task='hgs_run_cosia_1', grid='arb3', clim_mode='daily_transient',
                       lcheckComplete=False, bias_correction=None, resample='D', time_axis='datetime')
#     ds = loadHGS_StnTS(experiment='g-ctrl', domain=2, period=(1979,1994),
#                        station='water_balance', task='hgs_run_COSIA_1',
# #                        well='W424', z_aggregation=None, z_layers=None,
#                        clim_mode='transient', lpad=True, bias_correction='MyBC_CRU')
#     ds = loadHGS_StnTS(experiment='NRCan', task='hgs_run_v2', 
#                        clim_mode='periodic', lpad=True)
    print(ds)
    if 'discharge' in ds:
        print('\n')
        print((ds.discharge))
        print((ds.discharge.mean(), ds.discharge.max(), ds.discharge.min(),))
        print((ds.discharge.plot.units))
    if 'head' in ds:
        print('\n')
        print((ds.head))
        print((ds.head.mean(), ds.head.max(), ds.head.min(),))
        print((ds.head.plot.units))
    if ds.hasAxis('z'):
        print('\n')
        print((ds.z))
        print((ds.z[:]))
    
  elif test_mode == 'ensemble':
    
    print('')
    print(ensemble_list)
    print('')
    
    # load an esemble of datasets
    ens = loadHGS_StnEns(ensemble=['g-mean','t-mean'], domain=2, clim_mode='clim', 
                         well='W424', z_aggregation='max', z_layers='screen',
                         name='{EXP_NAME:s}_{RESOLUTION:s}', title='{EXP_NAME:s}_{RESOLUTION:s}',
                         period=[(1984,1994),(2050,2060),(2090,2100)], obs_period=(1974,2004),
                         lskipNaN=True, lcheckComplete=True, lold=False, bias_correction='AABC',
                         outer_list=['ensemble','period'], lensemble=True)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print((ens[0]))