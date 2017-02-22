'''
Created on Sep 4, 2016

This module contains a meta data for HGS simulations for the GRW and a wrapper to load them. 

@author: Andre R. Erler, GPL v3
'''
import numpy as np
import hgs.HGS as hgs # need to prevent name collisions here
import projects.WSC_basins as wsc
from projects.GreatLakes.WRF_experiments import WRF_exps
from collections import namedtuple
from geodata.misc import ArgumentError, DatasetError
from datasets.common import BatchLoad

project_name = 'GRW'
project_prefix = 'grw_omafra'
# some parameters
main_gage  = 'Brantford' # gage station(see station_list for proper HGS/WSC names)
main_basin = 'GRW' # main basin name
main_grid  = 'grw2' # climate data grid
main_task  = 'hgs_run' # HGS run folder
# project folders
project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project_name) # the dataset root folder
project_folder_pattern = '{PROJECT_FOLDER:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_DIR:s}/{TASK:s}/'

# plotting parameters for HGS simulations
hgs_plotargs = dict() 
hgs_plotargs['Observations'] = dict(color='#959595') # gray
hgs_plotargs['NRCan']        = dict(color='green')
hgs_plotargs['NRCan (old)']  = dict(color='green', linestyle='--')
# hgs_plotargs['Observations'] = dict(color='#68615E') # dark gray
hgs_plotargs['Transient']    = dict(color='#62A1C6') # blue
hgs_plotargs['Steady-State'] = dict(color='#AAA2D8') # purple
hgs_plotargs['Periodic']     = dict(color='#E24B34') # red
hgs_plotargs['WRF 90km']     = dict(color='#E24B34') # red
hgs_plotargs['WRF 30km']     = dict(color='#AAA2D8') # purple
hgs_plotargs['WRF 10km']     = dict(color='#62A1C6') # blue                                     
hgs_plotargs['1980']         = dict(color='#62A1C6') # blue
hgs_plotargs['2050']         = dict(color='#AAA2D8') # purple
hgs_plotargs['2100']         = dict(color='#E24B34') # red
hgs_plotargs['1979-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs['1984-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs['2045-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs['2050-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs['2085-2100']    = dict(color='#E24B34') # red
hgs_plotargs['2090-2100']    = dict(color='#E24B34') # red

# mapping of WSC station names to HGS hydrograph names
Station = namedtuple('Station', ('HGS','WSC'),)
station_list = dict( # short names for gages in this basin and their HGS/WSC names
    Brantford=Station(HGS='GR_Brantford',WSC='Grand River_Brantford') ) # this is the main gage of the GRW

# experiment aliases (for more systematic access)
exp_aliases = {'erai-g_d00':'erai-g3_d01','erai-t_d00':'erai-t3_d01',
               'g-ensemble_d00':'g3-ensemble_d01','t-ensemble_d00':'t3-ensemble_d01',
               'g-ensemble-2050_d00':'g3-ensemble-2050_d01','t-ensemble-2050_d00':'t3-ensemble-2050_d01',
               'g-ensemble-2100_d00':'g3-ensemble-2100_d01','t-ensemble-2100_d00':'t3-ensemble-2100_d01'}
obs_datasets = ['NRCan','CRU']
## wrapper functions to load HGS station timeseries with GRW parameters

# simple dataset loader
def loadHGS_StnTS(station=main_gage, varlist=None, varatts=None, name=None, title=None, period=None, 
                  experiment=None, domain=None, exp_aliases=exp_aliases, run_period=15, clim_mode=None,
                  folder=project_folder_pattern, project_folder=None, project=project_name, 
                  grid=main_grid, task=main_task, prefix=project_prefix, WSC_station=None, 
                  basin=main_basin, basin_list=None, lpad=True, **kwargs):
  ''' Get a properly formatted HGS dataset with a regular time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if experiment is None or clim_mode is None: raise ArgumentError
  # resolve station name
  if station in station_list: 
      if WSC_station is None: WSC_station = station_list[station].WSC
      station = station_list[station].HGS; 
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  # resolve folder arguments
  if project_folder is None and project is not None: 
      project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project)
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
  if experiment in exp_aliases:
      lWRF = True 
      # resolve aliases (always return full string format)
      experiment = exp_aliases[experiment] # resolve alias
      exp_name = experiment[:-4]; domain = int(experiment[-2:]) # always full string format
  if experiment in WRF_exps: # i.e. if domain is None
      lWRF = True
      exp_name = experiment
      domain = WRF_exps[experiment].domains # innermost domain
      experiment = '{:s}_d{:02d}'.format(experiment,domain) # used to determine aliases
  # add period extensions of necessary/possible
  if isinstance(run_period, (tuple,list)):
      start_year, end_year = run_period
      if isinstance(period, (int,np.integer)): period = (end_year-period,end_year)
  elif isinstance(period, (tuple,list)) and isinstance(run_period, (int,np.integer)):
      end_year = period[1]
      start_year = end_year - run_period
  else: start_year = end_year = None # not used
  # append conventional period name to name, if it appears to be missing
  if start_year is not None:
      if start_year > 2080:
        if exp_name[-5:] != '-2100': exp_name += '-2100'
        if old_name[-5:] != '-2100': old_name += '-2100'
      elif start_year > 2040:
        if exp_name[-5:] != '-2050': exp_name += '-2050'
        if old_name[-5:] != '-2050': old_name += '-2050'
  if 'exp_name' not in kwargs: kwargs['exp_name'] = old_name # save for later use in name
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
  if end_year is None: end_year = start_year + run_period # period length
  else: run_period = end_year - start_year
  assert isinstance(run_period,(np.integer,int)), run_period
  # construct period information  
  prdstr = '{:04d}-{:04d}'.format(start_year, end_year)
  if 'prdstr' not in kwargs: kwargs['prdstr'] = prdstr # for name expansion (will be capitalized)
  # append conventional period name to name, if it appears to be missing
  if start_year > 2080: prdext = '2100'
  elif start_year > 2040: prdext = '2050'
  elif start_year < 1990: prdext = '1980'
  else:  raise NotImplementedError("Unable to determine period extension for start date '{:d}'.".format(start_year))
  if 'prdext' not in kwargs: kwargs['prdext'] = prdext # for name expansion (will be capitalized)
  # resolve climate input mode (time aggregation; assuming period is just length, i.e. int;)
  if 'clim_mode' not in kwargs: kwargs['clim_mode'] = clim_mode # for name expansion
  # translate climate mode into path convention, while retaining clim_mode
  if clim_mode.lower() in ('timeseries','transient'): clim_dir = 'timeseries'
  elif clim_mode.lower() in ('clim','climatology','periodic'): clim_dir = 'clim_{:02d}'.format(run_period)
  elif clim_mode.lower() in ('mean','annual','steady-state'): clim_dir = 'annual_{:02d}'.format(run_period)
  # call load function from HGS module
  dataset = hgs.loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                              folder=folder, experiment=experiment, period=period, start_date=start_date,
                              project_folder=project_folder, project=project, grid=grid, task=task, 
                              prefix=prefix, clim_dir=clim_dir, WSC_station=WSC_station, basin=basin, 
                              basin_list=basin_list, filename=hgs.station_file, lpad=lpad, **kwargs)
  # slice time axis for period
  if start_year is not None and end_year is not None:
    dataset = dataset(years=(start_year,end_year))
  # add WRF attributes to dataset
  if lWRF:
    for key,value in exp.__dict__.items():
      dataset.atts['WRF_'+key] = value
    dataset.atts['WRF_resolution'] = resolution
  return dataset

# an enhanced ensemble loader
@BatchLoad
def loadHGS_StnEns(station=main_gage, varlist=None, varatts=None, name=None, title=None, period=None, 
                  experiment=None, domain=None, exp_aliases=exp_aliases, run_period=15, clim_mode=None,
                  folder=project_folder_pattern, project_folder=None, obs_period=None,
                  project=project_name, grid=main_grid, task=main_task, prefix=project_prefix, 
                  WSC_station=None, basin=main_basin, basin_list=None, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations instead '''
  if experiment.upper() == 'WSC':
    # translate parameters
    station = station if WSC_station is None else WSC_station
    period = period if obs_period is None else obs_period
    filetype = 'monthly'
    # load gage station with slightly altered parameters
    return loadGageStation_TS(station=station, name=name, title=title, basin=basin, basin_list=basin_list, 
                              varlist=varlist, varatts=varatts, period=period, filetype=filetype)
  else:
    # load HGS simulation
    return loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, period=period, 
                         experiment=experiment, domain=domain, exp_aliases=exp_aliases, run_period=run_period, 
                         clim_mode=clim_mode, folder=folder, project_folder=project_folder, project=project, 
                         grid=grid, task=task, prefix=prefix, WSC_station=WSC_station, basin=basin, 
                         basin_list=basin_list, **kwargs)

# custom gage station loader
def loadGageStation_TS(station=main_gage, name=None, title=None, basin=main_basin, basin_list=None,
                       varlist=None, varatts=None, period=None, filetype='monthly'):
  ''' a wrapper to load gage stations for the GRW with appropriate default values'''
  # load gage station
  if station in station_list: station = station_list[station].WSC
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  dataset = wsc.loadGageStation_TS(basin=basin, station=station, varlist=varlist, varatts=varatts, 
                                   filetype=filetype, folder=None, name=name, basin_list=basin_list)
  if period: dataset = dataset(years=period) # apply slice
  return dataset

# abuse for testing
if __name__ == '__main__':
    
#   test_mode = 'gage_station'
#   test_mode = 'dataset'
  test_mode = 'ensemble'

  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadGageStation_TS(period=(1974,2004), )
    print(ds)
    
  elif test_mode == 'dataset':

    # load single dataset
    ds = loadHGS_StnTS(experiment='NRCan', domain=None, period=(1984,1994), 
                       clim_mode='periodic', lpad=True)
    print(ds)
    print(ds.discharge[:])
    
#     # load single dataset
#     ds = loadHGS_StnTS(experiment='erai-g', domain=0, period=(1984,1994), clim_mode='periodic', )
#     print(ds)
#     assert ds.time[0] == 60, ds.time[:]
    
  elif test_mode == 'ensemble':
    
    # load an esemble of datasets
    ens = loadHGS_StnEns(experiment=['g-ensemble','t-ensemble'], domain=0, clim_mode='mean', 
                         name='{EXP_NAME:s}_{RESOLUTION:s}',
                         period=[(1984,1994),(2050,2060),(2090,2100)], obs_period=(1974,2004),
                         outer_list=['experiment','period'], lensemble=True)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print(ens[0])