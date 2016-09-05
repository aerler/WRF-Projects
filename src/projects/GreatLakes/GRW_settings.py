'''
Created on Sep 4, 2016

This module contains a meta data for HGS simulations for the GRW and a wrapper to load them. 

@author: Andre R. Erler, GPL v3
'''

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
project_folder_pattern = '{PROJECT_FOLDER:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_MODE:s}/{TASK:s}/'

# mapping of WSC station names to HGS hydrograph names
Station = namedtuple('Station', ('HGS','WSC'),)
station_list = dict( # short names for gages in this basin and their HGS/WSC names
    Brantford=Station(HGS='GR_Brantford',WSC='Grand River_Brantford') ) # this is the main gage of the GRW

# experiment aliases (for more systematic access)
exp_aliases = {'erai-g_d00':'erai-g3_d01','erai-t_d00':'erai-t3_d01',
               'g-ensemble_d00':'g3-ensemble_d01','t-ensemble_d00':'t3-ensemble_d01',
               'g-ensemble-2050_d00':'g3-ensemble-2050_d01','t-ensemble-2050_d00':'t3-ensemble-2050_d01',
               'g-ensemble-2100_d00':'g3-ensemble-2100_d01','t-ensemble-2100_d00':'t3-ensemble-2100_d01'}

## wrapper functions to load HGS station timeseries with GRW parameters

# simple dataset loader
def loadHGS_StnTS(station=main_gage, varlist=None, varatts=None, name=None, title=None, 
                  experiment=None, domain=None, exp_aliases=exp_aliases, period=15, clim_mode=None,
                  folder=project_folder_pattern, project_folder=None,
                  project=project_name, grid=main_grid, task=main_task, prefix=project_prefix, 
                  WSC_station=None, basin=main_basin, basin_list=None, resampling='1M', **kwargs):
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
  if domain is not None: 
    exp_name = experiment
    experiment = '{:s}_d{:02d}'.format(experiment,domain) # used to determine aliases
  elif experiment in WRF_exps: # i.e. if domain is None
    exp_name = experiment
    domain = WRF_exps[experiment].domains # innermost domain
    experiment = '{:s}_d{:02d}'.format(experiment,domain) # used to determine aliases
  else:
    exp_name = experiment[:-4] # cut off domain string
    domain = int(experiment[-2:]) # last two digits
  # extract resolution for later use
  if domain == 0: resolution = '90km'
  elif domain == 1: resolution = '30km'
  elif domain == 2: resolution = '10km'
  else: resolution = None
  # resolve aliases (always return full string format)
  if experiment in exp_aliases: 
    experiment = exp_aliases[experiment] # resolve alias
    exp_name = experiment[:-4]; domain = int(experiment[-2:]) # always full string format
  # add period extensions of necessary/possible
  if isinstance(period, (tuple,list)):
    start_year, end_year = period
    period = end_year - start_year # period length
    # append conventional period name to name, if it appears to be missing
    if start_year > 2080 and exp_name[-5:] != '-2100': exp_name += '-2100'
    elif start_year > 2040 and exp_name[-5:] != '-2050': exp_name += '-2050'
  else: start_year = None
  # reassemble experiment name with domain extension
  experiment = '{:s}_d{:02d}'.format(exp_name,domain)
  # validate WRF experiment
  if exp_name not in WRF_exps: raise DatasetError("Unknown WRF experiment: '{}'".format(exp_name))
  exp = WRF_exps[exp_name]
  if domain > exp.domains or domain < 1: raise DatasetError("Invalid domain number: {}".format(domain))
  # get start date from experiment (if not defined in period)
  if start_year is not None: start_date = start_year
  else: start_date = tuple(int(d) for d in exp.begindate.split('-')) # year,month,day tuple  
  # resolve climate input mode (time aggregation; assuming period is just length, i.e. int;)
  if clim_mode.lower() in ('timeseries','transient'): clim_mode = 'timeseries'
  elif clim_mode.lower() in ('clim','climatology','periodic'): clim_mode = 'clim_{:02d}'.format(period)
  elif clim_mode.lower() in ('mean','annual','steady-state'): clim_mode = 'annual_{:02d}'.format(period)
  # call load function from HGS module
  dataset = hgs.loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                              folder=folder, experiment=experiment, period=period, start_date=start_date,
                              project_folder=project_folder, project=project, grid=grid, task=task, 
                              prefix=prefix, clim_mode=clim_mode, WSC_station=WSC_station, basin=basin, 
                              basin_list=basin_list, filename=hgs.station_file, date_parser=hgs.date_parser, 
                              resampling=resampling, **kwargs)
  # add WRF attributes to dataset
  for key,value in exp.__dict__.items():
    dataset.atts['WRF_'+key] = value
  dataset.atts['WRF_resolution'] = resolution
  return dataset

# an enhanced ensemble loader
@BatchLoad
def loadHGS_StnEns(station=main_gage, varlist=None, varatts=None, name=None, title=None, 
                  experiment=None, domain=None, exp_aliases=exp_aliases, period=15, clim_mode=None,
                  folder=project_folder_pattern, project_folder=None, obs_period=None,
                  project=project_name, grid=main_grid, task=main_task, prefix=project_prefix, 
                  WSC_station=None, basin=main_basin, basin_list=None, resampling='1M', **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations instead '''
  if experiment.upper() == 'WSC':
    # translate parameters
    station = station if WSC_station is None else WSC_station
    period = period if obs_period is None else obs_period
    if resampling == '1M': filetype = 'monthly'
    else: raise NotImplementedError("Currently only monthly time-series are supported (not '{}').".format(resampling))
    # load gage station with slightly altered parameters
    return loadGageStation_TS(station=station, name=name, title=title, basin=basin, basin_list=basin_list, 
                              varlist=varlist, varatts=varatts, period=period, filetype=filetype)
  else:
    # load HGS simulation
    return loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                         experiment=experiment, domain=domain, exp_aliases=exp_aliases, period=period, 
                         clim_mode=clim_mode, folder=folder, project_folder=project_folder, project=project, 
                         grid=grid, task=task, prefix=prefix, WSC_station=WSC_station, basin=basin, 
                         basin_list=basin_list, resampling=resampling, **kwargs)

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
    
  test_mode = 'gage_station'
#   test_mode = 'dataset'
#   test_mode = 'ensemble'

  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadGageStation_TS(period=(1979,1994), )
    print(ds)
    
  elif test_mode == 'dataset':
    
    # load single dataset
    ds = loadHGS_StnTS(experiment='erai-g', domain=0, period=(1979,1994), clim_mode='periodic', )
    print(ds)
    
  elif test_mode == 'ensemble':
    
    # load an esemble of datasets
    ens = loadHGS_StnEns(experiment=['g-ensemble','t-ensemble'], domain=0, clim_mode='mean', 
                         name='{EXPERIMENT:s}',
                         period=[(1979,1994),(2045,2060),(2085,2100)], obs_period=(1979,1994),
                         outer_list=['experiment','period'], lensemble=True)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print(ens[0])