'''
Created on Sep 4, 2016

This module contains a meta data for HGS simulations for the GRW and a wrapper to load them. 

@author: Andre R. Erler, GPL v3
'''
import numpy as np
from collections import OrderedDict
import hgs.HGS as hgs # need to prevent name collisions here
import projects.WSC_basins as wsc
from projects.GreatLakes.WRF_experiments import WRF_exps
from collections import namedtuple
from geodata.misc import ArgumentError, DatasetError

project_name = 'GRW'
project_prefix = 'grw_omafra'
# some parameters
gage_scalefactors = 1e-4 # default scalefactor for discharge plots
main_gage  = 'Brantford' # gage station(see station_list for proper HGS/WSC names)
main_basin = 'GRW' # main basin name
main_grid  = 'grw2' # climate data grid
# project folders
project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project_name) # the dataset root folder
project_folder_pattern = '{PROJECT_FOLDER:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_DIR:s}/{TASK:s}/'
station_file_v1 = '{PREFIX:s}o.hydrograph.Station_{STATION:s}.dat' # general HGS naming convention
station_file_v2 = '{PREFIX:s}o.hydrograph.{WSC_ID0:s}.dat' # general HGS naming convention


# plotting parameters for HGS simulations
hgs_plotargs = dict() 
hgs_plotargs['Observations'] = dict(color='#959595') # gray
hgs_plotargs['WSC']          = dict(color='#959595') # gray
hgs_plotargs['HGS (V1)']     = dict(color='green') #, linewidth=3)
hgs_plotargs['NRCan']        = dict(color='green')
hgs_plotargs['NRCan (V1)']   = dict(color='gray', linestyle='--')
hgs_plotargs['NRCan (V2)']   = dict(color='gray')
hgs_plotargs['NRCan (V2k)']  = dict(color='#AAA2D8') # purple
hgs_plotargs['NRCan (V2f)']  = dict(color='red') # red
hgs_plotargs['NRCan (V3f)']  = dict(color='magenta')
hgs_plotargs['NRCan (V3s)']  = dict(color='black') 
hgs_plotargs['NRCan (V3w)']  = dict(color='green')
hgs_plotargs['NRCan (V3m2)'] = dict(color='blue')
hgs_plotargs['NRCan (V3m3)'] = dict(color='purple')
hgs_plotargs['NRCan (V3m4)'] = dict(color='red')
hgs_plotargs['NRCan (V3m5)'] = dict(color='green')
hgs_plotargs['NRCan (V3)']   = dict(color='green')
hgs_plotargs['V3 (Prairies)']  = dict(color='#AAA2D8')
hgs_plotargs['V3 (Maritime)']  = dict(color='#62A1C6')
hgs_plotargs['V3 (Ephemeral)'] = dict(color='#E24B34')
hgs_plotargs['NRCan (Prairies)']  = dict(color='cyan')
hgs_plotargs['NRCan (Ephemeral)'] = dict(color='coral')
hgs_plotargs['NRCan (hires)']     = dict(color='magenta')
# hgs_plotargs['Observations'] = dict(color='#68615E') # dark gray
hgs_plotargs['Steady-State (V1)'] = dict(color='#AAA2D8') # purple
hgs_plotargs['Periodic (V1)']     = dict(color='green') # red
hgs_plotargs['Steady-State (V2)'] = dict(color='#E24B34') # red
hgs_plotargs['Periodic (V2)']     = dict(color='coral') # red
hgs_plotargs['Transient']    = dict(color='#62A1C6') # blue
hgs_plotargs['Steady-State'] = dict(color='#AAA2D8') # purple
hgs_plotargs['Periodic']     = dict(color='#E24B34') # red
hgs_plotargs['WRF 90km']     = dict(color='#AAA2D8') # purple
hgs_plotargs['WRF 30km']     = dict(color='#E24B34') # red
hgs_plotargs['WRF 10km']     = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF (AABC)']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF (Delta)']  = dict(color='#E24B34') # red
hgs_plotargs['WRF T 10km']   = dict(color='#E24B34') # red
hgs_plotargs['WRF G 10km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF T 30km']   = dict(color='#E24B34') # red
hgs_plotargs['WRF G 30km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['T Mean']       = dict(color='#E24B34') # red
hgs_plotargs['G Mean']       = dict(color='black') # blue                                     
hgs_plotargs['Ensemble']     = dict(color='red')
hgs_plotargs['Mean']         = dict(color='black')
hgs_plotargs['1980']         = dict(color='#62A1C6') # blue
hgs_plotargs['2050']         = dict(color='#AAA2D8') # purple
hgs_plotargs['2100']         = dict(color='#E24B34') # red
hgs_plotargs['1979-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs['1984-1994']    = dict(color='#62A1C6') # blue
hgs_plotargs['2045-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs['2050-2060']    = dict(color='#AAA2D8') # purple
hgs_plotargs['2085-2100']    = dict(color='#E24B34') # red
hgs_plotargs['2090-2100']    = dict(color='#E24B34') # red
# adjust line thickness
for plotargs in hgs_plotargs.values(): 
  if 'linewidth' not in plotargs: plotargs['linewidth'] = 1.5
# extended color scheme for ensemble
color_args = {'Ctrl':'blue', 'Ens-A':'purple', 'Ens-B':'green','Ens-C':'coral','90km':'purple','30km':'red','10km':'blue'}
for key,value in color_args.items(): hgs_plotargs[key] = dict(color=value, linewidth=.75)
# add "old" versions
for expname,plotarg in hgs_plotargs.items():
    if expname[-1] != ')' and expname+' (old)' not in hgs_plotargs:
        oldarg = plotarg.copy(); 
        oldarg['linestyle'] = '--' 
        hgs_plotargs[expname+' (old)'] = oldarg

# mapping of WSC station names to HGS hydrograph names
Station = namedtuple('Station', ('HGS','WSC','ylim'),)
station_list = OrderedDict() # this is the main gage of the GRW
station_list['Grand River at Brantford']  = Station(HGS='GR_Brantford',WSC='Grand River_Brantford',ylim=15)
station_list['Grand River at Marsville']  = Station(HGS='GR_Marsville_(near_it)',WSC='Grand River_Marsville',ylim=3)
station_list['Conestogo at Glen Allan']   = Station(HGS='Conestogo_River_at_Glen_Allan',WSC='Conestogo River_Glen Allan',ylim=2)
station_list['Speed River at Guelph']     = Station(HGS='Speed_River_near_Guelph(moved_North)',WSC='Speed River_Guelph',ylim=2)
station_list['Whitemans at Mount Vernon'] = Station(HGS='Station_Whitemans_Creek_near_Mt_Vernon',WSC='Whitemans Creek_Mount Vernon',ylim=1.2)
station_list['Fairchild at Brantford']    = Station(HGS='Fairchild_Creek_near_Brantford',WSC='Fairchild Creek_Brantford',ylim=1.2 )
# look-up tables for WSC/HGS station name conversion                           
WSC_station_list = {stn.WSC:stn.HGS for stn in station_list.values()}
HGS_station_list = {stn.HGS:stn.WSC for stn in station_list.values()}
# short names for gages in this basin and their HGS/WSC names
station_list_etal = dict(**station_list) # not ordered
station_list_etal['Brantford']  = station_list['Grand River at Brantford'] # alias


# experiment aliases (for more systematic access)
exp_aliases = {'erai-g_d00':'erai-g3_d01','erai-t_d00':'erai-t3_d01',
               'g-ensemble_d00':'g3-ensemble_d01','t-ensemble_d00':'t3-ensemble_d01',
               'g-ensemble-2050_d00':'g3-ensemble-2050_d01','t-ensemble-2050_d00':'t3-ensemble-2050_d01',
               'g-ensemble-2100_d00':'g3-ensemble-2100_d01','t-ensemble-2100_d00':'t3-ensemble-2100_d01'}
obs_datasets = ('NRCan','CRU')
gage_datasets = ('wsc','obs','observations')
# ensemble definitions for GRW project 
ensemble_list = {'g-mean':('g-ctrl','g-ens-A','g-ens-B','g-ens-C'),
                 't-mean':('t-ctrl','t-ens-A','t-ens-B','t-ens-C')}
for name,members in ensemble_list.items():
    for prd in ('-2050','-2100'):
        ensemble_list[name+prd] = tuple(member+prd for member in members)
    

## wrapper functions to load HGS station timeseries with GRW parameters

#TODO: the function definitions should be moved into a separate module a la clim and eva, since they are pretty general

# simple dataset loader
def loadHGS_StnTS(experiment=None, domain=None, period=None, varlist=None, varatts=None, name=None, title=None, 
                  exp_aliases=exp_aliases, run_period=None, clim_mode=None, clim_period=None, station=main_gage, 
                  bias_correction=None, folder=project_folder_pattern, project_folder=None, project=project_name, 
                  grid=main_grid, task=None, prefix=project_prefix, WSC_station=None, basin=main_basin, 
                  basin_list=None, lpad=True, ENSEMBLE=None, scalefactors=gage_scalefactors, lWSCID=False, 
                  lold=False, **kwargs):
  ''' Get a properly formatted HGS dataset with a regular time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  experiment = experiment or ENSEMBLE # an alias to support the ensemble loader
  if ENSEMBLE is not None:
      if experiment is not None and experiment != ENSEMBLE: raise ArgumentError(ENSEMBLE,experiment)
      experiment = ENSEMBLE # use as experiment name, but also pass on as kwarg to HGS loader
  if experiment is None or clim_mode is None: raise ArgumentError(experiment,clim_mode)
  # resolve station name
  if station in station_list_etal: 
      if WSC_station is None: WSC_station = station_list_etal[station].WSC
      station = station_list_etal[station].HGS
  elif station is None and WSC_station in WSC_station_list:
      station = WSC_station_list[station]
  elif station in WSC_station_list and WSC_station is None:
      WSC_station = station; station = WSC_station_list[station]
  elif station in HGS_station_list and WSC_station is None: 
      WSC_station =  HGS_station_list[station]
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
  # set some defaults based on version
  if lold:
      station_file = station_file_v1
      if task is None: task = 'hgs_run_wrfpet' if bias_correction else 'hgs_run' 
      if run_period is None: run_period = 15 # in years
      if clim_period is None: clim_period = 15 # in years
  else:
      station_file = station_file_v2
      if task is None: task = 'hgs_run_v3_wrfpet' if lWRF else 'hgs_run_v3'
      if run_period is None: run_period = 10 # in years
      if clim_period is None: clim_period = 15 if lWRF else 30 # in years      
  # add period extensions of necessary/possible
  if isinstance(run_period, (tuple,list)):
      start_year, end_year = run_period
      if isinstance(period, (int,np.integer)): period = (end_year-period,end_year)
  elif isinstance(period, (tuple,list)) and isinstance(run_period, (int,np.integer)):
      end_year = period[1]
      start_year = end_year - run_period
  elif period is None and isinstance(run_period, (int,np.integer)):
      start_year = 1979; end_year = start_year + run_period
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
  elif clim_mode.lower() in ('clim','climatology','periodic'): clim_dir = 'clim_{:02d}'.format(clim_period)
  elif clim_mode.lower() in ('mean','annual','steady-state'): clim_dir = 'annual_{:02d}'.format(clim_period)
  else: raise ArgumentError(clim_mode)
  # select bias-correction method, if applicable
  if bias_correction is not None:
      clim_dir = '{:s}_{:s}'.format(bias_correction,clim_dir)
  # call load function from HGS module
  dataset = hgs.loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, ENSEMBLE=ENSEMBLE, 
                              folder=folder, experiment=experiment, period=period, start_date=start_date, grid=grid,
                              project_folder=project_folder, project=project, task=task, BIAS_CORRECTION=bias_correction,
                              prefix=prefix, clim_dir=clim_dir, WSC_station=WSC_station, basin=basin, lpad=lpad,  
                              basin_list=basin_list, filename=station_file, scalefactors=scalefactors, **kwargs)
  # slice time axis for period
  if start_year is not None and end_year is not None:
    dataset = dataset(years=(start_year,end_year))
  # add WRF attributes to dataset
  if lWRF:
    for key,value in exp.__dict__.items():
      dataset.atts['WRF_'+key] = value
    dataset.atts['WRF_resolution'] = resolution
  return dataset


# custom gage station loader
def loadWSC_StnTS(station=main_gage, name=None, title=None, basin=main_basin, basin_list=None, varlist=None, 
                  varatts=None, period=None, filetype='monthly', scalefactors=gage_scalefactors, **kwargs):
  ''' a wrapper to load gage stations for the GRW with appropriate default values'''
  # load gage station
  if station in station_list: station = station_list[station].WSC # get WSC name for gage station (different from HGS name...)
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  return wsc.loadWSC_StnTS(basin=basin, station=station, varlist=varlist, varatts=varatts, filetype=filetype,
                           name=name, basin_list=basin_list, period=period, scalefactors=scalefactors, **kwargs)
  

# wrapper to load HGS ensembles, otherwise the same
def loadHGS_StnEns(ensemble=None, station=main_gage, varlist=None, varatts=None, name=None, title=None, 
                   period=None, domain=None, exp_aliases=exp_aliases, run_period=None, clim_mode=None,
                   folder=project_folder_pattern, project_folder=None, obs_period=None, clim_period=None, 
                   ensemble_list=ensemble_list, ensemble_args=None, observation_list=gage_datasets, # ensemble and obs lists for project
                   project=project_name, grid=main_grid, task=None, prefix=project_prefix, 
                   WSC_station=None, basin=main_basin, basin_list=None, scalefactors=gage_scalefactors, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations and assemble ensembles '''  
  return hgs.loadHGS_StnEns(ensemble=ensemble, station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                            period=period, run_period=run_period, folder=folder, obs_period=obs_period, clim_period=clim_period, 
                            ensemble_list=ensemble_list, ensemble_args=ensemble_args, observation_list=observation_list, 
                            loadHGS_StnTS=loadHGS_StnTS, loadWSC_StnTS=loadWSC_StnTS, # use local versions of loaders
                            prefix=prefix, WSC_station=WSC_station, basin=basin, basin_list=basin_list, 
                            domain=domain, project_folder=project_folder, project=project, grid=grid, clim_mode=clim_mode, 
                            exp_aliases=exp_aliases, task=task, scalefactors=scalefactors, **kwargs)  


# abuse for testing
if __name__ == '__main__':
    
#   test_mode = 'gage_station'
#   test_mode = 'dataset'
  test_mode = 'ensemble'

  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadWSC_StnTS(period=(1974,2004), station='Speed River at Guelph')
    print(ds)
    
  elif test_mode == 'dataset':

    # load single dataset
    ds = loadHGS_StnTS(experiment='erai-g', domain=2, period=(1984,1994), 
                       clim_mode='periodic', lpad=True, bias_correction='AABC')
#     ds = loadHGS_StnTS(experiment='NRCan', task='hgs_run_v2', 
#                        clim_mode='periodic', lpad=True)
    print(ds)
    print(ds.discharge.mean(), ds.discharge.max(), ds.discharge.min(),)
    print(ds.discharge.plot.units)
#     # load single dataset
#     ds = loadHGS_StnTS(experiment='erai-g', domain=0, period=(1984,1994), clim_mode='periodic', )
#     print(ds)
#     assert ds.time[0] == 60, ds.time[:]
    
  elif test_mode == 'ensemble':
    
    print('')
    print(ensemble_list)
    print('')
    
    # load an esemble of datasets
    ens = loadHGS_StnEns(ensemble=['g-mean','t-mean'], domain=1, clim_mode='clim', 
                         name='{EXP_NAME:s}_{RESOLUTION:s}', title='{Name:s}',
                         period=[(1984,1994),(2050,2060),(2090,2100)], obs_period=(1974,2004),
                         lskipNaN=True, lcheckComplete=True,
                         outer_list=['ensemble','period'], lensemble=True)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print(ens[0])
    print(ens[0].discharge.plot.units)