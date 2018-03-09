'''
Created on Sep 4, 2016

This module contains a meta data for HGS simulations for the GRW and a wrapper to load them. 

@author: Andre R. Erler, GPL v3
'''
import numpy as np
import hgs.HGS as hgs # need to prevent name collisions here
import projects.WSC_basins as wsc
from collections import namedtuple
from geodata.misc import ArgumentError, DatasetError

# project folders
project_folder_pattern = hgs.root_folder+'/{PROJECT:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_DIR:s}/{TASK:s}/'
station_file = '{PREFIX:s}o.hydrograph.{WSC_ID0:s}.dat' # general HGS naming convention

# mapping of WSC station names to HGS hydrograph names
Station = namedtuple('Station', ('HGS','WSC','ylim'),)

# plotting parameters for HGS simulations
hgs_plotargs = dict() 
# stream observations
hgs_plotargs['Observations'] = dict(color='#959595') # gray
hgs_plotargs['Obs.']         = dict(color='#959595') # gray
hgs_plotargs['WSC Obs.']     = dict(color='#959595') # gray
hgs_plotargs['WSC']          = dict(color='#959595') # gray
# NRCan forcing, different model versions
hgs_plotargs['NRCan']        = dict(color='green')
hgs_plotargs['NRCan (Prairies)']  = dict(color='cyan')
hgs_plotargs['NRCan (Ephemeral)'] = dict(color='coral')
hgs_plotargs['NRCan (hires)']     = dict(color='magenta')
# Temporal Aggregation
hgs_plotargs['Transient']    = dict(color='#62A1C6') # blue
hgs_plotargs['Steady-State'] = dict(color='#AAA2D8') # purple
hgs_plotargs['Periodic']     = dict(color='#E24B34') # red
# WRF forcing
hgs_plotargs['WRF 10km']     = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF 30km']     = dict(color='#E24B34') # red
hgs_plotargs['WRF 90km']     = dict(color='#AAA2D8') # purple
hgs_plotargs['WRF (AABC)']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF (Delta)']  = dict(color='#E24B34') # red
hgs_plotargs['WRF G 10km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF G 30km']   = dict(color='#E24B34') # red                                      
hgs_plotargs['WRF G 90km']   = dict(color='#AAA2D8') # purple                                   
hgs_plotargs['WRF T 10km']   = dict(color='#62A1C6') # blue    
hgs_plotargs['WRF T 30km']   = dict(color='#E24B34') # red     
hgs_plotargs['WRF T 90km']   = dict(color='#AAA2D8') # purple  
# hgs_plotargs['WRF T 10km']   = dict(color='#E24B34') # red
# hgs_plotargs['WRF G 10km']   = dict(color='#62A1C6') # blue                                     
# hgs_plotargs['WRF T 30km']   = dict(color='#E24B34') # red
# hgs_plotargs['WRF G 30km']   = dict(color='#62A1C6') # blue                                     
# hgs_plotargs['WRF T 90km']   = dict(color='#E24B34') # red
# hgs_plotargs['WRF G 90km']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['WRF T']        = dict(color='#E24B34') # red
hgs_plotargs['WRF G']        = dict(color='#62A1C6') # blue                                     
hgs_plotargs['T Ensemble']   = dict(color='#E24B34') # red
hgs_plotargs['G Ensemble']   = dict(color='#62A1C6') # blue                                     
hgs_plotargs['ERAI-T']       = dict(color='#E24B34') # red
hgs_plotargs['ERAI-G']       = dict(color='#62A1C6') # blue                                     
hgs_plotargs['T Mean']       = dict(color='red') # red
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
  if 'linewidth' not in plotargs: plotargs['linewidth'] = 1.
# extended color scheme for ensemble
color_args = {'Ctrl':'blue', 'Ens-A':'purple', 'Ens-B':'green','Ens-C':'coral','90km':'purple','30km':'red','10km':'blue'}
for key,value in color_args.items(): hgs_plotargs[key] = dict(color=value, linewidth=.75)
# add "old" versions
for expname,plotarg in hgs_plotargs.items():
    if expname[-1] != ')' and expname+' (old)' not in hgs_plotargs:
        oldarg = plotarg.copy(); 
        oldarg['linestyle'] = '--' 
        hgs_plotargs[expname+' (old)'] = oldarg


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
    


#TODO: the function definitions should be moved into a separate module a la clim and eva, since they are pretty general

# helper function to determine experiment parameters
def experimentParameters(experiment=None, domain=None, clim_mode=None, bias_correction=None, 
                         clim_period=None, run_period=None, period=None, exp_aliases=exp_aliases,
                         task='hgs_run', WRF_exps=None, **kwargs):
    # figure out experiment name and domain
    lWRF = False # default: not a WRF experiment
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
    # resolve climate input mode (time aggregation; assuming period is just length, i.e. int;)
    if 'clim_mode' not in kwargs: kwargs['clim_mode'] = clim_mode # for name expansion
    # translate climate mode into path convention, while retaining clim_mode
    if clim_mode.lower() in ('ts','trans','timeseries','transient'): clim_dir = 'timeseries'; lts = True
    elif clim_mode.lower() in ('mn','norm','clim','peri','normals','climatology','periodic'): clim_dir = 'clim'; lts = False
    elif clim_mode.lower() in ('ss','const','mean','annual','steady-state'): clim_dir = 'annual'; lts = False
    else: raise ArgumentError(clim_mode)
    # set some defaults based on version
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
            start_year = 1979; end_year = start_year + run_period
            period = (end_year-period,end_year)
        else: raise ArgumentError(period,run_period)
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

            
## wrapper functions to load HGS station timeseries with GRW parameters

# simple dataset loader
def loadHGS_StnTS(experiment=None, domain=None, period=None, varlist=None, varatts=None, name=None, 
                  title=None, run_period=None, clim_mode=None, clim_period=None, 
                  station=None, well=None, bias_correction=None, project_folder=None, 
                  task=None, WSC_station=None, Obs_well=None, lpad=True, lWSCID=False,
                  project=None, folder=project_folder_pattern, station_file=station_file,
                  grid=None, conservation_authority=None, basin=None, basin_list=None, scalefactors=None, 
                  main_gage=None, station_list=None, experimentParameters=experimentParameters, **kwargs):
  ''' Get a properly formatted HGS dataset with a regular time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if experiment is None or clim_mode is None: raise ArgumentError(experiment,clim_mode)
  # resolve folder arguments
  if project_folder is None and project is not None: 
      project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project)
  
  # resolve station/well name
  if WSC_station and not station: station = WSC_station
  if not station and not well: station = main_gage
  if well and station: raise ArgumentError
  elif station:
      if station_list and station in station_list: 
          if WSC_station is None: WSC_station = station_list[station].WSC
          station = station_list[station].HGS if hasattr(station_list[station], 'HGS') else '{WSC_ID0:s}'
      if WSC_station is None:
          raise NotImplementedError
  elif well:
      if Obs_well is None:
          raise NotImplementedError
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  
  # determine various WRF/climatology related parameters
  params = experimentParameters(experiment=experiment, domain=domain, clim_mode=clim_mode, 
                                bias_correction=bias_correction, task=task, 
                                clim_period=clim_period, run_period=run_period, period=period, **kwargs)
  ( lWRF, task,experiment,exp, clim_dir, start_date,start_year,end_year, 
                                      clim_period,run_period,period, kwargs ) = params
    
  # call load function from HGS module
  dataset = hgs.loadHGS_StnTS(station=station, well=well, varlist=varlist, varatts=varatts, name=name, 
                              title=title, folder=folder, EXPERIMENT=experiment, 
                              period=period, run_period=run_period, start_date=start_date, 
                              PROJECT_FOLDER=project_folder, PROJECT=project, TASK=task, 
                              BIAS_CORRECTION=bias_correction, CLIM_DIR=clim_dir, 
                              WSC_station=WSC_station, Obs_well=Obs_well, basin=basin, lpad=lpad, 
                              GRID=grid, conservation_authority=conservation_authority,
                              basin_list=basin_list, filename=station_file, scalefactors=scalefactors, 
                              **kwargs)
  
  # slice time axis for period
  if start_year is not None and end_year is not None:
    dataset = dataset(years=(start_year,end_year))
  # add WRF attributes to dataset
  if lWRF:
    for key,value in exp.__dict__.items():
      dataset.atts['WRF_'+key] = value
    dataset.atts['WRF_resolution'] = kwargs['resolution']
  return dataset


# custom gage station loader
def loadWSC_StnTS(station=None, name=None, title=None, basin=None, basin_list=None, varlist=None, varatts=None,  
                  period=None, filetype='monthly', scalefactors=None, station_list=None, **kwargs):
  ''' a wrapper to load gage stations for the GRW with appropriate default values'''
  # load gage station
  if station_list and station in station_list: station = station_list[station].WSC # get WSC name for gage station (different from HGS name...)
  if basin_list is None: basin_list = wsc.basin_list # default basin list
  return wsc.loadWSC_StnTS(basin=basin, station=station, varlist=varlist, varatts=varatts, filetype=filetype,
                           name=name, basin_list=basin_list, period=period, scalefactors=scalefactors, **kwargs)
  

# wrapper to load HGS ensembles, otherwise the same
def loadHGS_StnEns(ensemble=None, station=None, varlist=None, varatts=None, name=None, title=None, 
                   period=None, domain=None, exp_aliases=exp_aliases, run_period=None, clim_mode=None,
                   folder=project_folder_pattern, project_folder=None, obs_period=None, clim_period=None, 
                   ensemble_list=ensemble_list, ensemble_args=None, observation_list=gage_datasets, # ensemble and obs lists for project
                   project=None, grid=None, task=None, bias_correction=None, conservation_authority=None,
                   WSC_station=None, basin=None, basin_list=None, scalefactors=None, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations and assemble ensembles '''  
  return hgs.loadHGS_StnEns(ensemble=ensemble, station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                            period=period, run_period=run_period, folder=folder, obs_period=obs_period, clim_period=clim_period, 
                            ensemble_list=ensemble_list, ensemble_args=ensemble_args, observation_list=observation_list, 
                            loadHGS_StnTS=loadHGS_StnTS, loadWSC_StnTS=loadWSC_StnTS, # use local versions of loaders
                            WSC_station=WSC_station, basin=basin, basin_list=basin_list, 
                            bias_correction=bias_correction, conservation_authority=conservation_authority,
                            domain=domain, project_folder=project_folder, project=project, grid=grid, clim_mode=clim_mode, 
                            exp_aliases=exp_aliases, task=task, scalefactors=scalefactors, **kwargs)  


# wrapper for HGS binary load function, which defines experiment folders based on WRF experiments
def loadHGS(experiment=None, varlist=None, name=None, title=None, basin=None, lkgs=False,  
            domain=1, clim_mode=None, clim_period=None, bias_correction=None, task='hgs_run', grid=None, 
            mode='climatology', file_mode='last_12', file_pattern='{PREFIX}o.head_olf.????', t_list=None, 
            varatts=None, constatts=None, project_folder=None, project=None, folder=project_folder_pattern,
            basin_list=None, metadata=None, conservation_authority=None, WRF_exps=None, 
            experimentParameters=experimentParameters, **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  
  if experiment is None or clim_mode is None: raise ArgumentError(experiment,clim_mode)
  # resolve folder arguments
  if project_folder is None and project is not None: 
      project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project)

  params = experimentParameters(experiment=experiment, domain=domain, clim_mode=clim_mode, 
                                clim_period=clim_period, run_period=None, period=None, 
                                bias_correction=bias_correction, exp_aliases=exp_aliases,
                                task=task, WRF_exps=WRF_exps, **kwargs)
  ( lWRF, task,experiment,exp, clim_dir, start_date,start_year,end_year, 
                                      clim_period,run_period,period, kwargs ) = params
  del start_date,start_year,end_year,run_period,period
  if basin_list is None: basin_list = wsc.basin_list # default basin list                                      
  # call load function from HGS module
  dataset = hgs.loadHGS(varlist=varlist, folder=folder, name=name, title=title, basin=basin, 
                        EXPERIMENT=experiment, PROJECT_FOLDER=project_folder, GRID=grid,
                        PROJECT=project, TASK=task, BIAS_CORRECTION=bias_correction, CLIM_DIR=clim_dir,
                        mode=mode, file_mode=file_mode, file_pattern=file_pattern, t_list=t_list, lkgs=lkgs, 
                        varatts=varatts, constatts=constatts, basin_list=basin_list, metadata=metadata, 
                        conservation_authority=conservation_authority)
   
  # add WRF attributes to dataset
  if lWRF:
    for key,value in exp.__dict__.items():
      dataset.atts['WRF_'+key] = value
    dataset.atts['WRF_resolution'] = kwargs['resolution']
  return dataset
                                  

# abuse for testing
if __name__ == '__main__':
    
#   test_mode = 'gage_station'
  test_mode = 'dataset'
#   test_mode = 'ensemble'

  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadWSC_StnTS(period=(1979,2009),) # station='Nith River at Canning')
    print(ds)
    
  elif test_mode == 'dataset':

    # load single dataset
    ds = loadHGS_StnTS(experiment='erai-g', domain=2, period=(1984,1994), 
                       well='W424', z_aggregation=None, z_layers=None,
                       clim_mode='periodic', lpad=True, bias_correction='AABC')
#     ds = loadHGS_StnTS(experiment='NRCan', task='hgs_run_v2', 
#                        clim_mode='periodic', lpad=True)
    print(ds)
    if 'discharge' in ds:
        print('\n')
        print(ds.discharge)
        print(ds.discharge.mean(), ds.discharge.max(), ds.discharge.min(),)
        print(ds.discharge.plot.units)
    if 'head' in ds:
        print('\n')
        print(ds.head)
        print(ds.head.mean(), ds.head.max(), ds.head.min(),)
        print(ds.head.plot.units)
    if ds.hasAxis('z'):
        print('\n')
        print(ds.z)
        print(ds.z[:])
    
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
    print(ens[0])