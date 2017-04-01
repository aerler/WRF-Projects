'''
Created on Feb 2, 2015

Utility functions related to loading basin-averaged data to support hydrological analysis.

@author: Andre R. Erler, GPL v3
'''

# internal imports
from geodata.base import Ensemble, Dataset
from utils.misc import defaultNamedtuple, reverse_enumerate, expandArgumentList
from datasets.common import loadEnsemble, loadEnsembles, BatchLoad, shp_params, stn_params
from datasets.common import observational_datasets, timeseries_datasets
from geodata.misc import ArgumentError, EmptyDatasetError, VariableError
from datasets.WSC import GageStationError, loadGageStation

# some definitions
VL = defaultNamedtuple('VarList', ('vars','files','label'))   
EX = defaultNamedtuple('Experiments', ('name','exps','styles','master','title','reference','target'),
                       defaults=dict(styles=['-','-.','--'], ))  

# wet-day thresholds
wetday_thresholds = [0.2,1,10,20]
wetday_extensions = ['_{:03.0f}'.format(threshold*10) for threshold in wetday_thresholds]

# dataset variables
obs_aliases = ['obs','Obs','observations','Observations']
CRU_vars = ('T2','Tmin','Tmax','dTd','Q2','pet','precip','cldfrc','wetfrq','frzfrq')
WSC_vars = ('runoff','sfroff','ugroff')

# internal method for do slicing for Ensembles and Obs
def _configSlices(slices=None, basins=None, provs=None, shapes=None, stations=None, period=None):
  ''' configure slicing based on basin/province/shape and period arguments '''
  if isinstance(slices,(list,tuple)):
      # handle slice and period lists recursively
      slices = [ _configSlices(slc, basins=basins, provs=provs, shapes=shapes, period=period) for slc in slices ]
  else:
      # the main function that augments a single slice dict
      if slices is None: slices = dict()
      else: slices = slices.copy()
      if shapes is not None:
          if not ( basins is None and provs  is None ): raise ArgumentError
          slices['shape_name'] = shapes
      if basins is not None: 
          if not ( shapes is None and provs  is None ): raise ArgumentError
          slices['shape_name'] = basins  
      if provs  is not None: 
          if not ( basins is None and shapes is None ): raise ArgumentError
          slices['shape_name'] = provs
      if stations is not None:
          slices['station_name'] = stations
      if period is not None:
          slices['years'] = period
  return slices

def _resolveVarlist(varlist=None, filetypes=None, params=None, variable_list=None, lforceList=True):
  # resolve variable list and filetype (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(params) # set required parameters
  filetypes = set() if filetypes is None else set(filetypes)
  for name in varlist: 
    if name in variable_list: 
      variables.update(variable_list[name].vars)
      filetypes.update(variable_list[name].files)
    elif lforceList: 
      raise VariableError("Variable list '{}' does not exist.".format(name))
    else: 
      variables.add(name) 
  # return variables and filetypes as list
  return list(variables), list(filetypes)

# define new load fct. for observations
@BatchLoad
def loadShapeObservations(obs=None, seasons=None, basins=None, provs=None, shapes=None, stations=None, varlist=None, slices=None,
                          aggregation='mean', dataset_mode='time-series', lWSC=True, WSC_period=None, shapetype=None, 
                          variable_list=None, basin_list=None, lforceList=True, obs_ts=None, obs_clim=None, 
                          name=None, title=None, obs_list=None, ensemble_list=None, ensemble_product='inner', **kwargs):
  ''' convenience function to load shape observations based on 'aggregation' and 'varlist' (mainly add WSC gage data) '''
  if obs_list is None: obs_list = observational_datasets
  if name is None: name = 'obs'
  if title is None: title = 'Observations'
  # variables for which ensemble expansion is not supported
  not_supported = ('season','seasons','varlist','mode','dataset_mode','provs','basins','shapes',) 
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
      if name in variable_list: variables.update(variable_list[name].vars)
      elif lforceList: raise VariableError("Variable list '{}' does not exist.".format(name))
      else: variables.add(name)
  variables = list(variables)
  # determine if we need gage dataset
  lWSC = isinstance(basins,basestring) and any([var in WSC_vars for var in variables]) and lWSC # doesn't work if multiple basins are loaded
  # default obs list
  if obs is None: obs = ['Observations',]
  elif isinstance(obs,basestring): obs = [obs]
  elif isinstance(obs,tuple): obs = list(obs)
  elif not isinstance(obs,list): raise TypeError(obs)
  # configure slicing (extract basin/province/shape and period)
  expand_vars = ('basins','stations','provs','shapes','slices') # variables that need to be added to slices (and expanded first)
  if ensemble_list: expand_list = [varname for varname in expand_vars if varname in ensemble_list]
  if ensemble_list and expand_list:
      local_vars = locals(); exp_args = dict()
      for varname in expand_vars: # copy variables to expand right away
          exp_args[varname] = local_vars[varname]
      for varname in expand_list: # remove entries from ensemble expansion
          if  varname != 'slices': ensemble_list.remove(varname) # only 'slices' will continue to be expanded
      if 'slices' not in ensemble_list: ensemble_list.append('slices')
      slices = [_configSlices(**arg_dict) for arg_dict in expandArgumentList(expand_list=expand_list, lproduct=ensemble_product, **exp_args)]
  else:
      slices = _configSlices(slices=slices, basins=basins, provs=provs, shapes=shapes, stations=stations, period=None)
  # substitute default observational dataset and seperate aggregation methods
  iobs = None; clim_ens = None
  for i,obs_name in reverse_enumerate(obs):
      # N.B.: we need to iterate in reverse order, so that deleting items does not interfere with the indexing
      if obs_name in obs_aliases or obs_name not in timeseries_datasets:
          if iobs is not None: raise ArgumentError("Can only resolve one default dataset: {}".format(obs))
          if aggregation == 'mean' and seasons is None and obs_clim is not None: 
              # remove dataset entry from list (and all the arguments)
              del obs[i]; iobs = i # remember position of default obs in ensemble              
              clim_args = kwargs.copy(); slc = slices; shp = shapetype
              # clean up variables for ensemble expansion, if necessary
              if ensemble_list and ensemble_product.lower() == 'inner':
                  if 'names' in ensemble_list:
                      obs_names = [obs_clim]
                      for arg in ensemble_list:
                          if arg in ('slices','shape'): pass # dealt with separately
                          elif arg in not_supported:
                              raise ArgumentError("Expansion of keyword '{:s}' is currently not supported in ensemble expansion.".format(arg))
                          elif arg in kwargs: 
                              clim_args[arg] = kwargs[arg][iobs]; del kwargs[arg][iobs]
                          else: 
                              raise ArgumentError("Keyword '{:s}' not found in keyword arguments.".format(arg))
                      if 'slices' in ensemble_list: slc = slices[iobs]; del slices[iobs]
                      if 'shape' in ensemble_list: shp = shapetype[iobs]; del shapetype[iobs]
                      clim_len = 1 # expect length of climatology ensemble
                  else: 
                      obs_names = obs_clim # no name expansion
                      clim_len = None # expect length of climatology ensemble
                      for arg in ensemble_list:
                          if arg in not_supported:
                              raise ArgumentError("Expansion of keyword '{:s}' is currently not supported in ensemble expansion.".format(arg))
                          elif 'slices' in ensemble_list: l = len(slc) 
                          elif 'shape' in ensemble_list: l = len(shp)
                          elif arg in clim_args: l = len(clim_args[arg])
                          else: raise ArgumentError("Keyword '{:s}' not found in keyword arguments.".format(arg))
                          if clim_len is None: clim_len = l
                          elif l != clim_len: raise ArgumentError(arg,l,clim_len)
              elif ensemble_list and ensemble_product.lower() == 'outer':
                  clim_len = 1
                  for arg in ensemble_list:
                      if arg != 'names':
                        assert isinstance(clim_args[arg],(list,tuple)), clim_args[arg] 
                        clim_len *= len(clim_args[arg])
                  obs_names = [obs_clim] if 'names' in ensemble_list else obs_clim
              else:
                  obs_names = [obs_clim]; clim_len = 1
              # now load climtology instead of time-series and skip aggregation
              try:
                  clim_ens = loadEnsemble(names=obs_names, season=seasons, aggregation=None, slices=slc, varlist=variables, 
                                          ldataset=False, dataset_mode='climatology', shape=shp,
                                          ensemble_list=ensemble_list, ensemble_product=ensemble_product, 
                                          obs_list=obs_list, basin_list=basin_list, **clim_args)
                  assert len(clim_ens) == clim_len, clim_ens
              except EmptyDatasetError: pass
          else: 
              obs[i] = obs_ts # trivial: just substitute default name and load time-series
  # prepare and load ensemble of observations
  if len(obs) > 0:
      if len(obs) == 1 and ensemble_list and 'names' not in ensemble_list: obs = obs[0]
      try:
          obsens = loadEnsemble(names=obs, season=seasons, aggregation=aggregation, slices=slices,
                                varlist=variables, ldataset=False, dataset_mode=dataset_mode, 
                                shape=shapetype, obs_list=obs_list, basin_list=basin_list, 
                                ensemble_list=ensemble_list, ensemble_product=ensemble_product, **kwargs)          
      except EmptyDatasetError:
          obsens = Ensemble(name=name, title=title, obs_list=obs_list, basetype=Dataset)
  else: 
      obsens = Ensemble(name=name, title=title, obs_list=obs_list, basetype=Dataset)
  # add default obs back in if they were removed earlier
  if clim_ens is not None:
      for clim_ds in clim_ens[::-1]: # add observations in correct order: adding backwards allows successive insertion ...
          obsens.insertMember(iobs,clim_ds) # ... at the point where the name block starts
  # load stream gage data from WSC; should not interfere with anything else; append to ensemble
  if lWSC: # another special case: river hydrographs
      try:
          if aggregation is not None and seasons is None: dataset_mode = 'climatology' # handled differently with gage data
          if WSC_period is None: WSC_period = kwargs.get('obs_period',kwargs.get('period',None))
          dataset = loadGageStation(basin=basins, varlist=['runoff'], aggregation=aggregation, period=WSC_period, 
                                    mode=dataset_mode, filetype='monthly', basin_list=basin_list, lfill=True, lexpand=True) # always load runoff/discharge
          if seasons:
              method = aggregation if aggregation.isupper() else aggregation.title() 
              if aggregation: dataset = getattr(dataset,'seasonal'+method)(season=seasons, taxis='time')
              else: dataset = dataset.seasonalSample(season=seasons)
          if slices is not None: dataset = dataset(**slices) # slice immediately
          obsens += dataset.load()
      except GageStationError: 
          pass # just ignore, if gage station data is missing 
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return obsens


# define new load fct. for experiments (not intended for observations)
@BatchLoad
def loadShapeEnsemble(names=None, seasons=None, basins=None, provs=None, shapes=None, varlist=None, 
                      aggregation='mean', slices=None, shapetype=None, filetypes=None, 
                      period=None, obs_period=None, WSC_period=None, name=None, title=None,
                      variable_list=None, WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, 
                      basin_list=None, lforceList=True, obs_list=None, obs_ts=None, obs_clim=None, 
                      ensemble_list=None, ensemble_product='inner', **kwargs):
  ''' convenience function to load shape ensembles (in Ensemble container) or observations; kwargs are passed to loadEnsembleTS '''
  names = list(names) # make a new list (copy)
  # separate observations
  if obs_list is None: obs_list = observational_datasets
  obs_names = []; iobs = []; ens_names = []; iens = []
  for i,name in enumerate(names):
      if name in obs_list or name in obs_aliases:
          obs_names.append(name); iobs.append(i)          
      else: 
          ens_names.append(name); iens.append(i)
  assert len(iens) == len(ens_names) and len(iobs) == len(obs_names) 
  if len(obs_names) > 0:       
      # assemble arguments
      obs_args = dict(obs=obs_names, seasons=seasons, basins=basins, provs=provs, shapes=shapes, varlist=varlist, 
                      slices=slices, aggregation=aggregation, shapetype=shapetype, 
                      period=period, obs_period=obs_period, obs_ts=obs_ts, obs_clim=obs_clim, 
                      variable_list=variable_list, basin_list=basin_list, WSC_period=WSC_period,
                      ensemble_list=ensemble_list, ensemble_product=ensemble_product, **kwargs)
      # check if we have to modify to preserve ensemble_list expansion
      if ensemble_list and ensemble_product == 'inner' and 'names' in ensemble_list and len(ensemble_list) > 1: 
          for key in ensemble_list:
              if key != 'names':
                  ens_list = obs_args[key]
                  obs_args[key] = [ens_list[i] for i in iobs]
      # observations for basins require special treatment to merge basin averages with gage values
      # load observations by redirecting to appropriate loader function
      obsens = loadShapeObservations(name=name, title=title, obs_list=obs_list, **obs_args)
  else: obsens = []
  if len(ens_names) > 0: # has to be a list
      # prepare arguments
      variables, filetypes = _resolveVarlist(varlist=varlist, filetypes=filetypes, 
                                             params=shp_params, variable_list=variable_list, lforceList=lforceList)
      # configure slicing (extract basin/province/shape and period)
      slices = _configSlices(slices=slices, basins=basins, provs=provs, shapes=shapes, period=period)
      # assemble arguments
      ens_args = dict(names=ens_names, season=seasons, slices=slices, varlist=variables, shape=shapetype, 
                      aggregation=aggregation, period=period, obs_period=obs_period, 
                      WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, filetypes=filetypes, 
                      ensemble_list=ensemble_list, ensemble_product=ensemble_product, **kwargs)
      # check if we have to remove obs datasets to preserve ensemble_list expansion
      if ensemble_list and ensemble_product == 'inner' and 'names' in ensemble_list and len(ensemble_list) > 1: 
          for key in ensemble_list:
              if key != 'names':
                  ens_list = ens_args[key]
                  ens_args[key] = [ens_list[i] for i in iens]
      # load ensemble (no iteration here)
      shpens = loadEnsemble(name=name, title=title, obs_list=obs_list, **ens_args)
  else: shpens = Ensemble(name=name, title=title, basetype='Dataset')
  # get resolution tag (will be added below)
  res = None
  for member in shpens:
      if 'resstr' in member.atts:
          if res is None: res = member.atts['resstr']
          elif res != member.atts['resstr']:
              res = None; break # no common resolution
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  if len(obsens) > 0 and len(shpens) > 0:
      for name,i in zip(obs_names,iobs): 
          shpens.insertMember(i,obsens[name]) # add known observations in correct order
          del obsens[name] # remove the ones we already know from list, so we can deal with the rest
      j = i + 1 # add remaining obs datasets after last one
      for i,obs in enumerate(obsens): shpens.insertMember(j+i,obs)
  elif len(obsens) > 0 and len(shpens) == 0:
      shpens = obsens
  shpens.resolution = res # ad resolution tag now, to make sure it is there 
  return shpens


# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
@BatchLoad
def loadStationEnsemble(names=None, seasons=None, provs=None, clusters=None, varlist=None, aggregation='mean', 
                        constraints=None, filetypes=None, cluster_name=None, stationtype=None, 
                        load_list=None, lproduct='outer', WRF_exps=None, CESM_exps=None, WRF_ens=None, 
                        CESM_ens=None, variable_list=None, default_constraints=None, lforceList=True, 
                        obs_list=None, obs_ts=None, master=None, **kwargs):
  ''' convenience function to load station data for ensembles (in Ensemble container); kwargs are passed to loadEnsembleTS '''
  if obs_list is None: obs_list = observational_datasets
  if load_list: load_list = load_list[:] # use a copy, since the list may be modified
  # figure out varlist  
  if isinstance(varlist,basestring) and not stationtype:
      if varlist.lower().find('prec') >= 0: 
        stationtype = 'ecprecip'
      elif varlist.lower().find('temp') >= 0: 
        stationtype = 'ectemp'
      else: raise ArgumentError, varlist
  if not isinstance(stationtype,basestring): raise ArgumentError, stationtype # not inferred
  if clusters and not cluster_name: raise ArgumentError
  params = stn_params  + [cluster_name] if cluster_name else stn_params # need to load cluster_name!
  variables, filetypes =  _resolveVarlist(varlist=varlist, filetypes=filetypes, 
                                          params=params, variable_list=variable_list, lforceList=lforceList)
  # replace default observations
  names = [obs_ts if name in obs_aliases else name for name in names]
  if master in obs_aliases: master = obs_ts
  # prepare arguments
  if provs or clusters:
    constraints = default_constraints.copy() if constraints is None else constraints.copy()
    constraint_list = []
    if load_list and 'provs' in load_list and 'clusters' in load_list: 
      raise ArgumentError, "Cannot expand 'provs' and 'clusters' at the same time."
    # figure out proper handling of provinces
    if provs:
      if not load_list or 'prov' not in load_list: 
        constraints['prov'] = provs; provs = None
      else:  
        if len(constraint_list) > 0: raise ArgumentError, "Cannot expand multiple keyword-constraints at once."
        for prov in provs:
          tmp = constraints.copy()
          tmp['prov'] = prov
          constraint_list.append(tmp)
        load_list[load_list.index('prov')] = 'constraints'
        constraints = constraint_list; provs = None
    # and analogously, handling of clusters!
    if clusters:
      if not load_list or 'cluster' not in load_list: 
        constraints['cluster'] = clusters; clusters = None
        if cluster_name: constraints['cluster_name'] = cluster_name
      else:  
        if len(constraint_list) > 0: raise ArgumentError, "Cannot expand multiple keyword-constraints at once."
        for cluster in clusters:
          tmp = constraints.copy()
          tmp['cluster'] = cluster
          if cluster_name: tmp['cluster_name'] = cluster_name # will be expanded next to cluster index
          constraint_list.append(tmp)
        load_list[load_list.index('cluster')] = 'constraints'
        constraints = constraint_list; clusters = None
  # load ensemble (no iteration here)
  stnens = loadEnsembles(names=names, season=seasons, prov=provs, station=stationtype, varlist=variables, 
                         aggregation=aggregation, constraints=constraints, filetypes=filetypes, 
                         WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, 
                         load_list=load_list, lproduct=lproduct, lcheckVar=False, master=master, **kwargs)
  # get resolution tag (will be added below)
  res = None
  for member in stnens:
      if 'resstr' in member.atts:
          if res is None: res = member.atts['resstr']
          elif res != member.atts['resstr']:
              res = None; break # no common resolution
  stnens.resolution = res # ad resolution tag now, to make sure it is there 
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return stnens


## abuse main section for testing
if __name__ == '__main__':
  
#   from projects.WesternCanada.analysis_settings import exps_rc, variables_rc, loadShapeObservations
#   from projects.WesternCanada.analysis_settings import loadShapeEnsemble, loadStationEnsemble
  from projects.GreatLakes.analysis_settings import exps_rc, variables_rc, loadShapeObservations
  from projects.GreatLakes.analysis_settings import loadShapeEnsemble, loadStationEnsemble
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

#   test = 'obs_timeseries'
  test = 'basin_timeseries'
#   test = 'station_timeseries'
#   test = 'province_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'obs_timeseries':
    
    # some settings for tests
    basins = ['GRW','GLB'] #; period = (1979,1994)
    stations = ['Brantford','Whitemans Creek_Mount Vernon']
    names = [stn.split('_')[-1] for stn in stations]
    varlist = ['precip','runoff',]

#     shpens = loadShapeObservations(obs='WSC', lWSC=False, name_tags=names, basins=basins, varlist=varlist, 
#                                    period=(1970,2000), ensemble_list=['stations','name_tags',], stations=stations, 
#                                    #WSC_period=(1970,2000), lWSC=True, #obs_clim='NRCan', obs_ts='CRU',
#                                    aggregation='mean', load_list=['varlist'],)
#     assert len(shpens) == len(varlist)
    shpens = loadShapeObservations(obs='Obs', lWSC=True, name_tags=None, basins=basins, varlist=varlist, 
                                   obs_period=(1970,2000), WSC_period=(1970,2000), 
                                   aggregation='mean', load_list=['varlist','basins'],)
    assert len(shpens) == len(varlist)*len(basins)
    # print diagnostics
    print(shpens[0]); print('')
    print(shpens[0].name); print('')
    print(shpens[0][0])
#     for i,basin in enumerate(basins): 
      #for ds in shpens[i]: print ds.atts.shape_name
#       for ds in shpens[i]: print ds.atts.shape_name
#       assert all(ds.atts.shape_name == basin for ds in shpens[i])
  
  
  # test load function for basin ensemble time-series
  elif test == 'basin_timeseries':
    
    # some settings for tests
#     exp = 'ctrl-obs'; basins = ['SSR'] 
    exp = 'erai'; basins = ['GLB','GRW',] 
    exps = exps_rc[exp].exps; #exps = ['Unity']
    aggregation = 'mean'
    grid = 'grw2'; period = (1979,1994); bias_correction = 'AABC'; aggregation = None; dataset_mode = 'climatology'
    varlists = ['precip','runoff','pet']; red = dict(i_s='mean'); domain = 1
    
    shpens = loadShapeEnsemble(names=exps, basins=basins, varlist=varlists, reduction=red, 
                               grid=grid, period=period, bias_correction=bias_correction, 
                               dataset_mode=dataset_mode, aggregation=aggregation,
                               load_list=['basins','varlist'], lproduct='outer', domain=domain,)
    assert shpens[0].resolution == ( '10km' if domain == 2 else '30km' ), shpens.resolution
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(varlists)
    print shpens[0][1]
      

  # test load function for station ensemble
  elif test == 'station_timeseries':
      
    # station selection criteria
    constraints_rc = dict()
    constraints_rc['min_len'] = 15 # for valid climatology
    constraints_rc['lat'] = (40,50) 
    constraints_rc['max_zerr'] = 100 # reduce sample size
    constraints_rc['prov'] = ('ON')
    constraints_rc['end_after'] = 1980
  
    # some settings for tests
    provs = None; clusters = None; lensembleAxis = False; sample_axis = None; lflatten = False
#     exp = 'val'; exps = ['EC', 'erai-max', 'max-ctrl']; provs = ('BC','AB')
#     exp = 'max-all'; exps = exps_rc[exp]; provs = ('BC','AB')
    exps = ['erai-g']; provs = ['ON']
    seasons = ['summer']; lfit = True; lrescale = True; lbootstrap = False
    lflatten = False; lensembleAxis = True
    varlist = ['MaxPrecip_1d', 'MaxPrecip_5d','MaxPreccu_1d']; filetypes = ['hydro']
#     varlist = ['precip', 'pet']; filetypes = ['hydro','aux']
    stnens = loadStationEnsemble(names=exps, provs=provs, clusters=clusters, varlist=varlist,  
                                 seasons=seasons, master=None, stationtype='ecprecip',
                                 domain=2, lensembleAxis=lensembleAxis, filetypes=filetypes,                                 
                                 variable_list=variables_rc, default_constraints=constraints_rc,
                                 load_list=['seasons','provs'], lproduct='outer', lforceList=False)
    # print diagnostics
    print stnens[0]; print ''
    assert len(stnens) == len(seasons)
    print stnens[0][0]
    print stnens[0][0].MaxPrecip_1d.mean()

      
  # test load function for province ensemble climatology
  if test == 'province_climatology':
    
    # some settings for tests
    exp = 'val'; exps = exps_rc[exp].exps
    provs = ['ON',] 
    varlists = ['precip','runoff']; aggregation = 'mean'

    shpens = loadShapeEnsemble(names=exps, provs=provs, varlist=varlists, aggregation=aggregation,
                               period=(1979,1994), obs_period=(1970,2000), # this does not work properly with just a number...
                               load_list=['provs',], lproduct='outer', filetypes=['srfc'],
                               variable_list=variables_rc)
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(provs)
    assert shpens[0][1].time.coord[0] == 1
    for i,basin in enumerate(provs):
      i0 = i*len(varlists); ie = len(varlists)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])    
    print ''; print shpens[0][0]
    