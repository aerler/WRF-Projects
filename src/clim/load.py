'''
Created on Feb 2, 2015

Utility functions related to loading basin-averaged data to support hydrological analysis.

@author: Andre R. Erler, GPL v3
'''

# internal imports
from geodata.base import Ensemble, Dataset
from utils.misc import defaultNamedtuple
from datasets.common import loadEnsembleTS, BatchLoad, loadDataset, shp_params, stn_params
from geodata.misc import ArgumentError, EmptyDatasetError
from datasets.WSC import GageStationError, loadGageStation

# some definitions
VL = defaultNamedtuple('VarList', ('vars','files','label'))   
EX = defaultNamedtuple('Experiments', ('name','exps','styles','master','title','reference','target'),
                       defaults=dict(styles=['-','-.','--'], ))  

# wet-day thresholds
wetday_thresholds = [0.2,1,10,20]
wetday_extensions = ['_{:03.0f}'.format(threshold*10) for threshold in wetday_thresholds]

# dataset variables
CRU_vars = ('T2','Tmin','Tmax','dTd','Q2','pet','precip','cldfrc','wetfrq','frzfrq')
WSC_vars = ('runoff','sfroff','ugroff')

# default shape type
shapetype = 'glbshp'

# internal method for do slicing for Ensembles and Obs
def _configSlices(slices=None, basins=None, provs=None, shapes=None, period=None):
  ''' configure slicing based on basin/province/shape and period arguments '''
  if slices is None: slices = dict()
  if shapes is not None:
    if not ( basins is None and provs  is None ): raise ArgumentError
    slices['shape_name'] = shapes
  if basins is not None: 
    if not ( shapes is None and provs  is None ): raise ArgumentError
    slices['shape_name'] = basins  
  if provs  is not None: 
    if not ( basins is None and shapes is None ): raise ArgumentError
    slices['shape_name'] = provs
  if period is not None:
    if slices is None: slices = dict()
    slices['years'] = period
  return slices

def _resolveVarlist(varlist=None, filetypes=None, params=None, variable_list=None):
  # resolve variable list and filetype (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(params) # set required parameters
  filetypes = set() if filetypes is None else set(filetypes)
  for name in varlist: 
    if name in variable_list: 
      variables.update(variable_list[name].vars)
      filetypes.update(variable_list[name].files)
    else: variables.add(name) 
  # return variables and filetypes as list
  return list(variables), list(filetypes)

# define new load fct. for observations
@BatchLoad
def loadShapeObservations(obs=None, seasons=None, basins=None, provs=None, shapes=None, varlist=None, slices=None,
                          aggregation='mean', shapetype=shapetype, period=None, variable_list=None, basin_list=None, **kwargs):
  ''' convenience function to load shape observations; the main function is to select sensible defaults 
      based on 'varlist', if no 'obs' are specified '''
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
    if name in variable_list: variables.update(variable_list[name].vars)
    else: variables.add(name)
  variables = list(variables)
  # figure out default datasets
  if obs is None: obs = 'Observations'
  lUnity = lCRU = lWSC = False
  if obs[:3].lower() in ('obs','wsc'):    
    if any(var in CRU_vars for var in variables): 
      if aggregation == 'mean' and seasons is None: 
        lUnity = True; obs = []
    if basins and any([var in WSC_vars for var in variables]):
      if aggregation.lower() in ('mean','std','sem','min','max') and seasons is None: 
        lWSC = True; obs = []
  if not isinstance(obs,(list,tuple)): obs = (obs,)
  # configure slicing (extract basin/province/shape and period)
  slices = _configSlices(slices=slices, basins=basins, provs=provs, shapes=shapes, period=period)
  if slices is not None:
    noyears = slices.copy(); noyears.pop('years',None) # slices for climatologies
  # prepare and load ensemble of observations
  obsens = Ensemble(name='obs',title='Observations', basetype=Dataset)
  if len(obs) > 0: # regular operations with user-defined dataset
    try:
      ensemble = loadEnsembleTS(names=obs, season=seasons, aggregation=aggregation, slices=slices, 
                                varlist=variables, shape=shapetype, ldataset=False, **kwargs)
      for ens in ensemble: obsens += ens
    except EmptyDatasetError: pass 
  if lUnity: # load Unity data instead of averaging CRU data
    if period is None: period = (1979,1994)
    dataset = loadDataset(name='Unity', varlist=variables, mode='climatology', period=period, shape=shapetype)
    if slices is not None: dataset = dataset(**noyears) # slice immediately
    obsens += dataset.load() 
  if lCRU: # this is basically regular operations with CRU as default
    obsens += loadEnsembleTS(names='CRU', season=seasons, aggregation=aggregation, slices=slices, 
                             varlist=variables, shape=shapetype, ldataset=True, **kwargs)    
  if lWSC: # another special case: river hydrographs
#     from datasets.WSC import loadGageStation, GageStationError
    try:
      dataset = loadGageStation(basin=basins, varlist=['runoff'], aggregation=aggregation, mode='climatology', 
                                filetype='monthly', basin_list=basin_list)
      if slices is not None: dataset = dataset(**noyears) # slice immediately
      obsens += dataset.load()
    except GageStationError: 
      pass # just ignore, if gage station data is missing 
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return obsens


# define new load fct. for experiments (not intended for observations)
@BatchLoad
def loadShapeEnsemble(seasons=None, basins=None, provs=None, shapes=None, varlist=None, aggregation='mean', 
                      slices=None, shapetype=shapetype, filetypes=None, period=None, variable_list=None, 
                      WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, **kwargs):
  ''' convenience function to load shape ensembles (in Ensemble container); kwargs are passed to loadEnsembleTS '''
  # prepare arguments
  variables, filetypes =  _resolveVarlist(varlist=varlist, filetypes=filetypes, 
                                          params=shp_params, variable_list=variable_list)
  # configure slicing (extract basin/province/shape and period)
  slices = _configSlices(slices=slices, basins=basins, provs=provs, shapes=shapes, period=period)
  # load ensemble (no iteration here)
  shpens = loadEnsembleTS(season=seasons, slices=slices, varlist=variables, shape=shapetype, 
                          aggregation=aggregation, filetypes=filetypes, WRF_exps=WRF_exps, CESM_exps=CESM_exps, 
                          WRF_ens=WRF_ens, CESM_ens=CESM_ens, **kwargs)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return shpens


# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
@BatchLoad
def loadStationEnsemble(seasons=None, provs=None, clusters=None, varlist=None, aggregation='mean', constraints=None, 
                        filetypes=None, cluster_name=None, stationtype=shapetype, load_list=None, lproduct='outer',
                        WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, 
                        variable_list=None, default_constraints=None, **kwargs):
  ''' convenience function to load station data for ensembles (in Ensemble container); kwargs are passed to loadEnsembleTS '''
  if load_list: load_list = load_list[:] # use a copy, since the list may be modified
 
#XXX: move this into helper function with Batch-decorator to allow batch-loading of varlists

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
                                          params=params, variable_list=variable_list)
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
  stnens = loadEnsembleTS(season=seasons, prov=provs, station=stationtype, varlist=variables, 
                          aggregation=aggregation, constraints=constraints, filetypes=filetypes, 
                          WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, 
                          load_list=load_list, lproduct=lproduct, lcheckVar=False, **kwargs)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return stnens


## abuse main section for testing
if __name__ == '__main__':
  
#   from projects.WesternCanada.WRF_experiments import WRF_exps, ensembles
#   from projects.WesternCanada.analysis_settings import exps_rc, variables_rc, loadShapeObservations  
  from projects.CESM_experiments import CESM_exps, CESM_ens
  from projects.GreatLakes.WRF_experiments import WRF_exps, WRF_ens
  from projects.GreatLakes.analysis_settings import exps_rc, variables_rc, loadShapeObservations
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

#   test = 'obs_timeseries'
  test = 'basin_timeseries'
#   test = 'station_timeseries'
#   test = 'province_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'obs_timeseries':
    
    # some settings for tests
    basins = ['GRW'] #; period = (1979,1994)
    varlist = ['sfcflx',]; aggregation = 'mean'

    shpens = loadShapeObservations(obs=None, basins=basins, varlist=varlist,
                                   aggregation=aggregation, load_list=['basins','varlist'],)
#                                    variable_list=variables_rc)
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins) # len(seasons)
    print shpens[0][0]
    for i,basin in enumerate(basins): 
      #for ds in shpens[i]: print ds.atts.shape_name
      assert all(ds.atts.shape_name == basin for ds in shpens[i])
  
  
  # test load function for basin ensemble time-series
  elif test == 'basin_timeseries':
    
    # some settings for tests
    exp = 'g-prj'; exps = exps_rc[exp].exps; #exps = ['Unity']
    basins = ['GRW']; seasons = ['summer','winter']
    varlist = ['sfcflx']; aggregation = 'mean'; red = dict(s='mean')

    shpens = loadShapeEnsemble(names=exps, basins=basins, seasons=seasons, varlist=varlist, 
                               aggregation=aggregation, filetypes=None, reduction=red, 
                               load_list=['basins','seasons',], lproduct='outer', domain=2,
                               WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens,
                               variable_list=variables_rc)
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(seasons)
    print shpens[0][0]
    for i,basin in enumerate(basins):
      i0 = i*len(seasons); ie = len(seasons)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])
      

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
    varlist = ['MaxPrecip_1d', 'MaxPrecip_5d','MaxPreccu_1d'][:1]; filetypes = ['hydro']
    varlist = ['precip', 'pet']; filetypes = ['hydro','aux']
    stnens = loadStationEnsemble(names=exps, provs=provs, clusters=clusters, varlist=varlist,  
                                 seasons=seasons, master=None, stationtype='ecprecip',
                                 domain=2, lensembleAxis=lensembleAxis, filetypes=filetypes,                                 
                                 variable_list=variables_rc, default_constraints=constraints_rc,
                                 WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens,
                                 load_list=['season','provs'], lproduct='outer',)
    # print diagnostics
    print stnens[0][0]; print ''
    assert len(stnens) == len(seasons)
    print stnens[0][0]
    print stnens[0][1].MaxPrecip_1d.mean()

      
  # test load function for province ensemble climatology
  if test == 'province_climatology':
    
    # some settings for tests
    exp = 'g-ens'; exps = exps_rc[exp].exps
    basins = ['GLB','GRW'] 
    varlists = ['precip','T2']; aggregation = 'mean'

    shpens = loadShapeEnsemble(names=exps, basins=basins, varlist=varlists, aggregation=aggregation,
                               period=(1979,1994), # this does not work properly with just a number...
                               load_list=['basins','varlist'], lproduct='outer', filetypes=['srfc'],
                               WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens,
                               variable_list=variables_rc)
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(varlists)
    assert shpens[0][1].time.coord[0] == 1
    for i,basin in enumerate(basins):
      i0 = i*len(varlists); ie = len(varlists)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])    