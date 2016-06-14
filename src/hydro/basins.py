'''
Created on Feb 2, 2015

Utility functions related to loading basin-averaged data to support hydrological analysis.

@author: Andre R. Erler, GPL v3
'''

# internal imports
from geodata.base import Ensemble, Dataset
from utils.misc import defaultNamedtuple
from datasets.common import loadEnsembleTS, BatchLoad, loadDataset, shp_params
from geodata.misc import ArgumentError, EmptyDatasetError
from datasets.WSC import GageStationError, loadGageStation

# some definitions
VL = defaultNamedtuple('VarList', ('vars','files','label'))   
EX = defaultNamedtuple('Experiments', ('name','exps','styles','master','title','reference','target'))  

# wet-day thresholds
wetday_thresholds = [0.2,1,10,20]
wetday_extensions = ['_{:03.0f}'.format(threshold*10) for threshold in wetday_thresholds]

# dataset variables
CRU_vars = ('T2','Tmin','Tmax','dTd','Q2','pet','precip','cldfrc','wetfrq','frzfrq')
WSC_vars = ('runoff','sfroff','ugroff')

# define new load fct. for observations
@BatchLoad
def loadShapeObservations(obs=None, seasons=None, basins=None, provs=None, varlist=None, slices=None,
                          aggregation='mean', shapefile=None, period=None, variable_atts=None, **kwargs):
  ''' convenience function to load basin & province observations; if no 'obs' are specified,
      sensible defaults are selected, based on 'varlist' '''
  # prepare arguments
  if shapefile is None: shapefile = 'shpavg' # really only one in use  
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
    if name in variable_atts: variables.update(variable_atts[name].vars)
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
  # configure slicing (extract basin)
  if basins is not None and provs is not None: raise ArgumentError
  elif provs is not None or basins is not None:
    if slices is None: slices = dict()
    if provs is not None: slices['shape_name'] = provs
    if basins is not None: slices['shape_name'] = basins  
  if period is not None:
    if slices is None: slices = dict()
    slices['years'] = period
  if slices is not None:
    noyears = slices.copy(); noyears.pop('years',None) # slices for climatologies
  # prepare and load ensemble of observations
  obsens = Ensemble(name='obs',title='Observations', basetype=Dataset)
  if len(obs) > 0: # regular operations with user-defined dataset
    try:
      ensemble = loadEnsembleTS(names=obs, season=seasons, aggregation=aggregation, slices=slices, 
                                varlist=variables, shape=shapefile, ldataset=False, **kwargs)
      for ens in ensemble: obsens += ens
    except EmptyDatasetError: pass 
  if lUnity: # load Unity data instead of averaging CRU data
    if period is None: period = (1979,1994)
    dataset = loadDataset(name='Unity', varlist=variables, mode='climatology', period=period, shape=shapefile)
    if slices is not None: dataset = dataset(**noyears) # slice immediately
    obsens += dataset.load() 
  if lCRU: # this is basically regular operations with CRU as default
    obsens += loadEnsembleTS(names='CRU', season=seasons, aggregation=aggregation, slices=slices, 
                             varlist=variables, shape=shapefile, ldataset=True, **kwargs)    
  if lWSC: # another special case: river hydrographs
#     from datasets.WSC import loadGageStation, GageStationError
    try:
      dataset = loadGageStation(basin=basins, varlist=['runoff'], aggregation=aggregation, mode='climatology', filetype='monthly')
      if slices is not None: dataset = dataset(**noyears) # slice immediately
      obsens += dataset.load()
    except GageStationError: 
      pass # just ignore, if gage station if data is missing 
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return obsens


# define new load fct. for experiments (not intended for observations)
@BatchLoad
def loadShapeEnsemble(names=None, seasons=None, basins=None, provs=None, varlist=None, aggregation='mean', 
                         slices=None, shapefile=None, filetypes=None, period=None, variable_atts=None, 
                         WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, **kwargs):
  ''' convenience function to load basin & province ensembles '''
  # prepare arguments
  if shapefile is None: shapefile = 'shpavg' # really only one in use  
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
    if name in variable_atts: variables.update(variable_atts[name].vars)
    else: variables.add(name) 
  variables = list(variables)
  # resolve filetypes (and maintain order)
  if filetypes is None and name in variable_atts:
    filetypes = []
    for name in varlist: 
      for ft in variable_atts[name].files: 
        if ft not in filetypes: filetypes.append(ft)   
  # configure slicing (extract basin/province)
  if basins is not None and provs is not None: raise ArgumentError
  elif provs is not None or basins is not None:
    if slices is None: slices = dict()
    if provs is not None: slices['shape_name'] = provs
    if basins is not None: slices['shape_name'] = basins  
  if period is not None:
    if slices is None: slices = dict()
    slices['years'] = period
  # load ensemble (no iteration here)
  shpens = loadEnsembleTS(names=names, season=seasons, aggregation=aggregation, slices=slices,
                          varlist=variables, shape=shapefile, filetypes=filetypes, 
                          WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, 
                          **kwargs)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return shpens


# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
def loadStationEnsemble(names=None, provs=None, seasons=None, clusters=None, varlist=None, 
                        station=None, constraints=None, master=None, lall=True, 
                        aggregation='max', lminmax=False, slices=None, obsslices=None, 
                        lensembleAxis=False, years=None, reduction=None,
                        filetypes=None, domain=None, name_tags=None, cluster_name=None,
                        ensemble_list=None, ensemble_product='inner', load_list=None, lproduct='outer',
                        WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, 
                        variable_atts=None, default_constraints=None):
  ''' convenience function to load station data for ensembles (in Ensemble container) '''
  load_list = [] if load_list is None else load_list[:] # use a copy, since the list may be modified
  
  # figure out varlist
  if isinstance(varlist,basestring):
    if station is None:
      if varlist.lower().find('precip') >= 0: 
        station = 'ecprecip'
      elif varlist.lower().find('temp') >= 0: 
        station = 'ectemp'
    if filetypes is None: filetypes = variable_atts[varlist].files[:]
    varlist = variable_atts[varlist].vars[:]
  if cluster_name: varlist = varlist[:] + [cluster_name] # need to load this variable!
    
  # prepare arguments
  if provs or clusters:
    if constraints is None: constraints = default_constraints.copy()
    constraint_list = []
    if 'provs' in load_list and 'clusters' in load_list: 
      raise ArgumentError, "Cannot expand 'provs' and 'clusters' at the same time."
    # figure out proper handling of provinces
    if provs:
      if 'prov' not in load_list: 
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
      if 'cluster' not in load_list: 
        constraints['cluster'] = clusters; clusters = None
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
  stnens = loadEnsembleTS(names=names, season=seasons, prov=provs, aggregation=aggregation, lminmax=lminmax, 
                          years=years, station=station, slices=slices, obsslices=obsslices, lall=lall,
                          constraints=constraints, varlist=varlist, filetypes=filetypes, domain=domain,
                          ensemble_list=ensemble_list, ensemble_product=ensemble_product, master=master, 
                          load_list=load_list, lproduct=lproduct, name_tags=name_tags, lcheckVar=False,
                          WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens,
                          lensembleAxis=lensembleAxis, reduction=reduction)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return stnens


## abuse main section for testing
if __name__ == '__main__':
  
  #from projects.WesternCanada.WRF_experiments import Exp, WRF_exps, ensembles
  #from projects.WesternCanada.settings import exps_rc
  from projects.GreatLakes.WRF_experiments import WRF_exps, ensembles
  from projects.GreatLakes.settings import exps_rc, variables_rc, loadShapeObservations
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

#  test = 'obs_timeseries'
#   test = 'basin_timeseries'
  test = 'station_timeseries'
#   test = 'province_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'obs_timeseries':
    
    # some settings for tests
    basins = ['GLB'] #; period = (1979,1994)
    varlist = ['precip',]; aggregation = 'mean'

    shpens = loadShapeObservations(obs='GPCC', basins=basins, varlist=varlist,
                                   aggregation=aggregation, load_list=['basins','varlist'],)
#                                    variable_atts=variables_rc)
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
    exp = 'g-ens'; exps = exps_rc[exp].exps; #exps = ['Unity']
    basins = ['GLB']; seasons = ['summer','winter']
    varlist = ['precip']; aggregation = 'mean'; red = dict(s='mean')

    shpens = loadShapeEnsemble(names=exps, basins=basins, seasons=seasons, varlist=varlist, 
                               aggregation=aggregation, filetypes=None, reduction=red, 
                               load_list=['basins','seasons',], lproduct='outer', domain=2,
                               WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                               variable_atts=variables_rc)
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
    constraints_rc['lat'] = (45,55) 
    constraints_rc['max_zerr'] = 100 # reduce sample size
    constraints_rc['prov'] = ('ON')
    constraints_rc['end_after'] = 1980
  
    # some settings for tests
#     exp = 'val'; exps = ['EC', 'erai-max', 'max-ctrl'] # 
#     exp = 'max-all'; exps = exps_rc[exp]; provs = ('BC','AB'); clusters = None
#     exp = 'marc-prj'; exps = exps_rc[exp]; provs = 'ON'; clusters = None
#     provs = ('BC','AB'); clusters = None
    provs = None; clusters = None; lensembleAxis = False; sample_axis = None; lflatten = False
    exp = 'g-ens'; exps = exps_rc[exp]; provs = None; clusters = [1,3]
    seasons = ['summer']; lfit = True; lrescale = True; lbootstrap = False
    lflatten = False; lensembleAxis = True; sample_axis = ('station','year')
    varlist = ['MaxPrecip_1d', 'MaxPrecip_5d','MaxPreccu_1d'][:1]; filetypes = ['hydro']
    stnens = loadStationEnsemble(names=exps.exps, provs=provs, clusters=clusters, varlist=varlist,  
                                 seasons=seasons, master=None, station='ecprecip',
                                 domain=2, lensembleAxis=lensembleAxis, filetypes=filetypes,                                 
                                 variable_atts=variables_rc, default_constraints=constraints_rc,
                                 WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                                 load_list=['season',], lproduct='outer',)
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
                               WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                               variable_atts=variables_rc)
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(varlists)
    assert shpens[0][1].time.coord[0] == 1
    for i,basin in enumerate(basins):
      i0 = i*len(varlists); ie = len(varlists)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])    