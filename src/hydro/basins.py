'''
Created on Feb 2, 2015

Utility functions related to loading basin-averaged data to support hydrological analysis.

@author: Andre R. Erler, GPL v3
'''

# external imports
from collections import namedtuple
# internal imports
from geodata.base import Ensemble, Dataset
from datasets.common import loadEnsembleTS, BatchLoad, loadDataset, shp_params
from geodata.misc import ArgumentError, EmptyDatasetError
from datasets.WSC import GageStationError, loadGageStation

# wet-day thresholds
wetday_thresholds = [0.2,1,10,20][:3]
wetday_extensions = ['_{:03.0f}'.format(threshold*10) for threshold in wetday_thresholds]
# variable collections
variables_rc = dict(); VL = namedtuple('VarList', ('vars','files','label'))     
variables_rc['temp']            = VL(vars=('T2', 'Tmax', 'Tmin'), files=('srfc','xtrm',), label='2m Temperature')
# variables_rc['temp']          = VL(vars=('T2',), files=('srfc',), label='2m Temperature')
variables_rc['precip_obs']      = VL(vars=('precip', 'solprec', 'wetfrq_010'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm_obs'] = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'wetprec_010'), files=('hydro',), label='Precipitation')
variables_rc['precip']          = VL(vars=('precip', 'preccu', 'solprec'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm']     = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'MaxPreccu_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_cesm']     = VL(vars=('MaxPrecip_1d', 'MaxPreccu_1d', ), files=('hydro',), label='Precipitation')
variables_rc['precip_alt']      = VL(vars=('MaxPrecip_1d', 'MaxPrecnc_1d', 'MaxSolprec_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_types']    = VL(vars=('precip','preccu','precnc'), files=('hydro',), label='Water Flux')
variables_rc['flux']            = VL(vars=('precip','snwmlt','p-et'), files=('hydro',), label='Water Flux')
variables_rc['flux_alt']        = VL(vars=('precip','snwmlt','solprec'), files=('hydro',), label='Water Flux')
variables_rc['flux_days']       = VL(vars=('wetfrq_010','snwmlt','p-et'), files=('hydro',), label='Water Flux')
variables_rc['wetprec']         = VL(vars=['wetprec'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Precip.')
variables_rc['wetdays']         = VL(vars=['wetfrq'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Ratio')
variables_rc['CWD']             = VL(vars=['CWD'+ext for ext in wetday_extensions]+['CNWD'], files=('hydro',), label='Continuous Wet-days')
variables_rc['CDD']             = VL(vars=['CDD'+ext for ext in wetday_extensions[:-1]]+['CNDD'], files=('hydro',), label='Continuous Dry-days')
variables_rc['sfcflx']          = VL(vars=('p-et','snwmlt','waterflx',), files=('hydro',), label='Surface Flux')
variables_rc['roffflx']         = VL(vars=('runoff','snwmlt','p-et'), files=('lsm','hydro'), label='Water Flux')
variables_rc['runoff']          = VL(vars=('waterflx','sfroff','runoff'), files=('lsm','hydro'), label='Runoff')
variables_rc['heat']            = VL(vars=('hfx','lhfx','rSM'),files=('srfc','lsm'), label='Energy Flux')
variables_rc['evap']            = VL(vars=('p-et','evap','pet',), files=('hydro',), label='Water Flux')
variables_rc['spei']            = VL(vars=('precip','evap','pet',), files=('hydro',), label='Water Flux')
variables_rc['Q2']              = VL(vars=('Q2',),files=('srfc',), label='2m Humidity')
variables_rc['aSM']             = VL(vars=('aSM',),files=('lsm',), label='Soil Moisture') 
variables_rc['rSM']             = VL(vars=('rSM',),files=('lsm',), label='Relative Soil Moisture')
# N.B.: Noah-MP does not have relative soil moisture and there is not much difference anyway

# corresponding filtype collections          
                                
# dataset collections
exps_rc = dict(); EX = namedtuple('Experiments', ('name','exps','styles','master','title'))
exps_rc['obs']     = EX(name='obs',exps=['CRU','WSC'], styles=['-','-.'], master='CRU', title='Observations')
exps_rc['val']     = EX(name='val',exps=['max-ens','erai-max','max-ens_d01'], master='max-ens',
                        styles=['-','--','-.'], title='WRF (IC Ens. Avg., ERA-I, 30km)')
exps_rc['prj']     = EX(name='prj',exps=['max-ens','max-ens-2050','max-ens-2100'], master='max-ens',
                        styles=['-','-.','--'], title='IC Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['max']     = EX(name='max',exps=['max-ens','max-ens-2050','max-ens-2100'], master='max-ens',
                        styles=['-','-.','--'], title='IC Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['ctrl']    = EX(name='ctrl',exps=['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], master='ctrl-ens',
                        styles=['-','-.','--'], title='Alt. IC Ens. Average (Hist., Mid-, End-Century)')
exps_rc['g-ens']   = EX(name='g-ens',exps=['g-ens','g-ens-2050','g-ens-2100'], master='g-ens',
                        styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['t-ens']   = EX(name='t-ens',exps=['t-ens','t-ens-2050','t-ens-2100'], master='t-ens',
                        styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['mex']     = EX(name='mex',exps=['mex-ens','mex-ens-2050','mex-ens-2100'], master='mex-ens',
                        styles=['-','-.','--'], title='WRF Ext. Ens. (Hist., Mid-, End-Century)')
exps_rc['phys']    = EX(name='phys',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], master='phys-ens',
                        styles=['-','-.','--'], title='WRF Phys. Ens. (Hist., Mid-, End-Century)')
exps_rc['nmp']     = EX(name='nmp',exps=['max-ens','max-ens-2050','max-nmp','max-nmp-2050'], master='max-nmp',
                        styles=['-.','.','--','-'], title='IC Ensemble Average & Noah-MP (Hist., Mid-Century)')
exps_rc['test']    = EX(name='test',exps=['max-ens','max-ctrl','max-hilev'], master='max-ens',
                        styles=['-','--','-.'], title='WRF (Ens. Avg., Ctrl-1, Hi-lev)')
exps_rc['1deg']    = EX(name='1deg',exps=['max-ens','max-ens_d01','max-1deg'], master='max-ens',
                        styles=['--','-.','-'], title='WRF (Ens. Avg., 30km, 1deg)')
exps_rc['seaice']  = EX(name='seaice',exps=['max-ens','max-ens-2050','max-seaice-2050','max-ens-2100','max-seaice-2100'],
                        styles=['-.','--','-','.--','.-'], master='max-ens', 
                        title='WRF (Hist., Mid-, End-Century; Ens. Avg., Sea-ice)')
exps_rc['max-Z']   = EX(name='max-Z',exps=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], master='max-ctrl',
                        styles=['-','-.','--'], title='WRF-Z (Hist., Mid-, End-Century)')
exps_rc['max-A']   = EX(name='max-A',exps=['max-ens-A','max-ens-A-2050','max-ens-A-2100'], master='max-ens-A',
                        styles=['-','-.','--'], title='WRF-A (Hist., Mid-, End-Century)')
exps_rc['max-B']   = EX(name='max-B',exps=['max-ens-B','max-ens-B-2050','max-ens-B-2100'], master='max-ens-B',
                        styles=['-','-.','--'], title='WRF-B (Hist., Mid-, End-Century)')
exps_rc['max-C']   = EX(name='max-C',exps=['max-ens-C','max-ens-C-2050','max-ens-C-2100'], master='max-ens-C',
                        styles=['-','-.','--'], title='WRF-C (Hist., Mid-, End-Century)')
exps_rc['si25']  = EX(name='si25',exps=['max-ens','max-ens-2050','max-seaice-2050'], master='max-ens-2050',
                        styles=[':','--','-'], title='WRF (Hist., Mid-Century; Ens. Avg., Sea-ice)')
exps_rc['si21']  = EX(name='si21',exps=['max-ens','max-ens-2100','max-seaice-2100'], master='max-ens-2100',
                        styles=[':','--','-'], title='WRF (Hist., End-Century; Ens. Avg., Sea-ice)')

# dataset variables
CRU_vars = ('T2','Tmin','Tmax','dTd','Q2','pet','precip','cldfrc','wetfrq','frzfrq')
WSC_vars = ('runoff','sfroff','ugroff')

# define new load fct. for observations
@BatchLoad
def loadShapeObservations(obs=None, seasons=None, basins=None, provs=None, varlist=None, slices=None,
                          aggregation='mean', shapefile=None, period=None, **kwargs):
  ''' convenience function to load basin & province observations; if no 'obs' are specified,
      sensible defaults are selected, based on 'varlist' '''
  # prepare arguments
  if shapefile is None: shapefile = 'shpavg' # really only one in use  
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
    if name in variables_rc: variables.update(variables_rc[name].vars)
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
def loadShapeEnsemble(names=None, seasons=None, basins=None, provs=None, varlist=None, aggregation='mean', slices=None,
                      shapefile=None, filetypes=None, period=None, **kwargs):
  ''' convenience function to load basin & province ensembles '''
  # prepare arguments
  if shapefile is None: shapefile = 'shpavg' # really only one in use  
  # resolve variable list (no need to maintain order)
  if isinstance(varlist,basestring): varlist = [varlist]
  variables = set(shp_params)
  for name in varlist: 
    if name in variables_rc: variables.update(variables_rc[name].vars)
    else: variables.add(name) 
  variables = list(variables)
  # resolve filetypes (and maintain order)
  if filetypes is None and name in variables_rc:
    filetypes = []
    for name in varlist: 
      for ft in variables_rc[name].files: 
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
                          varlist=variables, shape=shapefile, filetypes=filetypes, **kwargs)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  return shpens

## abuse main section for testing
if __name__ == '__main__':
  
  test = 'obs_timeseries'
#   test = 'basin_timeseries'
#   test = 'province_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'obs_timeseries':
    
    # some settings for tests
    basins = ['FRB','ARB'] #; period = (1979,1994)
    varlist = ['precip']; aggregation = 'std'

    shpens = loadShapeObservations(obs=None, seasons=None, basins=basins, varlist=varlist, period=None, 
                                   aggregation=aggregation, load_list=['basins',], lproduct='outer')
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
    exp = 'val'; exps = exps_rc[exp].exps; exps = ['Unity']
    basins = ['FRB','ARB']; seasons = ['summer','winter']
    varlist = ['wetfrq_010']; aggregation = 'mean'; red = dict(s='mean')

    shpens = loadShapeEnsemble(names=exps, basin=basins, season=seasons, varlist=varlist, 
                               aggregation=aggregation, filetypes=['lsm'], reduction=red, 
                               load_list=['basin','season',], lproduct='outer')
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(seasons)
    print shpens[0][0]
    for i,basin in enumerate(basins):
      i0 = i*len(seasons); ie = len(seasons)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])
      
      
  # test load function for province ensemble climatology
  if test == 'province_climatology':
    
    # some settings for tests
    exp = 'prj'; exps = exps_rc[exp].exps
    basins = ['FRB','ARB'] 
    varlists = ['precip','runoff']; aggregation = 'std'

    shpens = loadShapeEnsemble(names=exps, basin=basins, varlist=varlists, aggregation=aggregation, 
                               load_list=['basin','varlist'], lproduct='outer')
    # print diagnostics
    print shpens[0]; print ''
    assert len(shpens) == len(basins)*len(varlists)
    assert shpens[0][1].time.coord[0] == 1
    for i,basin in enumerate(basins):
      i0 = i*len(varlists); ie = len(varlists)*(i+1)
      assert all(all(ds.atts.shape_name == basin for ds in ens) for ens in shpens[i0:ie])    