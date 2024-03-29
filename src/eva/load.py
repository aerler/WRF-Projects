'''
Created on Jan 30, 2015

Utility functions related to loading station  data to support extreme value analysis.

@author: Andre R. Erler, GPL v3
'''
# external imports
import numpy as np
# internal imports
from geodata.base import Dataset, Ensemble
from geodata.misc import ArgumentError, AxisError
from geodata.stats import VarRV
from datasets.common import expandArgumentList
# imports from clim
from clim.load import loadShapeEnsemble, loadStationEnsemble
from geodata.netcdf import DatasetNetCDF


# convenience function to extract a station (or skip, if flat)
def extractLocation(*datasets, **kwargs):
    lflatten = kwargs.pop('lflatten',None)
    if len(kwargs) == 0: lflatten = True
    elif len(kwargs) > 1: raise KeyError
    # extract data
    if lflatten:
        # just rename...
        newsets = [None]*2+list(datasets)
    else:
        # interpret arguments
        var,val = list(kwargs.items())[0] 
        # select coordinate
        ds0 = getattr(datasets[0],var)
        if isinstance(ds0,Ensemble): ds0 = ds0[0]
        coord = ds0.findValue(val)
        idx = ds0.findValue(val,lidx=True)
        cargs = {ds0.axes[0].name:coord} # may not have the same index everywhere
        newsets = [coord, idx] # first return arguments
        for dataset in datasets:
            newsets.append(dataset(**cargs))
    return newsets


# helper method for generateStatistics
def _rescaleSample(smpl, loc, bs_axis=None):
  ''' helper method to rescale samples'''
  #print smpl.shape, np.nanmean(smpl), loc
  if isinstance(loc, np.ndarray):
    if np.any(loc != 1): 
      if bs_axis is not None: loc = loc.take([0], axis=bs_axis).squeeze() # untie bootstraps
      if loc.ndim == smpl.ndim: 
        smpl = smpl*loc
      elif loc.ndim < smpl.ndim:
        smpl = smpl*loc.reshape(loc.shape+(1,)*(smpl.ndim-loc.ndim)) # broadcast
      else: raise ValueError(loc.shape)
  elif loc is not None and loc != 1: smpl = smpl*loc
  return smpl

# function to generate a rescaled dataset or ensemble of datasets
def rescaleDistributions(fits, samples=None, reference=None, target=None, lscale=False, lsourceScale=False, 
                         suffixes=None, lglobal=False):
  ''' Rescale fits (datasets), so that the mean of each variable matches the corresponding variable in the
      reference dataset; if a target is specified, the target scale factors are applied to all
      datasets, if target is None, each dataset is rescaled individually. '''
  if not isinstance(fits, (list,tuple,Ensemble)): raise TypeError(fits)
  if lsourceScale: 
    if samples is None: raise ArgumentError('Need source samples for source rescaling!')
    elif not isinstance(samples, (list,tuple,Ensemble)): raise TypeError(samples)
  elif samples is None: samples = [None]*len(fits) # avoid errors
  if isinstance(fits,Ensemble) and isinstance(reference,str):
    reference = fits[reference]
  elif not isinstance(reference,Dataset): raise TypeError
  if target is None or target == 'auto': pass # every dataset is scaled individually or based on suffixes
  elif isinstance(fits,Ensemble) and isinstance(target,str):
    target = fits[target]
  elif not isinstance(target,Dataset): raise TypeError(target)
  if suffixes is None: suffixes = ('-2050','2100') # suffixes for scaling heuristic
  # determine scale factor
  def scaleFactor(reference, target, lscale=False, lglobal=False):
    ''' internal function to compute rescaling factors for common variables '''
    scalefactors = dict() # return dict with scalefactors for all applicable variables 
    for varname,refvar in reference.variables.items():
      if varname in target and isinstance(refvar,VarRV): # only varaibles that appear in both sets
        tgtvar = target.variables[varname]
        iloc = 1 if refvar.shape[-1] == 3 else 0
        # insert dummy ensemble axis, if necessary
        refvar = refvar.insertAxes(new_axes=tgtvar.axes, lcopy=True, asVar=True, linplace=False)  
        if refvar.axes[-1].name.startswith('params'): refdata = refvar.data_array.take(iloc, axis=-1)
        else: raise AxisError(refvar.axes[-1])
        if refvar.ndim < tgtvar.ndim:
          # N.B.: this is necessary, because WRF (target) can have an extra ensemble dimension that obs
          #       typically don't have; then we just replicate the obs for each ensemble element
          from warnings import warn
          if lglobal: warn("Scalefactors are being averaged over extra target dimensions (e.g. 'ensemble' axis)")
          dimdiff = tgtvar.ndim-refvar.ndim
          if refvar.shape != tgtvar.shape[dimdiff:]: raise AxisError("{:s} != {:s}".format(tgtvar, refvar))
          refdata = refdata.reshape((1,)*dimdiff+refvar.shape[:-1])
        elif refvar.shape != tgtvar.shape: raise AxisError("{:s} != {:s}".format(tgtvar, refvar))
        tgtdata = tgtvar.data_array.take(iloc, axis=-1)
        if lglobal: loc = np.mean(refdata) / np.mean(tgtdata)
        else: loc = refdata / tgtdata
        if lscale:
          iscale = 2 if refvar.shape[-1] == 3 else 1  
          if lglobal:
            scale = np.mean(refvar.data_array.take(iscale, axis=-1)) / np.mean(tgtvar.data_array.take(iscale, axis=-1))
          else: scale = refvar.data_array.take(iscale, axis=-1) / tgtvar.data_array.take(iscale, axis=-1)
          scalefactors[varname] = loc, (scale/loc)
        else: scalefactors[varname] = loc
    return scalefactors # return dict with scale factors for variables    
  # compute general scalefactors
  if target == 'auto': 
    scalefactor_collection = dict()
  elif target is not None: 
    scalefactors = scaleFactor(reference, target, lscale=lscale, lglobal=lglobal) 
  # loop over datasets
  rescaled_fits = []
  for dataset,sample in zip(fits,samples):
    if dataset == reference:
      # determine variables that can be scaled (VarRV's)
      varlist = [varname for varname,var in dataset.variables.items() if isinstance(var,VarRV)]
      rescaled_dataset = dataset.copy(varlist=varlist)
      # add mock scale factors for consistency
      for var in rescaled_dataset.variables.values():
        var.atts['loc_factor'] = 1
        var.atts['scale_factor'] = 1
        var.atts['shape_factor'] = 1
    else:
      # generate new dataset (without variables, and in-memory)
      if isinstance(dataset, DatasetNetCDF):
        rescaled_dataset = dataset.copy(varlist=[], asNC=False)
      else: rescaled_dataset = dataset.copy(varlist=[]) 
      # individual scaling      
      if target is None or target == 'auto':
        parent = None
        if target == 'auto' and dataset.name.endswith(suffixes):
          for suffix in suffixes: 
            if dataset.name.endswith(suffix): # check, which suffix, and remove it
              parent = dataset.name[:-(len(suffix)+1)]; break        
          if parent and '-' not in parent: parent += '-1' # convention for WRF names
        if parent and parent in scalefactor_collection:
          scalefactors = scalefactor_collection[parent] # use scale factors from parent
        else: # scale individually
          scalefactors = scaleFactor(reference, dataset, lscale=lscale, lglobal=lglobal)
          if target == 'auto': scalefactor_collection[dataset.name] = scalefactors # for later use 
      # loop over variables
      for varname,scalefactor in scalefactors.items():
        if varname in dataset:
          # rescale and add variable to new dataset
          var = dataset.variables[varname]
          if lscale: rsvar = var.rescale(loc=scalefactor[0], scale=scalefactor[1])
          else: rsvar = var.rescale(loc=scalefactor)
          rescaled_dataset.addVariable(rsvar)
          # if desired, also rescale the source/sample
          if lsourceScale:
            assert sample is not None
            if lscale: raise NotImplementedError
            svar = sample.variables[varname]
            bsi = var.axisIndex(var.bootstrap_axis) if var.bootstrap_axis else None # tie bootstraps
            svar.load(_rescaleSample(svar.getArray(), scalefactor, bs_axis=bsi), lrecast=True)
#             svar *= scalefactor # in-place scaling
            svar.atts['rescaled'] = True
    # add dataset to list
    rescaled_fits.append(rescaled_dataset)
  # put everythign into Ensemble, if input was Ensemble
  if isinstance(fits,Ensemble): 
    rescaled_fits = Ensemble(*rescaled_fits, name=fits.ens_name, title=fits.ens_title)
  # return fits/ensemble
  return rescaled_fits

  
def addDistFit(ensemble=None, lfit=True, lflatten=None, lrescale=False, reference=None, target=None,  
               lbootstrap=False, nbs=30, sample_axis=None, lglobalScale=False, lsourceScale=False, 
               lcrossval=False, ncv=0.2, dist=None, dist_args=None, load_list=None, lproduct='outer', **kwargs): 
  ''' add distribution fits to ensemble; optionally also rescale; kwargs are necessary for correct list expansion '''
  
  # find appropriate sample axis
  if lflatten: 
    if sample_axis is not None: raise ArgumentError(sample_axis)
  elif sample_axis is None: # auto-detect
    for saxis in ('time','year'):
      if all([all(ens.hasAxis(saxis)) for ens in ensemble]): 
        sample_axis = saxis; break
    if sample_axis is None: raise AxisError("No sample axis detected") 
  else:
    if isinstance(sample_axis,str):
      if not all([all(ens.hasAxis(sample_axis)) for ens in ensemble]): raise AxisError(sample_axis) 
    elif isinstance(sample_axis,(tuple,list)):
      # check that axes are there
      for ax in sample_axis: 
        if not all([all(ens.hasAxis(ax)) for ens in ensemble]): raise AxisError(ax)
      # merge axes
      ensemble = [ens.mergeAxes(axes=sample_axis, new_axis='sample', asVar=True, linplace=False, \
                              lcheckAxis=False) for ens in ensemble]
      sample_axis = 'sample'
    else: raise AxisError(sample_axis)
  # perform fit or return dummy
  if dist_args is None: dist_args = dict()
  if lfit: fitens = [ens.fitDist(lflatten=lflatten, axis=sample_axis, lcrossval=lcrossval, ncv=ncv,
                                 lignoreParams=True, lbootstrap=lbootstrap, nbs=nbs,
                                 dist=dist, **dist_args) for ens in ensemble]
  else: fitens = [None]*len(ensemble)
  
  # rescale fitted distribution (according to certain rules)
  if lrescale:
    if not reference: raise ArgumentError(str(reference))
    # expand target list
    if isinstance(target, (list,tuple)):  
      expand_list = load_list[:]
      if 'names' in expand_list: expand_list[expand_list.index('names')] = 'target' 
      kwarg_list = expandArgumentList(target=target, expand_list=expand_list, lproduct=lproduct, **kwargs)
      targets = [kwarg['target'] for kwarg in kwarg_list]
    else: targets = [target]*len(fitens)
    if isinstance(reference, (list,tuple)): raise NotImplementedError # don't expand reference list
    # use global reference, if necessary
    if isinstance(reference,str) and not all(reference in fit for fit in fitens):
      i = 0 
      while i < len(fitens) and reference not in fitens[i]: i += 1
      if i >= len(fitens): raise ArgumentError("Reference {:s} not found in any dataset!".format(reference))
      reference = fitens[i][reference]
    sclens = [rescaleDistributions(fit, samples=ens, reference=reference, target=tgt, 
                                   lglobal=lglobalScale, lsourceScale=lsourceScale) for fit,ens,tgt in zip(fitens,ensemble,targets)]
  else: sclens = [None]*len(ensemble)
  
  # return results
  return fitens, sclens


## function to load ensembles of time-series data and compute associated distribution ensembles
# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
def loadEnsembleFit(lfit=True, dist=None, dist_args=None, reference=None, target=None,
                    lglobalScale=False, lrescale=False, lsourceScale=False, lflatten=True, 
                    sample_axis=None, lcrossval=False, ncv=0.2, lbootstrap=False, nbs=100,
                    WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, variable_list=None, 
                    load_list=None, lproduct='outer', lshort=True, datatype=None, 
                    obs_list=None, obs_ts=None, lforceList=False, **kwargs):
  ''' convenience function to load ensemble time-series data and compute associated distribution ensembles '''
  if lrescale and not lfit: raise ArgumentError
  load_list = [] if load_list is None else load_list[:] # use a copy, since the list may be modified

  # load shape ensemble
  if datatype.lower() == 'shape':
    ensemble = loadShapeEnsemble(variable_list=variable_list, WRF_exps=WRF_exps, CESM_exps=CESM_exps, 
                                 WRF_ens=WRF_ens, CESM_ens=CESM_ens, load_list=load_list, lforceList=lforceList,
                                 lproduct=lproduct, obs_list=obs_list, obs_ts=obs_ts, **kwargs)
  if datatype.lower() == 'station':
    ensemble = loadStationEnsemble(variable_list=variable_list, WRF_exps=WRF_exps, CESM_exps=CESM_exps, 
                                   WRF_ens=WRF_ens, CESM_ens=CESM_ens, load_list=load_list, lforceList=lforceList,
                                   lproduct=lproduct, obs_list=obs_list, obs_ts=obs_ts, **kwargs)
  # N.B.: kwargs are first passed on to loadShapeEnsemble/loadStationEnsemble and then to loadEnsembleTS

  # generate matching datasets with fitted distributions
  fitens, sclens = addDistFit(ensemble=ensemble, lfit=lfit, dist=dist, dist_args=dist_args, lflatten=lflatten, 
                              lrescale=lrescale, reference=reference, target=target, lglobalScale=lglobalScale,   
                              lsourceScale=lsourceScale, lbootstrap=lbootstrap, nbs=nbs, lcrossval=lcrossval, 
                              ncv=ncv, sample_axis=sample_axis, load_list=load_list, lproduct=lproduct, **kwargs)
  # N.B.: kwargs are mainly needed to infer the expanded shape of the ensemble list
  
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  if lshort:
    if lfit and lrescale: 
      if lsourceScale: return ensemble, sclens
      else: return ensemble, fitens, sclens
    elif lfit: return ensemble, fitens
    else: return ensemble
  else:
    return ensemble, fitens, sclens

# specialized ensemble-fit function for shape data
def loadShapeFit(**kwargs):
  ''' specialized ensemble-fit function for shape data; arguments are passed in this order to: 
      loadEnsembleFit, loadShapeEnsemble, loadEnsembleTS, loadDataset, and the dataset module '''
  if 'datatype' in kwargs: raise ArgumentError
  return loadEnsembleFit(datatype='shape', **kwargs)

# specialized ensemble-fit function for station data
def loadStationFit(**kwargs):
  ''' specialized ensemble-fit function for shape data; arguments are passed in this order to: 
      loadEnsembleFit, loadShapeEnsemble, loadEnsembleTS, loadDataset, and the dataset module '''
  if 'datatype' in kwargs: raise ArgumentError
  return loadEnsembleFit(datatype='station', **kwargs)


## abuse main section for testing
if __name__ == '__main__':
  
  from projects.WesternCanada.WRF_experiments import Exp, WRF_exps, ensembles
  from projects.WesternCanada import exps_rc, variables_rc
  from projects.WesternCanada import loadShapeFit, loadStationFit
  #from projects.GreatLakes.WRF_experiments import WRF_exps, ensembles
  #from projects.GreatLakes.settings import exps_rc
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

#   test = 'shape_ensemble'
  test = 'station_ensemble'
#   test = 'rescaling'

  # station selection criteria
  constraints_rc = dict()
#   constraints_rc['min_len'] = 15 # for valid climatology
#   constraints_rc['lat'] = (45,55) 
#   constraints_rc['max_zerr'] = 300 # reduce sample size
#   constraints_rc['prov'] = ('BC','AB')
#   constraints_rc['end_after'] = 1980
  
  # test load function for station ensemble
  if test == 'shape_ensemble':
    
    # some settings for tests
    exp = 'max-prj'; exps = exps_rc[exp]
    basins = ['FRB','ARB']; seasons = 'jas'
    load_list = ['basins']

    bsnens, fitens = loadShapeFit(names=exps.exps, basins=basins, seasons=seasons, varlist=['aSM'], 
                                       filetypes=['lsm'], lfit=True, aggregation='mean', dist='norm',
                                       load_list=load_list, lproduct='outer',
#                                        WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                                       variable_list=variables_rc,)
    # print diagnostics
    print(bsnens[0]); print('')
    print(fitens[0][0])
    #print fitens[0] if fitens is not None else fitens ; print ''
    assert len(bsnens) == len(basins)
    assert fitens[0][0].atts.shape_name == basins[0]
    
  # test load function for station ensemble
  elif test == 'station_ensemble':
    
    # some settings for tests
    exp = 'val'; exps = ['EC', 'max-ctrl'] # 
    exp = 'max-obs'; exps = exps_rc[exp]; slices = None; obsslice = dict(years=(1952,2012)) # 60 years 
#     exp = 'marc-prj'; exps = exps_rc[exp]; provs = 'ON'; clusters = None
    provs = ('AB',); clusters = None; cluster_name = None; load_list = ['seasons','provs']
#     provs = None; clusters = None; lensembleAxis = False; sample_axis = None; lflatten = False
#     exp = 'max-obs'; exps = exps_rc[exp]; slices = None; 
#     provs = None; clusters = [1,8,]; cluster_name = 'cluster_projection2'; load_list = ['seasons','clusters']
    seasons = ['summer','winter']; lfit = True; lrescale = True; lbootstrap = False
#     lflatten = False; lensembleAxis = True; sample_axis = ('station','year')
    lflatten = True; lensembleAxis = False; sample_axis = None
    varlist = ['MaxPrecip_1d', 'MaxPrecip_5d','MaxPreccu_1d'][:1]; filetypes = ['hydro']
#     slices = [dict(years=(1979,1995))]*3+[None]
    stnens, fitens, sclens  = loadStationFit(names=exps.exps, provs=provs, clusters=clusters, 
                                             cluster_name=cluster_name, lforceList=False,
                                             seasons=seasons, lfit=lfit, master=None, stationtype='ecprecip',
                                             lrescale=lrescale, reference=exps.reference, target=exps.target,
                                             lflatten=lflatten, domain=None, lbootstrap=lbootstrap, nbs=10,
                                             lensembleAxis=lensembleAxis, sample_axis=sample_axis,
                                             varlist=varlist, filetypes=filetypes, slices=slices,
                                             obsslices=obsslice, lminmax=True, aggregation='max',
                                             ensemble_list=['names','slices'] if slices else None,
                                             variable_list=variables_rc, default_constraints=constraints_rc,
                                             #WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                                             load_list=load_list, lproduct='outer', lcrossval=None,)
    # print diagnostics
    print(stnens[0][0]); print('')
#     print(fitens[0][0].MaxPrecip_1d.atts.sample_axis); print('')
#     print(fitens[0][1] if fitens is not None else fitens) ; print('')
#     print(sclens[0] if sclens is not None else sclens) ; print('')
#     assert len(stnens) == len(seasons) * (len(clusters) if clusters else len(provs))
    print(stnens[0][0].MaxPrecip_1d.max())
  
  elif test == 'rescaling':
    
    # some settings for tests
    exps = ['max-val','max-all']; exps = [exps_rc[exp].exps for exp in exps]
    provs = ('BC',); seasons = 'summer'
    lflatten = True; lfit = True

    # rescale    
    stnens,fitens,scalens = loadStationFit(names=exps, provs=provs, seasons=seasons, varlist='precip', 
                                           lrescale=True, reference='EC', target='max-ens',
                                           filetypes='hydro', domain=2, lfit=lfit, lflatten=lflatten,
                                           lbootstrap=False, nbs=10, name='test_ensemble',
                                           variable_list=variables_rc, default_constraints=constraints_rc,
                                           #WRF_exps=WRF_exps, CESM_exps=None, WRF_ens=ensembles, CESM_ens=None,
                                           load_list=['names'], lproduct='outer')
    # print diagnostics
    print(stnens[0]); print('')
    print(fitens[0]); print('')
    
    # special rescaling for projection plots
    scalens = [rescaleDistributions(ens, reference=fitens[0]['EC'], target='max-ens') for ens in fitens]

    # print diagnostics
    print(scalens[0]); print('')
    print(scalens[1][1])
