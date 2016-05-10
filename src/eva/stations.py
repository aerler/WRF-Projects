'''
Created on Jan 30, 2015

Utility functions related to loading station data to support extreme value analysis.

@author: Andre R. Erler, GPL v3
'''
# external imports
from collections import namedtuple
import numpy as np
# internal imports
from geodata.base import Dataset, Ensemble
from geodata.misc import ArgumentError, AxisError
from geodata.stats import VarRV
from datasets.common import loadEnsembleTS, expandArgumentList
# imports from hydro
from hydro.basins import loadShapeEnsemble
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
        var,val = kwargs.items()[0] 
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


# function to generate a rescaled dataset or ensemble of datasets
def rescaleDistributions(datasets, reference=None, target=None, lscale=False, suffixes=None, lglobal=False):
  ''' Rescale datasets, so that the mean of each variable matches the corresponding variable in the
      reference dataset; if a target is specified, the target scale factors are applied to all
      datasets, if target is None, each dataset is rescaled individually. '''
  if not isinstance(datasets, (list,tuple,Ensemble)): raise TypeError
  if isinstance(datasets,Ensemble) and isinstance(reference,basestring):
    reference = datasets[reference]
  elif not isinstance(reference,Dataset): raise TypeError
  if target is None or target == 'auto': pass # every dataset is scaled individually or based on suffixes
  elif isinstance(datasets,Ensemble) and isinstance(target,basestring):
    target = datasets[target]
  elif not isinstance(target,Dataset): raise TypeError, target
  if suffixes is None: suffixes = ('-2050','2100') # suffixes for scaling heuristic
  # determine scale factor
  def scaleFactor(reference, target, lscale=False, lglobal=False):
    ''' internal function to compute rescaling factors for common variables '''
    scalefactors = dict() # return dict with scalefactors for all applicable variables 
    for varname,refvar in reference.variables.iteritems():
      if varname in target and isinstance(refvar,VarRV): # only varaibles that appear in both sets
        tgtvar = target.variables[varname]
        iloc = 1 if refvar.shape[-1] == 3 else 0
        # insert dummy ensemble axis, if necessary
        refvar = refvar.insertAxes(new_axes=tgtvar.axes, lcopy=True, asVar=True, linplace=False)  
        if refvar.axes[-1].name.startswith('params'): refdata = refvar.data_array.take(iloc, axis=-1)
        else: raise AxisError, refvar.axes[-1]
        if refvar.ndim < tgtvar.ndim:
          # N.B.: this is necessary, because WRF (target) can have an extra ensemble dimension that obs
          #       typically don't have; then we just replicate the obs for each ensemble element
          from warnings import warn
          if lglobal: warn("Scalefactors are being averaged over extra target dimensions (e.g. 'ensemble' axis)")
          dimdiff = tgtvar.ndim-refvar.ndim
          if refvar.shape != tgtvar.shape[dimdiff:]: raise AxisError, "{:s} != {:s}".format(tgtvar, refvar)
          refdata = refdata.reshape((1,)*dimdiff+refvar.shape[:-1])
        elif refvar.shape != tgtvar.shape: raise AxisError, "{:s} != {:s}".format(tgtvar, refvar)
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
  rescaled_datasets = []
  for dataset in datasets:
    if dataset == reference:
      # determine variables that can be scaled (VarRV's)
      varlist = [varname for varname,var in dataset.variables.iteritems() if isinstance(var,VarRV)]
      rescaled_dataset = dataset.copy(varlist=varlist)
      # add mock scale factors for consistency
      for var in rescaled_dataset.variables.itervalues():
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
      for varname,scalefactor in scalefactors.iteritems():
        if varname in dataset:
          # rescale and add variable to new dataset
          var = dataset.variables[varname]
          if lscale: rsvar = var.rescale(loc=scalefactor[0], scale=scalefactor[1])
          else: rsvar = var.rescale(loc=scalefactor)
          rescaled_dataset.addVariable(rsvar)
    # add dataset to list
    rescaled_datasets.append(rescaled_dataset)
  # put everythign into Ensemble, if input was Ensemble
  if isinstance(datasets,Ensemble): 
    rescaled_datasets = Ensemble(*rescaled_datasets, name=datasets.ens_name, title=datasets.ens_title)
  # return datasets/ensemble
  return rescaled_datasets

# wet-day thresholds
wetday_thresholds = [0.2,1,10,20]
wetday_extensions = ['_{:03.0f}'.format(threshold*10) for threshold in wetday_thresholds]
# variable collections
variables_rc = dict(); VL = namedtuple('VarList', ('vars','files','label'))
variables_rc['runoff'] = VL(vars=['precip', 'p-et', 'snwmlt', 'waterflx','sfroff','runoff','aSM'],
                            files=['hydro','lsm'], label='')
variables_rc['soil'] = VL(vars=['aSM'], files=['lsm'], label='Soil Moisture')
#                          ['MaxPrecip_7d', 'MaxWaterflx_7d', 'wetfrq', 'CDD', 'CWD']
variables_rc['precip_all'] = VL(vars=['precip', 'MaxPrecip_1h', 'MaxPrecip_6h', 'MaxPrecip_1d', 'MaxPrecip_5d',
                                  'MaxWaterflx_5d', 'wetfrq', 'CDD', 'CWD'], files=['hydro','srfc'], label='Precipitation')
variables_rc['temp']   = VL(vars=['T2', 'MaxT2', 'MaxTmax', 'MinTmin', 'MaxTmax_7d', 'MinTmin_7d', 
                                  'frzfrq', 'CFD', 'CDFD', 'CNFD'], files=['srfc','xtrm'], label='Temperature')                          
variables_rc['precip_short'] = VL(vars=['MaxPreccu_1h', 'MaxPrecnc_1h', 'MaxPrecip_6h', 'MaxPreccu_6h'], files=['xtrm','srfc'], label='MaxPrecip')
variables_rc['precip_long'] = VL(vars=['MaxPrecip_1d', 'MaxPreccu_1d', 'MaxPrecip_5d',], files=['hydro'], label='MaxPrecip')
variables_rc['precip_CDD'] = VL(vars=['CNDD','CNWD']+[var+threshold for threshold in wetday_extensions for var in 'CDD','CWD'], 
                                 files=['hydro'], label='CDD')

# station selection criteria
constraints_rc = dict()
constraints_rc['min_len'] = 15 # for valid climatology
constraints_rc['lat'] = (45,55) 
constraints_rc['max_zerr'] = 300 # can't use this, because we are loading EC data separately from WRF
constraints_rc['prov'] = ('BC','AB')
constraints_rc['end_after'] = 1980

# dataset collections
exps_rc = dict(); EX = namedtuple('Experiments', ('name','exps','colors','master','reference','target','title'))
exps_rc['val-res'] = EX(name='val-res', exps=['Observations', 'max-ens_d01', 'max-ens', 'Ens'], # ,'erai-max','ctrl-1','max-1deg' 
                        colors=None, master='max-ens', reference='Observations', target=None, title='Resolution')
exps_rc['val-all'] = EX(name='val-all', exps=['Observations', 'Ens', 'max-ens', 'ctrl-ens'], # ,'erai-max','ctrl-1','max-1deg' 
                        colors=None, master='max-ens', reference='Observations', target=None, title='Resolution')
exps_rc['max-all'] = EX(name='max-all', exps=['Observations', 'max-ens','max-ens_d01','erai-max','max-ens-2050','max-ens-2100'], # ,'erai-max','ctrl-1','max-1deg' 
                        colors=None, master='max-ens', reference='Observations', target='auto', title='Validation & Projection')
exps_rc['cesm-all'] = EX(name='cesm-all', exps=['Observations', 'Ens','Ens-2050','Ens-2100'],  
                        colors=None, master='Ens', reference='Observations', target='Ens', title='CESM Validation & Projection')
exps_rc['cesm-mal'] = EX(name='cesm-ensa', exps=['Observations', 'MEns','MEns-2050','MEns-2100'],  
                        colors=None, master='MEns', reference='Observations', target='MEns', title='CESM Validation & Projection')
exps_rc['cesm-ens'] = EX(name='cesm-ens', exps=['MEns','MEns-2050','MEns-2100'],  
                        colors=None, master='MEns', reference='MEns', target='MEns', title='CESM Projection')
exps_rc['max-val'] = EX(name='max-val', exps=['Observations','max-ens_d01','max-ens','erai-max'], # ,'erai-max','ctrl-1','max-1deg' 
                        colors=None, master='max-ens', reference='Observations', target=None, title='Validation')
exps_rc['max-prj'] = EX(name='prj', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                        colors=None, master='max-ens', reference=None, target=None, title='Projection')
exps_rc['sens-res']  = EX(name='sens-res', exps=['Observations', 'Ens', 'max-ctrl','max-ctrl_d01','max-1deg','max-lowres_d01'], 
                          colors=None, master='max-ctrl', reference='Observations', target=None, title='Sensitivity to Resolution')
exps_rc['sens-max']  = EX(name='sens-max',exps=['max-ctrl','max-nosub','max-kf','max-nmp','max-noflake'], 
                          colors=None, master='max-ctrl', reference=None, target=None, title='Sensitivity Tests (Max)')
exps_rc['sens-phys'] = EX(name='sens-phys',exps=['Observations','max-ctrl','old-ctrl','ctrl-1','new-ctrl','new-v361'], 
                          colors=None, master='max-ctrl', reference='Observations', target=None, title='Sensitivity (Phys. Ens.)')
exps_rc['sens-phys-2100'] = EX(name='sens-phys-2100',exps=['max-ctrl-2100','old-ctrl-2100','ctrl-2100','new-ctrl-2100','new-v361-2100'], 
                               colors=None, master='max-ctrl-2100', reference=None, target=None, title='Sensitivity (Phys. Ens. 2100)')
exps_rc['sens-ens']  = EX(name='sens-ens',exps=['Observations','max-ens', 'max-ctrl','max-ens-A','max-ens-B','max-ens-C','erai-max'], 
                          colors=None, master='max-ens', reference='Observations', target=None, title='Sensitivity (Phys. Ens.)')
# exps_rc['sens-cu']     = EX(name='sens-cu', exps=['grell-ens','kf-ens','grell-ens-2050','kf-ens-2050','grell-ens-2100','kf-ens-2100'], 
#                             colors=None, master=['grell-ens','grell-ens-2050','grell-ens-2100'], reference=None, target=None, title='Grell vs. KF')
exps_rc['sens-cu']     = EX(name='sens-cu', exps=['grell-ens','kf-ens','grell-ens-2100','kf-ens-2100'], 
                            colors=None, master=['grell-ens','grell-ens-2100'], reference=None, target=None, title='Grell vs. KF')
exps_rc['kf-ens']     = EX(name='kf-ens', exps=['kf-ens','kf-ens-2050','kf-ens-2100'], 
                            colors=None, master='kf-ens', reference=None, target=None, title='KF Ens.')
exps_rc['max']     = EX(name='max', exps=['Observations','max-ens','max-ens-2050','max-ens-2100'], 
                        colors=None, master='max-ens', reference='Observations', target='max-ens', title='IC Ensemble')
exps_rc['maxs']     = EX(name='maxs', exps=['Observations','max-ens','max-ens-2100'], 
                        colors=None, master='max-ens', reference='Observations', target='max-ens', title='IC Ensemble')
exps_rc['ctrl']     = EX(name='ctrl', exps=['Observations','ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         colors=None, master='ctrl-ens', reference='Observations', target='ctrl-ens', title='DE Ensemble')
exps_rc['ctrl-all'] = EX(name='ctrl-all', exps=['Observations','ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         colors=None, master='ctrl-ens', reference='Observations', target='ctrl-ens', title='DE Ensemble')
exps_rc['ctrl-prj'] = EX(name='ctrl-prj', exps=['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         colors=None, master='ctrl-ens', reference=None, target=None, title='DE Ensemble')
# exps_rc['max-phys']     = EX(name='max-phys', exps=['max-ens','phys-ens','max-ens-2050','phys-ens-2050','max-ens-2100','phys-ens-2100'], 
#                         colors=None, master=['max-ens','max-ens-2050','max-ens-2100'], reference=None, target=None, title='IC vs. Phys.')
exps_rc['max-phys']     = EX(name='max-phys', exps=['max-ens','phys-ens','max-ens-2100','phys-ens-2100'], 
                        colors=None, master=['max-ens','max-ens-2100'], reference=None, target=None, title='IC vs. Phys.')
exps_rc['max-ctrl']  = EX(name='max-ctrl', exps=['max-ens','max-ctrl','max-ens-2050','max-ctrl-2050','max-ens-2100','max-ctrl-2100'], 
                        colors=None, master=['max-ens','max-ens-2050','max-ens-2100'], reference=None, target=None, title='IC & Max')
exps_rc['max-shape'] = EX(name='shape', exps=['Observations','max-ens','max-ens-2050','max-ens-2100'], 
                          colors=None, master='max-ens', reference='Observations', target=None, title='Projected Shape Change')
exps_rc['max-shape-cu'] = EX(name='shape-cu', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                          colors=None, master='max-ens', reference='max-ens', target=None, title='Projected Shape Change')
exps_rc['new']     = EX(name='new', exps=['max-ctrl','new-v361','max-ctrl-2050','new-v361-2050','max-ctrl-2100','new-v361-2100'], 
                        colors=None, master=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], reference=None, target=None, title='New V3.6.1')
exps_rc['mex']     = EX(name='mex', exps=['Observations','mex-ens','mex-ens-2050','mex-ens-2100'], 
                        colors=None, master='mex-ens', reference='Observations', target='mex-ens', title='IC Ensemble (Ext.)')
exps_rc['phys']    = EX(name='phys',exps=['Observations','phys-ens','phys-ens-2050','phys-ens-2100'], 
                        colors=None, master='phys-ens', reference='Observations', target='phys-ens', title='Physics Ensemble')
exps_rc['phys-val']    = EX(name='phys-val',exps=['Observations','phys-ens','phys-ens_d01'], 
                        colors=None, master='phys-ens', reference='Observations', target=None, title='Physics Validation')
exps_rc['phys-prj']    = EX(name='phys-prj',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], 
                            colors=None, master='phys-ens', reference=None, target=None, title='Physics Projection')
exps_rc['phys-shape'] = EX(name='shape', exps=['Observations','phys-ens','phys-ens-2050','phys-ens-2100'], 
                           colors=None, master='phys-ens', reference='Observations', target=None, title='Projected Shape Change')
# exps_rc['seaice']  = EX(name='seaice', exps=['Observations','max-ctrl','max-seaice-2050','max-seaice-2100'], 
#                         colors=None, master='max-ctrl', reference=None, target=None, title='Seaice Exp.')
exps_rc['seaice']  = EX(name='seaice', exps=['max-ens-2050','max-seaice-2050','max-ens-2100','max-seaice-2100'], 
                        colors=None, master=['max-ens-2050','max-ens-2100'], reference=None, target=None, title='Seaice Exp.')
exps_rc['ens-all']  = EX(name='ens-all', exps=['Observations','max-ens','max-ctrl','erai-max','phys-ens','max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050','max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                        colors=None, master=['max-ens','max-ens-2050','max-ens-2100'], reference='Observations', target=None, title='Ensembles')
exps_rc['ens-1980']  = EX(name='ens-1980', exps=['Observations','max-ens','max-ctrl','erai-max','phys-ens'], 
                        colors=None, master='max-ens', reference='Observations', target=None, title='Ensembles 1980')
exps_rc['ens-2050']  = EX(name='ens-2050', exps=['max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050'], 
                        colors=None, master='max-ens-2050', reference=None, target=None, title='Ensembles 2050')
exps_rc['ens-2100']  = EX(name='ens-2100', exps=['max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                        colors=None, master='max-ens-2100', reference=None, target=None, title='Ensembles 2100')
# Marc's simulations for Ontario
exps_rc['marc-prj'] = EX(name='marc-prj',exps=['marc-m','marc-mm', 'marc-m-2050','marc-mm-2050'], 
                         colors=None, master=None, reference='Observations', target=None, title='')
exps_rc['marc-val'] = EX(name='marc-val',exps=['Observations', 'marc-g','marc-gg', 'marc-m','marc-mm', 'marc-t'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')
exps_rc['marc-g']   = EX(name='marc-g',  exps=['Observations', 'marc-g', 'marc-g-2050'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')
exps_rc['marc-gg']  = EX(name='marc-gg', exps=['Observations', 'marc-g','marc-gg', 'marc-g-2050','marc-gg-2050'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')
exps_rc['marc-m']   = EX(name='marc-m',  exps=['Observations', 'marc-m', 'marc-m-2050'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')
exps_rc['marc-mm']  = EX(name='marc-mm', exps=['Observations', 'marc-mm', 'marc-mm-2050'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')
exps_rc['marc-t']   = EX(name='marc-t',  exps=['Observations', 'marc-t', 'marc-t-2050'], 
                         colors=None, master='Observations', reference='Observations', target=None, title='')


# function to load basin-averaged time-series and compute distribution
# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
def loadBasinEnsemble(names=None, seasons=None, varlist=None, basins=None, provs=None,
                      lfit=True, dist=None, dist_args=None, reference=None, target=None,
                      lrescale=False, lflatten=True, sample_axis=None,
                      shapefile=None, lglobalScale=False, lcrossval=False, ncv=0.2, lshort=True,
                      aggregation=None, filetypes=None, domain=None, lbootstrap=False, nbs=100,
                      load_list=None, lproduct='outer', **kwargs):
  ''' convenience function to load basin-averaged ensemble '''
  if load_list is None: load_list = []  
  # use varlist from here or hydro-module
  if isinstance(varlist,basestring) and varlist in variables_rc: 
    if filetypes is None: filetypes = variables_rc[varlist].files    
    varlist = variables_rc[varlist].vars
  # load shape ensemble
  bsnens = loadShapeEnsemble(names=names, seasons=seasons, basins=basins, provs=provs, shapefile=shapefile, 
                             varlist=varlist, aggregation=aggregation, filetypes=filetypes, domain=domain, 
                             load_list=load_list, lproduct=lproduct, **kwargs)
#   # generate argument list
#   args_list = expandArgumentList(names=names, seasons=seasons, aggregation=aggregation, basin=basins, 
#                                  varlist=varlist, filetypes=filetypes, domain=domain, 
#                                  expand_list=load_list, lproduct=lproduct, **kwargs)
#   # loop over argument pairs
#   bsnens = []
#   for args in args_list:
#     print args
#     # load ensemble (no iteration here)
#     tmpens = loadShapeEnsemble(load_list=None, lproduct=lproduct, **args)
#     bsnens.append(tmpens)
#     # generate matching datasets with fitted distributions
#     dist = args.pop('dist',None)
#     if lfit:
#       if dist is None: 
#         season = args.get('season',None)
#         aggregation = args.get('aggregation',None)
#         if ( ( season is not None and season.lower() == 'annual' ) and 
#              ( aggregation is not None and aggregation.lower() == 'mean' ) ): 
#           # N.B.: annual averages are usually normally distributed
#           dist = 'norm'
#       # compute fit
#       fitens.append(tmpens.fitDist(dist=dist, axis='time', lbootstrap=lbootstrap, nbs=nbs))

  # generate matching datasets with fitted distributions
  fitens, sclens = addDistFit(ensemble=bsnens, provs=provs, seasons=seasons, varlist=varlist, 
                              basins=basins, lfit=lfit, lflatten=lflatten, 
                              lrescale=lrescale, reference=reference, target=target, aggregation=aggregation,   
                              lbootstrap=lbootstrap, nbs=nbs, sample_axis=sample_axis, lglobalScale=lglobalScale,
                              lcrossval=lcrossval, ncv=ncv, dist=dist, dist_args=dist_args, filetypes=filetypes, 
                              domain=domain, load_list=load_list, lproduct=lproduct)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  if lshort:
    if lfit and lrescale: return bsnens, fitens, sclens
    elif lfit: return bsnens, fitens
    else: return bsnens
  else:
    return bsnens, fitens, sclens

# define new load fct. (batch args: load_list=['season','prov',], lproduct='outer')
def loadStationEnsemble(names=None, provs=None, seasons=None, clusters=None, basins=None, varlist=None,
                        station=None, constraints=None, lfit=True, lflatten=None, master=None, lall=True, 
                        lrescale=False, reference=None, target=None, aggregation='max', lminmax=False,
                        lbootstrap=False, nbs=30, sample_axis=None, slices=None, obsslices=None, 
                        lglobalScale=False, ensemble_axis='ensemble', lensembleAxis=False, 
                        lcrossval=False, ncv=0.2, dist=None, dist_args=None, lshort=True, reduction=None,
                        filetypes=None, domain=None, name_tags=None, cluster_name=None, years=None,
                        ensemble_list=None, ensemble_product='inner', load_list=None, lproduct='outer'):
  ''' convenience function to load station and basin data from ensembles (concatenated) '''
  if lrescale and not lfit: raise ArgumentError
  load_list = [] if load_list is None else load_list[:] # use a copy, since the list may be modified
  
  # figure out varlist
  if isinstance(varlist,basestring):
    if station is None:
      if varlist.lower().find('precip') >= 0: 
        station = 'ecprecip'
      elif varlist.lower().find('temp') >= 0: 
        station = 'ectemp'
    if filetypes is None: filetypes = variables_rc[varlist].files[:]
    varlist = variables_rc[varlist].vars[:]
  if cluster_name: varlist = varlist[:] + [cluster_name] # need to load this variable!
    
  # prepare arguments
  if provs or clusters:
    if constraints is None: constraints = constraints_rc.copy()
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
                          lensembleAxis=lensembleAxis, reduction=reduction)
  # generate matching datasets with fitted distributions
  fitens, sclens = addDistFit(ensemble=stnens, provs=provs, seasons=seasons, varlist=varlist, station=station, 
                              constraints=constraints, basins=basins, lfit=lfit, lflatten=lflatten, 
                              lrescale=lrescale, reference=reference, target=target, aggregation=aggregation,   
                              lbootstrap=lbootstrap, nbs=nbs, sample_axis=sample_axis, lglobalScale=lglobalScale,
                              lcrossval=lcrossval, ncv=ncv, dist=dist, dist_args=dist_args, filetypes=filetypes, 
                              domain=domain, load_list=load_list, lproduct=lproduct)
  # return ensembles (will be wrapped in a list, if BatchLoad is used)
  if lshort:
    if lfit and lrescale: return stnens, fitens, sclens
    elif lfit: return stnens, fitens
    else: return stnens
  else:
    return stnens, fitens, sclens
  
def addDistFit(ensemble=None, provs=None, seasons=None, varlist=None, station=None, constraints=None, basins=None, 
               lfit=True, lflatten=None, lrescale=False, reference=None, target=None, aggregation='max',   
               lbootstrap=False, nbs=30, sample_axis=None, lglobalScale=False, lcrossval=False, ncv=0.2, 
               dist=None, dist_args=None, filetypes=None, domain=None, load_list=None, lproduct='outer'): 
  ''' add distribution fits to ensemble; optionally also rescale'''
  
  # find appropriate sample axis
  if lflatten: 
    if sample_axis is not None: raise ArgumentError, sample_axis
  elif sample_axis is None: # auto-detect
    for saxis in ('time','year'):
      if all([all(ens.hasAxis(saxis)) for ens in ensemble]): 
        sample_axis = saxis; break
    if sample_axis is None: raise AxisError, "No sample axis detected" 
  else:
    if isinstance(sample_axis,basestring):
      if not all([all(ens.hasAxis(sample_axis)) for ens in ensemble]): raise AxisError, sample_axis 
    elif isinstance(sample_axis,(tuple,list)):
      # check that axes are there
      for ax in sample_axis: 
        if not all([all(ens.hasAxis(ax)) for ens in ensemble]): raise AxisError, ax
      # merge axes
      ensemble = [ens.mergeAxes(axes=sample_axis, new_axis='sample', asVar=True, linplace=False, \
                              lcheckAxis=False) for ens in ensemble]
      sample_axis = 'sample'
    else: raise AxisError, sample_axis
  # perform fit or return dummy
  if dist_args is None: dist_args = dict()
  if lfit: fitens = [ens.fitDist(lflatten=lflatten, axis=sample_axis, lcrossval=lcrossval, ncv=ncv,
                                 lignoreParams=True, lbootstrap=lbootstrap, nbs=nbs,
                                 dist=dist, **dist_args) for ens in ensemble]
  else: fitens = [None]*len(ensemble)
  
  # rescale fitted distribution (according to certain rules)
  if lrescale:
    # expand target list
    if isinstance(target, (list,tuple)):  
      expand_list = load_list[:]
      if 'names' in expand_list: expand_list[expand_list.index('names')] = 'target' 
      kwarg_list = expandArgumentList(target=target, season=seasons, prov=provs, aggregation=aggregation, 
                                      basins=basins, station=station, constraints=constraints, 
                                      varlist=varlist, filetypes=filetypes, domain=domain,
                                      expand_list=expand_list, lproduct=lproduct)
      targets = [kwarg['target'] for kwarg in kwarg_list]
    else: targets = [target]*len(fitens)
    if isinstance(reference, (list,tuple)): raise NotImplementedError # don't expand reference list
    # use global reference, if necessary
    if isinstance(reference,basestring) and not all(reference in fit for fit in fitens):
      i = 0 
      while i < len(fitens) and reference not in fitens[i]: i += 1
      if i >= len(fitens): raise ArgumentError, "Reference {:s} not found in any dataset!".format(reference)
      reference = fitens[i][reference]
    sclens = [rescaleDistributions(fit, reference=reference, target=tgt, lglobal=lglobalScale) for fit,tgt in zip(fitens,targets)]
  else: sclens = [None]*len(ensemble)
  
  # return results
  return fitens, sclens



## abuse main section for testing
if __name__ == '__main__':
  
#   test = 'basin_ensemble'
  test = 'station_ensemble'
#   test = 'rescaling'

  # reduce sample size
  constraints_rc['min_len'] = 0 
  constraints_rc['lat'] = (45,55) 
  constraints_rc['max_zerr'] = 100
  
  # test load function for station ensemble
  if test == 'basin_ensemble':
    
    # some settings for tests
    exp = 'max-prj'; exps = exps_rc[exp]
    basins = ['FRB','ARB']; seasons = 'jas'
    load_list = ['basins']

    bsnens, fitens = loadBasinEnsemble(names=exps.exps, basins=basins, seasons=seasons, varlist=['aSM'], 
                                       filetypes=['lsm'], lfit=True, aggregation='mean', dist='norm',
                                       load_list=load_list, lproduct='outer')
    # print diagnostics
    print bsnens[0]; print ''
    print fitens[0][0]
    #print fitens[0] if fitens is not None else fitens ; print ''
    assert len(bsnens) == len(basins)
    assert fitens[0][0].atts.shape_name == basins[0]
    
  # test load function for station ensemble
  elif test == 'station_ensemble':
    
    # some settings for tests
#     exp = 'val'; exps = ['EC', 'erai-max', 'max-ctrl'] # 
#     exp = 'max-all'; exps = exps_rc[exp]; provs = ('BC','AB'); clusters = None
#     exp = 'marc-prj'; exps = exps_rc[exp]; provs = 'ON'; clusters = None
#     provs = ('BC','AB'); clusters = None
    provs = None; clusters = None; lensembleAxis = False; sample_axis = None; lflatten = False
    exp = 'cesm-all'; exps = exps_rc[exp]; provs = None; clusters = [1,3]
    seasons = ['summer']; lfit = True; lrescale = True; lbootstrap = False
    lflatten = False; lensembleAxis = True; sample_axis = ('station','year')
    varlist = ['MaxPrecip_1d', 'MaxPrecip_5d','MaxPreccu_1d'][:1]; filetypes = ['hydro']
    stnens, fitens, sclens = loadStationEnsemble(names=exps.exps, provs=provs, clusters=clusters, 
                                         seasons=seasons, lfit=lfit, master=None, station='ecprecip',
                                         lrescale=lrescale, reference=exps.reference, target=exps.target,
                                         lflatten=lflatten, domain=2, lbootstrap=lbootstrap, nbs=10,
                                         lensembleAxis=lensembleAxis, sample_axis=sample_axis,
                                         varlist=varlist, filetypes=filetypes,
                                         load_list=['season',], lproduct='outer', lcrossval=None)
    # print diagnostics
    print stnens[0][0]; print ''
    print fitens[0][0].MaxPrecip_1d.atts.sample_axis; print ''
    print fitens[0][1] if fitens is not None else fitens ; print ''
    print sclens[0] if sclens is not None else sclens ; print ''
    assert len(stnens) == len(seasons)
    print stnens[0][0]
    print stnens[0][1].MaxPrecip_1d.mean()
  
  elif test == 'rescaling':
    
    # some settings for tests
    exps = ['val','prj']; exps = [exps_rc[exp][:3] for exp in exps]
    provs = ('BC',); seasons = 'summer'
    lflatten = True; lfit = True

    # rescale    
    stnens,fitens,scalens = loadStationEnsemble(names=exps, provs=provs, seasons=seasons, varlist='precip', 
                                                lrescale=True, reference='EC', target='max-ens',
                                                filetypes='hydro', domain=None, lfit=lfit, lflatten=lflatten,
                                                lbootstrap=False, nbs=10,  
                                                load_list=['names'], lproduct='outer')
    
    # special rescaling for projection plots
    scalens = [rescaleDistributions(ens, reference=fitens[0]['EC'], target='max-ens') for ens in fitens]

    # print diagnostics
    print fitens[0]; print ''
    print scalens[0]; print ''
    print scalens[1][1]
