'''
Created on 2013-11-08

This module contains meta data for WRF simulations for the Great Lakes region. 

@author: Andre R. Erler, GPL v3
'''

from collections import OrderedDict
from datasets.common import addLoadFcts
from datasets.WRF import Exp as WRF_Exp

## EXP class with specific default values
class Exp(WRF_Exp): 
  parameters = WRF_Exp.parameters.copy()
  defaults = WRF_Exp.defaults.copy()
  # set some project specific defaults
  defaults['parent'] = 'Ctrl-1' # CESM simulations that is driving most of the WRF runs   
  defaults['project'] = 'GreatLakes' # most WRF runs so far are from this project
  defaults['domains'] = 2 # most WRF runs have two domains
  defaults['begindate'] = '1979-01-01' # most WRF runs start in 1979
  defaults['grid'] = 'glb1' # Marc's Great Lakes domain
  
## list of experiments
# N.B.: This is the reference list, with unambiguous, unique keys and no aliases/duplicate entries  
experiments = OrderedDict() # dictionary of experiments
## Great Lakes experiments
# new experiments with V3.9 and lake model
experiments['erai-conus']    = Exp(shortname='conus', name='erai-conus', title='Conus (ERA-I)', parent='ERAI', domains=1, begindate='2005-01-01',)
experiments['erai-c-glerl25']    = Exp(shortname='c-glerl25', name='erai-c-glerl25', title='Conus (small)', parent='ERAI', domains=1, begindate='2005-01-01', grid='glb2')
experiments['erai-t-glerl25']    = Exp(shortname='t-glerl25', name='erai-t-glerl25', title='T (ERA-I, GLERL 25l)', parent='ERAI', domains=1, begindate='2005-01-01')
experiments['erai-t-glerl']    = Exp(shortname='t-glerl', name='erai-t-glerl', title='T (ERA-I, GLERL Lake)', parent='ERAI', domains=1, begindate='2005-01-01')
experiments['erai-t-nolake06']  = Exp(shortname='t-nolake', name='erai-t-nolake06', title='T (ERA-I, No Lake)', parent='ERAI', domains=1, begindate='2005-01-01')
experiments['erai-t-lake06']    = Exp(shortname='t-lake06', name='erai-t-lake06', title='T (ERA-I, CLM Lake)', parent='ERAI', domains=1, begindate='2005-01-01')
experiments['erai-g-lake06']    = Exp(shortname='g-lake06', name='erai-g-lake06', title='G (ERA-I, CLM Lake)', parent='ERAI', domains=1, begindate='2005-01-01')
experiments['erai-t-lake']    = Exp(shortname='t-lake', name='erai-t-lake', title='T80 (ERA-I, CLM Lake)', parent='ERAI', domains=1,)
experiments['erai-g-lake']    = Exp(shortname='g-lake', name='erai-g-lake', title='G80 (ERA-I, CLM Lake)', parent='ERAI', domains=1,)
# new experiments with V3.6 and new physics: V-Ensemble
experiments['erai-v36'] = Exp(shortname='erai-v36', name='erai-v36', title='V3.6 (ERA-I)', parent='ERAI', domains=1)
# lake sensitivity experiments
experiments['erai-v36-shallow'] = Exp(shortname='v36-shallow', name='erai-v36-shallow', title='V3.6 (ERA-I, Shallow)', parent='ERAI', domains=1)
experiments['erai-v36-noflake'] = Exp(shortname='v36-noflake', name='erai-v36-noflake', title='V3.6 (ERA-I, No FLake)', parent='ERAI', domains=1)
experiments['erai-v36-lake']    = Exp(shortname='v36-lake', name='erai-v36-lake', title='V3.6 (ERA-I, CLM Lake)', parent='ERAI', domains=1)
experiments['g-ctrl-noflake']      = Exp(shortname='g-ctrl-noflake', name='g-noflake', title='G (No FLake)', parent='Ctrl-1',)
experiments['g-ctrl-2100-noflake'] = Exp(shortname='g-ctrl-2100-noflake', name='g-2100-noflake', title='G (2100, No FLake)', begindate='2085-01-01', parent='Ctrl-1-2100')
experiments['g-ctrl-shallow']      = Exp(shortname='g-ctrl-shallow', name='g-shallow', title='G (Shallow)', parent='Ctrl-1',)
experiments['g-ctrl-2100-shallow'] = Exp(shortname='g-ctrl-2100-shallow', name='g-2100-shallow', title='G (2100, Shallow)', begindate='2085-01-01', parent='Ctrl-1-2100')
# G-Ensemble
experiments['erai-max'] = Exp(shortname='erai-max', name='erai-max', title='Max (ERA-I)', parent='ERAI', project='WesternCanada', grid='arb2')
experiments['erai-g'] = Exp(shortname='erai-g', name='erai-g', title='G (ERA-I)', parent='ERAI')
experiments['g-ctrl'] = Exp(shortname='g-ctrl', name='g-ctrl', title='G-Ctrl', parent='Ctrl-1', ensemble='g-ensemble')
experiments['g-ctrl-2050'] = Exp(shortname='g-ctrl-2050', name='g-ctrl-2050', title='G-Ctrl (2050)', begindate='2045-01-01', parent='Ctrl-1-2050', ensemble='g-ensemble-2050')
experiments['g-ctrl-2100'] = Exp(shortname='g-ctrl-2100', name='g-ctrl-2100', title='G-Ctrl (2100)', begindate='2085-01-01', parent='Ctrl-1-2100', ensemble='g-ensemble-2100')
experiments['g-ens-A'] = Exp(shortname='g-ens-A', name='g-ens-A', title='G-Ens-A', parent='Ens-A', ensemble='g-ensemble')
experiments['g-ens-A-2050'] = Exp(shortname='g-ens-A-2050', name='g-ens-A-2050', title='G-Ens-A (2050)', begindate='2045-01-01', parent='Ens-A-2050', ensemble='g-ensemble-2050')
experiments['g-ens-A-2100'] = Exp(shortname='g-ens-A-2100', name='g-ens-A-2100', title='G-Ens-A (2100)', begindate='2085-01-01', parent='Ens-A-2100', ensemble='g-ensemble-2100')
experiments['g-ens-B'] = Exp(shortname='g-ens-B', name='g-ens-B', title='G-Ens-B', parent='Ens-B', ensemble='g-ensemble')
experiments['g-ens-B-2050'] = Exp(shortname='g-ens-B-2050', name='g-ens-B-2050', title='G-Ens-B (2050)', begindate='2045-01-01', parent='Ens-B-2050', ensemble='g-ensemble-2050')
experiments['g-ens-B-2100'] = Exp(shortname='g-ens-B-2100', name='g-ens-B-2100', title='G-Ens-B (2100)', begindate='2085-01-01', parent='Ens-B-2100', ensemble='g-ensemble-2100')
experiments['g-ens-C'] = Exp(shortname='g-ens-C', name='g-ens-C', title='G-Ens-C', parent='Ens-C', ensemble='g-ensemble')
experiments['g-ens-C-2050'] = Exp(shortname='g-ens-C-2050', name='g-ens-C-2050', title='G-Ens-C (2050)', begindate='2045-01-01', parent='Ens-C-2050', ensemble='g-ensemble-2050')
experiments['g-ens-C-2100'] = Exp(shortname='g-ens-C-2100', name='g-ens-C-2100', title='G-Ens-C (2100)', begindate='2085-01-01', parent='Ens-C-2100', ensemble='g-ensemble-2100')
experiments['g-ensemble'] = Exp(shortname='g-ens', name='g-ensemble', title='G-Ensemble', parent='Ens')
experiments['g-ensemble-2050'] = Exp(shortname='g-ens-2050', name='g-ensemble-2050', title='G-Ens. (2050)', begindate='2045-01-01', parent='Ens-2050')
experiments['g-ensemble-2100'] = Exp(shortname='g-ens-2100', name='g-ensemble-2100', title='G-Ens. (2100)', begindate='2085-01-01', parent='Ens-2100')
# low-res G-Ensemble (90km resolution)
experiments['erai-g3'] = Exp(shortname='erai-g3', name='erai-g3', title='G (ERA-I, 90km)', grid='glb1-90km', domains=1, parent='ERAI')
experiments['g3-ctrl'] = Exp(shortname='g3-ctrl', name='g3-ctrl', title='G-Ctrl (90km)', parent='Ctrl-1', ensemble='g3-ensemble', grid='glb1-90km', domains=1)
experiments['g3-ctrl-2050'] = Exp(shortname='g3-ctrl-2050', name='g3-ctrl-2050', title='G-Ctrl (2050, 90km)', begindate='2045-01-01', parent='Ctrl-1-2050', ensemble='g3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['g3-ctrl-2100'] = Exp(shortname='g3-ctrl-2100', name='g3-ctrl-2100', title='G-Ctrl (2100, 90km)', begindate='2085-01-01', parent='Ctrl-1-2100', ensemble='g3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['g3-ens-A'] = Exp(shortname='g3-ens-A', name='g3-ens-A', title='G-Ens-A', parent='Ens-A (90km)', ensemble='g3-ensemble', grid='glb1-90km', domains=1)
experiments['g3-ens-A-2050'] = Exp(shortname='g3-ens-A-2050', name='g3-ens-A-2050', title='G-Ens-A (2050, 90km)', begindate='2045-01-01', parent='Ens-A-2050', ensemble='g3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['g3-ens-A-2100'] = Exp(shortname='g3-ens-A-2100', name='g3-ens-A-2100', title='G-Ens-A (2100, 90km)', begindate='2085-01-01', parent='Ens-A-2100', ensemble='g3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['g3-ens-B'] = Exp(shortname='g3-ens-B', name='g3-ens-B', title='G-Ens-B', parent='Ens-B (90km)', ensemble='g3-ensemble', grid='glb1-90km', domains=1)
experiments['g3-ens-B-2050'] = Exp(shortname='g3-ens-B-2050', name='g3-ens-B-2050', title='G-Ens-B (2050, 90km)', begindate='2045-01-01', parent='Ens-B-2050', ensemble='g3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['g3-ens-B-2100'] = Exp(shortname='g3-ens-B-2100', name='g3-ens-B-2100', title='G-Ens-B (2100, 90km)', begindate='2085-01-01', parent='Ens-B-2100', ensemble='g3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['g3-ens-C'] = Exp(shortname='g3-ens-C', name='g3-ens-C', title='G-Ens-C', parent='Ens-C (90km)', ensemble='g3-ensemble', grid='glb1-90km', domains=1)
experiments['g3-ens-C-2050'] = Exp(shortname='g3-ens-C-2050', name='g3-ens-C-2050', title='G-Ens-C (2050, 90km)', begindate='2045-01-01', parent='Ens-C-2050', ensemble='g3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['g3-ens-C-2100'] = Exp(shortname='g3-ens-C-2100', name='g3-ens-C-2100', title='G-Ens-C (2100, 90km)', begindate='2085-01-01', parent='Ens-C-2100', ensemble='g3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['g3-ensemble'] = Exp(shortname='g3-ens', name='g3-ensemble', title='G-Ens. (90km)', parent='Ens', grid='glb1-90km', domains=1)
experiments['g3-ensemble-2050'] = Exp(shortname='g3-ens-2050', name='g3-ensemble-2050', title='G-Ens. (2050, 90km)', begindate='2045-01-01', parent='Ens-2050', grid='glb1-90km', domains=1)
experiments['g3-ensemble-2100'] = Exp(shortname='g3-ens-2100', name='g3-ensemble-2100', title='G-Ens. (2100, 90km)', begindate='2085-01-01', parent='Ens-2100', grid='glb1-90km', domains=1)
# T-Ensemble
experiments['erai-t'] = Exp(shortname='erai-t', name='erai-t', title='T (ERA-I)', parent='ERAI')
experiments['t-ctrl'] = Exp(shortname='t-ctrl', name='t-ctrl', title='T-Ctrl', parent='Ctrl-1', ensemble='t-ensemble-2100')
experiments['t-ctrl-2050'] = Exp(shortname='t-ctrl-2050', name='t-ctrl-2050', title='T-Ctrl (2050)', begindate='2045-01-01', parent='Ctrl-1-2050', ensemble='t-ensemble-2050')
experiments['t-ctrl-2100'] = Exp(shortname='t-ctrl-2100', name='t-ctrl-2100', title='T-Ctrl (2100)', begindate='2085-01-01', parent='Ctrl-1-2100', ensemble='t-ensemble-2100')
experiments['t-ens-A'] = Exp(shortname='t-ens-A', name='t-ens-A', title='T-Ens-A', parent='Ens-A', ensemble='t-ensemble')
experiments['t-ens-A-2050'] = Exp(shortname='t-ens-A-2050', name='t-ens-A-2050', title='T-Ens-A (2050)', parent='Ens-A-2050', ensemble='t-ensemble-2050', begindate='2045-01-01')
experiments['t-ens-A-2100'] = Exp(shortname='t-ens-A-2100', name='t-ens-A-2100', title='T-Ens-A (2100)', parent='Ens-A-2100', ensemble='t-ensemble-2100', begindate='2085-01-01')
experiments['t-ens-B'] = Exp(shortname='t-ens-B', name='t-ens-B', title='T-Ens-B', parent='Ens-B', ensemble='t-ensemble')
experiments['t-ens-B-2050'] = Exp(shortname='t-ens-B-2050', name='t-ens-B-2050', title='T-Ens-B (2050)', begindate='2045-01-01', parent='Ens-B-2050', ensemble='t-ensemble-2050')
experiments['t-ens-B-2100'] = Exp(shortname='t-ens-B-2100', name='t-ens-B-2100', title='T-Ens-B (2100)', begindate='2085-01-01', parent='Ens-B-2100', ensemble='t-ensemble-2100')
experiments['t-ens-C'] = Exp(shortname='t-ens-C', name='t-ens-C', title='T-Ens-C', parent='Ens-C', ensemble='t-ensemble')
experiments['t-ens-C-2050'] = Exp(shortname='t-ens-C-2050', name='t-ens-C-2050', title='T-Ens-C (2050)', begindate='2045-01-01', parent='Ens-C-2050', ensemble='t-ensemble-2050')
experiments['t-ens-C-2100'] = Exp(shortname='t-ens-C-2100', name='t-ens-C-2100', title='T-Ens-C (2100)', begindate='2085-01-01', parent='Ens-C-2100', ensemble='t-ensemble-2100')
experiments['t-ensemble'] = Exp(shortname='t-ens', name='t-ensemble', title="T-Ensemble", parent='Ens')
experiments['t-ensemble-2050'] = Exp(shortname='t-ens-2050', name='t-ensemble-2050', title="T-Ens. (2050)", begindate='2045-01-01', parent='Ens-2050')
experiments['t-ensemble-2100'] = Exp(shortname='t-ens-2100', name='t-ensemble-2100', title="T-Ens. (2100)", begindate='2085-01-01', parent='Ens-2100')
# low-res T-Ensemble (90km resolution)
experiments['erai-t3'] = Exp(shortname='erai-t3', name='erai-t3', title='T (ERA-I, 90km)', grid='glb1-90km', domains=1)
experiments['t3-ctrl'] = Exp(shortname='t3-ctrl', name='t3-ctrl', title='T-Ctrl (90km)', ensemble='t3-ensemble', grid='glb1-90km', domains=1)
experiments['t3-ctrl-2050'] = Exp(shortname='t3-ctrl-2050', name='t3-ctrl-2050', title='T-Ctrl (2050, 90km)', begindate='2045-01-01', ensemble='t3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['t3-ctrl-2100'] = Exp(shortname='t3-ctrl-2100', name='t3-ctrl-2100', title='T-Ctrl (2100, 90km)', begindate='2085-01-01', ensemble='t3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['t3-ens-A'] = Exp(shortname='t3-ens-A', name='t3-ens-A', title='T-Ens-A', parent='Ens-A (90km)', ensemble='t3-ensemble', grid='glb1-90km', domains=1)
experiments['t3-ens-A-2050'] = Exp(shortname='t3-ens-A-2050', name='t3-ens-A-2050', title='T-Ens-A (2050, 90km)', parent='Ens-A-2050', ensemble='t3-ensemble-2050', begindate='2045-01-01', grid='glb1-90km', domains=1)
experiments['t3-ens-A-2100'] = Exp(shortname='t3-ens-A-2100', name='t3-ens-A-2100', title='T-Ens-A (2100, 90km)', parent='Ens-A-2100', ensemble='t3-ensemble-2100', begindate='2085-01-01', grid='glb1-90km', domains=1)
experiments['t3-ens-B'] = Exp(shortname='t3-ens-B', name='t3-ens-B', title='T-Ens-B', parent='Ens-B (90km)', ensemble='t3-ensemble', grid='glb1-90km', domains=1)
experiments['t3-ens-B-2050'] = Exp(shortname='t3-ens-B-2050', name='t3-ens-B-2050', title='T-Ens-B (2050, 90km)', begindate='2045-01-01', parent='Ens-B-2050', ensemble='t3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['t3-ens-B-2100'] = Exp(shortname='t3-ens-B-2100', name='t3-ens-B-2100', title='T-Ens-B (2100, 90km)', begindate='2085-01-01', parent='Ens-B-2100', ensemble='t3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['t3-ens-C'] = Exp(shortname='t3-ens-C', name='t3-ens-C', title='T-Ens-C', parent='Ens-C (90km)', ensemble='t3-ensemble', grid='glb1-90km', domains=1)
experiments['t3-ens-C-2050'] = Exp(shortname='t3-ens-C-2050', name='t3-ens-C-2050', title='T-Ens-C (2050, 90km)', begindate='2045-01-01', parent='Ens-C-2050', ensemble='t3-ensemble-2050', grid='glb1-90km', domains=1)
experiments['t3-ens-C-2100'] = Exp(shortname='t3-ens-C-2100', name='t3-ens-C-2100', title='T-Ens-C (2100, 90km)', begindate='2085-01-01', parent='Ens-C-2100', ensemble='t3-ensemble-2100', grid='glb1-90km', domains=1)
experiments['t3-ensemble'] = Exp(shortname='t3-ens', name='t3-ensemble', title='T-Ens. (90km)', parent='Ens', grid='glb1-90km', domains=1)
experiments['t3-ensemble-2050'] = Exp(shortname='t3-ens-2050', name='t3-ensemble-2050', title='T-Ens. (2050, 90km)', begindate='2045-01-01', parent='Ens-2050', grid='glb1-90km', domains=1)
experiments['t3-ensemble-2100'] = Exp(shortname='t3-ens-2100', name='t3-ensemble-2100', title='T-Ens. (2100, 90km)', begindate='2085-01-01', parent='Ens-2100', grid='glb1-90km', domains=1)
# sensitivity experiments
experiments['gg-ctrl'] = Exp(shortname='gg-ctrl', name='gg-ctrl', title='G-Ctrl (N-MP)')
experiments['gg2-ctrl'] = Exp(shortname='gg-ctrl', name='gg-ctrl', title='G-Ctrl (N-MP)')
experiments['gg-ctrl-2050'] = Exp(shortname='gg-ctrl-2050', name='gg-ctrl-2050', title='G-Ctrl (2050, N-MP)', begindate='2045-01-01')
experiments['gg-ctrl-2100'] = Exp(shortname='gg-ctrl-2100', name='gg-ctrl-2100', title='G-Ctrl (2100, N-MP)', begindate='2085-01-01')
experiments['gg2-ctrl-2100'] = Exp(shortname='gg2-ctrl-2100', name='gg2-ctrl-2100', title='G-Ctrl (2100, N-MP, 2nd)', begindate='2085-01-01')
experiments['m-ctrl'] = Exp(shortname='m-ctrl', name='m-ctrl', title='M-Ctrl')
experiments['m-ctrl-2050'] = Exp(shortname='m-ctrl-2050', name='m-ctrl-2050', title='M-Ctrl (2050)', begindate='2045-01-01')
experiments['m-ctrl-2100'] = Exp(shortname='m-ctrl-2100', name='m-ctrl-2100', title='M-Ctrl (2100)', begindate='2085-01-01')
experiments['mm-ctrl'] = Exp(shortname='mm-ctrl', name='mm-ctrl', title='M-Ctrl (N-MP)')
experiments['mm-ctrl-2050'] = Exp(shortname='mm-ctrl-2050', name='mm-ctrl-2050', title='M-Ctrl (2050, N-MP)', begindate='2045-01-01')
experiments['mm-ctrl-2100'] = Exp(shortname='mm-ctrl-2100', name='mm-ctrl-2100', title='M-Ctrl (2100, N-MP)', begindate='2085-01-01')
# Marc's Physics Ensemble
experiments['physics-ensemble'] = Exp(shortname='phys-ens', name='physics-ensemble', title="Physics Ensemble", parent='Ctrl-1')
experiments['physics-ensemble-2050'] = Exp(shortname='phys-ens-2050', name='physics-ensemble-2050', title="Phys. Ens. (2050)", begindate='2045-01-01', parent='Ctrl-1-2050')
experiments['physics-ensemble-2100'] = Exp(shortname='phys-ens-2100', name='physics-ensemble-2100', title="Phys. Ens. (2100)", begindate='2085-01-01', parent='Ctrl-1-2100')

## an alternate dictionary using short names and aliases for referencing
exps = OrderedDict()
# use short names where available, normal names otherwise
for key,item in experiments.items():
  exps[key] = item # this prevents name collisions between regions
  if item.shortname is not None: 
    exps[item.shortname] = item
  # both, short and long name are added to list
# add aliases here
WRF_exps = exps # alias for whole dict (short and proper names)
WRF_experiments = experiments # alias for dict with proper names

## dict of ensembles
ensembles = OrderedDict()
# Great Lakes ensembles
ensembles['phys-ens'] = ('g-ctrl','gg-ctrl', 'm-ctrl','mm-ctrl', 't-ctrl')
ensembles['g-ens-1'] = ('g-ctrl', 'g-ens-A',)
ensembles['t-ens-1'] = ('t-ctrl', 't-ens-A',)
ensembles['g3-ens-1'] = ('g3-ctrl', 'g3-ens-A',)
ensembles['t3-ens-1'] = ('t3-ctrl', 't3-ens-A',)
ensembles['g-ens-2'] = ('g-ens-B', 'g-ens-C',)
ensembles['t-ens-2'] = ('t-ens-B', 't-ens-C',)
ensembles['g3-ens-2'] = ('g3-ens-B', 'g3-ens-C',)
ensembles['t3-ens-2'] = ('t3-ens-B', 't3-ens-C',)
ensembles['g-ens'] = ('g-ctrl', 'g-ens-A', 'g-ens-B', 'g-ens-C')
ensembles['t-ens'] = ('t-ctrl', 't-ens-A', 't-ens-B', 't-ens-C')
ensembles['g3-ens'] = ('g3-ctrl', 'g3-ens-A', 'g3-ens-B', 'g3-ens-C')
ensembles['t3-ens'] = ('t3-ctrl', 't3-ens-A', 't3-ens-B', 't3-ens-C')
# N.B.: static & meta data for the ensemble is copied from the first-listed member;
#       this includes station attributes, such as the elevation error 
# add future versions
for ensname,enslist in list(ensembles.items()):
  for suffix in '2050','2100':
    suffix = '-'+suffix
    ensembles[ensname+suffix] = tuple(expname[:-2]+suffix if expname[-2:] == '-1' else expname+suffix for expname in enslist)
# replace names with experiment instances
for ensname,enslist in ensembles.items():
  ensembles[ensname] = tuple(experiments[expname] for expname in enslist)
# make sorted copy
WRF_ens = OrderedDict()
name_list = list(ensembles.keys()); name_list.sort()
for key in name_list:
  WRF_ens[key] = ensembles[key]
ensembles = enss = WRF_ens


## generate loadWRF* versions with these experiments
from datasets.WRF import loadWRF, loadWRF_Shp, loadWRF_Stn, loadWRF_TS, loadWRF_ShpTS, loadWRF_StnTS, loadWRF_Ensemble, loadWRF_ShpEns, loadWRF_StnEns
addLoadFcts(locals(), locals(), exps=WRF_exps, enses=WRF_ens)


if __name__ == '__main__':
    
  ## view/test ensembles
  for name,members in WRF_ens.items():
    s = '  {:s}: '.format(name)
    for member in members: s += ' {:s},'.format(member.name)
    print(s)
