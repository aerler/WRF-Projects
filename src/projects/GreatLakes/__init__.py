'''
Created on 2016-04-13

A package that contains settings for the GreatLakes region projects for use with the geodata package.

@author: Andre R. Erler, GPL v3
'''

# import figure settings
from .figure_settings import getVariableSettings, getFigureSettings, figure_folder
from warnings import warn

# import map projection settings (basemap)
try: 
    from .map_settings import getSetup, map_folder
except IOError:
    warn("Error importing map settings - possibly due to missing shape data.")
except ImportError:
    warn("Error importing map settings - 'basemap' is likely not installed.")
# N.B.: apparently Basemap is not maintained anymore and dedent does not exist in matplotlib anymore...
#       but replacing dedent with cleandoc as shown below in mpl_toolkits/basemap/proj.py and
#       mpl_toolkits/basemap/__init__.py seems to do the trick
# try:
#     from matplotlib.cbook import dedent
# except ImportError:
#     from inspect import cleandoc as dedent

## import load functions with GreatLakes experiments into local namespace

try:
    # import relevant WRF experiments
    from .WRF_experiments import WRF_exps, WRF_ens
    # import WRF load functions
    from .WRF_experiments import loadWRF, loadWRF_Shp, loadWRF_Stn, loadWRF_TS, loadWRF_ShpTS, loadWRF_StnTS, loadWRF_Ensemble, loadWRF_ShpEns, loadWRF_StnEns
except (ImportError,IOError):
    WRF_exps = None; WRF_ens = None
    warn("Error importing WRF experiments.")

try:
    # also load CESM experiments and functions
    from projects.CESM_experiments import CESM_exps, CESM_ens
    # import CESM load functions
    from projects.CESM_experiments import loadCESM, loadCESM_Shp, loadCESM_Stn, loadCESM_TS, loadCESM_ShpTS, loadCESM_StnTS, loadCESM_Ensemble, loadCESM_ShpEns, loadCESM_StnEns
except (ImportError,IOError):
    CESM_exps = None; CESM_ens = None
    warn("Error importing CESM experiments.")

# add relevant experiments to general load functions
from datasets.common import loadDataset, loadEnsembleTS, addLoadFcts
try:
    from datasets.Unity import loadUnity, loadUnity_Shp, loadUnity_Stn, loadUnity_ShpTS, loadUnity_StnTS # loadUnity_TS doesn't exist
except (ImportError,IOError):
    warn("Error importing Unified Observational Dataset 'Unity'.")
unity_grid = 'glb1_d02' # Unified Dataset default grid
# N.B.: it is recommended to import Unity load fcts. from here
# modify functions (wont affect modified WRF/CESM functions)
addLoadFcts(locals(), locals(), unity_grid=unity_grid , WRF_exps=WRF_exps, WRF_ens=WRF_ens, 
            CESM_exps=CESM_exps, CESM_ens=CESM_ens)


## import shape dictionaries
# provinces, major basins and lakes etc.
try:
    from projects.WSC_basins import basins, provinces, great_lakes, gauges # import the dicts with unique entries
except (ImportError,IOError):
    warn("Error importing shape files and/or WSC module.")
# southern Ontario watersheds
try:
    from .SON_settings import son_ws_names, son_watersheds
except (ImportError,IOError):
    warn("Error importing shape files from SON module.")


# import figure with hydro settings
from .analysis_settings import loadStationEnsemble, loadShapeEnsemble, loadShapeObservations  # load datasets 
from .analysis_settings import loadStationFit, loadShapeFit
from .analysis_settings import exps_rc, variables_rc, constraints_rc
from .analysis_settings import climFigAx, climPlot, evaFigAx, distPlot, quantPlot # plotting
