'''
Created on Oct. 24, 2020

This module contains some settings for southern Ontario basins

@author: Andre R. Erler, GPL v3
'''

from collections import OrderedDict
# internal imports
from datasets.WSC import Basin, shape_root
import hgs.HGS as hgs # need to prevent name collisions here
# imports from HGS_settings
# import projects.HGS_settings as default
# from projects.HGS_settings import Station, subfolder_pattern

project_name = 'SON'
project_prefix = ''
# some parameters
gage_scalefactors = 1. # default scalefactor for discharge plots
main_basin  = 'SON' # main basin name
main_grid   = 'son2' # climate data grid
binary_grid = 'son2' # grid used to interpolated binary output
# project folders
project_folder = '{:s}/{:s}/'.format(hgs.root_folder,project_name) # the dataset root folder
project_folder_pattern = '{PROJECT_FOLDER:s}/{GRID:s}/{EXPERIMENT:s}/{CLIM_DIR:s}/{TASK:s}/'
station_file = '{PREFIX:s}o.hydrograph.{WSC_ID0:s}.dat' # general HGS naming convention

# mapping of WSC station names to HGS hydrograph names
station_list = OrderedDict() # this is the main gage of the GRW
# ... using Station named tuples

# look-up tables for WSC/HGS station name conversion                           
WSC_station_list = {stn.WSC:stn.HGS for stn in list(station_list.values())}
HGS_station_list = {stn.HGS:stn.WSC for stn in list(station_list.values())}

# list of groundwater observation wells
pgmn_wells = []

# plotting parameters for HGS simulations
hgs_plotargs = dict() 
# stream observations
hgs_plotargs['Observations'] = dict(color='#959595') # gray
hgs_plotargs['Obs.']         = dict(color='#959595') # gray
hgs_plotargs['WSC Obs.']     = dict(color='#959595') # gray
hgs_plotargs['WSC']          = dict(color='#959595') # gray
# NRCan forcing, different model versions
hgs_plotargs['NRCan']        = dict(color='green')
hgs_plotargs['Steady-State'] = dict(color='#AAA2D8') # purple
hgs_plotargs['Periodic']     = dict(color='#E24B34') # red
hgs_plotargs['Transient']    = dict(color='#62A1C6') # blue
# SnoDAS
hgs_plotargs['SnoDAS']  = dict(color='#62A1C6') # blue
hgs_plotargs['RF-BC']   = dict(color='#E24B34') # red
# ...


## southern Ontario subbasins/watersheds
son_ws_names = OrderedDict()
son_ws_names['AMW'] = 'Ausible_Maitland'
son_ws_names['Big'] = 'Big'
son_ws_names['BSW'] = 'BlackRiver_LakeSimcoe'
son_ws_names['Cat'] = 'Cataraqui'
son_ws_names['C16N'] = 'Credit_16Mile_Niagara'
son_ws_names['DHW'] = 'Don_Highland'
son_ws_names['EMW'] = 'Etobicoke_Mimco'
son_ws_names['GRW'] = 'GrandRiver'
son_ws_names['GKSOT'] = 'Gull_KawarthaLakes_Scugog_Otonabee_Trent-Crowe'
son_ws_names['HDG'] = 'Humber_Don_Ganaraska'
son_ws_names['HLW'] = 'Humber_Local'
son_ws_names['MNPE'] = 'Moira_Napanee_PrinceEdwardCounty'
son_ws_names['Not'] = 'Nottawasaga'
son_ws_names['RDW'] = 'Rouge_Duffins'
son_ws_names['SPBP'] = 'Saugeen_Penetangore_SWGeorgianBay_BrucePeninsula'
son_ws_names['Syd'] = 'Sydenham'
son_ws_names['TCR'] = 'Thames_Cedar_Rondeau'
# now construct proper basin
son_watersheds = OrderedDict()
for name,long_name in son_ws_names.items():
    folder = "{:s}/Basins/Southern Ontario/{:s}/".format(shape_root,long_name)
    shapefile = long_name + '.shp'
    son_watersheds[name] = Basin(name=name, long_name=long_name.replace('_',' '), 
                                 shapefile=shapefile, folder=folder, data_source='Aquanty',
                                 load=False, ldebug=False, subbasins=None, rivers=None, stations=None, shapetype='BSN')


# abuse for testing
if __name__ == '__main__':
     
  test_mode = 'create_grid'
#   test_mode = 'watersheds'
#   test_mode = 'gage_station'
#   test_mode = 'dataset_regrid'
#   test_mode = 'binary_dataset'
#   test_mode = 'dataset'
#   test_mode = 'ensemble'

  ## create a new grid
  if test_mode == 'watersheds':
    
    print("SON Watershed")
    for name,obj in son_watersheds.items():
        print(name,obj.long_name)
  
  ## create a new grid
  elif test_mode == 'create_grid':
    
    import os
    from geodata.gdal import GridDefinition, pickleGridDef, loadPickledGridDef, grid_folder
    
    convention='Proj4'
    ## parameters for UTM 17 Great Lakes grids
#     name = 'glb1' # 5km resolution
#     geotransform = [ -709489.58091004, 5.e3, 0, 4148523.7226861, 0, 5.e3]; size = (405,371)
#     projection = "+proj=utm +zone=17 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
    ## grid for Elisha
#     name = 'uph1' # 5km resolution
#     geotransform = [ 443000, 5.e3, 0, 4774000, 0, 5.e3]; size = (3,6)
#     projection = "+proj=utm +zone=17 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
    ## Southern Ontario grids (also in Great Lakes, UTM 17)
#     name = 'son1' # 5km resolution
#     geotransform = [320920.,5.e3,0,4624073.,0,5.e3]; size = (118,82)
#     name = 'son2' # 1km resolution
#     geotransform = [320920.,1.e3,0,4624073.,0,1.e3]; size = (590,410)
#     projection = "+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ## Hugo's grid for Quebec
#     name = 'hd1' # 5km resolution
#     geotransform = [-479184.769227,5.e3,0,68508.4877898,0,5.e3]; size = (70,49)
#     projection = "+proj=lcc +ellps=NAD83 +datum=NAD83 +lat_0=44.0 +lat_1=46.0 +lat_2=60.0 +lon_0=-68.5  +x_0=0 +y_0=0 +units=m +no_defs"
    ## UTM 18 parameters for South Nation grids
#     name = 'snw1' # 9km resolution
#     geotransform = [401826.125365249,9.e3,0,4851533.71730136,0,9.e3]; size = (22,29)
#     projection = "+proj=utm +zone=18 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
#     convention='Wkt'; projection = 'PROJCS["NAD_1983_UTM_Zone_14N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["false_easting",500000.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-99.0],PARAMETER["scale_factor",0.9996],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
    name = 'snw2' # 2km resolution
    geotransform = [438.e3,2.e3,0,4940.e3,0,2.e3]; size = (44,55)
    projection = "+proj=utm +zone=18 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    # N.B.: (x_0, dx, 0, y_0, 0, dy); (xl,yl)
    #       GT(0),GT(3) are the coordinates of the bottom left corner
    #       GT(1) & GT(5) are pixel width and height
    #       GT(2) & GT(4) are usually zero for North-up, non-rotated maps
    # create grid
    griddef = GridDefinition(name=name, projection=projection, geotransform=geotransform, size=size, 
                             xlon=None, ylat=None, lwrap360=False, geolocator=True, convention='Proj4')

    # save pickle to standard location
    filepath = pickleGridDef(griddef, folder=grid_folder, filename=None, lfeedback=True)
    assert os.path.exists(filepath)
    print('')
    
    # load pickle to make sure it is right
    del griddef
    griddef = loadPickledGridDef(grid=name, res=None, folder=grid_folder)
    print(griddef)

