'''
Created on Aug 27, 2016

This module contains meta data for all available WSC river basins in Canada. 

@author: Andre R. Erler, GPL v3
'''

from collections import OrderedDict
from datasets.common import addLoadFcts
from datasets.WSC import BasinSet, Basin, LakeSet, Lake, Nation, Province, Catchment

# expand province names
province_names = OrderedDict() # make sure they are sorted alphabetically
province_names['AB'] = 'Alberta'
province_names['BC'] = 'British Columbia'
province_names['MB'] = 'Manitoba'
province_names['NB'] = 'New Brunswick'
province_names['NL'] = 'Newfoundland and Labrador'
province_names['NS'] = 'Nova Scotia'
province_names['NT'] = 'Northwest Territories'
province_names['NU'] = 'Nunavut'
province_names['PE'] = 'Prince Edward Island'
province_names['ON'] = 'Ontario'
province_names['QC'] = 'Quebec'
province_names['SK'] = 'Saskatchewan'
province_names['YT'] = 'Yukon Territory'
province_names['CAN'] = 'Canada'
# generate province info objects
province_list = OrderedDict() # just provinces
provinces = OrderedDict() # include nation for shape averages
for key,val in province_names.items():
  shapetype = 'NAT' if key == 'CAN' else 'PRV' 
  prov = Province(name=key, long_name=val, shapefile=val, data_source='?', shapetype=shapetype)
  province_list[key] = prov
  provinces[key] = prov
# the whole nation
provinces['CAN'] = Nation(name='CAN', long_name=province_names['CAN'], 
                              provinces=list(province_list.values()), data_source='?')


# dictionary with basin meta data
basin_list = OrderedDict() # maintain order
# meta data for specific basins
basin_list['AY']  = BasinSet(name='AY', long_name='Alaska and Yukon', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeAY'])
basin_list['AO']  = BasinSet(name='AO', long_name='Arctic Ocean', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeAO'])
basin_list['ARB'] = BasinSet(name='ARB', long_name='Athabasca River Basin', rivers=['Athabasca'], data_source='WSC',
                             stations=dict(Athabasca=['EmbarrasAirport','FortMcMurray','Athabasca','Windfall','Hinton','Jasper']),
                             subbasins=['WholeARB','UpperARB','LowerARB'])
basin_list['CRB'] = BasinSet(name='CRB', long_name='Columbia River Basin', rivers=['Columbia'], data_source='WSC',
                             stations=dict(), subbasins=['WholeCRB'])
basin_list['MW']  = BasinSet(name='MW', long_name='Mica Watershed', rivers=['Mica'], data_source='UNBC',
                             stations=dict(), subbasins=['WholeMW'])
basin_list['FRB'] = BasinSet(name='FRB', long_name='Fraser River Basin', rivers=['Fraser'], data_source='WSC',
                             stations=dict(Fraser=['PortMann','Mission']),
                             subbasins=['WholeFRB','UpperFRB','LowerFRB'])
basin_list['GLB'] = BasinSet(name='GLB', long_name='Great Lakes Basin', rivers=['St. Lawrence'], data_source='WSC',
                             stations={'St. Lawrence':['Cornwall']}, subbasins=['WholeGLB','LandGLB'])
basin_list['GRW'] = BasinSet(name='GRW', long_name='Grand River Watershed', data_source='Aquanty',
                             rivers=['Grand River', 'Nith River', 'Conestogo River', 'Fairchild Creek', 'Speed River', 'Whitemans Creek'], 
                             stations={'Grand River':['Brantford','Marsville'], 'Conestogo River':['Glen Allan'],
                                       'Fairchild Creek':['Brantford'],'Speed River':['Guelph'],
                                       'Whitemans Creek':['Mount Vernon'], 'Nith River':['Canning']}, 
                             subbasins=['WholeGRW','UpperGRW','LowerGRW','NorthernGRW','SouthernGRW','WesternGRW'])
basin_list['PRW'] = BasinSet(name='PRW', long_name='Payne River Watershed', data_source='Aquanty',
                             rivers=['Payne River', ], stations={'Payne River':['Berwick',],}, subbasins=['WholePRW',])
basin_list['SNW'] = BasinSet(name='SNW', long_name='South Nation Watershed', rivers=['South Nation'], data_source='Aquanty',
                             stations={'South Nation':['',]}, subbasins=['WholeSNW'])
basin_list['SON'] = BasinSet(name='SON', long_name='Southern Ontario', rivers=['Grand River'], data_source='Aquanty',
                             stations={'Grand River':['Brantford',]}, subbasins=['WholeSON','WholeGRW'])
basin_list['GSL'] = BasinSet(name='GSL', long_name='Great Slave Lake', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeGSL'])
basin_list['LS']  = BasinSet(name='LS', long_name='Labrador Sea', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeLS'])
basin_list['MKB'] = BasinSet(name='MKB', long_name='MacKenzie Basin', rivers=['MacKenzie'], data_source='',
                             stations=dict(), subbasins=['WholeMKB'])
basin_list['MRB'] = BasinSet(name='MRB', long_name='Missouri River Basin', rivers=['Missouri'], data_source='WSC',
                             stations=dict(), subbasins=['WholeMRB'])
basin_list['NRB'] = BasinSet(name='NRB', long_name='Nelson River Basin', rivers=['Nelson'], data_source='WSC',
                             stations=dict(), subbasins=['WholeNRB'])
basin_list['NHB'] = BasinSet(name='NHB', long_name='Northern Hudson Bay', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeNHB'])
basin_list['NO']  = BasinSet(name='NO', long_name='Northern Ontario', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholeNO'])
basin_list['PO']  = BasinSet(name='PO', long_name='Pacific Ocean', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholePO'])
basin_list['PSB'] = BasinSet(name='PSB', long_name='Pacific Seaboard', rivers=[], data_source='WSC',
                             stations=dict(), subbasins=['WholePSB','NorthernPSB','SouthernPSB'])
basin_list['SLR'] = BasinSet(name='SLR', long_name='Saint Lawrence River', rivers=['Saint Lawrence'], data_source='WSC',
                             stations=dict(), subbasins=['WholeSLR'])
basin_list['SSR'] = BasinSet(name='SSR', long_name='South Saskatchewan River', rivers=['South Saskatchewan'], data_source='Aquanty',
                             stations={'South Saskatchewan':['Saskatoon',]}, subbasins=['WholeSSR'])
                             # N.B.: there seem to be some issues with the St. Louis station in the SSRB - the flow is consistently 
                             #       less than Saskatoon, eventhough the drainage area is larger...

basin_sets = basin_list.copy() # dict that only contains basin sets
# dictionary of basins
basins = OrderedDict() # just the basin, and every basin only once (for shape averaging)
for name,basin in list(basin_list.items()):
  # add plain Basin instance of main basin under proper name
  basins[name] = basin # main basin under basin_list key (i.e. no Whole*)
  basin_list[basin.long_name] = basin # also make available by long_name
  # add all subbasins (including main basin with 'Whole' prefix)
  for subname,subbasin in basin.subbasins.items():
    if subname != basin.outline: # don't list outline/main basin twice
      basins[subname] = subbasin # list with all Basin instances
      basin_list[subname] = subbasin # list with all basins, BasinSet and Basin instances

# N.B.: to add new gage stations add the name to the stations dict and download the CSV files for monthly values and meta data
#       from the WSC historical archive (no missing days): http://wateroffice.ec.gc.ca/search/historical_e.html
#       all shapefiles for basins were obtained from the Water Survey of Canada; shapefiles for the lakes were optained from
#       the Natural Earth website (http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_lakes.zip)


# dictionary with lake meta data
great_lakes = OrderedDict() # maintain order
# meta data for specific lakes
great_lakes['Ontario']   = Lake(name='LakeOntario', long_name='Lake Ontario', data_source='Natural Earth',)
great_lakes['Huron']     = Lake(name='LakeHuron', long_name='Lake Huron', data_source='Natural Earth',)
great_lakes['Erie']      = Lake(name='LakeErie', long_name='Lake Erie', data_source='Natural Earth',)
great_lakes['Michigan']  = Lake(name='LakeMichigan', long_name='Lake Michigan', data_source='Natural Earth',)
great_lakes['Superior']  = Lake(name='LakeSuperior', long_name='Lake Superior', data_source='Natural Earth',)
# the Great Lakes in groups
great_lakes['NorthernLakes'] = LakeSet(name='NorthernLakes', long_name='Northern Great Lakes', data_source='Natural Earth',
                                       lakes=['NorthernLakes']+[great_lakes[lake] for lake in ['Superior','Michigan','Huron']])
great_lakes['GreatLakes'] = LakeSet(name='GreatLakes', long_name='Great Lakes', data_source='Natural Earth',
                                    lakes=['GreatLakes',]+list(great_lakes.values()))
# N.B.: individual lakes have to be initialized seperately and added to the set post-hoc, 
#       because they don't share a folder

# add aliases
lake_list = OrderedDict()
for name,lake in list(great_lakes.items()): # do not iter, since dict is changes in for-loop
  lake_list[name] = lake
  lake_list[lake.name] = lake
  lake_list[lake.long_name] = lake


## a function to create a Catchment object for a WSC gauge stations based on WSC IDs
def getCatchmentList(WSC_ID):
    ''' create and return a Catchment object based on WSC gauge station ID '''
    return Catchment(name=WSC_ID)
gauges = getCatchmentList # alias for shpavg naming
  

## generate loadGageStation* versions with these basins
from datasets.WSC import loadGageStation, loadWSC_Shp, loadWSC_ShpTS, loadWSC_Stn, loadWSC_StnTS
addLoadFcts(locals(), locals(), basins=basins, basin_list=basin_list)


## some testing/verification
if __name__ == '__main__':
    
  ## print basins
  for name,basin in basin_list.items():
    s = '  {:3s} ({:s}): '.format(name,basin.shapetype)
    if hasattr(basin, 'subbasins'):
      for subbasin in basin.subbasins: s += ' {:9s}'.format('{:s},'.format(subbasin))
    print(s)
#     if basin.maingage: print(basin.maingage.monthly_file)
    
#   ## print Great Lakes
#   for name,lake in great_lakes.iteritems():
#     s = '  {:3s} ({:s}): '.format(name,lake.shapetype)
#     if hasattr(lake, 'lakes'):
#       for sublake in lake.lakes: s += ' {:9s}'.format('{:s},'.format(sublake))
#     print(s)
#     #print(lake.folder)
#     
#   ## print Canadian Provinces
#   for name,prov in provinces.iteritems():
#     s = '  {:3s} ({:s}): '.format(name,prov.shapetype)
#     if hasattr(prov, 'provinces'):
#       for province in prov.provinces: s += ' {:3s}'.format('{:s},'.format(province))
#     print(s)
#     #print(prov.folder)
