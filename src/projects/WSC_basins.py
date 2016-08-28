'''
Created on Aug 27, 2016

This module contains meta data for all available WSC river basins in Canada. 

@author: Andre R. Erler, GPL v3
'''

from collections import OrderedDict
from datasets.common import addLoadFcts
from datasets.WSC import BasinInfo, Basin

# dictionary with basin meta data
basins_info = OrderedDict() # maintain order
# meta data for specific basins

basins_info['AY']  = BasinInfo(name='AY', long_name='Alaska and Yukon', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeAY'])
basins_info['AO']  = BasinInfo(name='AO', long_name='Arctic Ocean', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeAO'])
basins_info['ARB'] = BasinInfo(name='ARB', long_name='Athabasca River Basin', rivers=['Athabasca'], data_source='WSC',
                               stations=dict(Athabasca=['Embarras','McMurray']),
                               subbasins=['WholeARB','UpperARB','LowerARB'])
basins_info['CRB'] = BasinInfo(name='CRB', long_name='Columbia River Basin', rivers=['Columbia'], data_source='WSC',
                               stations=dict(), subbasins=['WholeCRB'])
basins_info['FRB'] = BasinInfo(name='FRB', long_name='Fraser River Basin', rivers=['Fraser'], data_source='WSC',
                               stations=dict(Fraser=['PortMann','Mission']),
                               subbasins=['WholeFRB','UpperFRB','LowerFRB'])
basins_info['GLB'] = BasinInfo(name='GLB', long_name='Great Lakes Basin', rivers=['Upper Saint Lawrence'], data_source='WSC',
                              stations=dict(), subbasins=['WholeGLB'])
basins_info['GRW'] = BasinInfo(name='GRW', long_name='Grand River Watershed', rivers=['Grand River'], data_source='Aquanty',
                               stations={'Grand River':['Brantford']}, subbasins=['WholeGRW','UpperGRW','LowerGRW','NorthernGRW','SouthernGRW','WesternGRW'])
basins_info['GSL'] = BasinInfo(name='GSL', long_name='Great Slave Lake', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeGSL'])
basins_info['LS']  = BasinInfo(name='LS', long_name='Labrador Sea', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeLS'])
basins_info['MKB'] = BasinInfo(name='MKB', long_name='MacKenzie Basin', rivers=['MacKenzie'], data_source='',
                               stations=dict(), subbasins=['WholeMKB'])
basins_info['MRB'] = BasinInfo(name='MRB', long_name='Missouri River Basin', rivers=['Missouri'], data_source='WSC',
                               stations=dict(), subbasins=['WholeMRB'])
basins_info['NRB'] = BasinInfo(name='NRB', long_name='Nelson River Basin', rivers=['Nelson'], data_source='WSC',
                               stations=dict(), subbasins=['WholeNRB'])
basins_info['NHB'] = BasinInfo(name='NHB', long_name='Northern Hudson Bay', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeNHB'])
basins_info['NO']  = BasinInfo(name='NO', long_name='Northern Ontario', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholeNO'])
basins_info['PO']  = BasinInfo(name='PO', long_name='Pacific Ocean', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholePO'])
basins_info['PSB'] = BasinInfo(name='PSB', long_name='Pacific Seaboard', rivers=[], data_source='WSC',
                               stations=dict(), subbasins=['WholePSB','NorthernPSB','SouthernPSB'])
basins_info['SLR'] = BasinInfo(name='SLR', long_name='Saint Lawrence River', rivers=['Saint Lawrence'], data_source='WSC',
                               stations=dict(), subbasins=['WholeSLR'])
basins_info['SSR'] = BasinInfo(name='SSR', long_name='South Sasketchewan River', rivers=['South Sasketchewan River'], data_source='Aquanty',
                               stations=dict(), subbasins=['WholeSSR'])

# N.B.: to add new gage stations add the name to the statins-dict and download the CSV files for monthly values and meta data
#       from the WSC historical archive (no missing days): http://wateroffice.ec.gc.ca/search/search_e.html?sType=h2oArc

# N.B.: all shapefiles from Water Survey of Canada

# dictionary of basins
basins = OrderedDict() # maintain order
for name,basin in basins_info.iteritems():
  # add main basin
  basins[basin.name] = Basin(basin=basin, subbasin=None)
  if len(basin.subbasins) > 1 :
    # preserve grouping
    for subbasin in basin.subbasins[1:]: # skip first
      basins[subbasin] = Basin(basin=basin, subbasin=subbasin)
    
# get hydrographs from WSC here: https://wateroffice.ec.gc.ca/search/search_e.html?sType=h2oArc

## generate loadGageStation* versions with these basins
from datasets.WSC import loadGageStation
addLoadFcts(locals(), locals(), basins=basins, basins_info=basins_info)


if __name__ == '__main__':
    
  ## print basins
  for name,basin in basins_info.iteritems():
    s = '  {:3s}: '.format(name)
    for subbasin in basin.subbasins: s += ' {:9s}'.format('{:s},'.format(subbasin))
    print(s)
    #print(basin.folder)