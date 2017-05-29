'''
Created on 2014-02-13

A package that contains project specific data for use with the geodata package.

@author: Andre R. Erler, GPL v3
'''
from warnings import warn

## import shape dictionaries
try:
    from projects.WSC_basins import basins, provinces, great_lakes # import the dicts with unique entries
except (ImportError,IOError):
    warn("Error importing shape files and/or WSC module.")
