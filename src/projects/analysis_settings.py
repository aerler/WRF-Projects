'''
Created on Jun 4, 2016

A module that provides various project specific definitions, mainly related to experiments and plotting. 

@author:  Andre R. Erler, GPL v3
'''

# external imports
import matplotlib as mpl
# internal imports
import clim.load as clim_load
import clim.plot as clim_plot
import eva.load as eva_load
import eva.plot as eva_plot
from geodata.misc import AttrDict
from plotting import figure
# from WRF_experiments import WRF_ens, WRF_exps
# from projects.CESM_experiments import CESM_ens, CESM_exps
# import projects.WSC_basins as wsc_basins
# N.B.: importing WRF_experments here causes a circular import, because WRF_experiments loads the projects.x.__init__ modules
#       which themselves import this module via the local analysis_settings (importing without __init__ is not possible)


# default shape type
default_shapetype = 'shpavg'
# default station type
default_stationtype = 'ecprecip'

# method to add defaults to specifics
def mergeAnnotation(specifics, defaults):
  annotation = AttrDict()
  for name,specs in specifics.iteritems():
    atts = defaults.copy()
    atts.update(specs)
    annotation[name] = atts 
  return annotation

# variable collections
wetday_extensions = clim_load.wetday_extensions[:3]
variables_rc = dict(); VL = clim_load.VL
# mostly for hydrological analysis
variables_rc['temp']            = VL(vars=('T2', 'Tmax', 'Tmin'), files=('srfc','xtrm',), label='2m Temperature')
variables_rc['temp_mean']       = VL(vars=('T2',), files=('srfc',), label='2m Temperature')
variables_rc['tmpx']            = VL(vars=('Tmax', 'Tmin'), files=('xtrm',), label='2m Temperature')
# variables_rc['temp']          = VL(vars=('T2',), files=('srfc',), label='2m Temperature')
variables_rc['precip_obs']      = VL(vars=('precip', 'solprec', 'wetfrq_010'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm_obs'] = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'wetprec_010'), files=('hydro',), label='Precipitation')
variables_rc['precip']          = VL(vars=('precip', 'preccu', 'solprec'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm']     = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'MaxPreccu_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_cesm']     = VL(vars=('MaxPrecip_1d', 'MaxPreccu_1d', ), files=('hydro',), label='Precipitation')
variables_rc['precip_alt']      = VL(vars=('MaxPrecip_1d', 'MaxPrecnc_1d', 'MaxSolprec_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_types']    = VL(vars=('precip','preccu','precnc'), files=('hydro',), label='Water Flux')
variables_rc['precip_conv']     = VL(vars=('precip','preccu','precnc'), files=('hydro',), label='Water Flux')
variables_rc['precip_net']      = VL(vars=('precip','solprec','p-et'), files=('srfc',), label='Water Flux')
variables_rc['precip_snow']     = VL(vars=('precip','snwmlt','solprec'), files=('hydro',), label='Water Flux')
variables_rc['snow']            = VL(vars=('snow',), files=('hydro',), label='Water Equivalent')
variables_rc['flux_snow']       = VL(vars=('precip','snwmlt','solprec'), files=('hydro',), label='Water Flux')
variables_rc['flux_days']       = VL(vars=('wetfrq_010','snwmlt','p-et'), files=('hydro',), label='Water Flux')
variables_rc['wetprec']         = VL(vars=['wetprec'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Precip.')
variables_rc['dryprec']         = VL(vars=['precip','dryprec_002','dryprec_010'], files=('hydro',), label='Wet-day Precip.')
variables_rc['wetdays']         = VL(vars=['wetfrq'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Ratio')
variables_rc['CWD']             = VL(vars=['CWD'+ext for ext in wetday_extensions]+['CNWD'], files=('hydro',), label='Continuous Wet-days')
variables_rc['CDD']             = VL(vars=['CDD'+ext for ext in wetday_extensions[:-1]]+['CNDD'], files=('hydro',), label='Continuous Dry-days')
variables_rc['sfcflx']          = VL(vars=('p-et','snwmlt','waterflx',), files=('hydro',), label='Surface Flux')
variables_rc['runoff']          = VL(vars=('runoff','sfroff','ugroff'), files=('lsm',), label='Runoff')
variables_rc['runoff_wflx']     = VL(vars=('waterflx','sfroff','runoff'), files=('lsm','hydro'), label='Runoff')
variables_rc['runoff_flux']     = VL(vars=('runoff','snwmlt','p-et'), files=('lsm','hydro'), label='Water Flux')
variables_rc['runoff_snow']     = VL(vars=('runoff','snwmlt','waterflx'), files=('lsm','hydro'), label='Water Flux')
variables_rc['hgs_forcing_wrf'] = VL(vars=('liqwatflx','pet_wrf',), files=('aux',), label='HGS Forcing')
variables_rc['hgs_precip_wrf']  = VL(vars=('liqprec','snwmlt','pet_wrf'), files=('aux',), label='HGS Forcing')
variables_rc['hgs_snow_wrf']    = VL(vars=('solprec','liqprec','pet_wrf'), files=('aux',), label='HGS Forcing')
variables_rc['hgs_forcing']     = VL(vars=('liqwatflx','pet',), files=('aux',), label='HGS Forcing')
variables_rc['hgs_precip']      = VL(vars=('liqprec','snwmlt','pet'), files=('aux',), label='HGS Forcing')
variables_rc['hgs_snow']        = VL(vars=('solprec','liqprec','pet'), files=('aux',), label='HGS Forcing')
variables_rc['heat']            = VL(vars=('hfx','lhfx','rSM'),files=('srfc','lsm'), label='Energy Flux')
variables_rc['evap']            = VL(vars=('p-et','evap','pet',), files=('hydro',), label='Water Flux')
variables_rc['spei']            = VL(vars=('precip','evap','pet',), files=('aux','hydro',), label='Water Flux')
variables_rc['wrfpet']          = VL(vars=('pet','pet_wrf','evap',), files=('aux','hydro',), label='Water Flux')
variables_rc['pet']             = VL(vars=('pet','petrad','petwnd'), files=('aux',), label='Water Flux')
variables_rc['rad']             = VL(vars=('SWDNB','netrad',), files=('aux','rad'), label='Radiative Flux')
variables_rc['vap']             = VL(vars=('Q2','vapdef',), files=('srfc','aux'), label='Vapor Pressure')
# variables_rc['rad']             = VL(vars=('SWD','netrad',), files=('aux','srfc'), label='Radiative Flux')
variables_rc['Q2']              = VL(vars=('Q2',),files=('srfc',), label='2m Humidity')
variables_rc['aSM']             = VL(vars=('aSM',),files=('lsm',), label='Soil Moisture') 
variables_rc['rSM']             = VL(vars=('rSM',),files=('lsm',), label='Relative Soil Moisture')
variables_rc['annual_flux']     = VL(vars=('waterflx',),files=('hydro',), label='Annual Net Water Forcing')
# N.B.: Noah-MP does not have relative soil moisture and there is not much difference anyway
# mostly for extreme value analysis
variables_rc['hydro_eva']    = VL(vars=['precip', 'p-et', 'snwmlt', 'waterflx','sfroff','runoff','aSM'],
                                  files=['hydro','lsm'], label='')
variables_rc['precip_eva']   = VL(vars=['precip', 'MaxPrecip_1h', 'MaxPrecip_6h', 'MaxPrecip_1d', 'MaxPrecip_5d',
                                        'MaxWaterflx_5d', 'wetfrq', 'CDD', 'CWD'], files=['hydro','srfc'], label='Precipitation')
variables_rc['temp_eva']     = VL(vars=['T2', 'MaxT2', 'MaxTmax', 'MinTmin', 'MaxTmax_7d', 'MinTmin_7d', 
                                        'frzfrq', 'CFD', 'CDFD', 'CNFD'], files=['srfc','xtrm'], label='Temperature')                          
variables_rc['precip_short'] = VL(vars=['MaxPreccu_1h', 'MaxPrecnc_1h', 'MaxPrecip_6h', 'MaxPreccu_6h'], 
                                  files=['xtrm','srfc'], label='MaxPrecip')
variables_rc['precip_long']  = VL(vars=['MaxPrecip_1d', 'MaxPreccu_1d', 'MaxPrecip_5d',], 
                                  files=['hydro'], label='MaxPrecip')
variables_rc['precip_CDD']   = VL(vars=['CNDD','CNWD']+[var+threshold for threshold in wetday_extensions for var in 'CDD','CWD'], 
                                  files=['hydro'], label='CDD')


# station selection criteria
constraints_rc = dict()
constraints_rc['min_len'] = 15 # for valid climatology
constraints_rc['lat'] = (45,55) 
constraints_rc['max_zerr'] = 100 # can't use this, because we are loading EC data separately from WRF
constraints_rc['end_after'] = 1990
                        
# dataset collections
exps_rc = dict(); EX = clim_load.EX
exps_rc['obs']       = EX(name='obs',exps=['CRU','WSC'], styles=['-','-.'], master='CRU', title='Observations')
exps_rc['cesm-obs']  = EX(name='cesm-all', exps=['Observations', 'Ens','Ens-2050','Ens-2100'],  
                         master='Ens', reference='Observations', target='Ens', title='CESM Validation & Projection')
exps_rc['cesm-mal']  = EX(name='cesm-ensa', exps=['Observations', 'MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='Observations', target='MEns', title='CESM Validation & Projection')
exps_rc['cesm-ens']  = EX(name='cesm-ens', exps=['MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='MEns', target='MEns', title='CESM Projection')

# set default variable atts for load functions from clim_load
def loadShapeObservations(variable_list=variables_rc, shapetype=default_shapetype, basin_list=None, **kwargs):
  ''' wrapper for clim.load.loadShapeObservations that sets variable lists '''
  return clim_load.loadShapeObservations(variable_list=variable_list, shapetype=shapetype, basin_list=basin_list, **kwargs)
def loadShapeEnsemble(variable_list=variables_rc, shapetype=default_shapetype, basin_list=None, **kwargs):
  ''' wrapper for clim.load.loadShapeEnsemble that sets experiment and variable lists '''
  return clim_load.loadShapeEnsemble(variable_list=variable_list, shapetype=shapetype, WRF_exps=None, 
                                     CESM_exps=None, WRF_ens=None, CESM_ens=None, basin_list=None, **kwargs)
def loadStationEnsemble(variable_list=variables_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for clim.load.loadStationEnsemble that sets experiment and variable lists '''
  return clim_load.loadStationEnsemble(variable_list=variable_list, stationtype=stationtype, WRF_exps=None, 
                                       CESM_exps=None, WRF_ens=None, CESM_ens=None, basin_list=None, **kwargs)
def loadShapeFit(variable_list=variables_rc, shapetype=default_shapetype, basin_list=None, **kwargs):
  ''' wrapper for eva.load.loadShapeEnsemble that sets experiment and variable lists '''
  return eva_load.loadShapeFit(variable_list=variable_list, shapetype=shapetype, WRF_exps=None, 
                               CESM_exps=None, WRF_ens=None, CESM_ens=None, basin_list=None, **kwargs)
def loadStationFit(variable_list=variables_rc, default_constraints=constraints_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for eva.load.loadStationEnsemble that sets experiment and variable lists etc. '''
  return eva_load.loadStationFit(variable_list=variable_list, default_constraints=default_constraints, stationtype=stationtype,
                                 WRF_exps=None, CESM_exps=None, WRF_ens=None, CESM_ens=None, **kwargs)
# N.B.: for a specific project the WRF/CESM_exps/ens variables and basin_list should be hardcoded with the appropriate lists


## settings for plotting 

# plot labels: translate internal names to something more presentable
plot_labels_rc = dict()
# datasets
plot_labels_rc['Unity']           = 'Uni. Obs.'
plot_labels_rc['Observations']    = 'EC Obs.'
plot_labels_rc['EC_1935']         = 'EC (1935)'
plot_labels_rc['EC_1965']         = 'EC (1965)'
plot_labels_rc['EC_1995']         = 'EC (1995)'
plot_labels_rc['Ens']             = 'CESM'  
plot_labels_rc['Ens-2050']        = 'CESM 2050' 
plot_labels_rc['Ens-2100']        = 'CESM 2100' 
plot_labels_rc['MEns']            = 'CESM*'  
plot_labels_rc['MEns-2050']       = 'CESM* 2050' 
plot_labels_rc['MEns-2100']       = 'CESM* 2100' 
# variables
plot_labels_rc['MaxPrecip_1d']    = 'Max Precip. (1d)'
plot_labels_rc['MaxPrecip_5d']    = 'Max Precip. (5d)'
plot_labels_rc['MaxPreccu_1d']    = 'Max Conv. (1d)'
plot_labels_rc['MaxSolprec_1d']   = 'Max Snow (1d)'
plot_labels_rc['MaxWaterFlx_1d']  = 'Max Flux (1d)'
plot_labels_rc['wetfrq_010']      = 'Wet-days'
plot_labels_rc['wetprec_010']     = 'Precip. Intensity'
plot_labels_rc['T2']              = 'T (2m)'   
plot_labels_rc['precip']          = 'Precip.'  
plot_labels_rc['liqprec']         = 'Liquid' 
plot_labels_rc['solprec']         = 'Snow' 
plot_labels_rc['dryprec']         = 'dryprec' 
plot_labels_rc['preccu']          = 'Conv.'  
plot_labels_rc['precnc']          = 'NC'  
plot_labels_rc['evap']            = 'ET'    
plot_labels_rc['p-et']            = 'Net Precip.'    
plot_labels_rc['pet']             = 'PET' 
plot_labels_rc['pet_wrf']         = 'PET (WRF)' 
plot_labels_rc['waterflx']        = 'Water Flux'
plot_labels_rc['liqwatflx']       = 'Water Forcing'
plot_labels_rc['snwmlt']          = 'Snow Melt' 
plot_labels_rc['runoff']          = 'Total Runoff' 
plot_labels_rc['ugroff']          = 'Undergr. R\'off' 
plot_labels_rc['sfroff']          = 'Surface Runoff' 
plot_labels_rc['Tmax']            = 'T (max)'    
plot_labels_rc['Tmin']            = 'T (min)'    
plot_labels_rc['hfx']             = 'Sens. Heat'     
plot_labels_rc['lhfx']            = 'Latent Heat'    
plot_labels_rc['Q2']              = 'Q (2m)'      
plot_labels_rc['aSM']             = 'aSM'     
plot_labels_rc['rSM']             = 'Soil Moist.'
plot_labels_rc['snow']             = 'SWE'
# HGS variables
plot_labels_rc['infil']           = 'Infiltration'
plot_labels_rc['exfil']           = 'Exfiltration'
plot_labels_rc['delta_storage']   = 'Storage'
plot_labels_rc['tot_et']          = 'Actual ET'
plot_labels_rc['tot_pet']         = 'PET'
plot_labels_rc['can_et']          = 'Canopy ET'
plot_labels_rc['tot_precip']      = 'Water Flx.'
plot_labels_rc['discharge']       = 'Discharge'
plot_labels_rc['outflow']         = 'Outflow'
plot_labels_rc['sfroff']          = 'Surface Runoff'
plot_labels_rc['totroff']         = 'Total Runoff'
plot_labels_rc['baseflow']        = 'Baseflow'
plot_labels_rc['recharge']        = 'Recharge'
plot_labels_rc['recharge_gwt']    = 'Recharge (GWT-based)'
plot_labels_rc['recharge_et']     = 'Recharge (ET-based)'
plot_labels_rc['d_tot']           = 'Total Storage'
plot_labels_rc['d_olf']           = 'Surface'
plot_labels_rc['d_pm']            = 'Subsurface'
plot_labels_rc['zs']              = 'Surface Elevation'
plot_labels_rc['zs_elm']          = 'Elevation (Elemental)'
plot_labels_rc['z_gw']            = 'Groundwater Table'


## plot styles for climatology plots
climds_plotargs_rc = dict()
# observational datasets
obs_args = AttrDict(marker='o', linestyle=' ', markeredgecolor='k', markeredgewidth=0.5) # 5*mpl.rcParams['lines.linewidth']
climds_plotargs_rc['Observations'] =  obs_args
climds_plotargs_rc['WSC']          =  obs_args
climds_plotargs_rc['Unity']        =  obs_args
climds_plotargs_rc['CRU']          =  obs_args
climds_plotargs_rc['GPCC']         =  obs_args
# reanalysis datasets
rea_args = AttrDict(marker='^', markersize=mpl.rcParams['lines.markersize'], linestyle=' ') # 'small'
climds_plotargs_rc['CFSR']         =  rea_args
climds_plotargs_rc['NARR']         =  rea_args
# variable settings
variable_plotargs_rc = dict()
variable_plotargs_rc['MaxPrecip_1d']   = AttrDict(color = 'green')
variable_plotargs_rc['MaxPrecip_5d']   = AttrDict(color = 'sienna')
variable_plotargs_rc['MaxPreccu_1d']   = AttrDict(color = 'magenta')
variable_plotargs_rc['MaxPrecnc_1d']   = AttrDict(color = 'grey')
variable_plotargs_rc['MaxSolprec_1d']  = AttrDict(color = 'blue')
variable_plotargs_rc['MaxWaterFlx_1d'] = AttrDict(color = 'blue')
variable_plotargs_rc['T2']             = AttrDict(color = 'green')
variable_plotargs_rc['precip']         = AttrDict(color = 'green')
variable_plotargs_rc['liqprec']        = AttrDict(color = 'green')
variable_plotargs_rc['solprec']        = AttrDict(color = 'dodgerblue') # 'blue'
variable_plotargs_rc['dryprec']        = AttrDict(color = 'green')
variable_plotargs_rc['preccu']         = AttrDict(color = 'magenta')
variable_plotargs_rc['precnc']         = AttrDict(color = 'coral')
variable_plotargs_rc['evap']           = AttrDict(color = 'red')
variable_plotargs_rc['p-et']           = AttrDict(color = 'red')
variable_plotargs_rc['pet']            = AttrDict(color = 'purple')
variable_plotargs_rc['pet_wrf']        = AttrDict(color = 'purple')
variable_plotargs_rc['waterflx']       = AttrDict(color = 'dodgerblue')
variable_plotargs_rc['liqwatflx']      = AttrDict(color = 'blue')
variable_plotargs_rc['snwmlt']         = AttrDict(color = 'orange')
variable_plotargs_rc['runoff']         = AttrDict(color = 'purple')
variable_plotargs_rc['ugroff']         = AttrDict(color = 'coral')
variable_plotargs_rc['sfroff']         = AttrDict(color = 'green')
variable_plotargs_rc['Tmax']           = AttrDict(color = 'red')
variable_plotargs_rc['Tmin']           = AttrDict(color = 'blue')
variable_plotargs_rc['hfx']            = AttrDict(color = 'red')
variable_plotargs_rc['lhfx']           = AttrDict(color = 'purple')
variable_plotargs_rc['Q2']             = AttrDict(color = 'blue')
variable_plotargs_rc['aSM']            = AttrDict(color = 'coral')
variable_plotargs_rc['rSM']            = AttrDict(color = 'green')
variable_plotargs_rc['CNWD']           = AttrDict(color = 'green')
variable_plotargs_rc['CNDD']           = AttrDict(color = 'green')
variable_plotargs_rc['petrad']         = AttrDict(color = 'red')
variable_plotargs_rc['petwnd']         = AttrDict(color = 'dodgerblue')
variable_plotargs_rc['vapdef']         = AttrDict(color = 'magenta')
variable_plotargs_rc['netrad']         = AttrDict(color = 'crimson')
variable_plotargs_rc['SWDNB']          = AttrDict(color = 'orange')
variable_plotargs_rc['snow']           = AttrDict(color = 'blue')
# HGS variables
# variable_plotargs_rc['baseflow']       = AttrDict(color = '#E24B34')
# variable_plotargs_rc['exfil']          = AttrDict(color = '#AAA2D8')
# variable_plotargs_rc['outflow']        = AttrDict(color = '#62A1C6')
variable_plotargs_rc['baseflow']       = AttrDict(color = 'red')
variable_plotargs_rc['recharge']       = AttrDict(color = 'purple')
variable_plotargs_rc['recharge_gwt']   = AttrDict(color = 'purple')
variable_plotargs_rc['recharge_et']    = AttrDict(color = 'purple')
variable_plotargs_rc['infil']          = AttrDict(color = 'green')
variable_plotargs_rc['exfil']          = AttrDict(color = 'purple')
variable_plotargs_rc['outflow']        = AttrDict(color = 'red')
variable_plotargs_rc['discharge']      = AttrDict(color = 'blue')
variable_plotargs_rc['delta_storage']  = AttrDict(color = 'gold')
variable_plotargs_rc['tot_et']         = AttrDict(color = 'red')
variable_plotargs_rc['tot_pet']        = AttrDict(color = 'purple')
variable_plotargs_rc['can_et']         = AttrDict(color = 'magenta')
variable_plotargs_rc['tot_precip']     = AttrDict(color = 'blue')
variable_plotargs_rc['sfroff']         = AttrDict(color = 'green')
variable_plotargs_rc['totroff']        = AttrDict(color = 'blue')
variable_plotargs_rc['d_tot']          = AttrDict(color = 'red')
variable_plotargs_rc['d_olf']          = AttrDict(color = 'green')
variable_plotargs_rc['d_pm']           = AttrDict(color = 'blue')


# add wet-day threshold dependent variables    
wetday_colors = ['steelblue', 'purple', 'crimson', 'orange']   
for wdext,color in zip(clim_load.wetday_extensions,wetday_colors):
  variable_plotargs_rc['wetprec'+wdext] = AttrDict(color = color)
  variable_plotargs_rc['wetfrq' +wdext] = AttrDict(color = color)
  variable_plotargs_rc['CWD'+wdext]     = AttrDict(color = color)
  variable_plotargs_rc['CDD' +wdext]    = AttrDict(color = color)

# wrapper with custom defaults to figure creator (plotargs and label positions)
def climFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, xtop=None, yright=None, **kwargs):
  if dataset_plotargs is None: dataset_plotargs = climds_plotargs_rc 
  if variable_plotargs is None: variable_plotargs = variable_plotargs_rc
  if plot_labels is None: plot_labels = plot_labels_rc
  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,
                         plot_labels=plot_labels, xtop=xtop, yright=yright, **kwargs)

## annotation for climatology plot
# defaults
clim_defaults = AttrDict(heat=(-30,130), Q2=(0,20), aSM=(0.1,0.5), temp=(245,305), wetprec=(0,25), wetfrq=(0,100))
# specific settings
clim_specifics = dict()
# add defaults to specifics
clim_annotation = mergeAnnotation(clim_specifics, clim_defaults)

# wrapper with annotation defaults for climPlot
def climPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = annotation
  if defaults is None: defaults = defaults
  if variable_list is None: variable_list = variables_rc
  return clim_plot.climPlot(annotation=clim_annotation, defaults=clim_defaults, 
                            variable_list=variable_list, **kwargs)


## plot styles for EVA plots
# dataset settings
evads_plotargs_rc = dict()
evads_plotargs_rc['EC']              = AttrDict(color='black')
evads_plotargs_rc['EC (1935)']       = AttrDict(color='blue')evads_plotargs_rc['EC (1940)']       = AttrDict(color='blue')evads_plotargs_rc['EC (1965)']       = AttrDict(color='purple')evads_plotargs_rc['EC (1995)']       = AttrDict(color='red')evads_plotargs_rc['EC (1990)']       = AttrDict(color='red')evads_plotargs_rc['Unity']           = AttrDict(color='black', marker='.')evads_plotargs_rc['Observations']    = AttrDict(color='black')evads_plotargs_rc['WSC']             = AttrDict(color='green', marker='.')evads_plotargs_rc['CRU']             = AttrDict(color='green', marker='.')evads_plotargs_rc['GPCC']            = AttrDict(color='purple', marker='.')evads_plotargs_rc['Ens']             = AttrDict(color='crimson',)evads_plotargs_rc['Ens-2050']        = AttrDict(color='darkorchid',)evads_plotargs_rc['Ens-2100']        = AttrDict(color='royalblue',)evads_plotargs_rc['MEns']            = AttrDict(color='royalblue',)    evads_plotargs_rc['MEns-2050']       = AttrDict(color='darkorchid',)   evads_plotargs_rc['MEns-2100']       = AttrDict(color='crimson',)     

# wrapper with custom defaults to figure creator (plotargs and label positions)
def evaFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, **kwargs):
  if dataset_plotargs is None: dataset_plotargs = evads_plotargs_rc 
  if variable_plotargs is None: variable_plotargs = None
  if plot_labels is None: plot_labels = plot_labels_rc
  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs, 
                         plot_labels=plot_labels, **kwargs)
## annotation for climatology plot
# distribution plot defaultsdist_defaults = AttrDict(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.35), waterflx=(-2,12), runoff=(-2,12),                         T2=(245,305), Tmax=(250,310), Tmin=(240,300), precip=(0,20), 
                         MaxWaterflx_5d=(5,35), CWD=(0,30), CDD=(0,50), CNWD=(0,30), CNDD=(0,80),                          CWD_002=(0,30), CDD_002=(0,50), CWD_010=(0,30), CDD_010=(0,80), 
                         CWD_100=(0,15), CDD_100=(0,100), CWD_200=(0,10), CDD_200=(0,100),
                         MaxPrecip_1d=(0,120), MaxPreccu_1d=(0, 40), MaxPrecnc_1d=(0,120),                          MaxSolprec_1d=(0, 40), MaxPrecip_6h=(0,150), MaxPreccu_1h=(0,500), 
                         MaxPrecnc_1h=(0,500), MaxPrecip_5d=(0,20),)
# specific settings
dist_specifics = dict() 
# add defaults to specifics
dist_annotation = mergeAnnotation(dist_specifics, dist_defaults)

# wrapper with custom annotation defaults for distPlot
def distPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = dist_annotation
  if defaults is None: defaults = dist_defaults
  return eva_plot.distPlot(annotation=annotation, defaults=defaults, **kwargs)

# quantile plot defaultsquant_defaults = AttrDict(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.25), waterflx=(-2,12),                           T2=(245,305), precip = (10,30),  runoff=(-2,12),                          MaxPrecip_6h  = (120,260), MaxPrecip_1d = (50,200), MaxPrecip_7d = (0,30),                           MaxWaterflx_5d = (0,50), CDD=(17.5,42.5), CWD=(20,100))# specific settings
quant_specifics = dict() 
# add defaults to specifics
quant_annotation = mergeAnnotation(quant_specifics, quant_defaults)

# wrapper with custom annotation defaults for distPlot
def quantPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = quant_annotation
  if defaults is None: defaults = quant_defaults
  return eva_plot.quantPlot(annotation=annotation, defaults=defaults, **kwargs)
