'''
Created on Jun 4, 2016

A module that provides various project specific definitions, mainly related to experiments and plotting. 

@author:  Andre R. Erler, GPL v3
'''

# internal imports
import clim.load as clim_load
import clim.plot as clim_plot
import eva.load as eva_load
import eva.plot as eva_plot
from geodata.misc import AttrDict
from plotting import figure
from WRF_experiments import WRF_ens, WRF_exps
from projects.CESM_experiments import CESM_ens, CESM_exps
import projects.WSC_basins as wsc_basins
# import defaults for analysis settings
import projects.analysis_settings as default

# default shape type
default_shapetype = 'glbshp'
# default station type
default_stationtype = 'ecprecip'


# variable collections
wetday_extensions = default.wetday_extensions[:] # make a copy to prevent interference
variables_rc = default.variables_rc.copy(); VL = clim_load.VL


# station selectioncriteria
constraints_rc = dict()
constraints_rc['min_len'] = 15 # for valid climatology
constraints_rc['lat'] = (40,50) 
constraints_rc['max_zerr'] = 300 # can't use this, because we are loading EC data separately from WRF
constraints_rc['prov'] = ('ON','QC')
constraints_rc['end_after'] = 1980


# dataset collections
exps_rc = dict(); EX = clim_load.EX
exps_rc['obs']       = EX(name='obs',exps=['CRU','WSC'], styles=['-','-.'], title='Observations', 
                          master='CRU', reference='CRU', target=None)
exps_rc['erai']      = EX(name='erai',exps=['erai-g','erai-t'], styles=['--','-',], 
                          title='ERA-I G & T', master='erai-g', reference=None, target=None)
exps_rc['erai-all']  = EX(name='erai-all',exps=['erai-g','erai-t','erai-g3','erai-t3'], styles=['--','-','-.',':'], 
                          title='G & T (10km & 90km, ERA-I)', master='erai-g', reference=None, target=None)
exps_rc['erai-3']    = EX(name='erai-3',exps=['erai-g3','erai-t3'], styles=['--','-',], 
                          title='G & T (90km, ERA-I)', master='erai-g3', reference=None, target=None)
exps_rc['hist']      = EX(name='hist',exps=['g-ens','erai-g','erai-t'], styles=['-','--','-.',], 
                          title='Historical WRF', master='erai-g', reference=None, target=None)
exps_rc['val-mod']   = EX(name='val-mod', exps=['g-ens', 'g-ens_d01', 'Ens','t-ens'],
                          master='g-ens', reference=None, target=None, title='Model Comparison')
exps_rc['cesm-all']  = EX(name='cesm-all', exps=['Observations', 'Ens','Ens-2050','Ens-2100'],  
                         master='Ens', reference='Observations', target='Ens', title='CESM Validation & Projection')
exps_rc['cesm-mal']  = EX(name='cesm-ensa', exps=['Observations', 'MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='Observations', target='MEns', title='CESM Validation & Projection')
exps_rc['cesm-prj']  = EX(name='cesm-ens', exps=['MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='MEns', target='MEns', title='CESM Projection')
# initial condition ensembles (including observations)
exps_rc['g-all']     = EX(name='g-all',exps=['Observations', 'g-ens','g-ens-2050','g-ens-2100'], 
                          styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)',
                          master='g-ens', reference='Observations', target='g-ens')
exps_rc['t-all']     = EX(name='t-all',exps=['Observations', 't-ens','t-ens-2050','t-ens-2100'], 
                          styles=['-','-.','--'], title='T Ensemble Average (Hist., Mid-, End-Century)',
                          master='t-ens', reference='Observations', target='g-ens')
exps_rc['g3-all']    = EX(name='g3-all',exps=['Observations', 'g3-ens','g3-ens-2050','g3-ens-2100'], 
                         styles=['-','-.','--'], title='G Ens. Avg. (90km, Hist., Mid-, End-Century)',
                         master='g3-ens', reference='Observations', target='g-ens')
exps_rc['t3-all']    = EX(name='t3-all',exps=['Observations', 't3-ens','t3-ens-2050','t3-ens-2100'], 
                         styles=['-','-.','--'], title='T Ens. Avg. (90km, Hist., Mid-, End-Century)',
                         master='t3-ens', reference='Observations', target='g-ens')
# initial condition ensembles (projection)
exps_rc['g-prj']     = EX(name='g-prj',exps=['g-ens','g-ens-2050','g-ens-2100'], 
                          styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)',
                          master='g-ens', reference=None, target=None)
exps_rc['t-prj']     = EX(name='t-prj',exps=['t-ens','t-ens-2050','t-ens-2100'], 
                          styles=['-','-.','--'], title='T Ensemble Average (Hist., Mid-, End-Century)',
                          master='t-ens', reference=None, target=None)
exps_rc['g3-prj']    = EX(name='g3-prj',exps=['g3-ens','g3-ens-2050','g3-ens-2100'], 
                         styles=['-','-.','--'], title='G Ens. Avg. (90km, Hist., Mid-, End-Century)',
                         master='g3-ens', reference=None, target=None)
exps_rc['t3-prj']    = EX(name='t3-prj',exps=['t3-ens','t3-ens-2050','t3-ens-2100'], 
                         styles=['-','-.','--'], title='T Ens. Avg. (90km, Hist., Mid-, End-Century)',
                         master='t3-ens', reference=None, target=None)
# G ensemble members
exps_rc['g-Z']       = EX(name='g-Z',exps=['g-ctrl','g-ctrl-2050','g-ctrl-2100'], 
                          styles=['-','-.','--'], title='WRF-Z (Hist., Mid-, End-Century)',
                          master='g-ctrl', reference=None, target=None)
exps_rc['g-A']       = EX(name='g-A',exps=['g-ens-A','g-ens-A-2050','g-ens-A-2100'], 
                          styles=['-','-.','--'], title='WRF-A (Hist., Mid-, End-Century)',
                          master='g-ens-A', reference=None, target=None)
exps_rc['g-B']       = EX(name='g-B',exps=['g-ens-B','g-ens-B-2050','g-ens-B-2100'], 
                          styles=['-','-.','--'], title='WRF-B (Hist., Mid-, End-Century)',
                          master='g-ens-B', reference=None, target=None)
exps_rc['g-C']       = EX(name='g-C',exps=['g-ens-C','g-ens-C-2050','g-ens-C-2100'], 
                          styles=['-','-.','--'], title='WRF-C (Hist., Mid-, End-Century)',
                          master='g-ens-C', reference=None, target=None)
# T ensemble members
exps_rc['t-Z']       = EX(name='t-Z',exps=['t-ctrl','t-ctrl-2050','t-ctrl-2100'], 
                          styles=['-','-.','--'], title='WRF-Z (Hist., Mid-, End-Century)',
                          master='t-ctrl', reference=None, target=None)
exps_rc['t-A']       = EX(name='t-A',exps=['t-ens-A','t-ens-A-2050','t-ens-A-2100'], 
                          styles=['-','-.','--'], title='WRF-A (Hist., Mid-, End-Century)',
                          master='t-ens-A', reference=None, target=None)
exps_rc['t-B']       = EX(name='t-B',exps=['t-ens-B','t-ens-B-2050','t-ens-B-2100'], 
                          styles=['-','-.','--'], title='WRF-B (Hist., Mid-, End-Century)',
                          master='t-ens-B', reference=None, target=None)
exps_rc['t-C']       = EX(name='t-C',exps=['t-ens-C','t-ens-C-2050','t-ens-C-2100'], 
                          styles=['-','-.','--'], title='WRF-C (Hist., Mid-, End-Century)',
                          master='t-ens-C', reference=None, target=None)
# physics ensemble  (withand without Observations)
exps_rc['phys-all']  = EX(name='phys',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], 
                        styles=['-','-.','--'], title='WRF Phys. Ens. (Hist., Mid-, End-Century)',
                        master='phys-ens', reference='Observations', target='phys-ens')
exps_rc['phys-prj']  = EX(name='phys',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], 
                        styles=['-','-.','--'], title='WRF Phys. Ens. (Hist., Mid-, End-Century)',
                        master='phys-ens', reference='phys-ens', target='phys-ens')
# physics ensemble members (with Observations)
exps_rc['gg-ctrl'] = EX(name='gg-ctrl', exps=['Observations', 'g-ctrl','gg-ctrl', 'g-ctrl-2050','gg-ctrl-2050'],
                        styles=['-','-.','--'], master='Observations', reference='Observations', target=None, title='')  
exps_rc['m-ctrl']  = EX(name='m-ctrl',  exps=['Observations', 'm-ctrl', 'm-ctrl-2050'],   
                        styles=['-','-.','--'], master='Observations', reference='Observations', target=None, title='')
exps_rc['mm-ctrl'] = EX(name='mm-ctrl', exps=['Observations', 'mm-ctrl', 'mm-ctrl-2050'], 
                        styles=['-','-.','--'], master='Observations', reference='Observations', target=None, title='')
exps_rc['t-ctrl']  = EX(name='t-ctrl',  exps=['Observations', 't-ctrl', 't-ctrl-2050'],
                        styles=['-','-.','--'], master='Observations', reference='Observations', target=None, title='')
# physics ensemble members (projection only)
exps_rc['gg-ctrl'] = EX(name='gg-ctrl', exps=['g-ctrl','gg-ctrl', 'g-ctrl-2050','gg-ctrl-2050'],
                        styles=['-','-.','--'], master='g-ctrl', reference=None, target=None, title='')  
exps_rc['m-ctrl']  = EX(name='m-ctrl',  exps=['m-ctrl', 'm-ctrl-2050'],   
                        styles=['-','-.','--'], master='m-ctrl', reference=None, target=None, title='')
exps_rc['mm-ctrl'] = EX(name='mm-ctrl', exps=['mm-ctrl', 'mm-ctrl-2050'], 
                        styles=['-','-.','--'], master='mm-ctrl', reference=None, target=None, title='')
exps_rc['t-ctrl']  = EX(name='t-ctrl',  exps=['t-ctrl', 't-ctrl-2050'],
                        styles=['-','-.','--'], master='t-ctrl', reference=None, target=None, title='')
# sea-ice experiments
# exps_rc['seaice']  = EX(name='seaice',exps=['g-ens-2050','g-seaice-2050','g-ens-2100','g-seaice-2100'],
#                         styles=['-.','--','-','.--','.-'], title='WRF (Mid-, End-Century; Ens. Avg., Sea-ice)',
#                         master=['g-ens-2050','g-ens-2100'], reference=None, target=None)
# exps_rc['si25']    = EX(name='si25',exps=['g-ens','g-ens-2050','g-seaice-2050'], 
#                         styles=[':','--','-'], title='WRF (Hist., Mid-Century; Ens. Avg., Sea-ice)',
#                         master='g-ens-2050', reference=None, target=None)
# exps_rc['si21']    = EX(name='si21',exps=['g-ens','g-ens-2100','g-seaice-2100'], 
#                         styles=[':','--','-'], title='WRF (Hist., End-Century; Ens. Avg., Sea-ice)',
#                         master='g-ens-2100', reference=None, target=None)


# set default variable atts for load functions from clim_load
def loadShapeObservations(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for clim.load.loadShapeObservations that sets variable lists '''
  return clim_load.loadShapeObservations(variable_list=variable_list, shapetype=shapetype, basin_list=basin_list, **kwargs)
def loadShapeEnsemble(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for clim.load.loadShapeEnsemble that sets experiment and variable lists '''
  return clim_load.loadShapeEnsemble(variable_list=variable_list, shapetype=shapetype, WRF_exps=WRF_exps, 
                                     CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, basin_list=basin_list, **kwargs)
def loadStationEnsemble(variable_list=variables_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for clim.load.loadStationEnsemble that sets experiment and variable lists '''
  return clim_load.loadStationEnsemble(variable_list=variable_list, stationtype=stationtype, WRF_exps=WRF_exps, 
                                       CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, **kwargs)
def loadShapeFit(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for eva.load.loadShapeEnsemble that sets experiment and variable lists '''
  return eva_load.loadShapeFit(variable_list=variable_list, shapetype=shapetype, WRF_exps=WRF_exps, 
                               CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, basin_list=basin_list, **kwargs)
def loadStationFit(variable_list=variables_rc, default_constraints=constraints_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for eva.load.loadStationEnsemble that sets experiment and variable lists etc. '''
  return eva_load.loadStationFit(variable_list=variable_list, default_constraints=default_constraints, stationtype=stationtype,
                                 WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, **kwargs)


## settings for plotting 

# plot labels: translate internal names to something more presentable
plot_labels_rc = dict()
# datasets
# plot_labels_rc['Unity']           = 'Uni. Obs.'
plot_labels_rc['Unity']           = 'Observations'
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
plot_labels_rc['phys-ens']        = 'Phy Ens.' 
plot_labels_rc['phys-ens-2050']   = 'Phy 2050' 
plot_labels_rc['phys-ens-2100']   = 'Phy 2100' 
plot_labels_rc['phys-ens_d01']    = 'Phy (D1)'
plot_labels_rc['g-ens']           = 'G Ens.'  
plot_labels_rc['g-ens-2050']      = 'G 2050' 
plot_labels_rc['g-ens-2100']      = 'G 2100' 
plot_labels_rc['g3-ens']          = 'G Ens. (90km)'  
plot_labels_rc['g3-ens-2050']     = 'G 2050 (90km)' 
plot_labels_rc['g3-ens-2100']     = 'G 2100 (90km)' 
# plot_labels_rc['g-ens']          = 'WRF Ens.'  
# plot_labels_rc['g-ens-2050']     = 'WRF 2050' 
# plot_labels_rc['g-ens-2100']     = 'WRF 2100' 
# plot_labels_rc['erai-g']          = 'ERA-I G'  
# plot_labels_rc['erai-t']          = 'ERA-I T'   
plot_labels_rc['erai-g']          = 'WRF G (ERA-I)'  
plot_labels_rc['erai-t']          = 'WRF T (ERA-I)'   
plot_labels_rc['erai-g3']         = 'WRF G (90km, ERA-I)'  
plot_labels_rc['erai-t3']         = 'WRF T (90km, ERA-I)'   
plot_labels_rc['t-ens']           = 'T Ens.'
plot_labels_rc['t-ens-2050']      = 'T 2050'
plot_labels_rc['t-ens-2100']      = 'T 2100'
plot_labels_rc['t3-ens']          = 'T Ens. (90km)'
plot_labels_rc['t3-ens-2050']     = 'T 2050 (90km)'
plot_labels_rc['t3-ens-2100']     = 'T 2100 (90km)'
# variables
plot_labels_rc['MaxPrecip_1d']    = 'Max Precip. (1d)'
plot_labels_rc['MaxPrecip_5d']    = 'Max Precip. (5d)'
plot_labels_rc['MaxPreccu_1d']    = 'Max Conv. (1d)'
plot_labels_rc['MaxSolprec_1d']   = 'Max Snow (1d)'
plot_labels_rc['MaxWaterFlx_1d']  = 'Max Flux (1d)'
plot_labels_rc['wetfrq_010']      = 'Wet-days'
plot_labels_rc['wetprec_010']     = 'Precip. Intensity'
plot_labels_rc['T2']              = 'T (2m)'   
plot_labels_rc['Q2']              = 'Humidity (2m)'   
plot_labels_rc['precip']          = 'Precip.'  
plot_labels_rc['liqprec']         = 'Liquid' 
plot_labels_rc['solprec']         = 'Snow' 
plot_labels_rc['dryprec']         = 'dryprec' 
plot_labels_rc['preccu']          = 'Conv.'  
plot_labels_rc['precnc']          = 'NC'  
plot_labels_rc['evap']            = 'ET'    
plot_labels_rc['p-et']            = 'Net Precip.'    
plot_labels_rc['pet']             = 'Pot. ET' 
plot_labels_rc['pet_wrf']         = 'WRF PET' 
plot_labels_rc['petrad']          = 'Rad. Term' 
plot_labels_rc['petwnd']          = 'Wind Term' 
plot_labels_rc['vapdef']          = 'Vapor Deficit' 
plot_labels_rc['netrad']          = 'Net Rad.' 
plot_labels_rc['SWDNB']           = 'SW Rad.' 
plot_labels_rc['waterflx']        = 'Water Flux'
plot_labels_rc['snwmlt']          = 'Snow Melt' 
plot_labels_rc['runoff']          = 'Total Runoff' 
plot_labels_rc['ugroff']          = 'Undergr. R\'off' 
plot_labels_rc['sfroff']          = 'Surface Runoff' 
plot_labels_rc['Tmax']            = 'T (max)'    
plot_labels_rc['Tmin']            = 'T (min)'    
plot_labels_rc['hfx']             = 'Sens. Heat'     
plot_labels_rc['lhfx']            = 'Latent Heat'    
plot_labels_rc['aSM']             = 'aSM'     
plot_labels_rc['rSM']             = 'Soil Moist.'       

## plot styles for climatology plots
climds_plotargs_rc = dict()
# observational datasets
obs_args = default.obs_args.copy()
climds_plotargs_rc['Observations'] =  obs_args
climds_plotargs_rc['Elora']        =  obs_args
climds_plotargs_rc['NRCan']        =  obs_args
climds_plotargs_rc['WSC']          =  obs_args
climds_plotargs_rc['Unity']        =  obs_args
climds_plotargs_rc['CRU']          =  obs_args
climds_plotargs_rc['GPCC']         =  obs_args
# reanalysis datasets
rea_args = default.rea_args.copy()
climds_plotargs_rc['CFSR']         =  rea_args
climds_plotargs_rc['NARR']         =  rea_args
# variable settings
variable_plotargs_rc = default.variable_plotargs_rc.copy()

# wrapper with custom defaults to figure creator (plotargs and label positions)
def climFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, xtop=None, yright=None, **kwargs):
  if dataset_plotargs is None: dataset_plotargs = climds_plotargs_rc 
  if variable_plotargs is None: variable_plotargs = variable_plotargs_rc
  if plot_labels is None: plot_labels = plot_labels_rc
  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,
                         plot_labels=plot_labels, xtop=xtop, yright=yright, **kwargs)

## annotation for climatology plot
# defaults
clim_defaults = default.clim_defaults.copy() # copy defaults to prevent interference, then update with regional defaults 
clim_defaults.update(heat=(-10,150), heat_flux=(-10,150), Q2=(0,20), aSM=(0.1,0.5), temp=(245,305), temp_mean=(250,300), wetprec=(0,25), wetfrq=(0,100), rad=(-50,350), vap=(-3.,22.), pet=(-0.5,5.), spei=(-1,6), wrfpet=(-1,6))
# specific settings
clim_specifics = dict()
clim_specifics['GLB']         = AttrDict(temp=(248,304), temp_mean=(255,300), water=(-1.,5.), precip_net=(-1.,5.5), precip_types=(-0.5,5.), precip_xtrm=(-1.,29.),
                                         runoff=(-1.,4.), rofflx=(-1.,5.), flux=(-1.,5.), flux_snow=(-1.,5.), spei=(-0.5,5.5), evap=(-1.5,4.5))
clim_specifics['GreatLakes']  = AttrDict(temp=(250,300), temp_mean=(255,300), water=(-1.,5.), precip_net=(-1.,5.5))
clim_specifics['LandGLB']     = clim_specifics['GreatLakes']
clim_specifics['GRW']         = AttrDict(temp=(250,308), water=(-1.,6.), precip_net=(-0.5,5.5), precip_types=(-0.5,5.5), precip_xtrm=(-1.,35.),
                                         runoff=(-0.5,4.5), rofflx=(-1.,4.), flux=(-1.,6.), flux_snow=(-1.,5.), spei=(-0.5,5.5), evap=(-1.5,5.5))
clim_specifics['NRB']         = AttrDict(temp=(248,304), water=(-1.4,2.2), precip=(-0.4,3.4), runoff=(-2.,2.), flux=(-2,5.5))
clim_specifics['Prairies']    = AttrDict(temp=(248,304), water=(-1.5,1.5), precip=(-0.5,3.5), precip_xtrm=(0,30), precip_cesm=(0,30), precip_alt=(0,30), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
# add defaults to specifics
clim_annotation = default.mergeAnnotation(clim_specifics, clim_defaults)
  
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
evads_plotargs_rc['EC (1935)']       = AttrDict(color='blue')evads_plotargs_rc['EC (1940)']       = AttrDict(color='blue')evads_plotargs_rc['EC (1965)']       = AttrDict(color='purple')evads_plotargs_rc['EC (1995)']       = AttrDict(color='red')evads_plotargs_rc['EC (1990)']       = AttrDict(color='red')evads_plotargs_rc['Unity']           = AttrDict(color='black', marker='.')evads_plotargs_rc['Observations']    = AttrDict(color='black')evads_plotargs_rc['WSC']             = AttrDict(color='green', marker='.')evads_plotargs_rc['CRU']             = AttrDict(color='green', marker='.')evads_plotargs_rc['GPCC']            = AttrDict(color='purple', marker='.')evads_plotargs_rc['Ens']             = AttrDict(color='crimson',)evads_plotargs_rc['Ens-2050']        = AttrDict(color='darkorchid',)evads_plotargs_rc['Ens-2100']        = AttrDict(color='royalblue',)evads_plotargs_rc['MEns']            = AttrDict(color='royalblue',)    evads_plotargs_rc['MEns-2050']       = AttrDict(color='darkorchid',)   evads_plotargs_rc['MEns-2100']       = AttrDict(color='crimson',)      evads_plotargs_rc['Ctrl-A']          = AttrDict(color='royalblue',)evads_plotargs_rc['Ctrl-A-2050']     = AttrDict(color='darkorchid',)evads_plotargs_rc['Ctrl-A-2100']     = AttrDict(color='crimson',)evads_plotargs_rc['max-ens']         = AttrDict(color='blue',  )evads_plotargs_rc['max-ens-2050']    = AttrDict(color='purple',)evads_plotargs_rc['max-ens-2100']    = AttrDict(color='red',   )evads_plotargs_rc['max-ens_d01']     = AttrDict(color='green',)evads_plotargs_rc['ctrl-ens']        = AttrDict(color='royalblue',  )evads_plotargs_rc['ctrl-ens-2050']   = AttrDict(color='darkorchid',)evads_plotargs_rc['ctrl-ens-2100']   = AttrDict(color='crimson',   )evads_plotargs_rc['ctrl-ens_d01']    = AttrDict(color='forestgreen',)evads_plotargs_rc['mex-ens']         = AttrDict(color='blue')evads_plotargs_rc['mex-ens-2050']    = AttrDict(color='purple')evads_plotargs_rc['mex-ens-2100']    = AttrDict(color='red')evads_plotargs_rc['mex-ens_d01']     = AttrDict(color='purple', ) evads_plotargs_rc['phys-ens']        = AttrDict(color='forestgreen',  ) evads_plotargs_rc['phys-ens-2050']   = AttrDict(color='brown',) evads_plotargs_rc['phys-ens-2100']   = AttrDict(color='coral',    ) evads_plotargs_rc['phys-ens_d01']    = AttrDict(color='sandybrown',)evads_plotargs_rc['grell-ens']       = AttrDict(color='blue',  )evads_plotargs_rc['grell-ens-2050']  = AttrDict(color='purple',)evads_plotargs_rc['grell-ens-2100']  = AttrDict(color='red',   )evads_plotargs_rc['kf-ens']          = AttrDict(color='seagreen',  ) evads_plotargs_rc['kf-ens-2050']     = AttrDict(color='brown',) evads_plotargs_rc['kf-ens-2100']     = AttrDict(color='coral',    ) evads_plotargs_rc['max-ctrl']        = AttrDict(color='gold')evads_plotargs_rc['max-ctrl-2050']   = AttrDict(color='greenyellow',)evads_plotargs_rc['max-ctrl-2100']   = AttrDict(color='magenta',   )evads_plotargs_rc['max-ctrl_d01']    = AttrDict(color='greenyellow',)evads_plotargs_rc['max-ens-A']       = AttrDict(color='orange')evads_plotargs_rc['max-ens-B']       = AttrDict(color='red')evads_plotargs_rc['max-ens-C']       = AttrDict(color='brown')evads_plotargs_rc['max-lowres_d01']  = AttrDict(color='sandybrown',)evads_plotargs_rc['max-lowres']      = AttrDict(color='sandybrown',)evads_plotargs_rc['max-kf']          = AttrDict(color='red')evads_plotargs_rc['max-nosub']       = AttrDict(color='purple')evads_plotargs_rc['max-nmp']         = AttrDict(color='green')evads_plotargs_rc['max-noflake']     = AttrDict(color='dodgerblue')evads_plotargs_rc['old-ctrl']        = AttrDict(color='green')evads_plotargs_rc['old-ctrl-2100']   = AttrDict(color='green')evads_plotargs_rc['ctrl-1']          = AttrDict(color='dodgerblue')evads_plotargs_rc['ctrl-2100']       = AttrDict(color='dodgerblue')evads_plotargs_rc['new-ctrl']        = AttrDict(color='purple')evads_plotargs_rc['new-ctrl-2100']   = AttrDict(color='magenta')evads_plotargs_rc['erai-max']        = AttrDict(color='green')evads_plotargs_rc['cfsr-max']        = AttrDict(color='yellow')evads_plotargs_rc['new-v361']        = AttrDict(color='crimson')evads_plotargs_rc['new-v361-2050']   = AttrDict(color='royalblue')evads_plotargs_rc['new-v361-2100']   = AttrDict(color='darkorchid')evads_plotargs_rc['max-seaice-2050'] = AttrDict(color='darkorchid')evads_plotargs_rc['max-seaice-2100'] = AttrDict(color='crimson')evads_plotargs_rc['max-1deg']        = AttrDict(color='green')evads_plotargs_rc['max-hilev']       = AttrDict(color='purple') # Marc's experiments for the Great Lakesevads_plotargs_rc['marc-g'] = AttrDict(color='blue')evads_plotargs_rc['marc-gg'] = AttrDict(color='darkcyan')evads_plotargs_rc['marc-m'] = AttrDict(color='red')evads_plotargs_rc['marc-mm'] = AttrDict(color='brown')evads_plotargs_rc['marc-t'] = AttrDict(color='gold')evads_plotargs_rc['marc-g-2050'] = AttrDict(color='purple')evads_plotargs_rc['marc-gg-2050'] = AttrDict(color='limegreen')evads_plotargs_rc['marc-m-2050'] = AttrDict(color='magenta')evads_plotargs_rc['marc-mm-2050'] = AttrDict(color='darkorange')evads_plotargs_rc['marc-t-2050'] = AttrDict(color='coral')# wrapper with custom defaults to figure creator (plotargs and label positions)def evaFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, **kwargs):  if dataset_plotargs is None: dataset_plotargs = evads_plotargs_rc   if variable_plotargs is None: variable_plotargs = None  if plot_labels is None: plot_labels = plot_labels_rc  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,                          plot_labels=plot_labels, **kwargs)## annotation for climatology plot
# distribution plot defaultsdist_defaults = default.dist_defaults.copy() # global defaults, followed by local defaults
dist_defaults.update(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.35), waterflx=(-2,12),                      T2=(245,305), Tmax=(250,310), Tmin=(240,300), precip=(0,20), 
                     MaxWaterflx_5d=(5,35), CWD=(0,30), CDD=(0,50), CNWD=(0,30), CNDD=(0,80),                      CWD_002=(0,30), CDD_002=(0,50), CWD_010=(0,30), CDD_010=(0,80), 
                     CWD_100=(0,15), CDD_100=(0,100), CWD_200=(0,10), CDD_200=(0,100),
                     MaxPrecip_1d=(0,120), MaxPreccu_1d=(0, 40), MaxPrecnc_1d=(0,120),                      MaxSolprec_1d=(0, 40), MaxPrecip_6h=(0,150), MaxPreccu_1h=(0,500), 
                     MaxPrecnc_1h=(0,500), MaxPrecip_5d=(0,20),)
dist_specifics = dict() 
# basins (annual)dist_specifics['GLB'] = AttrDict(MaxWaterflx_5d=(2,14), precip=(0,6), waterflx=(0,6))dist_specifics['GRW'] = AttrDict(water=(-2.,4.), runoff=(-2.,2.), flux=(-2,5.5), waterflx=(0,7))# clusters (seasonal)
dist_specifics['maritime_summer']  = AttrDict(MaxPreccu_1h=(0,600), MaxPrecnc_1h=(0,600), MaxPrecip_6h=(0,180), MaxPrecnc_6h=(0,180), MaxPreccu_6h=(0,100), MaxPrecip_1d=(0, 60), MaxSolprec_1d=(0, 60), MaxPrecnc_1d=(0, 60), MaxPreccu_1d=(0,30), MaxPrecip_5d=(0,30), MaxSolprec_5d=(0,30), MaxWaterflx_5d=(-4,10), CWD_002=(0,30), CDD_002=(0,50), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,20), CNDD=(0,120))dist_specifics['prairies_summer'] = AttrDict(MaxPreccu_1h=(0,800), MaxPrecnc_1h=(0,800), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,120), MaxPrecip_1d=(0,100), MaxSolprec_1d=(0,100), MaxPrecnc_1d=(0, 90), MaxPreccu_1d=(0,35), MaxPrecip_5d=(0,45), MaxSolprec_5d=(0,45), MaxWaterflx_5d=(-5,15), CWD_002=(0,40), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,60), CNWD=(0,20), CNDD=(0, 90))dist_specifics['maritime_winter']  = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,400), MaxPrecip_6h=(0,120), MaxPrecnc_6h=(0,120), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0, 40), MaxSolprec_1d=(0, 40), MaxPrecnc_1d=(0, 40), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,22), MaxSolprec_5d=(0,22), MaxWaterflx_5d=(-5,15), CWD_002=(0,40), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,60), CNDD=(0, 45))dist_specifics['prairies_winter'] = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,250), MaxPrecip_6h=(0, 90), MaxPrecnc_6h=(0, 90), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0, 40), MaxSolprec_1d=(0, 40), MaxPrecnc_1d=(0, 40), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,22), MaxSolprec_5d=(0,22), MaxWaterflx_5d=(-4,10), CWD_002=(0,30), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,60), CNDD=(0, 45))# provinces (seasonal)
dist_specifics['ON_summer']       = dist_specifics['maritime_summer'] 
dist_specifics['ON_winter']       = dist_specifics['maritime_winter']dist_annotation = default.mergeAnnotation(dist_specifics, dist_defaults)
# wrapper with custom annotation defaults for distPlot
def distPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = dist_annotation
  if defaults is None: defaults = dist_defaults
  return eva_plot.distPlot(annotation=annotation, defaults=defaults, **kwargs)

# quantile plot defaults
quant_defaults = default.quant_defaults.copy() # global defaults, followed by local defaults
quant_defaults.update(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.25), waterflx=(-2,12),                       T2=(245,305), precip = (10,30),                       MaxPrecip_6h  = (120,260), MaxPrecip_1d = (50,200), MaxPrecip_7d = (0,30),                       MaxWaterflx_5d = (0,50), CDD=(17.5,42.5), CWD=(20,100))quant_specifics = dict()# basins (annual)
quant_specifics['GLB'] = AttrDict(waterflx=(1.5,4.5), MaxWaterflx_5d = (5,15), aSM=(0.16,0.21))quant_specifics['GRW'] = AttrDict(waterflx=(1,4), aSM=(0.16,0.21))
# provinces (seasonal)quant_specifics['ON'] = AttrDict(MaxPrecip_1d=(10,80))
# clusters (seasonal)quant_specifics['maritime_summer']  = AttrDict(MaxPreccu_1h=(300,1200), MaxPrecnc_1h=(300,1200), MaxPrecip_6h=( 90,280), MaxPreccu_6h=( 80,220), MaxPrecip_1d=(30,120), MaxPrecnc_1d=(15, 90), MaxPreccu_1d=(15, 60), MaxPrecip_5d = (10,30), MaxWaterflx_5d = (40,120)) # plateau_summer quant_specifics['prairies_summer'] = AttrDict(MaxPreccu_1h=(400,1800), MaxPrecnc_1h=(400,1600), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(60,160), MaxPrecnc_1d=(30,120), MaxPreccu_1d=(20, 80), MaxPrecip_5d = (10,40), MaxWaterflx_5d = (40,120)) # prairies_summer quant_specifics['maritime_winter']  = AttrDict(MaxPreccu_1h=(200,1200), MaxPrecnc_1h=(100, 500), MaxPrecip_6h=( 60,180), MaxPreccu_6h=( 30,150), MaxPrecip_1d=(15, 75), MaxPrecnc_1d=(15, 75), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = ( 3,22), MaxWaterflx_5d = (40,120)) # plateau_winter quant_specifics['prairies_winter'] = AttrDict(MaxPreccu_1h=(200,1600), MaxPrecnc_1h=( 50, 300), MaxPrecip_6h=( 50,120), MaxPreccu_6h=( 20,100), MaxPrecip_1d=(15, 75), MaxPrecnc_1d=(15, 75), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = ( 3,22), MaxWaterflx_5d = (40,120)) # prairies_winter quant_annotation = default.mergeAnnotation(quant_specifics, quant_defaults)

# wrapper with custom annotation defaults for distPlot
def quantPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = quant_annotation
  if defaults is None: defaults = quant_defaults
  return eva_plot.quantPlot(annotation=annotation, defaults=defaults, **kwargs)
