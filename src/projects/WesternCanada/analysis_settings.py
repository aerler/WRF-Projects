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
from projects.WesternCanada.WRF_experiments import WRF_ens, WRF_exps
from projects.CESM_experiments import CESM_ens, CESM_exps
import projects.WSC_basins as wsc_basins
# import defaults for analysis settings
import projects.analysis_settings as default

# default shape type
default_shapetype = 'wcshp'
# default station type
default_stationtype = 'ecprecip'


# variable collections
wetday_extensions = default.wetday_extensions[:] # make a copy to prevent interference
variables_rc = default.variables_rc.copy(); VL = clim_load.VL


# station selection criteria
constraints_rc = default.constraints_rc.copy()
constraints_rc['min_len'] = 15 # for valid climatology
constraints_rc['lat'] = (45,55) 
constraints_rc['max_zerr'] = 100 # can't use this, because we are loading EC data separately from WRF
constraints_rc['prov'] = ('BC','AB')
constraints_rc['end_after'] = 1990
                        
                        
# dataset collections
exps_rc = default.exps_rc.copy(); EX = clim_load.EX
exps_rc['obs']       = EX(name='obs',exps=['CRU','WSC'], styles=['-','-.'], master='CRU', title='Observations')
exps_rc['erai-3km']  = EX(name='erai-3km', exps=['EC', 'erai-3km_d01', 'erai-3km_d02', 'erai-3km_d03'], # ,'erai-max','ctrl-1','max-1deg' 
                          master='erai-3km_d03', reference='EC', title='WRF 3km Validation',
                          styles=['--',':','-']) # obs are extra... only WRF
exps_rc['3km-prj']   = EX(name='3km-prj', exps=['max-3km', 'max-3km-2100'], styles=['--','-'], 
                          master='max-3km', reference='max-3km', title='WRF 3km Projections')
exps_rc['max-3km']   = EX(name='max-3km', exps=['EC','max-3km', 'max-3km-2100'], styles=['--','-'], # obs are extra...
                          master='max-3km', reference='EC', target='max-3km', title='WRF 3km Projections')
exps_rc['erai']      = EX(name='erai', exps=['EC', 'erai-max_d01', 'erai-max', 'erai-3km'], # ,'erai-max','ctrl-1','max-1deg' 
                          master='erai-3km', reference='EC', title='ERA-I Resolution')
exps_rc['val-res']   = EX(name='val-res', exps=['EC', 'max-ens_d01', 'max-ens', 'Ens'], # ,'erai-max','ctrl-1','max-1deg' 
                          master='EC', reference='EC', title='Resolution')
exps_rc['val-all']   = EX(name='val-all', exps=['EC', 'Ens', 'max-ens', 'ctrl-ens'], # ,'erai-max','ctrl-1','max-1deg' 
                          master='max-ens', reference='EC', target=None, title='Resolution')
exps_rc['val-mod']   = EX(name='val-mod', exps=['max-ens', 'max-ens_d01', 'Ens','ctrl-ens'],
                          master='max-ens', reference=None, target=None, title='Model Comparison')
exps_rc['cesm-obs']  = EX(name='cesm-all', exps=['EC', 'Ens','Ens-2050','Ens-2100'],  
                         master='Ens', reference='EC', target='Ens', title='CESM Validation & Projection')
exps_rc['cesm-mal']  = EX(name='cesm-ensa', exps=['EC', 'MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='EC', target='MEns', title='CESM Validation & Projection')
exps_rc['cesm-ens']  = EX(name='cesm-ens', exps=['MEns','MEns-2050','MEns-2100'],  
                         master='MEns', reference='MEns', target='MEns', title='CESM Projection')
exps_rc['max-val']   = EX(name='max-val', exps=['EC','max-ens',], # ,'erai-max','ctrl-1','max-1deg' 
                          master='max-ens', reference='EC', target=None, title='WRF G Val.')
exps_rc['max-prj']   = EX(name='prj', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                          master='max-ens', reference=None, target=None, title='WRF G Proj.')
exps_rc['max-all']   = EX(name='max-all', exps=['EC', 'max-ens','max-ens_d01','erai-max','max-ens-2050','max-ens-2100'], # ,'erai-max','ctrl-1','max-1deg' 
                          master='max-ens', reference='EC', target='auto', title='Validation & Projection')
exps_rc['max-obs']   = EX(name='max', exps=['EC','max-ens','max-ens-2050','max-ens-2100'], styles=['--',':','-'], 
                          master='max-ens', reference='EC', target='max-ens', title='WRF G Ensemble')
exps_rc['min-obs']   = EX(name='max', exps=['EC','max-ens','max-ens-2100'], styles=['-','--',], 
                          master='max-ens', reference='EC', target='max-ens', title='WRF G Ensemble')
exps_rc['max-1-obs'] = EX(name='max-1', exps=['EC','max-ctrl', 'max-ctrl-2100'], styles=['--', '-'], 
                          master='max-ctrl', reference='EC', target='max-ctrl', title='WRF-1 (G3)')
exps_rc['max-sum']   = EX(name='maxs', exps=['EC','max-ens','max-ens-2100'], 
                          master='max-ens', reference='EC', target='max-ens', title='IC Ensemble')
exps_rc['ctrl-obs']  = EX(name='ctrl', exps=['EC','ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], styles=['--',':','-'],
                          master='ctrl-ens', reference='EC', target='ctrl-ens', title='WRF T Ensemble')
exps_rc['ctrl-val']  = EX(name='ctrl', exps=['EC','ctrl-ens'], styles=['-'],
                          master='ctrl-ens', reference='EC', target=None, title='WRF T Val.')
exps_rc['ctrl-1-obs'] = EX(name='ctrl-1', exps=['EC','ctrl-1', 'ctrl-2100'], styles=['--', '-'], 
                           master='ctrl-1', reference='EC', target='ctrl-1', title='WRF-1 (KF)')
exps_rc['ctrl-sum']  = EX(name='ctrl-sum', exps=['EC','ctrl-ens','ctrl-ens-2100'], 
                          master='ctrl-ens', reference='EC', target='ctrl-ens', title='Alt. Ens.')
exps_rc['ctrl-prj']  = EX(name='ctrl-prj', exps=['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                          master='ctrl-ens', reference=None, target=None, title='WRF T Proj.')
exps_rc['mex-prj']   = EX(name='mex',exps=['mex-ens','mex-ens-2050','mex-ens-2100'], master='mex-ens',
                          styles=['-','-.','--'], title='WRF Ext. Ens. (Hist., Mid-, End-Century)')
exps_rc['phys-prj']  = EX(name='phys',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], master='phys-ens',
                          styles=['-','-.','--'], title='WRF Phys. Ens. (Hist., Mid-, End-Century)')
exps_rc['max-nmp']   = EX(name='nmp',exps=['max-ens','max-ens-2050','max-nmp','max-nmp-2050'], master='max-nmp',
                          styles=['-.','.','--','-'], title='IC Ensemble Average & Noah-MP (Hist., Mid-Century)')
exps_rc['max-hilev'] = EX(name='test',exps=['max-ens','max-ctrl','max-hilev'], master='max-ens',
                          styles=['-','--','-.'], title='WRF (Ens. Avg., Ctrl-1, Hi-lev)')
exps_rc['max-1deg']  = EX(name='1deg',exps=['max-ens','max-ens_d01','max-1deg'], master='max-ens',
                          styles=['--','-.','-'], title='WRF (Ens. Avg., 30km, 1deg)')
exps_rc['max-Z-prj'] = EX(name='max-Z',exps=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], master='max-ctrl',
                          styles=['-','-.','--'], title='WRF-Z (Hist., Mid-, End-Century)')
exps_rc['max-A-prj'] = EX(name='max-A',exps=['max-ens-A','max-ens-A-2050','max-ens-A-2100'], master='max-ens-A',
                          styles=['-','-.','--'], title='WRF-A (Hist., Mid-, End-Century)')
exps_rc['max-B-prj'] = EX(name='max-B',exps=['max-ens-B','max-ens-B-2050','max-ens-B-2100'], master='max-ens-B',
                          styles=['-','-.','--'], title='WRF-B (Hist., Mid-, End-Century)')
exps_rc['max-C-prj'] = EX(name='max-C',exps=['max-ens-C','max-ens-C-2050','max-ens-C-2100'], master='max-ens-C',
                          styles=['-','-.','--'], title='WRF-C (Hist., Mid-, End-Century)')
exps_rc['max-cice']  = EX(name='seaice',exps=['max-ens','max-ens-2050','max-seaice-2050','max-ens-2100','max-seaice-2100'],
                          styles=['-.','--','-','.--','.-'], master='max-ens', 
                          title='WRF (Hist., Mid-, End-Century; Ens. Avg., Sea-ice)')
exps_rc['cice-mid']  = EX(name='si25',exps=['max-ens','max-ens-2050','max-seaice-2050'], master='max-ens-2050',
                          styles=[':','--','-'], title='WRF (Hist., Mid-Century; Ens. Avg., Sea-ice)')
exps_rc['cice-end']  = EX(name='si21',exps=['max-ens','max-ens-2100','max-seaice-2100'], master='max-ens-2100',
                          styles=[':','--','-'], title='WRF (Hist., End-Century; Ens. Avg., Sea-ice)')
exps_rc['sens-res']  = EX(name='sens-res', exps=['EC', 'Ens', 'max-ctrl','max-ctrl_d01','max-1deg','max-lowres_d01'], 
                          master='max-ctrl', reference='EC', target=None, title='Sensitivity to Resolution')
exps_rc['sens-max']  = EX(name='sens-max',exps=['max-ctrl','max-nosub','max-kf','max-nmp','max-noflake'], 
                          master='max-ctrl', reference=None, target=None, title='Sensitivity Tests (Max)')
exps_rc['sens-phy']  = EX(name='sens-phys',exps=['EC','max-ctrl','old-ctrl','ctrl-1','new-ctrl','new-v361'], 
                           master='max-ctrl', reference='EC', target=None, title='Sensitivity (Phys. Ens.)')
exps_rc['sens-phy-2100'] = EX(name='sens-phys-2100',exps=['max-ctrl-2100','old-ctrl-2100','ctrl-2100','new-ctrl-2100','new-v361-2100'], 
                               master='max-ctrl-2100', reference=None, target=None, title='Sensitivity (Phys. Ens. 2100)')
exps_rc['sens-ens']  = EX(name='sens-ens',exps=['EC','max-ens', 'max-ctrl','max-ens-A','max-ens-B','max-ens-C','erai-max'], 
                          master='max-ens', reference='EC', target=None, title='Sensitivity (Phys. Ens.)')
exps_rc['sens-cu']   = EX(name='sens-cu', exps=['grell-ens','kf-ens','grell-ens-2100','kf-ens-2100'], 
                          master=['grell-ens','grell-ens-2100'], reference=None, target=None, title='Grell vs. KF')
exps_rc['kf-ens']    = EX(name='kf-ens', exps=['kf-ens','kf-ens-2050','kf-ens-2100'], 
                          master='kf-ens', reference=None, target=None, title='KF Ens.')
exps_rc['max-phys']  = EX(name='max-phys', exps=['max-ens','phys-ens','max-ens-2100','phys-ens-2100'], 
                        master=['max-ens','max-ens-2100'], reference=None, target=None, title='IC vs. Phys.')
exps_rc['max-ctrl']  = EX(name='max-ctrl', exps=['max-ens','max-ctrl','max-ens-2050','max-ctrl-2050','max-ens-2100','max-ctrl-2100'], 
                        master=['max-ens','max-ens-2050','max-ens-2100'], reference=None, target=None, title='IC & Max')
exps_rc['max-shape'] = EX(name='shape', exps=['EC','max-ens','max-ens-2050','max-ens-2100'], 
                          master='max-ens', reference='EC', target=None, title='Projected Shape Change')
exps_rc['max-shape-cu'] = EX(name='shape-cu', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                             master='max-ens', reference='max-ens', target=None, title='Projected Shape Change')
exps_rc['new-val']       = EX(name='new', exps=['max-ctrl','new-v361','max-ctrl-2050','new-v361-2050','max-ctrl-2100','new-v361-2100'], 
                          master=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], reference=None, target=None, title='New V3.6.1')
exps_rc['mex']       = EX(name='mex', exps=['EC','mex-ens','mex-ens-2050','mex-ens-2100'], 
                          master='mex-ens', reference='EC', target='mex-ens', title='IC Ensemble (Ext.)')
exps_rc['phys-obs']      = EX(name='phys',exps=['EC','phys-ens','phys-ens-2050','phys-ens-2100'], 
                          master='phys-ens', reference='EC', target='phys-ens', title='Physics Ensemble')
exps_rc['phys-val']  = EX(name='phys-val',exps=['EC','phys-ens','phys-ens_d01'], 
                      master='phys-ens', reference='EC', target=None, title='Physics Validation')
exps_rc['phys-prj']  = EX(name='phys-prj',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], 
                          master='phys-ens', reference=None, target=None, title='Physics Projection')
exps_rc['phys-shape']  = EX(name='shape', exps=['EC','phys-ens','phys-ens-2050','phys-ens-2100'], 
                            master='phys-ens', reference='EC', target=None, title='Projected Shape Change')
exps_rc['ens-all']  = EX(name='ens-all', exps=['EC','max-ens','max-ctrl','erai-max','phys-ens','max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050','max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                        master=['max-ens','max-ens-2050','max-ens-2100'], reference='EC', target=None, title='Ensembles')
exps_rc['ens-1980']  = EX(name='ens-1980', exps=['EC','max-ens','max-ctrl','erai-max','phys-ens'], 
                        master='max-ens', reference='EC', target=None, title='Ensembles 1980')
exps_rc['ens-2050']  = EX(name='ens-2050', exps=['max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050'], 
                        master='max-ens-2050', reference=None, target=None, title='Ensembles 2050')
exps_rc['ens-2100']  = EX(name='ens-2100', exps=['max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                          master='max-ens-2100', reference=None, target=None, title='Ensembles 2100')


# set default variable atts for load functions from clim_load
def loadShapeObservations(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for clim.load.loadShapeObservations that sets variable lists '''
  return clim_load.loadShapeObservations(variable_list=variable_list, shapetype=shapetype, basin_list=basin_list, 
                                         obs_clim='Unity', obs_ts='Unity', **kwargs)
def loadShapeEnsemble(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for clim.load.loadShapeEnsemble that sets experiment and variable lists '''
  return clim_load.loadShapeEnsemble(variable_list=variable_list, shapetype=shapetype, WRF_exps=WRF_exps, 
                                     CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, basin_list=basin_list, 
                                     obs_clim='Unity', obs_ts='Unity', **kwargs)
def loadStationEnsemble(variable_list=variables_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for clim.load.loadStationEnsemble that sets experiment and variable lists '''
  return clim_load.loadStationEnsemble(variable_list=variable_list, stationtype=stationtype, WRF_exps=WRF_exps, 
                                       CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, obs_ts='EC', **kwargs)
def loadShapeFit(variable_list=variables_rc, shapetype=default_shapetype, basin_list=wsc_basins.basin_list, **kwargs):
  ''' wrapper for eva.load.loadShapeEnsemble that sets experiment and variable lists '''
  return eva_load.loadShapeFit(variable_list=variable_list, shapetype=shapetype, WRF_exps=WRF_exps, 
                               CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, 
                               basin_list=basin_list, obs_ts='Unity', **kwargs)
def loadStationFit(variable_list=variables_rc, default_constraints=constraints_rc, stationtype=default_stationtype, **kwargs):
  ''' wrapper for eva.load.loadStationEnsemble that sets experiment and variable lists etc. '''
  return eva_load.loadStationFit(variable_list=variable_list, default_constraints=default_constraints, stationtype=stationtype,
                                 WRF_exps=WRF_exps, CESM_exps=CESM_exps, WRF_ens=WRF_ens, CESM_ens=CESM_ens, 
                                 obs_ts='EC', **kwargs)


## settings for plotting 

# plot labels: translate internal names to something more presentable
plot_labels_rc = default.plot_labels_rc.copy()
# datasets
plot_labels_rc['Unity']           = 'Merged Obs.'
plot_labels_rc['Observations']    = 'EC Obs.'
plot_labels_rc['EC']              = 'EC Obs.'
plot_labels_rc['EC_1935']         = 'EC (1935)'
plot_labels_rc['EC_1965']         = 'EC (1965)'
plot_labels_rc['EC_1995']         = 'EC (1995)'
plot_labels_rc['Ens']             = 'CESM'  
plot_labels_rc['Ens-2050']        = 'CESM 2050' 
plot_labels_rc['Ens-2100']        = 'CESM 2100' 
plot_labels_rc['MEns']            = 'CESM*'  
plot_labels_rc['MEns-2050']       = 'CESM* 2050' 
plot_labels_rc['MEns-2100']       = 'CESM* 2100' 
plot_labels_rc['max-ens']         = 'IC Ens.'  
plot_labels_rc['max-ens-2050']    = 'IC 2050' 
plot_labels_rc['max-ens-2100']    = 'IC 2100' 
# plot_labels_rc['max-ens']         = 'WRF G Ens.'  
# plot_labels_rc['max-ens-2050']    = 'WRF G 2050' 
# plot_labels_rc['max-ens-2100']    = 'WRF G 2100' 
plot_labels_rc['max-ctrl']        = 'WRF-A  (G3)'  
plot_labels_rc['max-ctrl-2050']   = 'A 2050 (G3)' 
plot_labels_rc['max-ctrl-2100']   = 'A 2100 (G3)' 
plot_labels_rc['max-ens_d01']     = 'IC (D1)'
plot_labels_rc['erai-max_d01']    = 'ERAI 30km'  
plot_labels_rc['erai-max']        = 'ERAI 10km'  
plot_labels_rc['erai-3km']        = 'ERAI 3km'  
plot_labels_rc['erai-3km_d01']    = 'ERAI 30km'  
plot_labels_rc['erai-3km_d02']    = 'ERAI 10km'  
plot_labels_rc['erai-3km_d03']    = 'ERAI  3km'  
plot_labels_rc['max-3km']         = 'WRF 3km'  
plot_labels_rc['max-3km_d01']     = 'WRF 30km'  
plot_labels_rc['max-3km_d02']     = 'WRF 10km'  
plot_labels_rc['max-3km_d03']     = 'WRF  3km'  
plot_labels_rc['max-3km-2100']     = '2100 (3km)'  
plot_labels_rc['max-3km-2100_d01'] = 'WRF 30km (2100)'  
plot_labels_rc['max-3km-2100_d02'] = 'WRF 10km (2100)'  
plot_labels_rc['max-3km-2100_d03'] = 'WRF  3km (2100)'  
# plot_labels_rc['ctrl-ens']        = 'Alt. Ens.'  
# plot_labels_rc['ctrl-ens-2050']   = 'AE 2050' 
# plot_labels_rc['ctrl-ens-2100']   = 'AE 2100'
plot_labels_rc['ctrl-ens']        = 'IC Ens.'  
plot_labels_rc['ctrl-ens-2050']   = 'IC 2050' 
plot_labels_rc['ctrl-ens-2100']   = 'IC 2100'
# plot_labels_rc['ctrl-ens']        = 'WRF T Ens.'  
# plot_labels_rc['ctrl-ens-2050']   = 'WRF T 2050' 
# plot_labels_rc['ctrl-ens-2100']   = 'WRF T 2100'
plot_labels_rc['ctrl-1']          = 'WRF-A  (KF)'
plot_labels_rc['ctrl-2050']       = 'A 2050 (KF)'
plot_labels_rc['ctrl-2100']       = 'A 2100 (KF)'
plot_labels_rc['ctrl-ens_d01']    = 'AE (D1)'
plot_labels_rc['max-seaice-2050'] = 'SI 2050' 
plot_labels_rc['max-seaice-2100'] = 'SI 2100' 
plot_labels_rc['phys-ens']        = 'Phy Ens.' 
plot_labels_rc['phys-ens-2050']   = 'Phy 2050' 
plot_labels_rc['phys-ens-2100']   = 'Phy 2100' 
plot_labels_rc['phys-ens_d01']    = 'Phy (D1)'
plot_labels_rc['grell-ens']       = 'G3 Ens.' 
plot_labels_rc['grell-ens-2050']  = 'G3 2050' 
plot_labels_rc['grell-ens-2100']  = 'G3 2100' 
plot_labels_rc['kf-ens']          = 'KF Ens.' 
plot_labels_rc['kf-ens-2050']     = 'KF 2050' 
plot_labels_rc['kf-ens-2100']     = 'KF 2100' 
plot_labels_rc['g-ens']           = 'G Ens.'  
plot_labels_rc['g-ens-2050']      = 'G 2050' 
plot_labels_rc['g-ens-2100']      = 'G 2100' 
plot_labels_rc['erai-g']          = 'ERA-I'  
plot_labels_rc['t-ens']           = 'T Ens.'
plot_labels_rc['t-ens-2050']      = 'T 2050'
plot_labels_rc['t-ens-2100']      = 'T 2100'
plot_labels_rc['erai-t']          = 'ERA-I'   
# variables
plot_labels_rc['MaxPrecip_1d']    = 'Max Precip. (1d)'
plot_labels_rc['MaxPrecip_5d']    = 'Max Precip. (5d)'
plot_labels_rc['MaxPreccu_1d']    = 'Max Conv. (1d)'
plot_labels_rc['MaxSolprec_1d']   = 'Max Snow (1d)'
plot_labels_rc['MaxSolprec_5d']   = 'Max Snow (5d)'
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
plot_labels_rc['waterflx']        = 'Water Flux'
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

## plot styles for climatology plots
climds_plotargs_rc = default.climds_plotargs_rc.copy()
# observational datasets
obs_args = default.obs_args.copy()
climds_plotargs_rc['Observations'] =  obs_args
climds_plotargs_rc['EC']           =  obs_args
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
clim_defaults.update(heat=(-30,130), Q2=(0,20), aSM=(0.1,0.5), temp=(245,305), wetprec=(0,25), wetfrq=(0,100))
# specific settings
clim_specifics = dict()
clim_specifics['ARB']         = AttrDict(temp=(245,300), water=(-1.5,2.5), precip=(-0.5,3.5), precip_types=(-0.5,3.5), precip_net=(-1.5,3.5), spei=(-1.5,5.5), runoff=(-1.2,2), flux_net=(-1.5,3.5), flux_snow=(-0.5,3.5), evap=(-1.5,5.5), precip_xtrm=(0,18), precip_cesm=(0,18))
clim_specifics['SSR']         = clim_specifics['ARB'].copy()
clim_specifics['CRB']         = AttrDict(temp=(255,305), water=(-3,6), precip=(-0.5,7.))
clim_specifics['FRB']         = AttrDict(temp=(255,300), water=(-2,7.), precip=(-1,7.), precip_types=(-1,7.), runoff=(-2,7), flux_alt=(-1,7.), spei=(-1,7.))
clim_specifics['NRB']         = AttrDict(temp=(245,305), water=(-1.4,2.2), precip=(-0.4,3.4), runoff=(-2.,2.), flux=(-2,5.5))
clim_specifics['PSB']         = AttrDict(temp=(255,295), water=(-2.,16.))
clim_specifics['NorthernPSB'] = AttrDict(temp=(255,295), water=(-2.,14.), precip=(-2,14))
clim_specifics['SouthernPSB'] = AttrDict(temp=(255,295), water=(-2.,16.), precip=(-2,16))
clim_specifics['Pacific']     = AttrDict(temp=(265,300), water=(-2,10)   , precip=(-1,12)    , precip_xtrm=(0,60), precip_cesm=(0,60), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,25), CDD=(0,25))
clim_specifics['Coast']       = AttrDict(temp=(265,300), water=(-4,6)    , precip=(-1,10)    , precip_xtrm=(0,50), precip_cesm=(0,40), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,20), CDD=(0,31))
clim_specifics['Plateau']     = AttrDict(temp=(255,305), water=(-2.,2.)  , precip=(-0.5,2.5), precip_xtrm=(0,20), precip_cesm=(0,20), precip_alt=(0,20), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
clim_specifics['Prairies']    = AttrDict(temp=(250,305), water=(-1.5,1.5), precip=(-0.5,4.), precip_xtrm=(0,28), precip_cesm=(0,28), precip_alt=(0,30), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
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
evads_plotargs_rc = default.evads_plotargs_rc.copy()
evads_plotargs_rc['EC']              = AttrDict(color='black')
evads_plotargs_rc['EC (1935)']       = AttrDict(color='blue')
evads_plotargs_rc['EC (1940)']       = AttrDict(color='blue')
evads_plotargs_rc['EC (1965)']       = AttrDict(color='purple')evads_plotargs_rc['EC (1995)']       = AttrDict(color='red')evads_plotargs_rc['EC (1990)']       = AttrDict(color='red')evads_plotargs_rc['Unity']           = AttrDict(color='black', marker='.')evads_plotargs_rc['Observations']    = AttrDict(color='black')evads_plotargs_rc['WSC']             = AttrDict(color='green', marker='.')evads_plotargs_rc['CRU']             = AttrDict(color='green', marker='.')evads_plotargs_rc['GPCC']            = AttrDict(color='purple', marker='.')evads_plotargs_rc['Ens']             = AttrDict(color='crimson',)evads_plotargs_rc['Ens-2050']        = AttrDict(color='darkorchid',)evads_plotargs_rc['Ens-2100']        = AttrDict(color='royalblue',)evads_plotargs_rc['MEns']            = AttrDict(color='royalblue',)    evads_plotargs_rc['MEns-2050']       = AttrDict(color='darkorchid',)   evads_plotargs_rc['MEns-2100']       = AttrDict(color='crimson',)      evads_plotargs_rc['Ctrl-A']          = AttrDict(color='royalblue',)evads_plotargs_rc['Ctrl-A-2050']     = AttrDict(color='darkorchid',)evads_plotargs_rc['Ctrl-A-2100']     = AttrDict(color='crimson',)evads_plotargs_rc['max-ens']         = AttrDict(color='blue',  )evads_plotargs_rc['max-ens-2050']    = AttrDict(color='purple',)evads_plotargs_rc['max-ens-2100']    = AttrDict(color='red',   )evads_plotargs_rc['max-ens_d01']     = AttrDict(color='green',)evads_plotargs_rc['ctrl-ens']        = AttrDict(color='royalblue',  )evads_plotargs_rc['ctrl-ens-2050']   = AttrDict(color='darkorchid',)evads_plotargs_rc['ctrl-ens-2100']   = AttrDict(color='crimson',   )evads_plotargs_rc['ctrl-ens_d01']    = AttrDict(color='forestgreen',)evads_plotargs_rc['mex-ens']         = AttrDict(color='blue')evads_plotargs_rc['mex-ens-2050']    = AttrDict(color='purple')evads_plotargs_rc['mex-ens-2100']    = AttrDict(color='red')evads_plotargs_rc['mex-ens_d01']     = AttrDict(color='purple', ) evads_plotargs_rc['phys-ens']        = AttrDict(color='forestgreen',  ) evads_plotargs_rc['phys-ens-2050']   = AttrDict(color='brown',) evads_plotargs_rc['phys-ens-2100']   = AttrDict(color='coral',    ) evads_plotargs_rc['phys-ens_d01']    = AttrDict(color='sandybrown',)evads_plotargs_rc['grell-ens']       = AttrDict(color='blue',  )evads_plotargs_rc['grell-ens-2050']  = AttrDict(color='purple',)evads_plotargs_rc['grell-ens-2100']  = AttrDict(color='red',   )evads_plotargs_rc['kf-ens']          = AttrDict(color='seagreen',  ) evads_plotargs_rc['kf-ens-2050']     = AttrDict(color='brown',) evads_plotargs_rc['kf-ens-2100']     = AttrDict(color='coral',    ) evads_plotargs_rc['max-ctrl']        = AttrDict(color='gold')evads_plotargs_rc['max-ctrl-2050']   = AttrDict(color='greenyellow',)evads_plotargs_rc['max-ctrl-2100']   = AttrDict(color='magenta',   )evads_plotargs_rc['max-ctrl_d01']    = AttrDict(color='greenyellow',)evads_plotargs_rc['max-3km']         = AttrDict(color='blue')
evads_plotargs_rc['max-3km-2100']    = AttrDict(color='red')
evads_plotargs_rc['max-ens-A']       = AttrDict(color='orange')evads_plotargs_rc['max-ens-B']       = AttrDict(color='red')evads_plotargs_rc['max-ens-C']       = AttrDict(color='brown')evads_plotargs_rc['max-lowres_d01']  = AttrDict(color='sandybrown',)evads_plotargs_rc['max-lowres']      = AttrDict(color='sandybrown',)evads_plotargs_rc['max-kf']          = AttrDict(color='red')evads_plotargs_rc['max-nosub']       = AttrDict(color='purple')evads_plotargs_rc['max-nmp']         = AttrDict(color='green')evads_plotargs_rc['max-noflake']     = AttrDict(color='dodgerblue')evads_plotargs_rc['old-ctrl']        = AttrDict(color='green')evads_plotargs_rc['old-ctrl-2100']   = AttrDict(color='green')evads_plotargs_rc['ctrl-1']          = AttrDict(color='dodgerblue')evads_plotargs_rc['ctrl-2100']       = AttrDict(color='dodgerblue')evads_plotargs_rc['new-ctrl']        = AttrDict(color='purple')evads_plotargs_rc['new-ctrl-2100']   = AttrDict(color='magenta')evads_plotargs_rc['erai-max']        = AttrDict(color='green')evads_plotargs_rc['erai-max_d01']    = AttrDict(color='blue')
evads_plotargs_rc['erai-3km']        = AttrDict(color='red')
evads_plotargs_rc['erai-3km_d01']    = AttrDict(color='blue')
evads_plotargs_rc['erai-3km_d02']    = AttrDict(color='green')
evads_plotargs_rc['erai-3km_d03']    = AttrDict(color='red')
evads_plotargs_rc['max-3km_d01']     = AttrDict(color='blue')
evads_plotargs_rc['max-3km_d02']     = AttrDict(color='green')
evads_plotargs_rc['max-3km_d03']     = AttrDict(color='red')
evads_plotargs_rc['cfsr-max']        = AttrDict(color='yellow')evads_plotargs_rc['new-v361']        = AttrDict(color='crimson')evads_plotargs_rc['new-v361-2050']   = AttrDict(color='royalblue')evads_plotargs_rc['new-v361-2100']   = AttrDict(color='darkorchid')evads_plotargs_rc['max-seaice-2050'] = AttrDict(color='darkorchid')evads_plotargs_rc['max-seaice-2100'] = AttrDict(color='crimson')evads_plotargs_rc['max-1deg']        = AttrDict(color='green')evads_plotargs_rc['max-hilev']       = AttrDict(color='purple') 
# wrapper with custom defaults to figure creator (plotargs and label positions)def evaFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, **kwargs):  if dataset_plotargs is None: dataset_plotargs = evads_plotargs_rc   if variable_plotargs is None: variable_plotargs = None  if plot_labels is None: plot_labels = plot_labels_rc  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,                          plot_labels=plot_labels, **kwargs)
## annotation for climatology plot
# distribution plot defaults
dist_defaults = default.dist_defaults.copy() # global defaults, followed by local defaults
dist_defaults.update(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.35), waterflx=(-2,12), runoff=(-2,12),                     T2=(245,305), Tmax=(250,310), Tmin=(240,300), precip=(0,20), 
                     MaxWaterflx_5d=(5,35), CWD=(0,30), CDD=(0,50), CNWD=(0,30), CNDD=(0,80),                      CWD_002=(0,30), CDD_002=(0,50), CWD_010=(0,30), CDD_010=(0,80), 
                     CWD_100=(0,15), CDD_100=(0,100), CWD_200=(0,10), CDD_200=(0,100),
                     MaxPrecip_1d=(0,120), MaxPreccu_1d=(0, 40), MaxPrecnc_1d=(0,120),                      MaxSolprec_1d=(0, 40), MaxPrecip_6h=(0,150), MaxPreccu_1h=(0,500), 
                     MaxPrecnc_1h=(0,500), MaxPrecip_5d=(0,20),)
dist_specifics = dict() 
# basins (annual)dist_specifics['ARB'] = AttrDict(MaxWaterflx_5d=(2,14), precip=(0,6), waterflx=(0,6))dist_specifics['CRB'] = AttrDict(water=(-2.5,8.5), waterflx=(2,16))dist_specifics['FRB'] = AttrDict(water=(-2.5,8.5), aSM=(0.2,0.4), precip=(0,12), waterflx=(2,14), runoff=(2,14))dist_specifics['NRB'] = AttrDict(water=(-2.,4.), runoff=(-2.,2.), flux=(-2,5.5), waterflx=(0,7))dist_specifics['PSB'] = AttrDict(water=(-2.,16.))
dist_specifics['SSR'] = AttrDict(MaxWaterflx_5d=(2,14), precip=(0,6), waterflx=(0,6), aSM=(0.05,0.4), runoff=(0,6),)
# clusters (seasonal)dist_specifics['coast_summer']    = AttrDict(MaxPreccu_1h=(0,500), MaxPrecnc_1h=(0,500), MaxPrecip_6h=(0,180), MaxPrecnc_6h=(0,180), MaxPreccu_6h=(0, 80), MaxPrecip_1d=(0, 80), MaxSolprec_1d=(0, 80), MaxPrecnc_1d=(0, 60), MaxPreccu_1d=(0,20), MaxPrecip_5d=(0,30), MaxSolprec_5d=(0,30), MaxWaterflx_5d=(-4,10), CWD_002=(0,40), CDD_002=(0,60), CWD_010=(0,30), CDD_010=(0,120), CNWD=(0,20), CNDD=(0,150))dist_specifics['pacific_summer']  = AttrDict(MaxPreccu_1h=(0,500), MaxPrecnc_1h=(0,500), MaxPrecip_6h=(0,180), MaxPrecnc_6h=(0,180), MaxPreccu_6h=(0, 80), MaxPrecip_1d=(0,100), MaxSolprec_1d=(0,100), MaxPrecnc_1d=(0, 90), MaxPreccu_1d=(0,25), MaxPrecip_5d=(0,50), MaxSolprec_5d=(0,30), MaxWaterflx_5d=(-4,10), CWD_002=(0,40), CDD_002=(0,60), CWD_010=(0,30), CDD_010=(0,120), CNWD=(0,20), CNDD=(0,150))dist_specifics['island_summer']   = dist_specifics['coast_summer']                                                    dist_specifics['plateau_summer']  = AttrDict(MaxPreccu_1h=(0,600), MaxPrecnc_1h=(0,600), MaxPrecip_6h=(0,180), MaxPrecnc_6h=(0,180), MaxPreccu_6h=(0,100), MaxPrecip_1d=(0, 60), MaxSolprec_1d=(0, 60), MaxPrecnc_1d=(0, 60), MaxPreccu_1d=(0,30), MaxPrecip_5d=(0,30), MaxSolprec_5d=(0,30), MaxWaterflx_5d=(-4,10), CWD_002=(0,30), CDD_002=(0,50), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,20), CNDD=(0,120))dist_specifics['north_summer']    = dist_specifics['plateau_summer']                                                  dist_specifics['prairies_summer'] = AttrDict(MaxPreccu_1h=(0,800), MaxPrecnc_1h=(0,800), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,120), MaxPrecip_1d=(0,100), MaxSolprec_1d=(0,100), MaxPrecnc_1d=(0, 90), MaxPreccu_1d=(0,35), MaxPrecip_5d=(0,45), MaxSolprec_5d=(0,45), MaxWaterflx_5d=(-5,15), CWD_002=(0,40), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,60), CNWD=(0,20), CNDD=(0, 90))dist_specifics['pacific_fall']    = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,700), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0,150), MaxSolprec_1d=(0,150), MaxPrecnc_1d=(0,150), MaxPreccu_1d=(0,60), MaxPrecip_5d=(0,110), MaxSolprec_5d=(0,60), MaxWaterflx_5d=(-5,30), CWD_002=(0,80), CDD_002=(0,40), CWD_010=(0,60), CDD_010=(0, 60), CNWD=(0,80), CNDD=(0, 45))
dist_specifics['coast_fall']      = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,700), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0,150), MaxSolprec_1d=(0,150), MaxPrecnc_1d=(0,150), MaxPreccu_1d=(0,60), MaxPrecip_5d=(0,110), MaxSolprec_5d=(0,60), MaxWaterflx_5d=(-5,30), CWD_002=(0,80), CDD_002=(0,40), CWD_010=(0,60), CDD_010=(0, 60), CNWD=(0,80), CNDD=(0, 45))dist_specifics['coast_winter']    = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,700), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0,120), MaxSolprec_1d=(0,120), MaxPrecnc_1d=(0,120), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,60), MaxSolprec_5d=(0,60), MaxWaterflx_5d=(-5,30), CWD_002=(0,80), CDD_002=(0,40), CWD_010=(0,60), CDD_010=(0, 60), CNWD=(0,80), CNDD=(0, 45))dist_specifics['pacific_winter']  = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,700), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0,150), MaxSolprec_1d=(0,150), MaxPrecnc_1d=(0,150), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,110), MaxSolprec_5d=(0,60), MaxWaterflx_5d=(-5,30), CWD_002=(0,80), CDD_002=(0,40), CWD_010=(0,60), CDD_010=(0, 60), CNWD=(0,80), CNDD=(0, 45))dist_specifics['island_winter']   = dist_specifics['coast_winter']                                                    dist_specifics['island_fall']     = dist_specifics['coast_fall']                                                        dist_specifics['plateau_winter']  = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,400), MaxPrecip_6h=(0,120), MaxPrecnc_6h=(0,120), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0, 40), MaxSolprec_1d=(0, 40), MaxPrecnc_1d=(0, 40), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,22), MaxSolprec_5d=(0,22), MaxWaterflx_5d=(-5,15), CWD_002=(0,40), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,60), CNDD=(0, 45))dist_specifics['north_winter']    = dist_specifics['plateau_winter']                                                  dist_specifics['prairies_winter'] = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,250), MaxPrecip_6h=(0, 90), MaxPrecnc_6h=(0, 90), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0, 40), MaxSolprec_1d=(0, 40), MaxPrecnc_1d=(0, 40), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,22), MaxSolprec_5d=(0,22), MaxWaterflx_5d=(-4,10), CWD_002=(0,30), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,60), CNDD=(0, 45))dist_specifics['BC_summer']       = dist_specifics['coast_summer']dist_specifics['AB_summer']       = dist_specifics['prairies_summer'] dist_specifics['ARB_summer'] = AttrDict(MaxPreccu_1h=(0,800), MaxPrecnc_1h=(0,800), MaxPrecip_6h=(0,240), MaxPrecnc_6h=(0,240), MaxPreccu_6h=(0,120), MaxPrecip_1d=(0, 90), MaxSolprec_1d=(0, 40), MaxPrecnc_1d=(0, 80), MaxPreccu_1d=(0,35), MaxPrecip_5d=(0,30), MaxSolprec_5d=(0,35), MaxWaterflx_5d=(-5,15), CWD_002=(0,40), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,60), CNWD=(0,20), CNDD=(0, 90))dist_specifics['BC_winter']       = dist_specifics['coast_winter']dist_specifics['AB_winter']       = dist_specifics['prairies_winter']dist_specifics['ARB_winter'] = AttrDict(MaxPreccu_1h=(0, 10), MaxPrecnc_1h=(0,250), MaxPrecip_6h=(0, 90), MaxPrecnc_6h=(0, 90), MaxPreccu_6h=(0,  0), MaxPrecip_1d=(0, 35), MaxSolprec_1d=(0, 35), MaxPrecnc_1d=(0, 35), MaxPreccu_1d=(0, 0), MaxPrecip_5d=(0,10), MaxSolprec_5d=(0,14), MaxWaterflx_5d=(-4,10), CWD_002=(0,30), CDD_002=(0,40), CWD_010=(0,45), CDD_010=(0,90), CNWD=(0,60), CNDD=(0, 45))
dist_annotation = default.mergeAnnotation(dist_specifics, dist_defaults)
# wrapper with custom annotation defaults for distPlot
def distPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = dist_annotation
  if defaults is None: defaults = dist_defaults
  return eva_plot.distPlot(annotation=annotation, defaults=defaults, **kwargs)

# quantile plot defaults
quant_defaults = default.quant_defaults.copy() # global defaults, followed by local defaults
quant_defaults.update(heat=(-30,120), Q2=(0,20), aSM=(0.15,0.25), waterflx=(-2,12),                       T2=(245,305), precip = (10,30),  runoff=(-2,12),                      MaxPrecip_6h  = (120,260), MaxPrecip_1d = (50,200), MaxPrecip_7d = (0,30),                       MaxWaterflx_5d = (0,50), CDD=(17.5,42.5), CWD=(20,100))quant_specifics = dict()quant_specifics['ARB'] = AttrDict(waterflx=(1.5,4.5), MaxWaterflx_5d = (5,15), aSM=(0.16,0.21))quant_specifics['CRB'] = AttrDict(waterflx=(5,12), aSM=(0.2,0.28))quant_specifics['FRB'] = AttrDict(waterflx=(5,13), runoff=(5,12), MaxWaterflx_5d = (10,40), aSM=(0.2,0.28))quant_specifics['NRB'] = AttrDict(waterflx=(1,4), aSM=(0.16,0.21))quant_specifics['PSB'] = AttrDict(waterflx=(6,14), MaxWaterflx_5d = (10,40), aSM=(0.2,0.28))
quant_specifics['SSR'] = AttrDict(waterflx=(1.5,4.5), runoff=(1,3), MaxWaterflx_5d = (5,15), aSM=(0.16,0.21))
quant_specifics['coast_summer']    = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(200,1000), MaxPrecip_6h=( 75,300), MaxPreccu_6h=( 50,200), MaxPrecip_1d=(30,120), MaxPrecnc_1d=(30,120), MaxPreccu_1d=( 0, 60), MaxPrecip_5d = (10,30), MaxWaterflx_5d = (40,120)) # coast_summer  quant_specifics['pacific_summer']  = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(200,1000), MaxPrecip_6h=( 75,300), MaxPreccu_6h=( 50,200), MaxPrecip_1d=(30,120), MaxPrecnc_1d=(30,120), MaxPreccu_1d=( 0, 60), MaxPrecip_5d = (10,30), MaxWaterflx_5d = (40,120)) # pacific_summer quant_specifics['island_summer']   = quant_specifics['coast_summer']                                                                                                                                                                                                    # island_summer quant_specifics['plateau_summer']  = AttrDict(MaxPreccu_1h=(300,1200), MaxPrecnc_1h=(300,1200), MaxPrecip_6h=( 90,280), MaxPreccu_6h=( 80,220), MaxPrecip_1d=(30,120), MaxPrecnc_1d=(15, 90), MaxPreccu_1d=(15, 60), MaxPrecip_5d = (10,30), MaxWaterflx_5d = (40,120)) # plateau_summer quant_specifics['north_summer']    = quant_specifics['plateau_summer']                                                                                                                                                                                                  # north_summer quant_specifics['prairies_summer'] = AttrDict(MaxPreccu_1h=(400,1800), MaxPrecnc_1h=(400,1600), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(60,160), MaxPrecnc_1d=(30,120), MaxPreccu_1d=(20, 80), MaxPrecip_5d = (10,40), MaxWaterflx_5d = (40,120)) # prairies_summer quant_specifics['pacific_fall']    = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(300, 900), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(90,290), MaxPrecnc_1d=(90,240), MaxPreccu_1d=(30,120), MaxPrecip_5d = (30,90), MaxWaterflx_5d = (40,120)) # pacific_fall  quant_specifics['coast_fall']    = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(300, 900), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(90,290), MaxPrecnc_1d=(90,240), MaxPreccu_1d=(30,120), MaxPrecip_5d = (30,90), MaxWaterflx_5d = (40,120)) # pacific_fall  quant_specifics['coast_winter']    = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(300, 900), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(60,200), MaxPrecnc_1d=(60,200), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = (20,50), MaxWaterflx_5d = (40,120)) # coast_winter  quant_specifics['pacific_winter']  = AttrDict(MaxPreccu_1h=(200,1000), MaxPrecnc_1h=(300, 900), MaxPrecip_6h=(150,320), MaxPreccu_6h=(120,260), MaxPrecip_1d=(90,240), MaxPrecnc_1d=(90,240), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = (30,90), MaxWaterflx_5d = (40,120)) # pacific_winter                                                  quant_specifics['island_winter']   = quant_specifics['coast_winter']                                                                                                                                                                                                    # island_winter quant_specifics['plateau_winter']  = AttrDict(MaxPreccu_1h=(200,1200), MaxPrecnc_1h=(100, 500), MaxPrecip_6h=( 60,180), MaxPreccu_6h=( 30,150), MaxPrecip_1d=(15, 75), MaxPrecnc_1d=(15, 75), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = ( 3,22), MaxWaterflx_5d = (40,120)) # plateau_winter quant_specifics['north_winter']    = quant_specifics['plateau_winter']                                                                                                                                                                                                  # north_winter quant_specifics['prairies_winter'] = AttrDict(MaxPreccu_1h=(200,1600), MaxPrecnc_1h=( 50, 300), MaxPrecip_6h=( 50,120), MaxPreccu_6h=( 20,100), MaxPrecip_1d=(15, 75), MaxPrecnc_1d=(15, 75), MaxPreccu_1d=( 0,  0), MaxPrecip_5d = ( 3,22), MaxWaterflx_5d = (40,120)) # prairies_winter quant_specifics['BC'] = AttrDict(MaxPreccu_1h=( 0,1500), MaxPrecip_1d=(60,180))quant_specifics['AB'] = AttrDict(MaxPrecip_1d=(15,75))quant_specifics['ON'] = AttrDict(MaxPrecip_1d=(10,80))quant_annotation = default.mergeAnnotation(quant_specifics, quant_defaults)

# wrapper with custom annotation defaults for distPlot
def quantPlot(annotation=None, defaults=None, variable_list=None, **kwargs):
  if annotation is None: annotation = quant_annotation
  if defaults is None: defaults = quant_defaults
  return eva_plot.quantPlot(annotation=annotation, defaults=defaults, **kwargs)
