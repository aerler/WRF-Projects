'''
Created on Jun 4, 2016

A module that provides various project specific definitions, mainly related to experiments and plotting. 

@author:  Andre R. Erler, GPL v3
'''

# external imports
import matplotlib as mpl
# internal imports
from geodata.misc import AttrDict
from plotting import figure
from clim import load_ens, plots
from eva import eva
from WRF_experiments import WRF_ens, WRF_exps
from projects.CESM_experiments import CESM_ens, CESM_exps


# variable collections
wetday_extensions = load_ens.wetday_extensions[:3]
variables_rc = dict(); VL = load_ens.VL
# mostly for hydrological analysis
variables_rc['temp']            = VL(vars=('T2', 'Tmax', 'Tmin'), files=('srfc','xtrm',), label='2m Temperature')
# variables_rc['temp']          = VL(vars=('T2',), files=('srfc',), label='2m Temperature')
variables_rc['precip_obs']      = VL(vars=('precip', 'solprec', 'wetfrq_010'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm_obs'] = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'wetprec_010'), files=('hydro',), label='Precipitation')
variables_rc['precip']          = VL(vars=('precip', 'preccu', 'solprec'), files=('hydro',), label='Precipitation')
variables_rc['precip_xtrm']     = VL(vars=('MaxPrecip_1d', 'MaxPrecip_5d', 'MaxPreccu_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_cesm']     = VL(vars=('MaxPrecip_1d', 'MaxPreccu_1d', ), files=('hydro',), label='Precipitation')
variables_rc['precip_alt']      = VL(vars=('MaxPrecip_1d', 'MaxPrecnc_1d', 'MaxSolprec_1d'), files=('hydro',), label='Precipitation')
variables_rc['precip_types']    = VL(vars=('precip','preccu','precnc'), files=('hydro',), label='Water Flux')
variables_rc['precip_net']      = VL(vars=('precip','solprec','p-et'), files=('hydro',), label='Water Flux')
variables_rc['flux_snow']       = VL(vars=('precip','snwmlt','solprec'), files=('hydro',), label='Water Flux')
variables_rc['flux_days']       = VL(vars=('wetfrq_010','snwmlt','p-et'), files=('hydro',), label='Water Flux')
variables_rc['wetprec']         = VL(vars=['wetprec'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Precip.')
variables_rc['wetdays']         = VL(vars=['wetfrq'+ext for ext in wetday_extensions], files=('hydro',), label='Wet-day Ratio')
variables_rc['CWD']             = VL(vars=['CWD'+ext for ext in wetday_extensions]+['CNWD'], files=('hydro',), label='Continuous Wet-days')
variables_rc['CDD']             = VL(vars=['CDD'+ext for ext in wetday_extensions[:-1]]+['CNDD'], files=('hydro',), label='Continuous Dry-days')
variables_rc['sfcflx']          = VL(vars=('p-et','snwmlt','waterflx',), files=('hydro',), label='Surface Flux')
variables_rc['runoff']          = VL(vars=('waterflx','sfroff','runoff'), files=('lsm','hydro'), label='Runoff')
variables_rc['runoff_flux']     = VL(vars=('runoff','snwmlt','p-et'), files=('lsm','hydro'), label='Water Flux')
variables_rc['heat']            = VL(vars=('hfx','lhfx','rSM'),files=('srfc','lsm'), label='Energy Flux')
variables_rc['evap']            = VL(vars=('p-et','evap','pet',), files=('hydro',), label='Water Flux')
variables_rc['spei']            = VL(vars=('precip','evap','pet',), files=('aux','hydro',), label='Water Flux')
variables_rc['Q2']              = VL(vars=('Q2',),files=('srfc',), label='2m Humidity')
variables_rc['aSM']             = VL(vars=('aSM',),files=('lsm',), label='Soil Moisture') 
variables_rc['rSM']             = VL(vars=('rSM',),files=('lsm',), label='Relative Soil Moisture')
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
constraints_rc['max_zerr'] = 300 # can't use this, because we are loading EC data separately from WRF
constraints_rc['prov'] = ('BC','AB')
constraints_rc['end_after'] = 1980
                        
# dataset collections
exps_rc = dict(); EX = load_ens.EX
exps_rc['obs']     = EX(name='obs',exps=['CRU','WSC'], styles=['-','-.'], master='CRU', title='Observations')
exps_rc['val']     = EX(name='val',exps=['max-ens','erai-max','max-ens_d01'], master='max-ens',
                        styles=['-','--','-.'], title='WRF (IC Ens. Avg., ERA-I, 30km)')
exps_rc['prj']     = EX(name='prj',exps=['max-ens','max-ens-2050','max-ens-2100'], master='max-ens',
                        styles=['-','-.','--'], title='IC Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['max']     = EX(name='max',exps=['max-ens','max-ens-2050','max-ens-2100'], master='max-ens',
                        styles=['-','-.','--'], title='IC Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['ctrl']    = EX(name='ctrl',exps=['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], master='ctrl-ens',
                        styles=['-','-.','--'], title='Alt. IC Ens. Average (Hist., Mid-, End-Century)')
exps_rc['g-ens']   = EX(name='g-ens',exps=['g-ens','g-ens-2050','g-ens-2100'], master='g-ens',
                        styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['t-ens']   = EX(name='t-ens',exps=['t-ens','t-ens-2050','t-ens-2100'], master='t-ens',
                        styles=['-','-.','--'], title='G Ensemble Average (Hist., Mid-, End-Century)')
exps_rc['mex']     = EX(name='mex',exps=['mex-ens','mex-ens-2050','mex-ens-2100'], master='mex-ens',
                        styles=['-','-.','--'], title='WRF Ext. Ens. (Hist., Mid-, End-Century)')
exps_rc['phys']    = EX(name='phys',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], master='phys-ens',
                        styles=['-','-.','--'], title='WRF Phys. Ens. (Hist., Mid-, End-Century)')
exps_rc['nmp']     = EX(name='nmp',exps=['max-ens','max-ens-2050','max-nmp','max-nmp-2050'], master='max-nmp',
                        styles=['-.','.','--','-'], title='IC Ensemble Average & Noah-MP (Hist., Mid-Century)')
exps_rc['test']    = EX(name='test',exps=['max-ens','max-ctrl','max-hilev'], master='max-ens',
                        styles=['-','--','-.'], title='WRF (Ens. Avg., Ctrl-1, Hi-lev)')
exps_rc['1deg']    = EX(name='1deg',exps=['max-ens','max-ens_d01','max-1deg'], master='max-ens',
                        styles=['--','-.','-'], title='WRF (Ens. Avg., 30km, 1deg)')
exps_rc['seaice']  = EX(name='seaice',exps=['max-ens','max-ens-2050','max-seaice-2050','max-ens-2100','max-seaice-2100'],
                        styles=['-.','--','-','.--','.-'], master='max-ens', 
                        title='WRF (Hist., Mid-, End-Century; Ens. Avg., Sea-ice)')
exps_rc['max-Z']   = EX(name='max-Z',exps=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], master='max-ctrl',
                        styles=['-','-.','--'], title='WRF-Z (Hist., Mid-, End-Century)')
exps_rc['max-A']   = EX(name='max-A',exps=['max-ens-A','max-ens-A-2050','max-ens-A-2100'], master='max-ens-A',
                        styles=['-','-.','--'], title='WRF-A (Hist., Mid-, End-Century)')
exps_rc['max-B']   = EX(name='max-B',exps=['max-ens-B','max-ens-B-2050','max-ens-B-2100'], master='max-ens-B',
                        styles=['-','-.','--'], title='WRF-B (Hist., Mid-, End-Century)')
exps_rc['max-C']   = EX(name='max-C',exps=['max-ens-C','max-ens-C-2050','max-ens-C-2100'], master='max-ens-C',
                        styles=['-','-.','--'], title='WRF-C (Hist., Mid-, End-Century)')
exps_rc['si25']    = EX(name='si25',exps=['max-ens','max-ens-2050','max-seaice-2050'], master='max-ens-2050',
                          styles=[':','--','-'], title='WRF (Hist., Mid-Century; Ens. Avg., Sea-ice)')
exps_rc['si21']    = EX(name='si21',exps=['max-ens','max-ens-2100','max-seaice-2100'], master='max-ens-2100',
                          styles=[':','--','-'], title='WRF (Hist., End-Century; Ens. Avg., Sea-ice)')
# dataset collections for EVA
exps_rc['val-res'] = EX(name='val-res', exps=['Observations', 'max-ens_d01', 'max-ens', 'Ens'], # ,'erai-max','ctrl-1','max-1deg' 
                        master='max-ens', reference='Observations', target=None, title='Resolution')
exps_rc['val-all'] = EX(name='val-all', exps=['Observations', 'Ens', 'max-ens', 'ctrl-ens'], # ,'erai-max','ctrl-1','max-1deg' 
                        master='max-ens', reference='Observations', target=None, title='Resolution')
exps_rc['max-all'] = EX(name='max-all', exps=['Observations', 'max-ens','max-ens_d01','erai-max','max-ens-2050','max-ens-2100'], # ,'erai-max','ctrl-1','max-1deg' 
                        master='max-ens', reference='Observations', target='auto', title='Validation & Projection')
exps_rc['cesm-all'] = EX(name='cesm-all', exps=['Observations', 'Ens','Ens-2050','Ens-2100'],  
                        master='Ens', reference='Observations', target='Ens', title='CESM Validation & Projection')
exps_rc['cesm-mal'] = EX(name='cesm-ensa', exps=['Observations', 'MEns','MEns-2050','MEns-2100'],  
                        master='MEns', reference='Observations', target='MEns', title='CESM Validation & Projection')
exps_rc['cesm-ens'] = EX(name='cesm-ens', exps=['MEns','MEns-2050','MEns-2100'],  
                        master='MEns', reference='MEns', target='MEns', title='CESM Projection')
exps_rc['max-val'] = EX(name='max-val', exps=['Observations','max-ens_d01','max-ens','erai-max'], # ,'erai-max','ctrl-1','max-1deg' 
                        master='max-ens', reference='Observations', target=None, title='Validation')
exps_rc['max-prj'] = EX(name='prj', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                        master='max-ens', reference=None, target=None, title='Projection')
exps_rc['sens-res']  = EX(name='sens-res', exps=['Observations', 'Ens', 'max-ctrl','max-ctrl_d01','max-1deg','max-lowres_d01'], 
                          master='max-ctrl', reference='Observations', target=None, title='Sensitivity to Resolution')
exps_rc['sens-max']  = EX(name='sens-max',exps=['max-ctrl','max-nosub','max-kf','max-nmp','max-noflake'], 
                          master='max-ctrl', reference=None, target=None, title='Sensitivity Tests (Max)')
exps_rc['sens-phys'] = EX(name='sens-phys',exps=['Observations','max-ctrl','old-ctrl','ctrl-1','new-ctrl','new-v361'], 
                          master='max-ctrl', reference='Observations', target=None, title='Sensitivity (Phys. Ens.)')
exps_rc['sens-phys-2100'] = EX(name='sens-phys-2100',exps=['max-ctrl-2100','old-ctrl-2100','ctrl-2100','new-ctrl-2100','new-v361-2100'], 
                               master='max-ctrl-2100', reference=None, target=None, title='Sensitivity (Phys. Ens. 2100)')
exps_rc['sens-ens']  = EX(name='sens-ens',exps=['Observations','max-ens', 'max-ctrl','max-ens-A','max-ens-B','max-ens-C','erai-max'], 
                          master='max-ens', reference='Observations', target=None, title='Sensitivity (Phys. Ens.)')
# exps_rc['sens-cu']     = EX(name='sens-cu', exps=['grell-ens','kf-ens','grell-ens-2050','kf-ens-2050','grell-ens-2100','kf-ens-2100'], 
#                             master=['grell-ens','grell-ens-2050','grell-ens-2100'], reference=None, target=None, title='Grell vs. KF')
exps_rc['sens-cu']     = EX(name='sens-cu', exps=['grell-ens','kf-ens','grell-ens-2100','kf-ens-2100'], 
                            master=['grell-ens','grell-ens-2100'], reference=None, target=None, title='Grell vs. KF')
exps_rc['kf-ens']     = EX(name='kf-ens', exps=['kf-ens','kf-ens-2050','kf-ens-2100'], 
                            master='kf-ens', reference=None, target=None, title='KF Ens.')
exps_rc['max']     = EX(name='max', exps=['Observations','max-ens','max-ens-2050','max-ens-2100'], 
                        master='max-ens', reference='Observations', target='max-ens', title='IC Ensemble')
exps_rc['maxs']     = EX(name='maxs', exps=['Observations','max-ens','max-ens-2100'], 
                        master='max-ens', reference='Observations', target='max-ens', title='IC Ensemble')
exps_rc['ctrl']     = EX(name='ctrl', exps=['Observations','ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         master='ctrl-ens', reference='Observations', target='ctrl-ens', title='DE Ensemble')
exps_rc['ctrl-all'] = EX(name='ctrl-all', exps=['Observations','ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         master='ctrl-ens', reference='Observations', target='ctrl-ens', title='DE Ensemble')
exps_rc['ctrl-prj'] = EX(name='ctrl-prj', exps=['ctrl-ens','ctrl-ens-2050','ctrl-ens-2100'], 
                         master='ctrl-ens', reference=None, target=None, title='DE Ensemble')
# exps_rc['max-phys']     = EX(name='max-phys', exps=['max-ens','phys-ens','max-ens-2050','phys-ens-2050','max-ens-2100','phys-ens-2100'], 
#                         master=['max-ens','max-ens-2050','max-ens-2100'], reference=None, target=None, title='IC vs. Phys.')
exps_rc['max-phys']     = EX(name='max-phys', exps=['max-ens','phys-ens','max-ens-2100','phys-ens-2100'], 
                        master=['max-ens','max-ens-2100'], reference=None, target=None, title='IC vs. Phys.')
exps_rc['max-ctrl']  = EX(name='max-ctrl', exps=['max-ens','max-ctrl','max-ens-2050','max-ctrl-2050','max-ens-2100','max-ctrl-2100'], 
                        master=['max-ens','max-ens-2050','max-ens-2100'], reference=None, target=None, title='IC & Max')
exps_rc['max-shape'] = EX(name='shape', exps=['Observations','max-ens','max-ens-2050','max-ens-2100'], 
                          master='max-ens', reference='Observations', target=None, title='Projected Shape Change')
exps_rc['max-shape-cu'] = EX(name='shape-cu', exps=['max-ens','max-ens-2050','max-ens-2100'], 
                          master='max-ens', reference='max-ens', target=None, title='Projected Shape Change')
exps_rc['new']     = EX(name='new', exps=['max-ctrl','new-v361','max-ctrl-2050','new-v361-2050','max-ctrl-2100','new-v361-2100'], 
                        master=['max-ctrl','max-ctrl-2050','max-ctrl-2100'], reference=None, target=None, title='New V3.6.1')
exps_rc['mex']     = EX(name='mex', exps=['Observations','mex-ens','mex-ens-2050','mex-ens-2100'], 
                        master='mex-ens', reference='Observations', target='mex-ens', title='IC Ensemble (Ext.)')
exps_rc['phys']    = EX(name='phys',exps=['Observations','phys-ens','phys-ens-2050','phys-ens-2100'], 
                        master='phys-ens', reference='Observations', target='phys-ens', title='Physics Ensemble')
exps_rc['phys-val']    = EX(name='phys-val',exps=['Observations','phys-ens','phys-ens_d01'], 
                        master='phys-ens', reference='Observations', target=None, title='Physics Validation')
exps_rc['phys-prj']    = EX(name='phys-prj',exps=['phys-ens','phys-ens-2050','phys-ens-2100'], 
                            master='phys-ens', reference=None, target=None, title='Physics Projection')
exps_rc['phys-shape'] = EX(name='shape', exps=['Observations','phys-ens','phys-ens-2050','phys-ens-2100'], 
                           master='phys-ens', reference='Observations', target=None, title='Projected Shape Change')
# exps_rc['seaice']  = EX(name='seaice', exps=['Observations','max-ctrl','max-seaice-2050','max-seaice-2100'], 
#                         master='max-ctrl', reference=None, target=None, title='Seaice Exp.')
exps_rc['seaice']  = EX(name='seaice', exps=['max-ens-2050','max-seaice-2050','max-ens-2100','max-seaice-2100'], 
                        master=['max-ens-2050','max-ens-2100'], reference=None, target=None, title='Seaice Exp.')
exps_rc['ens-all']  = EX(name='ens-all', exps=['Observations','max-ens','max-ctrl','erai-max','phys-ens','max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050','max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                        master=['max-ens','max-ens-2050','max-ens-2100'], reference='Observations', target=None, title='Ensembles')
exps_rc['ens-1980']  = EX(name='ens-1980', exps=['Observations','max-ens','max-ctrl','erai-max','phys-ens'], 
                        master='max-ens', reference='Observations', target=None, title='Ensembles 1980')
exps_rc['ens-2050']  = EX(name='ens-2050', exps=['max-ens-2050','max-ctrl-2050','max-seaice-2050','phys-ens-2050'], 
                        master='max-ens-2050', reference=None, target=None, title='Ensembles 2050')
exps_rc['ens-2100']  = EX(name='ens-2100', exps=['max-ens-2100','max-ctrl-2100','max-seaice-2100','phys-ens-2100'], 
                          master='max-ens-2100', reference=None, target=None, title='Ensembles 2100')


# set default variable atts for load functions from load_ens
def loadShapeObservations(variable_atts=None, **kwargs):
  ''' wrapper for clim.load_ens.loadShapeObservations that sets variable lists '''
  if variable_atts is None: variable_atts = variables_rc
  return load_ens.loadShapeObservations(variable_atts=variable_atts, **kwargs)
def loadShapeEnsemble(variable_atts=None, **kwargs):
  ''' wrapper for clim.load_ens.loadShapeEnsemble that sets experiment and variable lists '''
  if variable_atts is None: variable_atts = variables_rc  
  return load_ens.loadShapeEnsemble(variable_atts=variable_atts, WRF_exps=WRF_exps, CESM_exps=CESM_exps, 
                                    WRF_ens=WRF_ens, CESM_ens=CESM_ens, **kwargs)
def loadStationEnsemble(variable_atts=None, **kwargs):
  ''' wrapper for clim.load_ens.loadStationEnsemble that sets experiment and variable lists '''
  if variable_atts is None: variable_atts = variables_rc  
  return load_ens.loadStationEnsemble(variable_atts=variable_atts, WRF_exps=WRF_exps, CESM_exps=CESM_exps, 
                                      WRF_ens=WRF_ens, CESM_ens=CESM_ens, **kwargs)


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
plot_labels_rc['max-ens']         = 'IC Ens.'  
plot_labels_rc['max-ens-2050']    = 'IC 2050' 
plot_labels_rc['max-ens-2100']    = 'IC 2100' 
plot_labels_rc['max-ens_d01']     = 'IC (D1)'
plot_labels_rc['erai-max']        = 'ERA-I'  
plot_labels_rc['ctrl-ens']        = 'Alt. Ens.'  
plot_labels_rc['ctrl-ens-2050']   = 'AE 2050' 
plot_labels_rc['ctrl-ens-2100']   = 'AE 2100' 
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

## custom default plot styles
obs_args = AttrDict(marker='o', linestyle=' ') # 5*mpl.rcParams['lines.linewidth']
rea_args = AttrDict(marker='^', markersize=mpl.rcParams['lines.markersize'], linestyle=' ') # 'small'
# datasets (mainly obs)
dataset_plotargs_rc = dict()
dataset_plotargs_rc['Observations'] =  obs_args
dataset_plotargs_rc['WSC']          =  obs_args
dataset_plotargs_rc['Unity']        =  obs_args
dataset_plotargs_rc['CRU']          =  obs_args
dataset_plotargs_rc['GPCC']         =  obs_args
dataset_plotargs_rc['CFSR']         =  rea_args
dataset_plotargs_rc['NARR']         =  rea_args
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
variable_plotargs_rc['liqprec']        = AttrDict(color = 'cyan')
variable_plotargs_rc['solprec']        = AttrDict(color = 'blue')
variable_plotargs_rc['dryprec']        = AttrDict(color = 'green')
variable_plotargs_rc['preccu']         = AttrDict(color = 'magenta')
variable_plotargs_rc['precnc']         = AttrDict(color = 'coral')
variable_plotargs_rc['evap']           = AttrDict(color = 'red')
variable_plotargs_rc['p-et']           = AttrDict(color = 'red')
variable_plotargs_rc['pet']            = AttrDict(color = 'purple')
variable_plotargs_rc['waterflx']       = AttrDict(color = 'dodgerblue') # 'dodgerblue'
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
# add wet-day threshold dependent variables    
wetday_colors = ['steelblue', 'purple', 'crimson', 'orange']   
for wdext,color in zip(load_ens.wetday_extensions,wetday_colors):
  variable_plotargs_rc['wetprec'+wdext] = AttrDict(color = color)
  variable_plotargs_rc['wetfrq' +wdext] = AttrDict(color = color)
  variable_plotargs_rc['CWD'+wdext]     = AttrDict(color = color)
  variable_plotargs_rc['CDD' +wdext]    = AttrDict(color = color)

# wrapper with custom defaults to figure creator (plotargs and label positions)
def getFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, xtop=None, yright=None, **kwargs):
  if dataset_plotargs is None: dataset_plotargs = dataset_plotargs_rc 
  if variable_plotargs is None: variable_plotargs = variable_plotargs_rc
  if plot_labels is None: plot_labels = plot_labels_rc
  return figure.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,
                         plot_labels=plot_labels, xtop=xtop, yright=yright, **kwargs)
climFigAx = getFigAx # alias for direct project import

## climatology plot and shape annotation
# defaults
shape_defaults_rc = AttrDict(heat=(-30,130), Q2=(0,20), aSM=(0.1,0.5), temp=(245,305), wetprec=(0,25), wetfrq=(0,100))
# specific settings
shape_specifics = dict()
shape_specifics['ARB']         = AttrDict(temp=(245,300), water=(-1.5,2.5), precip=(-0.5,3.5), precip_types=(-0.5,3.5), spei=(-1.5,5.5),
                                          runoff=(-1.2,2), flux=(-1.5,3.5), flux_alt=(-0.5,3.5), evap=(-1.5,5.5))
shape_specifics['CRB']         = AttrDict(temp=(255,305), water=(-3,6), precip=(-0.5,7.))
shape_specifics['FRB']         = AttrDict(temp=(255,300), water=(-2,7.), precip=(-1,7.), precip_types=(-1,7.), runoff=(-2,7), flux_alt=(-1,7.), spei=(-1,7.))
shape_specifics['NRB']         = AttrDict(temp=(245,305), water=(-1.4,2.2), precip=(-0.4,3.4), runoff=(-2.,2.), flux=(-2,5.5))
shape_specifics['PSB']         = AttrDict(temp=(255,295), water=(-2.,16.))
shape_specifics['NorthernPSB'] = AttrDict(temp=(255,295), water=(-2.,14.), precip=(-2,14))
shape_specifics['SouthernPSB'] = AttrDict(temp=(255,295), water=(-2.,16.), precip=(-2,16))
shape_specifics['Pacific']     = AttrDict(temp=(265,300), water=(-2,10)   , precip=(0,12)    , precip_xtrm=(0,60), precip_cesm=(0,60), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,25), CDD=(0,25))
shape_specifics['Coast']       = AttrDict(temp=(265,300), water=(-4,6)    , precip=(-1,8)    , precip_xtrm=(0,40), precip_cesm=(0,40), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,20), CDD=(0,31))
shape_specifics['Plateau']     = AttrDict(temp=(255,305), water=(-2.,2.)  , precip=(-0.5,2.5), precip_xtrm=(0,20), precip_cesm=(0,20), precip_alt=(0,20), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
shape_specifics['Prairies']    = AttrDict(temp=(250,305), water=(-1.5,1.5), precip=(-0.5,3.5), precip_xtrm=(0,30), precip_cesm=(0,30), precip_alt=(0,30), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
# add defaults to specifics
shape_annotation_rc = dict()
for basin,specs in shape_specifics.iteritems():
  atts = shape_defaults_rc.copy()
  atts.update(specs)
  atts.update(zip([key+'_obs' for key in atts.iterkeys()],atts.itervalues())) # add obs-only versions
  shape_annotation_rc[basin] = atts
  
# wrapper with custom basin/sape/etc. defaults for climPlot
def climPlot(shape_annotation=None, shape_defaults=None, variable_atts=None, **kwargs):
  if shape_annotation is None: shape_annotation = shape_annotation_rc 
  if shape_defaults is None: shape_defaults = shape_defaults_rc
  if variable_atts is None: variable_atts = variables_rc
  return plots.climPlot(shape_annotation=shape_annotation, shape_defaults=shape_defaults, 
                        variable_atts=variable_atts, **kwargs)


if __name__ == '__main__':
    pass