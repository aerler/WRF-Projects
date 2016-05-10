'''
Created on Feb 10, 2015

Utility module to provide specialized plotting functionality for hydro-climatological analysis
on a river basin scale.

@author: Andre R. Erler, GPL v3
'''

# external imports
import matplotlib as mpl
import numpy as np
# internal imports
from datasets.common import name_of_month, BatchLoad
from geodata.misc import ArgumentError, AttrDict
import plotting.figure as pltfig
from plotting.figure import show # don't import getFigAx directly, to avoid recursion
from basins import loadShapeEnsemble, loadShapeObservations, variables_rc
from geodata.base import Dataset
from hydro.basins import wetday_extensions

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
plot_labels_rc['MEns']             = 'CESM*'  
plot_labels_rc['MEns-2050']        = 'CESM* 2050' 
plot_labels_rc['MEns-2100']        = 'CESM* 2100' 
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
plot_labels_rc['MaxPrecip_1d']   = 'Max Precip. (1d)'
plot_labels_rc['MaxPrecip_5d']   = 'Max Precip. (5d)'
plot_labels_rc['MaxPreccu_1d']   = 'Max Conv. (1d)'
plot_labels_rc['MaxSolprec_1d']  = 'Max Snow (1d)'
plot_labels_rc['MaxWaterFlx_1d'] = 'Max Flux (1d)'
plot_labels_rc['wetfrq_010']     = 'Wet-days'
plot_labels_rc['wetprec_010']    = 'Precip. Intensity'
plot_labels_rc['T2']       = 'T (2m)'   
plot_labels_rc['precip']   = 'Precip.'  
plot_labels_rc['liqprec']  = 'Liquid' 
plot_labels_rc['solprec']  = 'Snow' 
plot_labels_rc['dryprec']  = 'dryprec' 
plot_labels_rc['preccu']   = 'Conv.'  
plot_labels_rc['precnc']   = 'NC'  
plot_labels_rc['evap']     = 'ET'    
plot_labels_rc['p-et']     = 'Net Precip.'    
plot_labels_rc['pet']      = 'PET' 
plot_labels_rc['waterflx'] = 'Water Flux'
plot_labels_rc['snwmlt']   = 'Snow Melt' 
plot_labels_rc['runoff']   = 'Total Runoff' 
plot_labels_rc['ugroff']   = 'Undergr. R\'off' 
plot_labels_rc['sfroff']   = 'Surface Runoff' 
plot_labels_rc['Tmax']     = 'T (max)'    
plot_labels_rc['Tmin']     = 'T (min)'    
plot_labels_rc['hfx']      = 'Sens. Heat'     
plot_labels_rc['lhfx']     = 'Latent Heat'    
plot_labels_rc['Q2']       = 'Q (2m)'      
plot_labels_rc['aSM']      = 'aSM'     
plot_labels_rc['rSM']      = 'Soil Moist.'       

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
variable_plotargs_rc['MaxPrecip_1d'] = AttrDict(color = 'green')
variable_plotargs_rc['MaxPrecip_5d'] = AttrDict(color = 'sienna')
variable_plotargs_rc['MaxPreccu_1d'] = AttrDict(color = 'magenta')
variable_plotargs_rc['MaxPrecnc_1d'] = AttrDict(color = 'grey')
variable_plotargs_rc['MaxSolprec_1d']  = AttrDict(color = 'blue')
variable_plotargs_rc['MaxWaterFlx_1d'] = AttrDict(color = 'blue')
variable_plotargs_rc['T2']       = AttrDict(color = 'green')
variable_plotargs_rc['precip']   = AttrDict(color = 'green')
variable_plotargs_rc['liqprec']  = AttrDict(color = 'cyan')
variable_plotargs_rc['solprec']  = AttrDict(color = 'blue')
variable_plotargs_rc['dryprec']  = AttrDict(color = 'green')
variable_plotargs_rc['preccu']   = AttrDict(color = 'magenta')
variable_plotargs_rc['precnc']   = AttrDict(color = 'coral')
variable_plotargs_rc['evap']     = AttrDict(color = 'red')
variable_plotargs_rc['p-et']     = AttrDict(color = 'red')
variable_plotargs_rc['pet']      = AttrDict(color = 'purple')
variable_plotargs_rc['waterflx'] = AttrDict(color = 'dodgerblue') # 'dodgerblue'
variable_plotargs_rc['snwmlt']   = AttrDict(color = 'orange')
variable_plotargs_rc['runoff']   = AttrDict(color = 'purple')
variable_plotargs_rc['ugroff']   = AttrDict(color = 'coral')
variable_plotargs_rc['sfroff']   = AttrDict(color = 'green')
variable_plotargs_rc['Tmax']     = AttrDict(color = 'red')
variable_plotargs_rc['Tmin']     = AttrDict(color = 'blue')
variable_plotargs_rc['hfx']      = AttrDict(color = 'red')
variable_plotargs_rc['lhfx']     = AttrDict(color = 'purple')
variable_plotargs_rc['Q2']       = AttrDict(color = 'blue')
variable_plotargs_rc['aSM']      = AttrDict(color = 'coral')
variable_plotargs_rc['rSM']      = AttrDict(color = 'green')
variable_plotargs_rc['CNWD']     = AttrDict(color = 'green')
variable_plotargs_rc['CNDD']     = AttrDict(color = 'green')
# add wet-day threshold dependent variables    
wetday_colors = ['steelblue', 'purple', 'crimson', 'orange']   
for wdext,color in zip(wetday_extensions,wetday_colors):
  variable_plotargs_rc['wetprec'+wdext] = AttrDict(color = color)
  variable_plotargs_rc['wetfrq' +wdext] = AttrDict(color = color)
  variable_plotargs_rc['CWD'+wdext]     = AttrDict(color = color)
  variable_plotargs_rc['CDD' +wdext]    = AttrDict(color = color)
# wrapper with custom defaults to figure creator (plotargs and label positions)
def getFigAx(subplot, dataset_plotargs=None, variable_plotargs=None, plot_labels=None, xtop=None, yright=None, **kwargs):
  if dataset_plotargs is None: dataset_plotargs = dataset_plotargs_rc 
  if variable_plotargs is None: variable_plotargs = variable_plotargs_rc
  if plot_labels is None: plot_labels = plot_labels_rc
  return pltfig.getFigAx(subplot, dataset_plotargs=dataset_plotargs, variable_plotargs=variable_plotargs,
                         plot_labels=plot_labels, xtop=xtop, yright=yright, **kwargs)

# basin annotation
# defaults
basin_defaults = AttrDict(heat=(-30,130), Q2=(0,20), aSM=(0.1,0.5), temp=(245,305), wetprec=(0,25), wetfrq=(0,100))
# specific settings
basin_specifics = dict()
basin_specifics['ARB'] = AttrDict(temp=(245,300), water=(-1.5,2.5), precip=(-0.5,3.5), precip_types=(-0.5,3.5), spei=(-1.5,5.5),
                                  runoff=(-1.2,2), flux=(-1.5,3.5), flux_alt=(-0.5,3.5), evap=(-1.5,5.5))
basin_specifics['CRB'] = AttrDict(temp=(255,305), water=(-3,6), precip=(-0.5,7.))
basin_specifics['FRB'] = AttrDict(temp=(255,300), water=(-2,7.), precip=(-1,7.), precip_types=(-1,7.), runoff=(-2,7), flux_alt=(-1,7.), spei=(-1,7.))
basin_specifics['NRB'] = AttrDict(temp=(245,305), water=(-1.4,2.2), precip=(-0.4,3.4), runoff=(-2.,2.), flux=(-2,5.5))
basin_specifics['PSB'] = AttrDict(temp=(255,295), water=(-2.,16.))
basin_specifics['NorthernPSB'] = AttrDict(temp=(255,295), water=(-2.,14.), precip=(-2,14))
basin_specifics['SouthernPSB'] = AttrDict(temp=(255,295), water=(-2.,16.), precip=(-2,16))
basin_specifics['Pacific']  = AttrDict(temp=(265,300), water=(-2,10)   , precip=(0,12)    , precip_xtrm=(0,60), precip_cesm=(0,60), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,25), CDD=(0,25))
basin_specifics['Coast']    = AttrDict(temp=(265,300), water=(-4,6)    , precip=(-1,8)    , precip_xtrm=(0,40), precip_cesm=(0,40), precip_alt=(0,40), wetprec=(0,30), wetdays=(0,80), CWD=(0,20), CDD=(0,31))
basin_specifics['Plateau']  = AttrDict(temp=(255,305), water=(-2.,2.)  , precip=(-0.5,2.5), precip_xtrm=(0,20), precip_cesm=(0,20), precip_alt=(0,20), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
basin_specifics['Prairies'] = AttrDict(temp=(250,305), water=(-1.5,1.5), precip=(-0.5,3.5), precip_xtrm=(0,30), precip_cesm=(0,30), precip_alt=(0,30), wetprec=(0,30), wetdays=(0,80), CWD=(0,15), CDD=(0,25))
# add defaults to specifics
basin_annotation = dict()
for basin,specs in basin_specifics.iteritems():
  atts = basin_defaults.copy()
  atts.update(specs)
  atts.update(zip([key+'_obs' for key in atts.iterkeys()],atts.itervalues())) # add obs-only versions
  basin_annotation[basin] = atts



## custom plotting functions
@BatchLoad
def hydroPlot(axes=None, expens=None, obsens=None, experr=None, obserr=None, varlist=None, master=None, 
              scalevars='auto', legend=0, shape_name=None, stnset_name=None, ylim=None,
              lperi=True, lparasiteMeans=True, axtitle=None, ylabel=True, xlabel=True, lyint=True,
              lprint=True, dataset_legend=False, dataset_labels=None, **plotargs):
  ''' plot the seasonal cycle over a basin ''' # linestyles=None, lineformats=None, colors=None, markers=None
  if axes is None: raise ArgumentError
  # define some meta data
  if varlist is None and expens.basetype is Dataset: raise TypeError
  if isinstance(varlist,basestring):
    varlist_name = varlist; varlist = variables_rc[varlist_name].vars # look up variable set
  elif isinstance(varlist,(tuple,list)):
    varlist_name = varlist[0] # use first entry as name
  else: raise TypeError
  assert isinstance(varlist_name,basestring)
  # reference dataset (for meta data)
  if expens is not None and len(expens) > 0: refds = expens[0]
  elif obsens is not None and len(obsens) > 0: refds = obsens[0]
  else: raise ArgumentError
  # x-axis
  xticks = axes.xaxis.get_ticklabels() # determine if labels are appropriate
  if xlabel and len(xticks) > 0 and xticks[-1].get_visible(): xlabel = 'Seasonal Cycle [{1:s}]' 
  else: xlabel = False
  # len(xticks) > 0 is necessary to avoid errors with AxesGrid, which removes invisible tick labels      
  xlim = 0.5,12.5
  # basin info
  if shape_name or stnset_name or 'shape_name' in refds.atts:
    if not (shape_name or stnset_name): basin_name = refds.atts.shape_name
    else: basin_name = shape_name or stnset_name
    basin_info = basin_annotation[basin_name] if basin_name in basin_annotation else basin_defaults 
    # some annotation
    axes.xpad += 3
    if axtitle is not None and axtitle.lower() == 'basin': 
      axtitle = basin_info.get('title',refds.atts.shp_long_name) # basin name
    #basin_area = meands.atts.shp_area # scale by area
    # y-axis labeling
    if varlist_name.lower() == 'temp': 
      if not ylim: ylim = basin_info['temp'] if 'temp' in basin_info else None
      if ylabel is True: ylabel = '2m Temperature [{1:s}]'
    elif varlist_name in basin_info: 
      if not ylim: ylim = basin_info[varlist_name]
      if ylabel is True and varlist_name in variables_rc: 
        ylabel = variables_rc[varlist_name].label + ' [{1:s}]'
      if varlist_name in ('precip',): axes.ypad += 3
      else: axes.ypad -= 2
    else: # fallback 
      if not ylim: ylim = basin_info.get('water',None) if basin_info else None
      if ylabel is True: ylabel = 'Water Flux [{1:s}]'
      axes.ypad -= 3 
  else: ylim = ylim or None
  
  # some default kwargs
  if 'bandalpha' not in plotargs: plotargs['bandalpha'] = 0.35
  if 'errorscale' not in plotargs: plotargs['errorscale'] = 0.5
  plotargs['parasite_axes'] = dict(offset=-0.5/(len(varlist)+1.))
  #if 'linestyles' not in plotargs: plotargs['linestyles'] = ('-','--','-.')
  
  # separate fractional variables ("percentages")
  if scalevars:
    if isinstance(scalevars,basestring) and scalevars == 'auto': lauto = True
    elif isinstance(scalevars,(tuple,list)): lauto = False 
    else: raise TypeError, scalevars
    tmp_varlist = []; tmp_scalevars = []
    for var in varlist:
        if lauto:
          if expens[0][var].units: tmp_varlist.append(var)
          else: tmp_scalevars.append(var)
        else: 
          if var in scalevars: tmp_scalevars.append(var)
          else: tmp_varlist.append(var)
    varlist,scalevars = tmp_varlist,tmp_scalevars

  ## loop over datasets and generate plots
  def makePlots(ens, err, varlist=None, scalevars=None, lleg=False, lerrbar=False, lerrbnd=False, **plotargs):
    # not that some arguments are not passed explicitly but are used from the contianing scope (e.g. kwargs)
    if scalevars:
      if varlist: raise ArgumentError, varlist
      else: varlist = scalevars; lrescale = True
    else: lrescale = False
    plts = [] # meant to be added to list 
    if ens is not None and len(ens) > 0:
      if err is None: err = [None]*len(ens)      
      imaster = master if isinstance(master, (int,np.int)) else len(ens)-1 
      # prepare plotargs
      if plotargs:
        plotargs = axes._expandArgumentList(labels=[None]*len(ens), plotargs=plotargs, expand_list=[], lproduct='inner')
        assert len(plotargs) == len(ens)
      else: plotargs = [dict()]*len(ens)
      # loop over datasets
      for n,vards,errds,plotarg in zip(xrange(len(ens)),ens,err,plotargs):
        lmaster = vards.name == master if isinstance(master, basestring) else n == imaster
        plotarg.pop('label',None) # that was just a dummy to get the length right
        bards = errds if lerrbar else None # draw as many bars as datasets
        bndds = errds if lerrbnd and lmaster else None # only draw bands for last dataset
        if any(var in vards for var in varlist):
          plt = axes.linePlot(vards, varname=varlist, errorbar=bards, errorband=bndds, lrescale=lrescale, 
                              legend=legend if lmaster and lleg else None, 
                              lperi=lperi, lparasiteMeans=lparasiteMeans, title=axtitle, llabel=lmaster and lleg, 
                              xlabel=xlabel, xlim=xlim, ylabel=ylabel, ylim=ylim, lprint=lprint, **plotarg)
          plts.append(plt[0]) # use first line object for each dataset
        else: plts.append(None)
    return plts
  
  # loop over observations
  plts = [] # collect representative line objects to create secondary legend (usually for datasets)
  if ylim is not None: axes.set_ylim(ylim) 
  if varlist:
    plts += makePlots(obsens, obserr, varlist=varlist, lleg=False, lerrbar=True, lerrbnd=False)
    plts += makePlots(expens, experr, varlist=varlist, lleg=True, lerrbar=False, lerrbnd=True, **plotargs)
  if scalevars:
    plts += makePlots(obsens, obserr, scalevars=scalevars, lleg=False, lerrbar=True, lerrbnd=False)
    plts += makePlots(expens, experr, scalevars=scalevars, lleg=True, lerrbar=False, lerrbnd=True, **plotargs)  
  
  # add dataset legend (if desired)
  if dataset_legend is not False and dataset_legend is not None:
    # stow away existing legend
    if axes.legend_handle is not None: axes.add_artist(axes.legend_handle) 
    # general positioning
    if dataset_legend is True: dataset_legend = dict(loc=0) # legend at default/optimal location
    elif isinstance(dataset_legend,(int,np.integer,float,np.inexact)): dataset_legend = dict(loc=dataset_legend)
    elif not isinstance(dataset_legend,dict): raise TypeError, dataset_legend
    # dataset label list
    if dataset_labels is None:
      dataset_labels = []
      if obsens: dataset_labels += [ds.name for ds in obsens]
      if expens: dataset_labels += [ds.name for ds in expens] # not always ensembles...
      if axes.plot_labels is not None: 
        dataset_labels = [axes.plot_labels.get(label,label) for label in dataset_labels]
    elif not isinstance(dataset_labels,(list,tuple)): raise TypeError, dataset_labels
    axes.addLegend(handles=plts, labels=dataset_labels, **dataset_legend) # handles fontsize and passes kwargs to legend()
  
  # use month as labels
  axes.xaxis.set_ticklabels(['']+[name[:3] for name in name_of_month[1::2]])
  axes.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
  # only use integers for y-lables (saves space)
  if lyint:
    axes.parasite_axes.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    axes.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
  # add zero-line
  if varlist_name.lower() == 'temp': axes.addHline(273.15, alpha=0.5) # freezing point 
  else: axes.addHline(0, alpha=0.5) # actual zero (flux direction)


if __name__ == '__main__':
  
  # imports for tests
  from basins import exps_rc
  
  #   from projects.WesternCanada.WRF_experiments import Exp, WRF_exps, ensembles
  from projects.GreatLakes.WRF_experiments import Exp, WRF_exps, ensembles
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

  test = 'simple_climatology'
#   test = 'advanced_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'simple_climatology':
    
    varlist = 'runoff'; basin = 'FRB';
    exp = 'max'; exps = exps_rc[exp]; obs = 'Observations'
    period = 1979,1994
#     period = 1979,2009
    # some settings for tests
    expens = None; experr = None; obsens = None; obserr = None
    expens = loadShapeEnsemble(names=exps.exps, basins=basin, varlist=varlist, aggregation='mean',)
    experr = loadShapeEnsemble(names=exps.exps, basins=basin, varlist=varlist, aggregation='std',)
    obsens = loadShapeObservations(obs=obs, basins=basin, varlist=varlist, aggregation='mean', period=period)
    obserr = loadShapeObservations(obs=obs, basins=basin, varlist=varlist, aggregation='std', period=period)
    # print diagnostics
    print expens[0] if expens else obsens[0]; print ''
    
    # set up plot    
    fig,ax = getFigAx((1,1), title=test, sharex=True, sharey=False, stylesheet='default', lpresentation=False)    
    
    # make plots
    hydroPlot(axes=ax, expens=expens, obsens=obsens, experr=experr, obserr=obserr, varlist=varlist, 
              legend=1, dataset_legend=2, lprint=True,
              lperi=True, lparasiteMeans=True, master='max-ens', lineformats=exps.styles)
    # adjust margins
    fig.updateSubplots(left=0.015, right=0.015, top=0.0, bottom=-0.0)
    
  # test load function for basin ensemble time-series
  elif test == 'advanced_climatology':
    
    # some settings for tests
    exp = 'max'; exps = exps_rc[exp]; obs = None
    varlists = ['flux','runoff']; basins = ['FRB','ARB']
    period = 1979,1994
#     period = 1979,2009
    # some settings for tests
    expens = None; experr = None; obsens = None; obserr = None
    kwargs = dict(basins=basins, varlist=varlists, load_list=['basins','varlist'], lproduct='outer')
    expenses = loadShapeEnsemble(names=exps.exps, aggregation='mean', **kwargs)
    experres = loadShapeEnsemble(names=exps.exps, aggregation='std', **kwargs)
    obsenses = loadShapeObservations(obs=obs, aggregation='mean', period=None, **kwargs)
    obserres = loadShapeObservations(obs=obs, aggregation='SEM', period=None, **kwargs)
    # print diagnostics
    print expenses[0] if expenses else obsenses[0]; print ''
    
    # set up plot    
    fig,axes = getFigAx((len(varlists),len(basins)), title=test, sharex=True, sharey=False, stylesheet='default', lpresentation=False)    
    
    # make plots
    axes = axes.ravel()
    for n,ax,expens,experr,obsens,obserr,varlist in zip(xrange(len(axes)),axes,expenses,experres,obsenses,obserres,varlists*len(basins)):
      hydroPlot(axes=ax, expens=expens, obsens=obsens, experr=experr, obserr=obserr, varlist=varlist, master=exps.master,
                legend=0 if n>1 else None, linestyles=('--','-.','-'), lparasiteMeans=True, axtitle=None)
    
    # adjust margins
    fig.updateSubplots(left=0.01, right=-0.035, top=-0.035, bottom=-0.025)
  # show plots
  show()