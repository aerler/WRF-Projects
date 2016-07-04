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
from geodata.misc import ArgumentError
from plotting.figure import show # don't import getFigAx directly, to avoid recursion
from geodata.base import Dataset


## custom plotting functions
@BatchLoad
def climPlot(axes=None, expens=None, obsens=None, experr=None, obserr=None, varlist=None, master=None, 
              scalevars='auto', legend=0, shape_name=None, stnset_name=None, ylim=None,
              lperi=True, lparasiteMeans=True, axtitle=None, ylabel=True, xlabel=True, lyint=True,
              lprint=True, dataset_legend=False, dataset_labels=None, variable_list=None,
              annotation=None, defaults=None, **plotargs):
  ''' plot the seasonal cycle over a basin ''' # linestyles=None, lineformats=None, colors=None, markers=None
  if axes is None: raise ArgumentError
  # define some meta data
  if varlist is None and expens.basetype is Dataset: raise TypeError
  if isinstance(varlist,basestring):
    varlist_name = varlist; varlist = variable_list[varlist_name].vars # look up variable set
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
    if not (shape_name or stnset_name): shape_name = refds.atts.shape_name
    else: shape_name = shape_name or stnset_name
    shape_info = annotation[shape_name] if shape_name in annotation else defaults 
    # some annotation
    axes.xpad += 3
    if axtitle is not None and axtitle.lower() == 'basin': 
      axtitle = shape_info.get('title',refds.atts.shp_long_name) # basin name
    #shape_area = meands.atts.shp_area # scale by area
    # y-axis labeling
    if varlist_name.lower() == 'temp': 
      if not ylim: ylim = shape_info['temp'] if 'temp' in shape_info else None
      if ylabel is True: ylabel = '2m Temperature [{1:s}]'
    elif varlist_name in shape_info: 
      if not ylim: ylim = shape_info[varlist_name]
      if ylabel is True and varlist_name in variable_list: 
        ylabel = variable_list[varlist_name].label + ' [{1:s}]'
      if varlist_name in ('precip',): axes.ypad += 3
      else: axes.ypad -= 2
    else: # fallback 
      if not ylim: ylim = shape_info.get('water',None) if shape_info else None
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
  
  from projects.WesternCanada.analysis_settings import variables_rc, clim_annotation, clim_defaults
  from projects.WesternCanada.analysis_settings import climFigAx, exps_rc, loadShapeEnsemble, loadShapeObservations
#   from projects.GreatLakes.analysis_settings import variables_rc, clim_annotation, clim_defaults
#   from projects.GreatLakes.analysis_settings import climFigAx, exps_rc, loadShapeEnsemble, loadShapeObservations
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail

#   test = 'simple_climatology'
  test = 'advanced_climatology'
  
  
  # test load function for basin ensemble time-series
  if test == 'simple_climatology':
    
    varlist = 'precip'; obs = 'GPCC'
#     exp = 'g-ens'; exps = exps_rc[exp]; basin = 'GLB'
    exp = 'max-prj'; exps = exps_rc[exp]; basin = 'ARB'
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
    fig,ax = climFigAx((1,1), title=test, sharex=True, sharey=False, stylesheet='myggplot', lpresentation=False)    
    
    # make plots
    climPlot(axes=ax, expens=expens, obsens=obsens, experr=experr, obserr=obserr, varlist=varlist, 
              legend=2, dataset_legend=4, lprint=True, variable_list=variables_rc, 
              annotation=clim_annotation, defaults=clim_defaults,
              lperi=True, lparasiteMeans=True, master=exps.master, lineformats=exps.styles)
    # adjust margins
    fig.updateSubplots(left=0.02, right=0.015, top=0.0, bottom=-0.0)
    
  # test load function for basin ensemble time-series
  elif test == 'advanced_climatology':
    
    # some settings for tests
#     exp = 'g-ens'; exps = exps_rc[exp]; basins = ['GLB']
    exp = 'max-prj'; exps = exps_rc[exp]; basins = ['ARB']
    obs = 'GPCC'; varlists = ['flux_snow','precip']
    period = 1979,1994
#     period = 1979,2009
    # some settings for tests
    expens = None; experr = None; obsens = None; obserr = None
    kwargs = dict(basins=basins, varlist=varlists, load_list=['basins','varlist'], lproduct='outer')
    expenses = loadShapeEnsemble(names=exps.exps, aggregation='mean', **kwargs)
    experres = loadShapeEnsemble(names=exps.exps, aggregation='SEM', **kwargs)
    obsenses = loadShapeObservations(obs=obs, aggregation='mean', period=None, **kwargs)
    obserres = loadShapeObservations(obs=obs, aggregation='SEM', period=None, **kwargs)
    # print diagnostics
    print expenses[0] if expenses else obsenses[0]; print ''
    
    # set up plot    
    fig,axes = climFigAx((len(varlists),len(basins)), title=test, sharex=True, sharey=False, stylesheet='default', lpresentation=False)    
    
    # make plots
    axes = axes.ravel()
    for n,ax,expens,experr,obsens,obserr,varlist in zip(xrange(len(axes)),axes,expenses,experres,obsenses,obserres,varlists*len(basins)):
      climPlot(axes=ax, expens=expens, obsens=obsens, experr=experr, obserr=obserr, varlist=varlist, 
                master=exps.master, variable_list=variables_rc, 
                annotation=clim_annotation, defaults=clim_defaults,
                legend=0 if n>1 else None, linestyles=('--','-.','-'), lparasiteMeans=True, axtitle=None)
    
    # adjust margins
    fig.updateSubplots(left=0.01, right=-0.04, top=+0.06, bottom=-0.03)
  # show plots
  show()