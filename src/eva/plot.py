'''
Created on Feb 2, 2015

Utility module to provide specialized plotting functionality for Extreme Value Analysis.

@author: Andre R. Erler, GPL v3
'''

# external imports
import matplotlib as mpl
# internal imports
import numpy as np
from geodata.base import Ensemble, Dataset, Variable, Axis, concatVars
from geodata.misc import ArgumentError, AxisError, DatasetError
from geodata.stats import ks_2samp, VarRV
from plotting.figure import show # don't import getFigAx directly, to avoid recursion
from plotting.axes import checkVarlist
from eva.load import _rescaleSample

def stationInfo(stnds, varname, name, titlestr=None, alttitle=None, lflatten=False, lmon=False,):
  ''' helper to generate an axes title with station info '''
  if stnds.hasAxis('station'): nstn = len(stnds.axes['station']) # number of stations        
  else: nstn = 1 # single station
  if stnds.name[:3].lower() == 'obs' and varname in stnds:
      ec = stnds[varname] # some variables are not present everywhere
      if ec.hasAxis('time') and ec.time.units[:3].lower() == 'mon': units = 'mon.'
      elif ec.hasAxis('year') and ec.year.units.lower() == 'year': units = 'yrs.'
      else: units = 'mon.' if lmon else 'yrs.'
      mask = ec.data_array.mask if isinstance(ec.data_array,np.ma.MaskedArray) else np.isnan(ec.data_array) 
      if lflatten: rec_len = (ec.data_array.size - mask.sum()) # valid years in obs/EC
      else: rec_len = int(np.round(ec.data_array.shape[-1] - mask.sum(axis=-1).mean())) # valid years in obs/EC
      if titlestr: axtitle = titlestr.format(name,nstn,rec_len) # axes label
      else: axtitle = "{:s} (#{:d}, {:d} {:s})".format(name,nstn,rec_len,units) # axes label
  else:
      if alttitle: axtitle = alttitle.format(name,nstn) # axes label
      elif titlestr: axtitle = titlestr.format(name,nstn) # axes label
      else: axtitle = "{:s} (#{:d}, WRF only)".format(name,nstn) # axes label
  return axtitle

# function to compute some statistics and print them
def generateStatistics(varname, ens, fit, scl=None, reference=None, mode='Ratio', plot_labels=None, 
                       nsamples=None, bootstrap_axis='bootstrap', lflatten=False, sample_axis='time', 
                       lcrossval=True):
  ''' Perform K-S test and compute ratio of means; return results in formatted string. '''
  # some average diagnosics
  idkey = 'dataset_name' if ens.basetype is Dataset else 'name'  
  varlist = Ensemble(*[ds[varname] for ds in ens if ds is not None and varname in ds], idkey=idkey)
  if not all(varlist[0].ndim==ndim for ndim in varlist.ndim):
    new_axes = varlist[np.argmax(varlist.ndim)].axes
    varlist = varlist.insertAxes(new_axes=new_axes, lcheckAxis=False)    
  mvars = varlist.mean() # growth rate
  lratio = mode.lower() == 'ratio'
  lshift = mode.lower() == 'shift'
  if plot_labels is None: plot_labels = dict()
  # figure out fillValue
  if np.issubdtype(varlist[0].dtype, np.floating): fillValue = np.NaN
  elif np.issubdtype(varlist[0].dtype, np.integer): fillValue = 0
  else: raise TypeError(varlist[0].dtype)
  # define reference
  if isinstance(reference,(list,tuple)): 
    reflist0 = list(reference); reference = reference[0]
  else: reflist0 = [] # dummy list
  if reference is None: iref0 = 0
  elif isinstance(reference,(int,np.integer)): iref0 = reference 
  elif isinstance(reference,str): iref0 = varlist.idkeys.index(reference)
  else: raise ArgumentError  
  # goodness of fit, reported on plot panels
  if fit:
    fitlist = Ensemble(*[ds[varname] for ds in fit if ds is not None and varname in ds], idkey=idkey)
    if any(fitlist.hasAxis(bootstrap_axis)): fitlist = fitlist(**{bootstrap_axis:0, 'lcheckAxis':False})
    if not all(fitlist[0].ndim==ndim for ndim in fitlist.ndim):
      new_axes = fitlist[np.argmax(fitlist.ndim)].axes
      fitlist = fitlist.insertAxes(new_axes=new_axes, lcheckAxis=False) 
#       for var in fitlist: 
#         print [ax.name for ax in var.axes], var.shape
#       assert  np.all(fitlist[0][1,:] == fitlist[0][2,:])
    assert not isinstance(reference,str) or iref0 == fitlist.idkeys.index(reference), reference
    if any([isinstance(dist,VarRV) for dist in fitlist]) or not scl:
      names = [plot_labels.get(getattr(dist,idkey),getattr(dist,idkey)) for dist in fitlist]  
      lnames = max([len(name) for name in names]) # allocate line space
      headline = 'Sample'; lhead = len(headline) # sample/exp header
      headline += ' '*max(lnames-lhead,0) # 'Exp.'+' '*max(lnames-4,0) if lnames < 8 else 'Experiment'
      string = '{:s}  Fit  {:s}\n'.format(headline,mode.title())
      namestr = '{{:>{:d}s}}  {{:s}}  '.format(max(lhead,lnames))
      iref = iref0; reflist = reflist0[:] # copy list
      for i,dist,var,name,mvar in zip(range(len(fitlist)),fitlist,varlist,names,mvars):
        if isinstance(dist,VarRV) or not scl:
          if isinstance(dist,VarRV):
            pval = dist.fittest(var, nsamples=nsamples, asVar=False, lcrossval=lcrossval) #lflatten=lflatten, axis_idx=var.axisIndex(sample_axis, lcheck=False))
#             print var.name, pval, pval.mean().__class__.__name__, '{:s}'.format(pval.mean())
#             pval = '{:3.2f}'.format(float(pval.mean())) # mean is only necessary to convert to scalar
            pval = '{:3.2f}'.format(float(np.median(pval))) # mean is only necessary to convert to scalar
            # for some reason masked array scalars appear string-type, rather than numbers... 
          else: pval = '  - '
          if len(reflist) > 0 and name == reflist[0]: # assign new reference 
            iref = i; del reflist[0] # pop element 
          if isinstance(mvar,np.ma.core.MaskedConstant) or isinstance(mvars[iref],np.ma.core.MaskedConstant): 
            string += namestr.format(name,' N/A\n')
          elif lratio: string += (namestr+'{:3.2f}\n').format(name,pval,(mvar/mvars[iref]).mean())
          elif lshift: string += (namestr+'{:+2.1f}\n').format(name,pval,(mvar-mvars[iref]).mean())
    else: string = ''
  else: raise NotImplementedError
  if scl:
    scllist = Ensemble(*[ds[varname] for ds in scl if ds is not None and varname in ds], idkey=idkey)
    bs_axes = scllist.axisIndex(bootstrap_axis, lcheck=False) # return None, if not present
    if bs_axes is None: bs_axes = [None]*len(scllist)
    scllist = scllist(**{bootstrap_axis:0, 'lcheckAxis':False})
    if not all(scllist[0].ndim==ndim for ndim in scllist.ndim):
      new_axes = scllist[np.argmax(scllist.ndim)].axes
      scllist = scllist.insertAxes(new_axes=new_axes, lcheckAxis=False) 
    assert not isinstance(reference,str) or iref0 == scllist.idkeys.index(reference), reference
    if len(scllist) != len(varlist): raise AxisError(scllist)
    # compute means
    mvars = []
    for svr,var in zip(scllist,varlist):
      if isinstance(svr,VarRV): mvar = svr.stats(moments='mv', asVar=False)[...,0] # only first moment
      else: mvar = var.mean()*svr.atts.get('loc_factor',1.)
      mvars.append(mvar)        
    # figure out label width and prepare header
    if len(varlist) > 1: # otherwise no comparison...
      names = [plot_labels.get(getattr(dist,idkey),getattr(dist,idkey)) for dist in scllist]  
      lnames = max([len(name) for name in names]) # allocate line space
      namestr = '{{:>{:d}s}}  {{:s}}  '.format(max(lhead,lnames))
      tmphead = 'Fit to {:s}:' if scl == fit else 'Rescaled to {:s}:' # new heading
      tmphead += ' '*(max(lnames-len(names[iref0]),0)+5)+'\n'
      string += tmphead.format(names[iref0])
      # prepare first reference sample for K-S test
      scale,shape = scllist[iref0].atts.get('scale_factor', 1),scllist[iref0].atts.get('shape_factor', 1)
      if not (scale is None or scale == 1) and not (shape is None or shape == 1): 
        raise NotImplementedError("Cannot rescale scale/variance and shape parameters of reference sample!")
      refsmpl = varlist[iref0].getArray(unmask=True, fillValue=fillValue) # only once
      loc0 = scllist[iref0].atts.get('loc_factor', 1)     
      refsmpl = _rescaleSample(refsmpl, loc0, bs_axis=bs_axes[iref0]) # apply rescaling (varies, dependign on loc-type)
  #     print varlist[iref0].dataset_name, [ax.name for ax in varlist[iref0].axes], refsmpl.shape, 
      # start loop
      iref = iref0; reflist = reflist0[:] # copy list
      for i,dist,varsmpl,mvar,bs_axis in zip(range(len(varlist)),scllist,varlist,mvars,bs_axes):
        name = getattr(dist,idkey)
        if len(reflist) > 0 and name == reflist[0]: # assign new reference 
          iref = i; del reflist[0] # pop element       
          # prepare subsequent reference sample for K-S test
          scale,shape = dist.atts.get('scale_factor', 1),dist.atts.get('shape_factor', 1)
          if not (scale is None or scale == 1) and not (shape is None or shape == 1): 
            raise NotImplementedError("Cannot rescale scale/variance and shape parameters of reference sample!")
          refsmpl = varsmpl.getArray(unmask=True, fillValue=fillValue) # only once
          if not varsmpl.atts.get('rescaled',False):
            refsmpl = _rescaleSample(refsmpl, dist.atts.get('loc_factor', 1), bs_axis=bs_axis) # apply rescaling (varies, dependign on loc-type)
        elif i != iref:
          scale,shape = dist.atts.get('scale_factor', 1),dist.atts.get('shape_factor', 1) 
          # perform K-S test
          if (scale is None or scale == 1) and (shape is None or shape == 1):
            # K-S test between actual samples is more realistic, and rescaling of mean is simple
            smpl = varsmpl.getArray(unmask=True, fillValue=fillValue) # only once
            if not varsmpl.atts.get('rescaled',False):
              smpl = _rescaleSample(smpl, dist.atts.get('loc_factor', 1), bs_axis=bs_axis) # apply rescaling (varies, dependign on loc-type)
  #           print varsmpl.dataset_name, [ax.name for ax in varsmpl.axes], smpl.shape
  #           print smpl.shape, np.nanmean(smpl), refsmpl.shape, np.nanmean(refsmpl)
  #           print lflatten, sample_axis
            pval = ks_2samp(refsmpl, smpl, asVar=False, lflatten=lflatten, 
                            axis_idx=varsmpl.axisIndex(sample_axis, lcheck=False))
  #           print dist.name, pval
  #           pval = '{:3.2f}'.format(float(pval.mean()))
            pval = '{:3.2f}'.format(float(np.median(pval)))
          else:
            # no straight-forward way to rescale samples, so have to compare distribution with 
            # reference sample, which means more noise (since the distribution will be randomly sampled)
            if isinstance(dist,VarRV): pval = '{:3.2f}'.format(float(dist.kstest(refsmpl).mean()))
            else: pval = '  - '
          # add column with ratio/difference of means after rescaling
          if name in plot_labels: name = plot_labels[name]  
          if isinstance(mvar,np.ma.core.MaskedConstant) or isinstance(mvars[iref],np.ma.core.MaskedConstant):
            string += namestr.format(name,' N/A\n') 
          elif lratio: string += (namestr+'{:3.2f}\n').format(name,pval,(mvar/mvars[iref]).mean())
          elif lshift: string += (namestr+'{:+2.1f}\n').format(name,pval,(mvar-mvars[iref]).mean())
  # return formatted table in string
  return string        

## custom default plot styles

  
## plotting functions

def selectDataset(ens, datasets):
  ''' helper function to extract specific datasets '''
  if ens is None:
    return None
  if datasets is not None:
    members = [] 
    for name in datasets:
      if name in ens:
        member = ens[name]
        members.append(member)
  else: members = ens[:]
  enscp = Ensemble(members=members, idkey=ens.idkey, basetype=ens.basetype, name=ens.ens_name, title=ens.ens_title)
  return enscp

# simple function to generate a single annotated plot
def distPlot(axes=None, varname=None, datasets=None, ens=None, fit=None, kde=None, scl=None, legend=True,
             shape_name=None, stnset_name=None, lsmall=False, axtitle=None, ylabel=False, xlabel=True,
             lbootstrap=False, bootstrap_axis='bootstrap', lsample=False, sample_axis='station',
             lrescale=None, reference=None, master=None, rescale_line='--', lrescaleBand=False, band_vars=None,
             xlim=None, ylim=None, lanno='right', xfactor=None, lthinx=True, lcrossval=False,
             percentiles=(0.025,0.975), lmedian=None, median_fmt='', lmean=None, mean_fmt='', 
             lvar=False, lvarBand=False, lsmooth=None, lprint=False, mode='ratio', 
             annotation=None, defaults=None, **kwargs):
  ''' Generate an EVA plot with annotations in a single axes; the EVA plot can contain histograms of the
      data, best fit curves (EVA), (re-)scaled curves, empiraically estimated fits (KDE), and text
      annotation indicating Goodness-ofFit (K-S test) and statistical moments. '''
  
  # figure out rescaling in a sensible manner
  lrescale = ( scl and any([varname in ds for ds in scl]) ) if lrescale is None else lrescale
  
  # get reasonable x-limits for variables based on region
  if xlim is None:
    if stnset_name is not None: xlim = annotation[stnset_name].get(varname,None)
    elif shape_name is not None: xlim = annotation[shape_name].get(varname,None)
    elif varname in defaults: xlim = defaults.get(varname,None)  
    else: xlim = axes.get_xlim()
    if stnset_name is not None and xfactor is None and (lrescale is False or not scl): 
      xfactor = 2./3. # empirical value...
  if xlim is None: 
    raise ArgumentError("Variable '{:s}' not found in station/shape settings '{:s}'.".format(varname,stnset_name or shape_name))
  # rescaling, after we found a good setting
  if xfactor is not None and xfactor is not False and xfactor != 1:
    xlim = [np.round(xl*xfactor,int(1-np.round(np.log10(xlim[1])))) for xl in xlim]
  
  # select datasets
  if datasets is not None:
    if not isinstance(datasets,(tuple,list)): raise TypeError(datasets)
    ens = selectDataset(ens, datasets)  
    fit = selectDataset(fit, datasets)
    kde = selectDataset(kde, datasets)
    scl = selectDataset(scl, datasets)

  ## handle annotation before any further modification to variables
  # legend position
  if not lanno: ann = None; leg = 0
  elif lanno is True or lanno.lower() == 'right': 
    if not isinstance(fit[-1][varname], VarRV): leg = 5 # right
    elif len(fit) < 4 or ( len(fit) <= 4 and (varname not in fit[0]) ): leg = 5
    elif len(fit) < 6: leg = dict(bbox_to_anchor=(1,0.2), loc='lower right', borderaxespad=1)
    else: leg = 4    
    ann = 1
  elif lanno.lower() == 'left': 
    leg = 3 if lsmall and scl and len(ens)>3 else 6
    ann = 2
  else: ann = None; leg = 0
  if legend is True: legend = leg
  elif isinstance(legend,dict): 
    if isinstance(leg,dict):
      leg.update(legend); legend = leg
    elif 'loc' not in legend: legend['loc'] = leg 
  elif legend is False: legend = None
  if legend is not None and not isinstance(legend, (int,np.integer,dict)):
    raise ArgumentError  
  # create annotation string    
  if lanno and ens and fit:
    # find a dataset that actually has the variables and get it
    mvar = fit[np.argmax(fit.hasVariable(varname))][varname]
    # retrieve meta data from variables
    lflatstats = mvar.atts.sample_axis == 'flattened'
    stats_axis = None if lflatstats else mvar.atts.sample_axis
    # determine if merging is necessary
    if not lflatstats and not all(ens.hasAxis(stats_axis)): raise DatasetError
    # generate statistics string
    string = generateStatistics(varname, ens, fit, scl=scl if lrescale else fit, 
                                lflatten=lflatstats, sample_axis=stats_axis, lcrossval=lcrossval,
                                mode=mode, reference=master or reference, plot_labels=axes.plot_labels)
    axes.addLabel(string, loc=ann, prop=dict(size='small', fontname='monospace'), borderpad=0.5)
  else: 
    string = ''
  
  # isolate variable for meta data
  for ds in ens or fit:
    if ds is None: pass
    elif isinstance(ds,Dataset) and varname in ds: var = ds[varname]; break
    elif isinstance(ds,Variable) and varname in ds.name: var = ds; varname = None; break
  
  # make histogram
  dx = (xlim[1]-xlim[0])//10
  if ens:
    axes.histogram(ens, varname=varname, bins=(xlim[0],xlim[1]-dx,20,), alpha=0.5, rwidth=0.7, lprint=lprint,
                   llabel=not bool(fit), xlabel=False if bool(fit) else xlabel, xlim=xlim, linewidth=0,
                   ylabel=ylabel, yticks=bool(ylabel),ylim=ylim, title=axtitle, **kwargs)
  
  # plot distributions
  bins = (xlim[0]/var.plot.scalefactor,xlim[1]/var.plot.scalefactor,200)
  
  # empirical fit (KDE)
  if kde:
    axes.linePlot(kde, varname=varname, support=bins, linestyle='--', alpha=0.8, linewidth=.75, llabel=False)
  
  if lsample and lmean is None: lmean = True # show mean instead of "true fit" 
  if lsample and lmedian is None: lmedian = False # don't show median
  # distribution fit
  if fit:    
    if lsample and any(ds.hasAxis(sample_axis, lany=True) for ds in fit): 
      axes.samplePlot(fit, varname=varname, support=bins, yticks=False, xticks=True, ylabel=ylabel, xlabel=xlabel, 
                      legend=legend, sample_axis=sample_axis, bootstrap_axis=bootstrap_axis,
                      linewidth=1.25, xlim=xlim, ylim=ylim, percentiles=percentiles, llabel=True,   
                      lmedian=lmedian, median_fmt=median_fmt, lmean=lmean, mean_fmt=mean_fmt, 
                      lsmooth=lsmooth, bandalpha=0.25, band_vars=band_vars, **kwargs)
    elif lbootstrap and any(ds.hasAxis(bootstrap_axis) for ds in fit): 
      axes.bootPlot(fit, varname=varname, support=bins, yticks=False, ylabel=ylabel, xlabel=xlabel, legend=legend,
                    bootstrap_axis=bootstrap_axis, band_vars=band_vars, 
                    linewidth=1.25, xlim=xlim, ylim=ylim, percentiles=percentiles, llabel=True,   
                    lmedian=lmedian, median_fmt=median_fmt, lmean=lmean, mean_fmt=mean_fmt, 
                    lvar=lvar, lvarBand=lvarBand, lsmooth=lsmooth, bandalpha=0.25, **kwargs)
    else: axes.linePlot(fit, varname=varname, support=bins, yticks=False, ylabel=ylabel, xlabel=xlabel, 
                        linewidth=1.25, llabel=True, legend=legend, xlim=xlim, ylim=ylim, **kwargs)
  
  # plot scaled distributions
  if lrescale:
    if lsample and any(ds.hasAxis(sample_axis, lany=True) for ds in scl):
      axes.samplePlot(scl, varname=varname, support=bins, linestyle=rescale_line, bandalpha=0.25, linewidth=1.5,
                      llabel=False, percentiles=percentiles if lrescaleBand else None, band_vars=band_vars, 
                      lmedian=lmedian, median_fmt=median_fmt, bootstrap_axis=bootstrap_axis, 
                      lmean=lmean, xlabel=xlabel, ylabel=False, xticks=True, yticks=False,
                      sample_axis=sample_axis, mean_fmt=mean_fmt, lsmooth=lsmooth, **kwargs)
    elif lbootstrap and lrescaleBand and any(ds.hasAxis(bootstrap_axis) for ds in scl):
#       bandalf = kwargs.pop('bandalhpa',0.5) * 0.7
      axes.bootPlot(scl, varname=varname, support=bins, linestyle=rescale_line, bandalpha=0.25, linewidth=1.5,
                    bootstrap_axis=bootstrap_axis, band_vars=band_vars, 
                    llabel=False, percentiles=percentiles, lmedian=lmedian, median_fmt=median_fmt, 
                    lmean=lmean, xlabel=xlabel, ylabel=ylabel, yticks=False,
                    mean_fmt=mean_fmt, lvar=lvar, lvarBand=lvarBand, lsmooth=lsmooth, **kwargs) 
    else: 
      axes.linePlot(scl, varname=varname, support=bins, linestyle=rescale_line, alpha=1, linewidth=1.5, llabel=False, 
                    xlabel=xlabel, ylabel=ylabel, yticks=False, **kwargs)                
  
  # make sure we don't waste any space
  axes.set_ylim(-0.1*axes.get_ylim()[1],axes.get_ylim()[1])  
  axes.addHline(0, color=mpl.rcParams['grid.color'], alpha=0.5) # mark zero line
  
  # thinning of grid labels
  if lthinx:
    xtls = axes.xaxis.get_ticklabels(); s = len(xtls)//10 +1; r = len(xtls)%s
    for z,xtl in enumerate(xtls): 
        if (z-r+1)%s != 0: xtl.set_visible(False) 
        #else: xtl.set_visible(True) 

  # return annotation string
  return string

  
# simple function to generate a single annotated plot
def quantPlot(axes=None, varname=None, datasets=None, fit=None, scl=None, legend=True, quantiles=None, 
              shape_name=None, stnset_name=None, ylim=None, xlim=None, ylabel=True, xlabel=True, 
              xticks=True, yticks=True, axtitle=None, lbootstrap=False, percentiles=(0.025,0.975), 
              band_vars=None,
              lmedian=None, median_fmt='', lrescaleBand=False, sampling_period=1, lmean=None, mean_fmt='', 
              lvar=False, lvarBand=False, lrescale=False, reference=0, bootstrap_axis='bootstrap', 
              lsample=False, sample_axis='station', annotation=None, defaults=None, **kwargs):
  ''' Generate an EVA quantile plot with annotations in a single axes, including (re-)scaled curves. 
      The y-axis is labeled with return periods and specific quantiles can be marked by vertical lines. '''
  
  # cast quantiles as array
  quantiles = np.asarray(quantiles) if quantiles else None
  # y-limits for upper and lower percentiles
  if ylim is None:
    if quantiles is not None:
      if np.min(quantiles) > 0.5: ylim = (0.95,1.01)
      else: ylim = (-0.01,0.065)
    else: ylim = (0.95,1.01)
  bins = ylim+(200,)
  if ylim[1] >= 1: axes.addHline(1, color='k', alpha=0.5) # mark one line
  if ylim[0] <= 0: axes.addHline(0, color='k', alpha=0.5) # mark zero line
  # get reasonable x-limits for variables based on region
  if xlim is None:
    if stnset_name is not None: xlim = annotation[stnset_name].get(varname,None)
    elif shape_name is not None: xlim = annotation[shape_name].get(varname,None)
    elif varname in defaults: xlim = defaults.get(varname,None)    
    else: xlim = axes.get_xlim()
  # select datasets
  if datasets is not None:
    if not isinstance(datasets,(tuple,list)): raise TypeError(datasets)
    fit = selectDataset(fit, datasets)
    scl = selectDataset(scl, datasets)
    
  if lsample and lmean is None: lmean = True # show mean instead of "true fit" 
  if lsample and lmedian is None: lmedian = False # don't show median  
  # generate quantile plot
  if lsample and any(ds.hasAxis(sample_axis, lany=True) for ds in fit):    
    plts = axes.samplePlot(fit,varname, support=bins, method='ppf', linewidth=1.25, title=axtitle,
                         flipxy=True, llabel=True, legend=legend, percentiles=percentiles,
                         sample_axis=sample_axis, bootstrap_axis=bootstrap_axis, band_vars=band_vars,
                         xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks, 
                         lmedian=lmedian, median_fmt=median_fmt, lmean=lmean, mean_fmt=mean_fmt,
                         lsmooth=False, bandalpha=0.25, **kwargs)    
  elif lbootstrap and any(ds.hasAxis(bootstrap_axis) for ds in fit):    
    plts = axes.bootPlot(fit,varname, support=bins, method='ppf', linewidth=1.25, title=axtitle,
                       flipxy=True, llabel=True, legend=legend, percentiles=percentiles, band_vars=band_vars,
                       xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks, 
                       lmedian=lmedian, median_fmt=median_fmt, lmean=lmean, mean_fmt=mean_fmt, 
                       lvar=lvar, lvarBand=lvarBand, lsmooth=False, bandalpha=0.25, **kwargs)
    # N.B.: smoothing causes (resolution-dependent) artifacts near the limit (1 or 0)
  else:
    plts = axes.linePlot(fit, varname, support=bins, method='ppf', linewidth=1.25, title=axtitle, 
                       xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks,
                       flipxy=True, llabel=True, legend=legend, **kwargs)
  # add percentile lines
  if quantiles is not None:    
    iref = fit.idkeys.index(reference) if isinstance(reference,str) else reference
    if lsample: 
      quants = checkVarlist(fit[iref], varname=varname, ndim=2, support=quantiles, method='ppf')[0].mean(asVar=True)      
    else: quants = checkVarlist(fit[iref], varname=varname, ndim=1, support=quantiles, method='ppf')[0]
    qvals = quants.data_array * quants.plot.scalefactor + quants.plot.offset
    axes.addVline(qvals, color=plts[iref].get_color(), alpha=0.5) # linestyle=plts[iref].get_linestyle()

  if lrescale and scl:
    # plot scaled quantile distributions
    if lsample and any(ds.hasAxis(sample_axis, lany=True) for ds in scl):    
      plts = axes.samplePlot(scl,varname, support=bins, method='ppf', linewidth=1., mean_fmt='--',
                           title=None, flipxy=True, llabel=False, legend=legend, 
                           percentiles=percentiles if lrescaleBand else None, band_vars=band_vars,
                           sample_axis=sample_axis, bootstrap_axis=bootstrap_axis, 
                           xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks, 
                           lmedian=lmedian, median_fmt='--', lmean=lmean, #mean_fmt=mean_fmt,
                           lsmooth=False, bandalpha=0.25, **kwargs)    
    elif lbootstrap and any(ds.hasAxis(bootstrap_axis) for ds in scl):    
#       bandalf = kwargs.pop('bandalhpa',0.5) * 0.7
      plts = axes.bootPlot(scl,varname, support=bins, method='ppf', percentiles=percentiles if lrescaleBand else None, 
                         linewidth=1., lineformat='--', bandalpha=0.25, flipxy=True, lmedian=lmedian, median_fmt='--',
                         xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks, 
                         lmean=lmean, mean_fmt='--', lvar=lvar, lvarBand=lvarBand, band_vars=band_vars,
                         llabel=False, lsmooth=False, **kwargs)
      # N.B.: smoothing causes (resolution-dependent) artifacts near the limit (1 or 0)
    else:
      plts = axes.linePlot(scl, varname, support=bins, method='ppf', linewidth=1., linestyle='--',
                         xlim=xlim, ylim=ylim, ylabel=ylabel, xlabel=xlabel, xticks=xticks, yticks=yticks,
                         flipxy=True, llabel=False, legend=None, **kwargs)
    # add percentile lines for scaled distributions
    if quantiles is not None:
      iref = scl.idkeys.index(reference) if isinstance(reference,str) else reference    
      if lsample: 
        quants = checkVarlist(scl[iref], varname=varname, ndim=2, support=quantiles, method='ppf')[0].mean(asVar=True)
      else: quants = checkVarlist(scl[iref], varname=varname, ndim=1, support=quantiles, method='ppf')[0]
      qvals = quants.data_array * quants.plot.scalefactor + quants.plot.offset    
      axes.addVline(qvals, color=plts[iref].get_color(), alpha=0.35) # linestyle=plts[iref].get_linestyle()
      
  # convert y-tick labels to return periods
  yticks = ['{:3.0f}'.format(sampling_period/(1-yt)) if yt < 1 else '' for yt in axes.get_yticks().tolist()]
  axes.set_yticklabels(yticks)
   
  # return
  return plts

# helper function to convert quantiles
def convQuant(q, smpl_prd=None, lround=True):
  ''' convert quantiles to return periods if sampling period is specified, else return qunatile '''
  if smpl_prd:
    r = smpl_prd/(1.-q)
    if lround: r = np.round(r) 
  else: r = q
  return r 

# simple function to compute return periods relative to reference
def quantPeriod(varname=None, datasets=None, fit=None, quantiles=None, name=None, units='years', 
                sampling_period=1, asVar=True,
                lmean=None, lmedian=None, lbootstrap=False, percentiles=(0.025,0.975), lvar=False, 
                reference=0, bootstrap_axis='bootstrap', lsample=False, sample_axis='station', **kwargs):
  ''' Compute the return periods of events with intensities equal to the given quantiles in the 
      reference sample. Confidence intervals can be estimated from bootstrapping or sample 
      distributions. The data is returned as a variable object with axes for quantiles, datasets, 
      and distribution percentiles; quantiles are converted to return periods (years, by default). '''
  
  # input requirements
  if not isinstance(fit, Ensemble): raise TypeError(fit)
  if not all([isinstance(ds, Dataset) for ds in fit]): raise TypeError(fit)
  # cast quantiles as array
  quantiles = np.asarray(quantiles)
  
  # select datasets
  if datasets is not None:
    if not isinstance(datasets,(tuple,list)): raise TypeError(datasets)
    fit = selectDataset(fit, datasets)
    assert datasets == fit.name
  else: 
    datasets = fit.name
      
  if lsample and lmean is None: lmean = True # show mean instead of "true fit" 
  if lsample and lmedian is None: lmedian = False # don't show median 
  
  # compute references quantiles
  if isinstance(reference,(int,np.integer)): reference = fit[reference].name
  if not isinstance(reference,str): raise TypeError(reference)  
  if lsample: 
    refvals = checkVarlist(fit[reference], varname=varname, ndim=2, support=quantiles, method='ppf')[0]
    refvals = refvals.mean(axis=sample_axis, asVar=False)      
  else: 
    refvals = checkVarlist(fit[reference], varname=varname, ndim=1, support=quantiles, method='ppf')[0][:] # just get values

  # compute qunatiles corresponding to reference quantiles
  if ( lsample and any(ds.hasAxis(sample_axis, lany=True) for ds in fit) or 
       lbootstrap and any(ds.hasAxis(bootstrap_axis) for ds in fit) ):
    ciaxis = bootstrap_axis if lbootstrap else sample_axis
    # compute confidence itnervals    
    civars = checkVarlist(fit, varname=varname, ndim=2, support=refvals, method='cdf', bootstrap_axis=None)
    q = (percentiles[0],0.5,percentiles[1]) if lmedian and lsample else percentiles # median is a percentile
    if lmean: qvars = [civar.mean(axis=sample_axis, asVar=True,) for civar in civars]
    civars = [civar.percentile(q=q, axis=ciaxis, asVar=True,) for civar in civars]
    if lmedian: qvars = [civar(percentile=0.5,) for civar in civars]
    lovars = [civar(percentile=percentiles[0],) for civar in civars]
    hivars = [civar(percentile=percentiles[1],) for civar in civars]
  if lbootstrap or not lsample:
    # get the one true realization    
    qvars = checkVarlist(fit, varname=varname, ndim=1, support=refvals, method='cdf', bootstrap_axis=bootstrap_axis)
  # returns a list of variables containing the quantile values
  if lbootstrap or lsample:
    # prepend the true/mean value to CI intervals (percentile axis)
    pax = Axis(coord=np.asarray((-1,)+percentiles), name='percentile', units='')
    tmplist = []
    for qvar,lovar,hivar in zip(qvars,lovars,hivars):
      tmp = concatVars((qvar,lovar,hivar), axis=pax, asVar=True, name=qvar.name, units=qvar.units, 
                       varatts=qvar.atts, lcheckAxis=True, lensembleAxis=True)
      tmplist.append(tmp)
    qvars = tmplist # new qvar list with percentiles / CI's
    
  # concatenate experiments by creating new axis
  if asVar:
    # create experiment axis
    expax = Axis(coord=np.arange(len(datasets)), name='experiment', units='#', atts=dict(datasets=datasets))
    # concat variables
    name = name or varname+'_periods'
    quantvar = concatVars(qvars, axis=expax, asVar=True, name=name, units=units, 
                           lcheckAxis=True, lensembleAxis=True)
    quantvar.data_array = convQuant(quantvar.data_array, smpl_prd=sampling_period, lround=True)
    # rename "_bins" variable
    for axis in quantvar.axes:
      if axis.name[-5:] == '_bins':
        axis.coord = convQuant(quantiles, smpl_prd=sampling_period, lround=True)
        axis.name = 'quantile'; axis.units = units
  else:
    quantvar = {dataset:qvar for dataset,qvar in zip(datasets,qvars)}
    
  # return
  # N.B.: the true/mean/median value is the first one along the percentile axis (percentile = -1)
  return quantvar


## abuse main section for testing
if __name__ == '__main__':
    
  from projects.WesternCanada.analysis_settings import dist_annotation, dist_defaults, quant_annotation, quant_defaults
  from projects.WesternCanada.analysis_settings import evaFigAx, exps_rc, loadStationFit, variables_rc
  # N.B.: importing Exp through WRF_experiments is necessary, otherwise some isinstance() calls fail
  
  # reduce sample size
  constraints = dict()
  constraints['min_len'] = 0 
  constraints['lat'] = (45,60) 
  constraints['max_zerr'] = 300
  
  test = 'annotated_plot'
#   test = 'quantile_plot'
#   test = 'quantile_table'
  
  # test quantile table
  if test == 'quantile_table':

    # some settings for tests
    exps = ['EC','max-ens','max-ens-2050','max-ens-2100',][:]
    prov = 'AB'; season = 'winter'
    varname = 'MaxPrecip_1d'; filetypes = ['hydro']
    lflatten = True; lfit = True; lrescale = True; lbootstrap = True
    sample_axis = 'station'; lsample = not lflatten
    # load some data
    stnens,fitens,sclens = loadStationFit(names=exps, provs=prov, seasons=season, varlist=[varname], lshort=False,
                                          lrescale=lrescale, reference='EC', target='max-ens', lfit=lfit,  
                                          stationtype='ecprecip', lflatten=lflatten, filetypes=filetypes,
                                          lbootstrap=lbootstrap, nbs=30)
    
    # compute quantiles    
    qvar = quantPeriod(varname=varname, datrasets=exps, fit=sclens[0], quantiles=(0.98,0.99), sampling_period=1., 
                       reference='max-ens', percentiles=(0.025,0.975), lmedian=None, lmean=None, 
                       sample_axis=sample_axis, stnset_name=prov, lbootstrap=lbootstrap, lsample=lsample,)
    var = qvar(percentile=-1, lidx=False, )
    print(var)
    print(var.data_array)

  # test annotated distribution plots
  elif test == 'annotated_plot':

    # some settings for tests
    clusters = None; provs = None; seasons = None
#     exps = ['EC', 'ctrl-ens', 'ctrl-ens-2050', 'ctrl-ens-2100'][:]
#     exps = ['EC', 'max-3km_d01', 'max-3km_d02', 'max-3km_d03',]
    exps = ['EC', 'max-3km', 'max-3km-2100'][:]
#     exps = ['EC', 'Ens', 'Ens-2050', 'Ens-2100'][:2]
#     exps = ['MEns', 'MEns-2050', 'MEns-2100']    
#     provs = ['BC','AB']; seasons = ['summer','winter']; load_list=['season','prov']
#     prov = None; cluster = range(6); season = ['annual']; load_list=['season','cluster']
    clusters = [8,]; seasons = ['summer',]; load_list=['seasons','clusters']
    varlist = ['MaxPrecip_1d']; filetypes = ['hydro']; cluster_name = 'cluster_projection'
    lflatten = True; lfit = True; lrescale = True; lbootstrap = True
    # station criteria (we don't want too many stations...)
    constraints['min_len'] = 15
    constraints['lat'] = (45,55) 
    constraints['max_zerr'] = 500
    obsslices = dict(years=(1952,2012))
#     obsslices = [dict(years=(1950,1970)),dict(years=(1970,1990)),dict(years=(1990,2010))]
#     name_tags = ['_1','_2','_3']
#     cluster = [0,2,3]; constraints_rc['max_zerr'] = 1000; constraints_rc['lat'] = (45,60)  
    # load some data
    stnens,fitens,sclens = loadStationFit(names=exps, provs=provs, seasons=seasons, clusters=clusters, 
                                          cluster_name=cluster_name, varlist=varlist, lsourceScale=lrescale,
                                          stationtype='ecprecip', obsslices=obsslices, lshort=False,
                                          lensembleAxis=False, constraints=constraints, lall=True,
                                          sample_axis=None if lflatten else ('year','station'),
                                          lrescale=lrescale, reference=exps[0], target=exps[1],
                                          filetypes=filetypes, domain=1, lflatten=lflatten, lfit=lfit,
#                                           ensemble_list=['obsslices','name_tags'], name_tags=name_tags,
                                          lbootstrap=lbootstrap, nbs=3, load_list=load_list,)
    print(fitens[0][1])
    # set up plot    
    if len(seasons) == 1: subplot = len(clusters or prov)
    elif len(clusters or prov) == 1: subplot = len(seasons)
    else: subplot = (len(seasons),len(clusters or prov))
    fig,ax = evaFigAx(subplot, title=stnens[0][-1][varlist[0]].plot.title, sharex=True, sharey=False, 
                      stylesheet='myggplot', lpresentation=False, lreduce=False)
    # make plots
    for n in range(len(stnens)):
      distPlot(axes=ax.ravel()[n], varname=varlist[0], ens=stnens[n], fit=fitens[n], master=exps[1],
               sample_axis= None if lflatten else ('ensemble','station'), band_vars=None,
               scl=sclens[n] if lrescale else None, lrescale=lrescale, lsample=not lflatten, lanno=True,
               lbootstrap=lbootstrap, legend= bool(n+1==len(stnens)), reference=0,
               annotation=dist_annotation, defaults=dist_defaults)
    # adjust margins      
    if n == 0: fig.updateSubplots(left=-0.05, bottom=-0.0, top=0.01)
    if n > 3: fig.updateSubplots(left=-0.06, bottom=-0.04, top=0.05)
  
  # test quantile plots
  elif test == 'quantile_plot':

    # some settings for tests
    exps = ['EC','max-ens','max-ens-2050','max-ens-2100',][:]
    prov = 'AB'; season = 'summer'
#     varlist = ['MaxPreccu_1h']; filetypes = ['xtrm']
    varlist = ['MaxPrecip_1d']; filetypes = ['hydro']
#     lflatten = True; lfit = True; lrescale = True; lbootstrap = True
    lflatten = True; lfit = True; lrescale = True; lbootstrap = False
    sample_axis = 'station'; lsample = not lflatten
    # load some data
    stnens,fitens,sclens = loadStationFit(names=exps, provs=prov, seasons=season, varlist=varlist, lshort=False,
                                          lrescale=lrescale, reference='EC', target='max-ens', lfit=lfit,  
                                          stationtype='ecprecip', lflatten=lflatten, filetypes=filetypes,
                                          lbootstrap=lbootstrap, nbs=10)
    
    # set up plot    
    fig,ax = evaFigAx(1, title='{} {}'.format(prov,season), sharex=True, sharey=False, 
                      stylesheet='myggplot', lpresentation=False)    
    # make plots
    if lrescale: quantPlot(axes=ax, varname=varlist[0], fit=fitens[0], scl=sclens[0], legend=True, reference='max-ens',
                           percentiles=None, lmedian=True, lmean=False, sample_axis=sample_axis, 
                           lsample=lsample, stnset_name=prov, quantiles=(0.98,0.99), lbootstrap=lbootstrap,
                           lrescale=True, annotation=quant_annotation, defaults=quant_defaults)
    else: quantPlot(axes=ax, varname=varlist[0], fit=fitens[0], scl=None, legend=True, reference='max-ens', 
                    percentiles=None, lmedian=True, lmean=False, sample_axis=sample_axis,
                    stnset_name=prov, quantiles=(0.98,0.99), lbootstrap=lbootstrap, lsample=lsample,
                    annotation=quant_annotation, defaults=quant_defaults)

    # adjust margins
    fig.updateSubplots(left=+0.03, right=0)

  # show plots
  show()