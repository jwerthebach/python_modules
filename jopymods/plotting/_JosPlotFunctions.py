#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

#cmap = plt.get_cmap('plasma')
cmap = plt.get_cmap('viridis')


color=['#440154', '#472777', '#3e4989', '#30678d', '#25828e', '#1e9d88',#
		'#35b778', '#6dce58', '#b5dd2b', '#fde724']

alpha=0.6


def JoHist2d(x, y, bins, xlabel='x', ylabel='y', clabel='# Events',file_name='Test.pdf'):
	sns.set(style='whitegrid', color_codes=True)

	plt.hist2d(x, y,
				bins=bins,
				cmap=cmap,
				cmin=0.1)

	plt.xlabel(xlabel)
	plt.ylabel(ylabel)

	cb = plt.colorbar() 
	cb.set_label(clabel)
	cmin, cmax = cb.get_clim()
	plt.clim(1, np.ceil(cmax))
	plt.savefig(file_name)
	plt.clf()


def JoHist1d(x, bins=50, MCweight=1, label='x', log=True, file_name='Test.pdf', color=color[0]):
	sns.set(style='whitegrid', color_codes=True)

	if isinstance(MCweight, int):
		hist = plt.hist(x,
						bins=bins,
						alpha=alpha,
						label=label,
						color=color,
						log=log)
	else:
		hist = plt.hist(x,
						bins=bins,
						alpha=alpha,
						weights=MCweight,
						label=label,
						color=color,
						log=log)

	plt.legend(loc=0)
	plt.ylabel('# Events')
	plt.xlabel(label)
	plt.savefig(file_name)
	plt.clf()


def JoSubPlots(x1, x2, x3, label1, label2, label3, pdf, bins=100, pltcolor=color[0], divE=0):
	'''Function to print three subplots horizontal.

	Plot three histograms in one figure and save the figure in a pdf file.
	x-limits of the lower two plots are 0-1.
	'''
	sns.set(style='whitegrid', color_codes=True)


	plt.subplot(311)
	pltbins = np.linspace(x1.min(), x1.max(), bins)
	hist1 = plt.hist(x1,
					 bins=pltbins,
					 alpha=alpha,
					 color=pltcolor,
					 log=True)
	ymin, ymax = plt.ylim()
	plt.ylim(0.1, ymax)
	#plt.xlim(x1.min(), x1.max())
	plt.xlabel(label1)

	plt.subplot(312)
	pltbins = np.linspace(0, 1, bins)
	hist2 = plt.hist(x2,
					 bins=pltbins,
					 alpha=alpha,
					 color=pltcolor,
					 log=True)
	ymin, ymax = plt.ylim()
	plt.ylim(0.1, ymax)
	plt.xlim(0, 1)
	plt.xlabel(label2)
	plt.ylabel('# Events')
	
	plt.subplot(313)
	pltbins = np.linspace(0, 1, bins)
	hist3 = plt.hist(x3,
					 bins=pltbins,
					 alpha=alpha,
					 color=pltcolor,
					 log=True)
	ymin, ymax = plt.ylim()
	plt.ylim(0.1, ymax)
	plt.xlim(0, 1)
	plt.xlabel(label3)
	plt.text(0.6, ymax/7., 'divE/E_sum = '+str(divE))

	plt.tight_layout()

	pdf.savefig()
	plt.close()


def GewErr(W):
	if len(W)==0:
		return 0.,0.
	return sum(W*W),sum(W)


def ErrCompPlot(x1, w1, x2, w2, bins=50, title=None,Log=False,xmin=None,xmax=None,space=0.07):
	'''Makes a comparison plot with weighted Histogramms.
	
	Will calculate the binerror for the weighted x1 and the simple sqrt(n)
	error for the unweighted x2.
	'''
	
	# default values
	linewidth  = 1.
	markersize = 3.
	alpha      = 0.7
	fmt        = 'o'
	colorData  = color[5]
	colorMC    = color[9]
	colorRatio = color[1]
		
	# --- prepare the histogramms
	# calculate binning
	if not xmin:
		xmin = min(x1)
	if not xmax:
		xmax = max(x1)
	bin_seq = np.linspace(xmin, xmax, bins)
	if Log:
		bin_seq = np.logspace(xmin, xmax, bins)
	
	# fill x1 values in the binning to calculate errors
	DIGI_X = np.digitize(x1,bins=bin_seq)-1
	# calculate the error
	Yerr = np.sqrt(
			np.array(
			map(lambda x: GewErr(w1[DIGI_X==x]),range(len(bin_seq)-1) )
				).transpose()[0]
			)
	
	sns.set(style='whitegrid', color_codes=True)
	
	
	gs = GridSpec(2, 1, height_ratios=[3, 1])
	gs.update(hspace=space)
	mainplt = plt.subplot(gs[0])
	ratioplt = plt.subplot(gs[1], sharex=mainplt)

	# main comparisson plot
	x_1, _ = np.histogram(x1, bins=bin_seq, weights=w1)
	x_2, _ = np.histogram(x2, bins=bin_seq)
	x_2 = x_2/w2
	
	bin_center = 0.5 * (bin_seq[1:] + bin_seq[:-1])
	hw = 0.5 * (bin_seq[1:] - bin_seq[:-1])
	
	mainplt.errorbar(bin_center,
					 x_1,
					 yerr=Yerr,
					 fmt=fmt,
					 lw=linewidth,
					 ms=markersize,
					 color=colorMC)
	mainplt.errorbar(bin_center,
					 x_1,
					 drawstyle='steps-mid',
					 color=colorMC,
					 alpha=alpha,
					 label='MC')
	
	mainplt.errorbar(bin_center,
					 x_2,
					 #yerr=np.sqrt(x_2),
					 fmt=fmt,
					 lw=linewidth,
					 ms=markersize,
					 color=colorData)
	mainplt.errorbar(bin_center,
					 x_2,
					 drawstyle='steps-mid',
					 color=colorData,
					 alpha=alpha,
					 label='Data')
	# plot config
	mainplt.legend(loc='best', frameon=True)
	plt.setp(mainplt.get_xticklabels(), visible=False)
	mainplt.set_yscale('log', nopose=True)
	mainplt.set_ylabel('# Events / Hz')
	
	# ratio plot
	y = ((x_1 / x_2) - 1.) * 100.
	mask = np.invert(np.isnan(y))
	_yerr = 100*(Yerr/x_1) 
	
	# plot line at 0 to guide the eye
	ratioplt.axhline(0., color='#000000',
					 linestyle='dashed',
					 linewidth=1,
					 zorder=1)
	ratioplt.errorbar(bin_center[mask],
					  y[mask],
					  yerr=_yerr[mask],
					  fmt=fmt,
					  lw=linewidth,
					  ms=markersize,
					  color=colorRatio,
					  zorder=2)
	# plot config
	ratioplt.set_xlabel(title)
	ratioplt.set_ylabel('rel. dev. / %')
	ratioplt.set_ylim(-50., 50.)
	ratioplt.set_xlim(xmin, xmax)

