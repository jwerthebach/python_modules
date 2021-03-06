
#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#cmap = plt.get_cmap("plasma")
cmap = plt.get_cmap("viridis")


color=["#440154", "#472777", "#3e4989", "#30678d", "#25828e", "#1e9d88",#
		"#35b778", "#6dce58", "#b5dd2b", "#fde724"]

alpha=0.6


def Hist2d(x, y, bins=10, xlabel="x", ylabel="y", clabel="# Events",
			file_name="Test.pdf", print_mean=False, print_events="None",
			colormap="viridis", view=False, xlim="None", ylim="None"):
	"""Creates a 2D Histogram with different options.

	This function is a wrapper around the standard `matplotlib.pyplot hist2d`.
	It sets the limit on both axis to the global minimum and maximum, so that a
	diagonal line corresponds to a perfect correlation. It also plots the
	`colorbar` and saves the figure to a file. There are other options you can
	use.

	Parameters
	----------
	x,y : array_like, shape (n,)
		Input values

	bins : [None | int | [int, int] | array_like | [array, array]]

		The bin specification:

		- If int, the number of bins for the two dimensions (nx=ny=bins).
		- If [int, int], the number of bins in each dimension (nx, ny = bins).
		- If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
		- If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).

		The default value is 10.

	xlabel,ylable, clable : string, optional
		Lable for the axises

		The default is `x = 'x'`, `y = 'y'`, `c = '# Events'`

	file_name : string, optional
		file name to save the figure. You can also specify a path.

		The default is `'Test.png'`

	print_mean : boolean, optional
		If set to true the mean and standard deviation are plotted for each x
		bin. The mean is marked by a red dot and the std is plotted as a red band.

		The default is `False`

	print_events : list, shape (3,), optional
		A list containing the string in which the total number of entries should
		be plotted as well as the position in the plot.

		- [0] : string, with exactly one '{}'. This text will be plotted.
		- [1] : float, x position relative from the left bottom (0...1)
		- [2] : float, y position relative from the left bottom (0...1)

		The default is 'None'.

	colormap : string, optional
		Define the `colormap`.

		The default is `'viridis'

	view: boolean
		If True the plot will be shown.

		Default is False

	xlim:  [int, int], optional
		Limits of the x axis

		The default is 'None'

	ylim:  [int, int], optional
		Limits of the y axis

		The default is 'None'
	"""

	# get min/max value
	if xlim == "None":
		x_min = x.min()
		x_max = x.max()
	else:
		x_min = xlim[0]
		x_max = xlim[1]
	if ylim == "None":
		y_min = y.min()
		y_max = y.max()
	else:
		y_min = ylim[0]
		y_max = ylim[1]


	# color map to use
	cmap = plt.get_cmap(colormap)

	# --- plot the correlation
	_, xedges, yedges, __ = plt.hist2d(x, y,
				bins=bins,
				cmap=cmap,
				cmin=0.1)

	if print_mean:
		# --- print mean and std of each bin
		# get the bin number for each event in the histogram
		indices  = np.digitize(x,xedges)-1
		# calculate bin centers to plot errors in the center of the bin
		bin_mid  = (xedges[:-1]+xedges[1:])/2
		# look for all events in one bin and calculate the mean and std
		bin_mean = np.array([np.mean(y[indices == i]) for i in range(len(bin_mid))])
		bin_std  = np.array([np.std(y[indices == i]) for i in range(len(bin_mid))])
		# print mean
		plt.plot(bin_mid, bin_mean, "r.")
		# finally print the error band
		plt.fill_between(bin_mid, bin_mean-bin_std, bin_mean+bin_std,
						facecolor="red", alpha=0.5)

	# --- final tweaks and save the figure
	# print the total number of events in the plot
	if print_events != "None":
		s = print_events[0].format(len(x))
		pos_x = print_events[1]
		pos_y = print_events[2]
		plt.text(pos_x , pos_y, s, fontsize=12)
	# plot x and y lables
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	# set the limits of x/y to the same value
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	# plot the colorbar, set the lable and set the limits
	cb = plt.colorbar()
	cb.set_label(clabel)
	cmin, cmax = cb.get_clim()
	plt.clim(1, np.ceil(cmax))
	# save the figure
	plt.savefig(file_name)
	if view:
		plt.show()
	plt.clf()


def Hist1d(x, bins=50, weights=1, xlabel="x", ylabel="Number of Events",
			dlabel=None, xlim=None, log=False, file_name="Test.pdf",
			color=None, view=False, histtype='step'):
	"""Creates a 1D histogram

	If `x` is multidimensional it plots the different histograms in one plot
	and uses weights if given.

	Parameters:
	-----------
	x : array_like, shape (n,m)
		Input values, `m` specifies the number of datasets.

	bins : [int | array_like], optional shape (1,m)

		The bin specification:

		- If int, the number of bins for all datasets (nx=bins).
		- If array_like, the number of bins specified for each dataset

		The default value is 50.

	MCweight : array_like, optional, shape (n,1|m)
		Containing the weights for the single events. The default uses the
		weight of 1 for all events.

		The default is None

	xlabel, ylabel : string, optional
		Label for the axises.

		The default is `x = 'x'`, `y = 'Number of Events'`

	dlabel : array_like, optional, shape (1,m)
		Sequence of strings to match multiple datasets.

		The default is None

	xlim : optional
		Length 2 sequence of floats with limit of x axis.

		The default is `None`

	log : boolean, optional
		If `True`, the histogram axis will be set to a log scale. If log is
		`True` and `x` is a 1D array, empty bins will be filtered out and only
		the non-empty (`n, bins, patches`) will be returned.

	Default is `False`


	file_name : string, optional
		file name to save the figure. You can also specify a path.

		The default is `'Test.pdf'`

	color : array_like, optional, shape (1,m)
		Color spec or sequence of color specs, one per dataset. Default (None)
		uses the standard line color sequence.

		Default is `None`

	view: boolean
		If True the plot will be shown.

		Default is False
	"""

	# get the number of datasets
	n_iter = len(x)

	# --- Prepare the options so that we build a loop
	# check for bins
	if isinstance(bins, int):
		bins = [bins]*n_iter
	# check if weights are present
	if isinstance(weights, int):
		weights = [None]*n_iter
	# check lables for datasets
	if dlabel == None:
		dlabel = [None]*n_iter
	elif isinstance(dlabel, str):
		dlabel = [dlabel]*n_iter
	# check for color
	if color is None:
		color = [None]*n_iter

	# iterate over all datasets
	for i in range(n_iter):

		hist = plt.hist(x[i],
						bins=bins[i],
						alpha=alpha,
						weights=weights[i],
						label=dlabel[i],
						color=color[i],
						log=log,
						histtype=histtype)

	plt.legend(loc='best')
	# plot x and y lables
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	# set the limits of x
	plt.xlim(xlim)
	# save the figure
	plt.savefig(file_name)
	if view:
		plt.show()
	plt.clf()


def ROC(fpr, tpr, fig_size=(5,5), color=None, lw=[1,1],
					file_name="Test.pdf", label=None, view=False):
	"""Creates a 1D histogram

	If `x` is multidimensional it plots the different histograms in one plot
	and uses weights if given.

	Parameters:
	-----------
	fpr, tpr : array_like
		Input values

	fig_size : optional
		w,h tuple in inches with size of the figure.

		The default is (10:5)

	color : array_like, optional
		Sequence of color specs. Default (None) uses the standard line color
		sequence.

		Default is `None`

	linewidth : array_like, optional
		Sequence of line widths.

		Default is `[1,1]`

	file_name : string, optional
		file name to save the figure. You can also specify a path.

		The default is `'Test.pdf'`

	label : str
		Label name of ROC curve.

		The default is `None`

	view: boolean
		If True the plot will be shown.

		Default is False
	"""
	# Set size of figure
	plt.figure(figsize=fig_size)
	# plot ROC curve
	plt.plot(fpr, tpr, color=color[0], lw=lw[0], label=label)
	# Plot diagonal line
	plt.plot([0, 1], [0, 1], color=color[1], lw=lw[1], linestyle='--')
	# set the limits of x
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.0])
	# plot x and y lables
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	# plot the legend
	plt.legend(loc="best", title="Data/MC classification")
	# save the figure
	plt.savefig(file_name)
	if view:
		plt.show()
	plt.clf()


def ImportancePlot(imp, feat, bins=50, xlabel="Feature Importance",
					ylabel="Number of Events", n_print=10, xlim=None,
					height=10, fig_size=(5,3), color=None, text_color='k',
					file_name="Test.pdf", view=False):
	"""Creates a 1D histogram

	If `x` is multidimensional it plots the different histograms in one plot
	and uses weights if given.

	Parameters:
	-----------
	imp : array_like
		Input values, containing the importance of each feature in `feat`

	feat : array_like
		the names of the features given in `imp`. Make sure that it has the
		same order as imp so that imp[0] belongs to feat[0].

	bins : [int | array_like], optional
		If an integer is given, `bins + 1` bin edges are returned, consistently
		with `numpy.histogram()` for `numpy version >= 1.3`.

		Unequally spaced bins are supported if `bins` is a sequence.

		default is 50

	xlabel, ylabel : string, optional
		Label for the axises.

		The default is `x = 'Feature Importance'`, `y = 'Number of Events'`

	n_print : int, optional
		The features with the `n_print` highest features are marked in the plot.

		The default is 10.

	xlim : optional
		Length 2 sequence of floats with limit of x axis.

		The default is `None`

	height : int, optional
		Starting height for the text annotations.

		The default is 0.

	fig_size : optional
		w,h tuple in inches with size of the figure.

		The default is (10:5)

	color, text_color : str, optional
		Color spec. Default (None) uses the standard line color sequence.

		Default is `None`

	file_name : string, optional
		file name to save the figure. You can also specify a path.

		The default is `'Test.pdf'`

	view: boolean
		If True the plot will be shown.

		Default is False
	"""
	# Set size of figure
	plt.figure(figsize=fig_size)
	# Plot the histogram
	plt.hist(imp, bins=bins, alpha=0.6, color=color)
	# Plot markings for the features
	plt.plot(imp, np.zeros(len(imp)), linestyle='', marker=2, color=color,
			markersize=10)

	index = np.argsort(imp)
	for i in range(n_print):
		n = index[-(i+1)]
		plt.annotate(feat[n], xy=(imp[n], 6), xytext=(imp[n], height+i*5),
					ha='left', va='bottom', rotation=30, color=text_color,
					arrowprops=dict(facecolor=text_color, shrink=0.01, width=0.2,
									headwidth=4, headlength=4))

	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.xlim(xlim)
	plt.savefig(file_name)
	if view:
		plt.show()
	plt.clf()


def SubPlots(x1, x2, x3, label1, label2, label3, pdf, bins=100,
			pltcolor=color[0], divE=0):
	"""Function to print three subplots horizontal.

	Plot three histograms in one figure and save the figure in a pdf file.
	x-limits of the lower two plots are 0-1.
	"""


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
	plt.ylabel("# Events")

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
	plt.text(0.6, ymax/7., "divE/E_sum = "+str(divE))

	plt.tight_layout()

	pdf.savefig()
	plt.close()


def GewErr(W):
	if len(W)==0:
		return 0.,0.
	return sum(W*W),sum(W)


def ErrCompPlot(x1, w1, x2, w2, bins=50, title=None,Log=False,xmin=None,
				xmax=None,space=0.07):
	"""Makes a comparison plot with weighted Histogramms.

	Will calculate the binerror for the weighted x1 and the simple sqrt(n)
	error for the unweighted x2.
	"""

	# default values
	linewidth  = 1.
	markersize = 3.
	alpha      = 0.7
	fmt        = "o"
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
					 drawstyle="steps-mid",
					 color=colorMC,
					 alpha=alpha,
					 label="MC")

	mainplt.errorbar(bin_center,
					 x_2,
					 #yerr=np.sqrt(x_2),
					 fmt=fmt,
					 lw=linewidth,
					 ms=markersize,
					 color=colorData)
	mainplt.errorbar(bin_center,
					 x_2,
					 drawstyle="steps-mid",
					 color=colorData,
					 alpha=alpha,
					 label="Data")
	# plot config
	mainplt.legend(loc="best", frameon=True)
	plt.setp(mainplt.get_xticklabels(), visible=False)
	mainplt.set_yscale("log", nopose=True)
	mainplt.set_ylabel("# Events / Hz")

	# ratio plot
	y = ((x_1 / x_2) - 1.) * 100.
	mask = np.invert(np.isnan(y))
	_yerr = 100*(Yerr/x_1)

	# plot line at 0 to guide the eye
	ratioplt.axhline(0., color="#000000",
					 linestyle="dashed",
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
	ratioplt.set_ylabel("rel. dev. / %")
	ratioplt.set_ylim(-50., 50.)
	ratioplt.set_xlim(xmin, xmax)


