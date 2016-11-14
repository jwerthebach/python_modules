#!/usr/bin/env python
# coding: utf-8

import os

import numpy as np
import cPickle as pickle
from tqdm import tqdm

import _i3hdf_to_df as i3df



def extract_run_numbers(file_list):
	"""Extracts the run numbers from the file list.

	Takes a file list with runs and extracts the run numbers from the file names.

	Parameter:
	----------
	file_list : list
		A list of strings containing the files which should be looked for run numbers.

	Returns:
	--------
	runs : list
		A list of integers containing all unique run numbers.
	"""

	# Extract the run numbers from the file list
	tmp_runs = []
	for f in file_list:
		# retrieve the file name
		name = f.split("/")[-1]
		# retrieve the run number
		run = name.split("_")[-2]
		# strip the run number from the prefix "RUN", convert it to int and add it to runs
		tmp_runs.append(int(run[3:]))
	return np.unique(tmp_runs)


def get_livetime(runs, goodrunlist):
	"""Returns the combined livetime of all runs.

	Parameter:
	----------
	runs : list
		A list of strings containing the six digit run numbers.
	goodrunlist : string
		Path to the GoodRunInfolist.

	Returns:
	--------
	livetime : int
		The summed up livetime of all runs given by `runs` in seconds.
	"""

	# Load the good run list
	run, live = np.loadtxt(goodrunlist, usecols=(0, 3), unpack=True, skiprows=2)
	# livetime in s
	livetime = 0
	# Crawl trough all given runs and look up the livetime
	for r in runs:
		mask = np.where(run == int(r))[0]
		# Add livetime from run to total livetime
		livetime += live[mask]
	return int(livetime)

def get_n_files_total(file_list):
	"""Returns the total number of i3files.

	Parameters:
	-----------
	file_list : array_like
		List containing hdf5 files.

	Returns:
	--------
	nfiles : int
		Number of i3files which were used to create the hdf5 files.
	"""

	# Get the base path of the files
	s = file_list[0].split('/')
	path = os.path.join(*s[:-2])
	# Initialize n_files with 0
	nfiles = 0
	# Loop over all files in file_list
	for f in file_list:
		# Split the filename to retrieve the set number and running number
		splitted = f.split("/")[-1].split(".")[0].split("_")
		# Reconstruct the file name of the pickle file that was used to create
		# the hdf5 file
		file_name = "i3Files_MC_{}_{}.pickle".format(splitted[3], splitted[-1])
		# Combine base path and pickle name
		prep_path = os.path.join(path, "prep_files", file_name)
		# Open the pickle file and retrieve the number of i3files
		with open(prep_path) as p:
			d = pickle.load(p)
		# Add number of i3files to total number of i3 files
		nfiles += len(d["i3_list"])-1
	return float(nfiles)

def _read_data(file_list, atts, silent=True):
	"""Reads all hdf5 files and converts them to a pandas data frame.

	Parameters:
	-----------
	file_list : array_like
		List containing hdf5 files to load.
	atts : array_like
		List containing all features to load.
		Must be a List even if just one feature should be loaded.

	Returns:
	--------
	df : pandas data frame
		Data frame containing the values.
	"""

	# Create container
	hdf_container = i3df.HDFcontainer(file_list=file_list,
									  exists_col="exists",
									  silent=silent)
	# Read all attributes
	df = hdf_container.get_df(atts)
	return df

def read_data(file_list, atts, n_events='None', silent=True):
	"""Reads hdf5 files from a list in chunks.

	Due to a memory leak in the i3_to_hdf script, which causes an exponential
	grow of used memory if one load a lot of files at once, I will load all
	files in chunks of two and append them to one big dataframe.

	This function replaces the old read_data function.

	I should check if this behavior is caused by the script or some of my data
	i.e. multiple use of the same ID.

	Parameters:
	-----------
	file_list : array_like
		List containing hdf5 files to load.
	atts : array_like
		List containing all features to load.
		Must be a List even if just one feature should be loaded.
	n_events : int
		Specify the number of events that should be loaded. It will load
		the first files until the number of events are greater than n_events.
		In the default all files are loaded.
		Default is 'None'
	silent : boolean
		Set to false to show a progress bar.
		Default is True

	Returns:
	--------
	df : pandas data frame
		Data frame containing the values.
	"""

	# Files loaded at once
	# Use two, because loading more than two at once causes the memory leak
	step_size=2
	# Load the first two and initialize the dataframe.
	df = _read_data(file_list[0:0+step_size], atts)

	# Load all the other files and append them to the existing one.
	with tqdm(total=len(file_list), unit='Files', disable=silent) as pbar:
		pbar.update(step_size)
		for i in np.arange(step_size,len(file_list),step_size):
			print()
			if n_events is 'None' or len(df.index) < n_events:
				df = df.append(_read_data(file_list[i:i+step_size], atts,
					silent=silent))
			if i*step_size > len(file_list):
				pbar.update(1)
			else:
				pbar.update(step_size)

	return df