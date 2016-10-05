# ------------------------------------------------------------------------------
#    Author:  Johannes Wethebach
#    Mail:  johannes.werthebach@udo.edu
#    Version:  3.2
#
# ------------------------------------------------------------------------------


#Import the icecube packages into python :
from __future__ import division

from icecube import icetray, dataclasses, dataio, common_variables, simclasses, phys_services
from icecube import tableio
from icecube.common_variables import direct_hits, hit_multiplicity, hit_statistics, track_characteristics
from icecube.icetray import I3Units
from I3Tray import I3Tray


import os
import copy
import numpy as np
import pickle
import math as m

from JoPlotFunctions import *


color=['#440154', '#472777', '#3e4989', '#30678d', '#25828e', '#1e9d88',#
	   '#35b778', '#6dce58', '#b5dd2b', '#fde724']


def D4RToMap(frame):
	'''Saves values from D4R reconstruction to a I3_double_container.

	Using the values from Patricks I3Segment 'Eme100' and save it to 'D4R'.
	'''

	cascade = 'PatrickAnalysis_D4R100SplineMPECascadeParams'
	params  = 'PatrickAnalysis_D4R100SplineMPEParams'

	if frame.Stop==icetray.I3Frame.Physics:
		I3_double_container = dataclasses.I3MapStringDouble()
		if frame.Has(cascade) and frame.Has(params):
			I3_double_container['cascade_energy'] = frame[cascade].cascade_energy
			I3_double_container['cascade_energy_sigma'] = frame[cascade].cascade_energy_sigma
			I3_double_container['nDOMsCascade'] = frame[cascade].nDOMsCascade
			I3_double_container['mean'] = frame[params].mean
			I3_double_container['median'] = frame[params].median
			I3_double_container['peak_energy'] = frame[params].peak_energy
			I3_double_container['peak_sigma'] = frame[params].peak_sigma 
			I3_double_container['rllh'] = frame[params].rllh
			I3_double_container['nDOMs'] = frame[params].nDOMs
			I3_double_container['N'] = frame[params].N
			I3_double_container['N_err'] = frame[params].N_err
			I3_double_container['b'] = frame[params].b
			I3_double_container['b_err'] = frame[params].b_err
			I3_double_container['cascade_slant_depth'] = frame[cascade].cascade_slant_depth
		else:
			print 'D4R Failed!!!'
			I3_double_container['cascade_energy'] = 0
			I3_double_container['cascade_energy_sigma'] = 0
			I3_double_container['nDOMsCascade'] = 0
			I3_double_container['mean'] = 0
			I3_double_container['median'] = 0
			I3_double_container['peak_energy'] = 0
			I3_double_container['peak_sigma'] = 0
			I3_double_container['rllh'] = 0
			I3_double_container['nDOMs'] = 0
			I3_double_container['N'] = 0
			I3_double_container['N_err'] = 0
			I3_double_container['b'] = 0
			I3_double_container['b_err'] = 0
			I3_double_container['cascade_slant_depth'] = 0
		
		frame['D4R'] = I3_double_container

def GenKeys(i3_files):
	'''Generating keys to keep in further processing.

	Reads all keys from first physics frame and remove then all puls maps.
	'''

	# Read first i3 file, index 1 because first file in list is the GCD.
	FILE = dataio.I3File(i3_files[1])
	
	frame = FILE.pop_frame()	# First frame is DAQ
	frame = FILE.pop_frame()	# Second frame is PHYSICS

	# Save all keys from phys. frame
	keep_all = sorted(frame.keys())

	# Delete all keys containing 'Pulses'
	keep_booking = [ x for x in keep_all if 'Pulses' not in x ]

	# Remove other keys, because there are no converters
	#remove = ['BestTrackName', 'CalibratedWaveformRange', 'DSTTriggers',\
	#		'I3MCPESeriesMap']
	#
	#for key in remove:
	#	keep_booking.remove(key)
	#
	
	# Save all keys which should be deleted
	delete_booking = [ x for x in keep_all if 'Pulses' in x ]

	return keep_all,keep_booking,delete_booking

def KeysToDict(keys):
	'''Converts a given list of strings to a dict with tableio.default as value.'''

	values = [tableio.default]*len(keys)

	return dict(zip(keys, values))

def CheckEventType(frame, divE1, divE2, e_sum, q_rel):
	"""Applying cuts to set the event type.

	
	EventTypes:
	-1     -- Uncharacterized
	 0     -- Background
	 1     -- HEMuon
	 2     -- HEMuon with two leading muons
	 31/32 -- Balloon
	"""

	try:
		if divE1/e_sum > 0.3:
			if q_rel > 0.2:
				return 31
			else:
				return 1
		elif divE1/e_sum < 0.05 and divE2/e_sum > 0.3:
			if q_rel > 0.2:
				return 32
			else:
				return 2
		else:
			return 0
	
	except:
		return -1

def PolKoord(phi,theta,r=1):
	'''Converts pol coordinates to Cartesian coordinates with r=1'''
	return r*np.sin(theta)*np.sin(phi),r*np.sin(theta)*np.cos(phi),r*np.cos(theta)

def GetDirFromParticle(ptcl):
	'''Read the pol coordinates of the particle and converts them to Cartesian'''
	return PolKoord(ptcl.dir.azimuth,ptcl.dir.zenith)



class QCuts(icetray.I3ConditionalModule):
	'''Applying quality cuts

	Parameters:
	CutQtot    -- min. value of Qtotal
	CutLDirC   -- min. value of dir_track_length
	CutSpMPE   -- (bool) true if SplineMPE should not failed
	CutQmax    -- min. value of Qmax/Qtot
	'''
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter('CutQtot',                 # name
						  'Qtot > x',                # doc
						   1000)                     # default
		self.AddParameter('CutLDirC',                # name
						  'LDirC > x',               # doc
						   600)                      # default
		self.AddParameter('CutSpMPE',                # name
						  'SplineMPE not failed',    # doc
						   True)                     # default
		self.AddParameter('CutQmax',                 # name
						  'Qmax/Qtot < x',           # doc
						   0.2)                      # default
		self.AddParameter('CutZen',                  # name
						  'cos(zen) > x',            # doc
						   0.1)                      # default

	def Configure(self):
		self.CutQtot  = self.GetParameter('CutQtot')
		self.CutLDirC = self.GetParameter('CutLDirC')
		self.CutSpMPE = self.GetParameter('CutSpMPE')
		self.CutQmax  = self.GetParameter('CutQmax')
		self.CutZen   = self.GetParameter('CutZen')
		return True
	
	def Physics(self, frame):
		# Check if fit status should be used
		# if not use all frames
		if self.CutSpMPE:
			fit_status = (frame['SplineMPE'].fit_status == 0)
		else:
			fit_status == True

		# Do quality cuts, if all passed push frame
		if (frame['HitStatisticsValues'].q_tot_pulses > self.CutQtot
			and fit_status
			and frame['SplineMPEDirectHitsC'].dir_track_length > self.CutLDirC
			and np.cos(frame['SplineMPE'].dir.zenith) > self.CutZen):

			# Splitted the cuts up to avoid a 'divide by zero' error
			if frame['HitStatisticsValues'].q_max_doms/frame['HitStatisticsValues'].q_tot_pulses < self.CutQmax:
				self.PushFrame(frame)

class PlotWeigth(icetray.I3ConditionalModule):
	'''Save the energy and weight of MCPrimary1.

	Parameters:
	Cache  -- Array to save energies (MCPrimary1)
	Weight -- Array to save weigths from MC set (GaisserH3aWeight)
	'''
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("Cache", "Array to save energies", None)
		self.AddParameter("Weight", "Array to save weigths", None)
	
	def Configure(self):
		self.cache = self.GetParameter('Cache')
		self.weight = self.GetParameter('Weight')
	
	def Physics(self, frame):
		
		self.cache.append(frame['MCPrimary1'].energy)
		self.weight.append(frame['GaisserH3aWeight'].value)

		self.PushFrame(frame)

class MergeWeigth(icetray.I3ConditionalModule):
	'''Merges the two weigths into one.

	- GaisserH3aWeight
	- Hoerandel5
	- TIG1996
	- GlobalFit

	Saves all weigths into one I3_double_container.
	'''

	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("dataset", "dataset", None)
	
	def Configure(self):
		self.dataset = self.GetParameter('dataset')
		return True
	
	def Physics(self, frame):
		#print 'MergeWeigth'

		I3_container = dataclasses.I3MapStringDouble()

		CombWeight = frame['GaisserH3aWeight_'+str(self.dataset[0])].value
		for i,Dset in enumerate(self.dataset[1:]):
			CombWeight += frame['GaisserH3aWeight_'+str(Dset)].value
		I3_container['GaisserH3a'] = float(CombWeight)
		
		CombWeight = frame['Hoerandel5_'+str(self.dataset[0])].value
		for i,Dset in enumerate(self.dataset[1:]):
			CombWeight += frame['Hoerandel5_'+str(Dset)].value
		I3_container['Hoerandel5'] = float(CombWeight)
		
		CombWeight = frame['TIG1996_'+str(self.dataset[0])].value
		for i,Dset in enumerate(self.dataset[1:]):
			CombWeight += frame['TIG1996_'+str(Dset)].value
		I3_container['TIG1996'] = float(CombWeight)
		
		CombWeight = frame['GlobalFit_'+str(self.dataset[0])].value
		for i,Dset in enumerate(self.dataset[1:]):
			CombWeight += frame['GlobalFit_'+str(Dset)].value
		I3_container['GlobalFit'] = float(CombWeight)


		frame['Jo_MC_Weight'] = I3_container

		self.PushFrame(frame)

class WeigthToMap(icetray.I3ConditionalModule):
	'''Saves all weigths into one I3_double_container.

	- GaisserH3aWeight
	- Hoerandel5
	- TIG1996
	- GlobalFit
	'''

	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("dataset", "dataset", None)
	
	def Configure(self):
		self.dataset = self.GetParameter('dataset')
		return True
	
	def Physics(self, frame):
		#print 'MergeWeigth'

		I3_container = dataclasses.I3MapStringDouble()

		I3_container['GaisserH3a'] = float(frame['GaisserH3aWeight_'+str(self.dataset)].value)
		I3_container['Hoerandel5'] = float(frame['Hoerandel5_'+str(self.dataset)].value)
		I3_container['TIG1996']    = float(frame['TIG1996_'+str(self.dataset)].value)
		I3_container['GlobalFit']  = float(frame['GlobalFit_'+str(self.dataset)].value)


		frame['Jo_MC_Weight'] = I3_container

		self.PushFrame(frame)

class EventType(icetray.I3ConditionalModule):
	"""Applying cuts to set the event type.

	Generating attributes and event type are stored.
	EventTypes:
	-1     -- Uncharacterized
	 0     -- Background
	 1     -- HEMuon
	 2     -- HEMuon with two leading muons
	 31/32 -- Balloon
	"""
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter('RecoName',                 # name
						  'RecoName',                 # doc
						  'JoAnalysis')               # default

	def Configure(self):
		self.RecoName = self.GetParameter('RecoName')
		return True
	
	def Physics(self, frame):

		#print 'Still Alive'
		
		divE1 = frame['TomaszAnalysis']['MC_Bundle_Leading'] - frame['TomaszAnalysis']['MC_Bundle_SecondLeading']
		divE2 = frame['TomaszAnalysis']['MC_Bundle_SecondLeading'] - frame['TomaszAnalysis']['MC_Bundle_ThirdLeading']
		e_sum = frame['TomaszAnalysis']['MC_Bundle_EinSum']
		q_rel = frame['HitStatisticsValues'].q_max_doms/frame['HitStatisticsValues'].q_tot_pulses
		
		I3_double_container = dataclasses.I3MapStringDouble()

		I3_double_container['MC_divE1rel']       = float( divE1/e_sum )
		I3_double_container['MC_divE2rel']       = float( divE2/e_sum )
		I3_double_container['MC_Qrel']        = float( q_rel )
		I3_double_container['MC_Bundle_Type'] = float( CheckEventType(frame, divE1, divE2, e_sum, q_rel) )

		frame[self.RecoName] = I3_double_container
		
		self.PushFrame(frame)

class CountEventType(icetray.I3ConditionalModule):
	"""Counting different event types and type of primary particle.

	
	The values are saved in a dict with the following structure:
	{
	 'JoLeading'					: (int),
	 'JoBundle'						: (int),
	 'JoBalloon'					: (int),
	 'PrimParticleType'				: (list),
	 'PrimParticleCount'			: (list),
	 'PrimCountWeigthed_GH3a'		: (list),
	 'PrimCountWeigthed_Hoerandel5'	: (list),
	 'PrimCountWeigthed_TIG1996'	: (list),
	 'PrimCountWeigthed_GlobalFit'	: (list),
	 'TomaszLeading'				: (int),
	 'TomaszBundle'					: (int),
	 'TomaszBalloon'				: (int)
	}
	The weighted counts are for just one processed file and need to be divided by
	the number of total files.
	It is the type of the reconstructed particle. 
	"""
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("dict_name",
						  "Path+Name of the dictionary in which to save the values",
						  None)
		self.AddParameter("dataset",
						  "Number of the dataset",
						  None)

	def Configure(self):
		self.dict_name = self.GetParameter('dict_name')
		self.dataset   = self.GetParameter('dataset')
		self.temp_dict = {
		'JoLeading'						: 0,
		'JoBundle'						: 0,
		'JoBalloon'						: 0,
		'PrimParticleType'				: [],
		'PrimParticleCount'				: [],
		'PrimCountWeigthed_GH3a'		: [],
		'PrimCountWeigthed_Hoerandel5'	: [],
		'PrimCountWeigthed_TIG1996'		: [],
		'PrimCountWeigthed_GlobalFit'	: [],
		'TomaszLeading'					: 0,
		'TomaszBundle'					: 0,
		'TomaszBalloon'					: 0}
		return True
	   
	def Physics(self, frame):
		#print 'CountEventType'

		eventType  = frame['JoAnalysis']['MC_Bundle_Type']
		tomaszType = frame['TomaszAnalysis']['MC_Bundle_Type']


		# Check the event type and count one up
		if eventType == 1:
			self.temp_dict['JoLeading'] += 1
		elif eventType == 31:
			self.temp_dict['JoBalloon'] += 1
		else:
			self.temp_dict['JoBundle'] += 1


		# Check Tomasz event type and count one up
		if tomaszType == 1:
			self.temp_dict['TomaszLeading'] += 1
		elif tomaszType == 2:
			self.temp_dict['TomaszBalloon'] += 1
		else:
			self.temp_dict['TomaszBundle'] += 1
		

		# save all the primary particles
		prims = [frame[s] for s in frame.keys() if 'MCPrimary' in s]

		# if there is only one save its primary type
		if len(prims)==1:
			primType   = str(prims[0].type)
		
		# else find the primary particle which is closest to the 'SplineMPE' track
		else:
			DirSpe = GetDirFromParticle(frame['SplineMPE'])
			mindiff = 1000
			# calculate for each primary the angular difference to 'SplineMPE'
			for i,prim in enumerate(prims):
				DirPrim = GetDirFromParticle(prim)
				angdiff = np.arccos(np.dot(DirSpe,DirPrim))
				# 
				if  angdiff < mindiff:
					primType   = str(prims[i].type)

		# Add primary particle to the list
		try:
			# cast list in np.array to look if particle is already in list  
			temp_type = np.array(self.temp_dict['PrimParticleType'])
			# search for particle
			temp_index = np.where(temp_type == primType)
			# if particle is found this entry exists and contains the index of the particle 
			# if it don't exist an IndexError id thrown and will be caught below 
			index = temp_index[0][0]

			# add 1 to the corresponding particle count
			self.temp_dict['PrimParticleCount'][index] += 1
			# add the weight to the corresponding particle count
			GH3a = frame['GaisserH3aWeight_'+str(self.dataset)].value
			H5   = frame['Hoerandel5_'+str(self.dataset)].value
			TIG  = frame['TIG1996_'+str(self.dataset)].value
			GF   = frame['GlobalFit_'+str(self.dataset)].value

			self.temp_dict['PrimCountWeigthed_GH3a'][index]       += GH3a
			self.temp_dict['PrimCountWeigthed_Hoerandel5'][index] += H5
			self.temp_dict['PrimCountWeigthed_TIG1996'][index]    += TIG
			self.temp_dict['PrimCountWeigthed_GlobalFit'][index]  += GF

		# catch the error thrown above
		except:
			# add the new particle type to the list
			self.temp_dict['PrimParticleType'].append(primType)
			# add a new entry in the count list, initialized with one
			self.temp_dict['PrimParticleCount'].append(1)

			# add a new entry in the count list, initialized with the weight
			GH3a = frame['GaisserH3aWeight_'+str(self.dataset)].value
			H5   = frame['Hoerandel5_'+str(self.dataset)].value
			TIG  = frame['TIG1996_'+str(self.dataset)].value
			GF   = frame['GlobalFit_'+str(self.dataset)].value

			self.temp_dict['PrimCountWeigthed_GH3a'].append(GH3a)
			self.temp_dict['PrimCountWeigthed_Hoerandel5'].append(H5)
			self.temp_dict['PrimCountWeigthed_TIG1996'].append(TIG)
			self.temp_dict['PrimCountWeigthed_GlobalFit'].append(GF)

		
		
		self.PushFrame(frame)

	def Finish(self):
		f = open(self.dict_name, 'wb')
		pickle.dump(self.temp_dict, f, protocol=2)
		f.close()

class RunIDCorrector(icetray.I3ConditionalModule):
	'''IceTray Module providing functionality to generate correct I3Event Headers.

	The module checks the Event IDS to track when the next file is read and adjust the I3EventHeader.

	Icetray module providing functionality to generate correct I3Event
	Headers. This module was used for some broken simulations.
	This Module needs a list of # DAQ frames in each file. It must be in the same order as the list given the I3Reader.
	'''
	def __init__(self, context):
		'''Standard module __init__, two parameters are added:
			- "i3files" expecting the list of files same list as
				I3Reader
			- "nframes" expecting the list of # DAQ frames per file in the same
				order as the list given the I3Reader
			- "report" False/True  whether to print control output in the
				Finish method'''
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter('i3_files',                      # name
						  'i3 files read in',              # doc
						  None)                            # default
		self.AddParameter('nframes',
						  'Dict of DAQ frame.',
						  None)
		self.AddParameter('report',                        # name
						  'Report switching of IDs',       # doc
						  False)                           # default

	def Configure(self):
		'''The two parameters added in the __init__ method are fetched.
			and variables to count frames are initialized.'''
		i3_files = self.GetParameter('i3_files')
		self.nframes = self.GetParameter('nframes')
		self.report = self.GetParameter('report')
		
		i3_files = i3_files[1:]

		self.f_name, self.run_ids = self.extract_runids(i3_files)
		self.current_run = 0
		
		self.report_list = []
		self.q_frame_counter = 0
		self.p_frame_counter = 0
		self.frame = 0

	def DAQ(self, frame):
		'''Count frame and adjust I3EventHeader.'''
		self.q_frame_counter += 1
		self.correct_id(frame)
		self.PushFrame(frame)

	def Physics(self, frame):
		'''Same as DAQ. Probably the correct_id is obsolete because
		Q-Frames and P-Frames have the same event header.y.'''
		self.p_frame_counter += 1
		self.correct_id(frame)
		self.PushFrame(frame)

	def Finish(self):
		'''Prepare the report and print when "report" was True.'''
		last_run_id = int(self.run_ids[self.current_run])
		self.report_list.append([last_run_id, self.q_frame_counter,
								 self.p_frame_counter])
		if self.report:
			for run in self.report_list:
				print('%d Q-Frames and %d P-Frames with Run ID: %d' % (run[1],
																	   run[2],
																	   run[0]))

	def extract_runids(self, i3_files):
		'''Function to extract the Run IDs (file numbers) from the
		in_path.'''
		run_ids = []
		f_name  = []
		for f in i3_files:
			file_name = os.path.basename(f)
			splitted_file_name = file_name.split('.')
			dataset = splitted_file_name[-4]
			file_nr = splitted_file_name[-3]
			run_ids.append((int(dataset)*100000) + int(file_nr))
			f_name.append(file_name)
		return f_name, run_ids

	def correct_id(self, frame):
		'''Function that copies the old event header, adjusts the runID
		and adds the new header to the frame.'''
		header = copy.copy(frame['I3EventHeader'])
		del frame['I3EventHeader']

		if self.q_frame_counter > self.nframes[self.f_name[self.current_run]]:
			last_run_id = int(self.run_ids[self.current_run])
			self.report_list.append([last_run_id, self.q_frame_counter,
									 self.p_frame_counter])
			self.current_run += 1
			self.q_frame_counter = 1
			self.p_frame_counter = 1
			
		#print self.f_name[self.current_run], self.q_frame_counter, self.nframes[self.f_name[self.current_run]]

		header.run_id = self.run_ids[self.current_run]
		frame['I3EventHeader'] = header

class CountDAQFrames(icetray.I3ConditionalModule):
	"""Counting the DAQ frames."""
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("DAQframes",
						  "Variable to save the # of DAQ frames.",
						  None)

	def Configure(self):
		self.DAQframes = self.GetParameter('DAQframes')
		self.nframes = 0
		return True
	   
	def DAQ(self, frame):
		self.nframes += 1


	def Finish(self):
		self.DAQframes.append(self.nframes)

class GenAttrDeviation(icetray.I3ConditionalModule):
	"""The Attribute called "deviations" is calculated.

	Measuring the standard deviation of different fit directions in a frame.
	"""

	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
	
	def Configure(self):
		return True
	
	def Physics(self, frame):

		I3_container = dataclasses.I3MapStringDouble()

		n_fits   = 0
		zeniths  = []
		azimuths = []
		alphas   = []
		e        = np.array((1,1,1))
		for key in frame.keys():
			if ('Fit' in key and 'Params' not in key):
				try:
					d     = np.array((frame[key].dir.x, frame[key].dir.y, frame[key].dir.z))
					alpha = m.acos(sum(d*e)/(m.sqrt(sum(d**2)*sum(e**2))))
					alphas.append(alpha)
					zeniths.append(frame[key].dir.zenith)
					azimuths.append(frame[key].dir.azimuth)
					n_fits += 1
				except:
					pass

		std_dev_zen   = np.std(zeniths)
		std_dev_azi   = np.std(azimuths)
		std_dev_alpha = np.std(alphas)


		I3_container['std_dev_zen']   = float(std_dev_zen)
		I3_container['std_dev_azi']   = float(std_dev_azi)
		I3_container['std_dev_alpha'] = float(std_dev_alpha)
		I3_container['n_fits']        = float(n_fits)

		frame['Jo_Deviations'] = I3_container

		self.PushFrame(frame)

class PlotEvent(icetray.I3ConditionalModule):
	'''Plot all muons of the event.

	Color:
	dark purple  --  Background
	yellow       --  HE muon
	light green  --  HE muon in my and Patricks analysis
	ligth blue   --  HE muon in Patricks analysis

	Parameters:
	pdf  --  pdf file in which the plots are saved.
	'''
	def __init__(self, context):
		icetray.I3ConditionalModule.__init__(self, context)
		self.AddParameter("pdf", "pdf file", None)
	
	def Configure(self):
		self.pdf = self.GetParameter('pdf')
	
	def Physics(self, frame):
		energy = []
		e_sum = 0.
		
		MMCTrack = frame['MMCTrackList']

		for particle in MMCTrack:
			energy.append(particle.Ei)
			e_sum += particle.Ei

		energy = np.array(energy)


		divE = frame['TomaszAnalysis']['MC_Bundle_Leading'] - frame['TomaszAnalysis']['MC_Bundle_SecondLeading']



		pltcolor = color[0]
		if frame['PatrickAnalysis']['CUTS_ALL'] == 1.:
			pltcolor = color[4]
		if divE/e_sum > 0.3:
			pltcolor = color[9]
		if (frame['PatrickAnalysis']['CUTS_ALL'] == 1. and divE/e_sum > 0.3):
			pltcolor = color[7]


		JoSubPlots(energy, energy/energy.max(), energy/e_sum,
					'E_muon[MeV]', 'E_muon/E_max', 'E_muon/E_sum',
					self.pdf, pltcolor=pltcolor, divE=round(divE/e_sum, 4))

		self.PushFrame(frame)