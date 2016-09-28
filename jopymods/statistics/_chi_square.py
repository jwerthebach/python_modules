#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Module to perform weighted/unweighted chi2 tests.
https://root.cern.ch/doc/master/classTH1.html#a11153bd9c45ceac48bbfac56cb62ea74
'''
from __future__ import print_function
import numpy as np
from scipy.special import gammainc


class LengthError(Exception):
	pass

def Prob(chi2, ndf):
	'''Computation of the probability for a certain Chi-squared (chi2)
	and number of degrees of freedom (ndf).

	Calculations are based on the incomplete gamma function P(a,x),
	where a=ndf/2 and x=chi2/2.
	P(a,x) represents the probability that the observed Chi-squared
	for a correct model should be less than the value chi2.
	The returned probability corresponds to 1-P(a,x),
	which denotes the probability that an observed Chi-squared exceeds
	the value chi2 by chance, even for a correct model.
	'''
	
	if (ndf <= 0):
		return 0 # Set CL to zero in case ndf<=0

	if (chi2 <= 0):
		if (chi2 < 0):
			return 0
		else:
			return 1
	
	x0 = 0
	return 1-gammainc(0.5*ndf, 0.5*(chi2-x0))


def Chi2Test(h1, h2, w1=None, w2=None, option=None, res=None):
	'''chi2 test for comparing weighed and unweigthed histograms

	Parameters:
		h1		the first histogram
		h2		the second histogram

		w1		weights of the first histogram
		w2		weights of the second histogram (use this for "UW")
		
		weights must be a float, list or np.array
		w1=w2  = None: "UU" = experiment experiment comparison (unweighted-unweighted)
		w2    != None: "UW" = experiment MC comparison (unweighted-weighted).
						Note that the first histogram should be unweighted
		w1+w2 != None: "WW" = MC MC comparison (weighted-weighted)
		
		option	"NORM" = to be used when one or both of the histograms is scaled
						but the histogram originally was unweighted
				"P"     = print chi2, ndf, p_value, igood
				"CHI2" = returns chi2 instead of p-value
				"CHI2/NDF" = returns chi2/ndf

		res		computes normalized residuals and returns them in this array
	'''

	# Check the input
	if len(h1) != len(h2):
		raise LengthError(
			"Length 'h1' must be the same as length 'h2")
	if w1 is not None:
		if not isinstance(w1, (float, np.ndarray, list)):
			raise TypeError("'w1' has to be float, a list or a numpy array")
		if isinstance(w1, float):
			w1 = np.ones_like(h1)*w1
		if len(w1) != len(h1):
			raise LengthError(
				"Length 'w1' must be the same as length 'h1")
	if w2 is not None:
		if not isinstance(w2, (float, np.ndarray, list)):
			raise TypeError("'w2' has to be float, a list or a numpy array")
		if isinstance(w2, float):
			w2 = np.ones_like(h2)*w2
		if len(w2) != len(h2):
			raise LengthError(
				"Length 'w2' must be the same as length 'h2")


	# Calculate the chi2
	global chi2, ndf, igood
	chi2 = 0.
	ndf = 0
	igood = 0


	prob = Chi2(h1,h2,w1,w2,option,res);

	# Finishing: Check for output options
	if 'P' in option:
		if 'S' not in option:
			print("Chi2 = {}, Prob = {}, NDF = {}, igood = {}\n".format(chi2,prob,ndf,igood))
		return (chi2,prob,ndf,igood)
	if 'CHI2/NDF' in option:
		if (ndf == 0):
			return 0
		else:
			return chi2/ndf
	elif 'CHI2' in option:
		return chi2

	return prob
	


def Chi2(h1, h2, w1, w2, option, res):
	'''The computation routine of the Chisquare test.

	For the method description, see Chi2Test() function.
	Returns p-value
	parameters:
	- igood:
		igood=0 - no problems
		For unweighted unweighted  comparison
		igood=1'There is a bin in the 1st histogram with less than 1 event'
		igood=2'There is a bin in the 2nd histogram with less than 1 event'
		igood=3'when the conditions for igood=1 and igood=2 are satisfied'
		For  unweighted weighted  comparison
		igood=1'There is a bin in the 1st histogram with less then 1 event'
		igood=2'There is a bin in the 2nd histogram with less then 10 effective number of events'
		igood=3'when the conditions for igood=1 and igood=2 are satisfied'
		For  weighted weighted  comparison
		igood=1'There is a bin in the 1st  histogram with less then 10 effective number of events'
		igood=2'There is a bin in the 2nd  histogram with less then 10 effective number of events'
		igood=3'when the conditions for igood=1 and igood=2 are satisfied'

	- chi2 - chisquare of the test
	- ndf  - number of degrees of freedom (important, when both histograms have the same empty bins)
	- res  - normalized residuals for further analysis
	'''

	global chi2, ndf, igood

	nbinx1 = len(h1)
	nbinx2 = len(h2)
	
	i_start  = 1
	i_end    = nbinx1

	ndf = i_end - i_start

	comparisonUU = (w1 is None and w2 is None)
	comparisonUW = (w1 is None and w2 is not None)
	comparisonWW = (w1 is not None and w2 is not None)
	if (w1 is not None and w2 is None):
		print('First histogram: unweigthed and second: weighted.')
		print('But h1 has weights!')
		return -1
	
	scaledHistogram  = ("NORM" in option)
	if (scaledHistogram and not comparisonUU):
		print("NORM option should be used together with UU. It is ignored")




	#get number of events in histogram
	if (comparisonUU and scaledHistogram):
		# for unweighted histograms error=sqrt(entries in bin)
		w1 = np.sqrt(h1)
		w2 = np.sqrt(h2)
		# sum weights
		sumw1 = np.sum(w1)
		sumw2 = np.sum(w2)
		# raise error
		if (sumw1 <= 0 or sumw2 <= 0):
			print("Cannot use option NORM when one histogram has all zero errors");
			return -1

		# --- scale histograms according to errors
		# create mask to cover division by 0
		mask1 = np.where(w1<=0.)
		mask2 = np.where(w2<=0.)
		# divide squared histogram entries by squared weights
		n = h1*h1/w1
		m = h2*h2/w2
		# replace all 'div. by 0' (inf) with 0
		n[mask] = 0
		m[mask] = 0
		#avoid rounding errors
		n = np.floor(n+0.5)
		m = np.floor(m+0.5)
		# sum contents
		N  = np.sum(n)
		M  = np.sum(m)

	else:
		n = h1
		m = h2
		# sum contents
		N  = np.sum(n)
		M  = np.sum(m)
		# sum weights
		if comparisonWW:
			sumw1 = np.sum(w1*w1)
		if (comparisonWW or comparisonUW):
			sumw2 = np.sum(w2*w2)
		

	#checks that the histograms are not empty
	if (N == 0 or M == 0):
	  print("one histogram is empty")
	  return -1

	if ( comparisonWW  and ( sumw1 <= 0 and sumw2 <=0 ) ):
	  print("Hist1 and Hist2 have both all zero errors")
	  return -1

	#### --- THE CHI2 TEST ---
	k = l = 0

	#Experiment - experiment comparison
	if (comparisonUU):
		np1 = N*(n+m)/(N+M)

		# check for zero bin entries
		i_n  = np.in1d(n,0)
		i_m  = np.in1d(m,0)
		mask = np.invert(i_n&i_m)
		# no data means one degree of freedom less
		ndf -= len(np.where(mask == False)[0])
		# flag error only when of the two histogram is zero
		if len(np.where(i_n == True)[0]) > 0:
			igood += 1
			if 'S' not in option:
				print("There is a bin in h1 with less than 1 event.")
		if len(np.where(i_m == True)[0]) > 0:
			igood += 2
			if 'S' not in option:
				print("There is a bin in h2 with less than 1 event.")

		# calculate residuals
		if (res is not None):
			print('Yippy')
			#Habermann correction for residuals
			correc = (1-N/(N+M))*(1-(n+m)/(N+M))
			res = (n-np1)/np.sqrt(np1*correc)

			print(res)

		# remove the entries with both 0 from arrays n and m
		n = n[mask]
		m = m[mask]
		# calculate chi2
		chi2 = 1./(M*N)*np.sum(((M*n-N*m)**2)/(n + m))




	#unweighted - weighted  comparison
	# case of err2 = 0 and m not zero is treated without problems
	# by excluding second chi2 sum
	# and can be considered as a comparison data-theory
	if (comparisonUW):
		# n = unweighted histogram
		# m = weighted histogram
		s2 = w2*w2 		# sum of squares of weights of events in each bin

		# check for zero bin entries
		i_n  = np.in1d(n,0)
		i_m  = np.in1d(m,0)
		i_s  = np.in1d(s2,0)
		mask = np.invert(i_n&i_m)
		# no data means one degree of freedom less
		ndf -= len(np.where(mask == False)[0])
		# case weighted histogram has zero bin content and error
		# use as approximated  error as 1 scaled by a scaling ratio
		# estimated from the total sum weight and sum weight squared
		tmpmask = i_m&i_s
		if len(np.where(tmpmask == True)[0]) > 0:
			if sumw2 > 0:
				s2[tmpmask] = sumw2/M
		  	else:
				# return error because infinite discrepancy here:
				# bin1 != 0 and bin2 =0 in a histogram with all errors zero
				print("Hist2 has in bin {} zero content and all zero errors".format(np.where(tmpmask == True)[0]))
				return -1
		# flag error only when of the two histogram is zero
		if len(np.where(i_n == True)[0]) > 0:
			igood += 1
			if 'S' not in option:
				print("There is a bin in h1 with less than 1 event.")
		if len(np.where((m*m/s2) < 10)[0]) > 0:
			igood += 2
			if 'S' not in option:
				print("There is a bin in h2 with less than 10 effective events.")

		# calculate residuals
		if (res is not None):
			print('Yippy')
			p_tmp1 = M*m - N*s2
			p_tmp2 = 4.*M**2*s2*n
			p = p_tmp1 + np.sqrt(p_tmp1**2 + p_tmp2)
			p = p/(2*M**2)

			z_div  = np.sqrt((N*s2 - m*M)**2 + p_tmp2)
			z_tmp1 = (M*s2) / z_div
			z_tmp2 = 1 + (N*s2 - m*M) / z_div
			z = N*p*(1-N*p) * z_tmp1**2 + s2/4. * z_tmp2**2

			res = (m - M*p)/z

		# remove the entries with both 0 from arrays n and m

		# if n == 0, approximate by adding +1 to n
		n[i_n] += 1

		n = n[mask]
		m = m[mask]
		s2 = s2[mask]
		# recalculate the sum of n
		N = np.sum(n)

		#print(n)
		#print(m)
		#print(s2)
		# recalculate p with the corrected arrays
		p_tmp1 = M*m - N*s2
		p_tmp2 = 4.*M**2*s2*n
		p = p_tmp1 + np.sqrt(p_tmp1**2 + p_tmp2)
		p = p/(2*M**2)

		#print(p_tmp1)
		#print(p_tmp2)
		#print(p)

		# calculate chi2
		VALUES = (n-N*p)**2/(N*p) + (m-M*p)**2/s2
		#print(VALUES)

		#print(zip(p_tmp1,p_tmp2,p,VALUES))
		chi2 = np.sum(VALUES)




	# weighted - weighted  comparison
	if (comparisonWW):
		s1 = w1*w1
		s2 = w2*w2

		# check for zero bin entries
		i_n  = np.in1d(n,0)
		i_m  = np.in1d(m,0)
		i_s1  = np.in1d(s1,0)
		i_s2  = np.in1d(s2,0)
		mask = np.invert(i_n&i_m)
		# no data means one degree of freedom less
		ndf -= len(np.where(mask == False)[0])
		# cannot treat case of booth histogram have zero zero errors
		tmpmask = i_s1&i_s2
		if len(np.where(tmpmask == True)[0]) > 0:
			print("h1 and h2 both have bin {} all zero errors".format(np.where(tmpmask == True)[0]))
			return -1
		# flag error only when of the two histogram is zero
		tmp_m1 = np.where((n*n/s1) < 10)[0]
		tmp_m2 = np.where(s1[tmp_m1] > 0)[0]
		if (len(tmp_m1) > 0 and len(tmp_m2) >0):
			igood += 1
			if 'S' not in option:
				print("There is a bin in h1 with less than 10 effective events.")
		tmp_m1 = np.where((m*m/s2) < 10)[0]
		tmp_m2 = np.where(s2[tmp_m1] > 0)[0]
		if (len(tmp_m1) > 0 and len(tmp_m2) >0):
			igood += 2
			if 'S' not in option:
				print("There is a bin in h2 with less than 10 effective events.")

		# calculate residuals
		if (res is not None):
			p = (n*N/s1 + m*M/s2) / (N**2/s1 + M**2/s2)
			res = (n - N*p) / (s1*np.sqrt(1 - 1/(1+(M**2*s1)/(N**2*s2))))

		# remove the entries with both 0 from arrays
		n = n[mask]
		m = m[mask]
		s1 = s1[mask]
		s2 = s2[mask]
		# calculate chi2
		chi2 = np.sum((M*n - N*m)**2 / (N**2*s2 + M**2*s1))		




	prob = Prob(chi2,ndf)

	return prob