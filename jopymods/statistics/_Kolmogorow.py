import numpy as np
import scipy as sc


def TomaszWrapKolg(x=np.ones(1), y=np.ones(1), xweights=[], yweights=[]):
	dmax = TomaszKolg(TomaszDistr(x, weights=xweights),
					  TomaszDistr(y, weights=yweights))
	N = len(x)
	M = len(y)
	# R.J.Barlow 8.4.5
	K = np.sqrt((1. * N * M) / (N + M))
	return dmax, K, dmax * K, 1 - sc.stats.kstwobign.cdf(dmax * K)


def TomaszDistr(data, weights=[]):
	if len(weights) == 0:
		weights = np.ones(len(data))

	srtidx = np.argsort(data)
	sumarray = np.zeros(len(data))
	xvalues = np.zeros(len(data))

	for index, value in enumerate(srtidx):
		sumarray[index] = weights[value]
		xvalues[index] = data[value]

	sumarray = np.cumsum(sumarray) / np.sum(sumarray)

	return xvalues, sumarray


def TomaszKolg(dist1, dist2):
	if len(dist1[0]) > len(dist2[0]):
		return TomaszKolgInner(dist2, dist1)
	return TomaszKolgInner(dist1, dist2)


def TomaszKolgInner(dist1, dist2):
	minimum = 3
	maximum = -3
	k = 0
	l = 0
	value = 0
	DoL = False

	while k != len(dist1[0]) and l != len(dist2[0]):
		value = dist2[1][l] - dist1[1][k]
		maximum = max(value, maximum)
		minimum = min(value, minimum)
		if (dist1[0][k] < dist2[0][l]):
			k += 1
		elif (dist1[0][k] > dist2[0][l]):
			l += 1
		else:
			if DoL:
				DoL = False
				l += 1
			else:
				DoL = True
				k += 1

	return max(abs(minimum), abs(maximum))


def testkolg():
	testsrt = [3, 12, 1, 18]
	testsrt2 = [2, 11, 13, 22]

	a = TomaszDistr(testsrt)
	b = TomaszDistr(testsrt2)
	TomaszKolg(a, b)
	N = 100000
	R1 = np.random.random_integers(0, N, size=N)
	R2 = np.random.random_integers(0, N, size=N/100)

	DW1 = TomaszDistr(R1)
	DW2 = TomaszDistr(R2)

	dmax = TomaszKolg(DW1, DW2)
	print "dmax:"
	print dmax

	alpha = 0.2
	Ka = np.sqrt(np.log(2 / alpha) / 2)
	Shit = N * np.sqrt(1 / (2. * N))
	print Shit, Ka

	from sc import stats
	print stats.ks_2samp(R1, R2)
	print TomaszWrapKolg(x=R1, y=R2)
	from matplotlib import pyplot as plt
	plt.plot(DW1[0], DW1[1])
	plt.plot(DW2[0], DW2[1])
	plt.show()