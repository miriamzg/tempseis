
import os
import sys
import numpy as np
#import matplotlib.pylab as plt
from numpy import linalg as LA

j = np.complex(0,1)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def load_f(filename):
	lines = open(filename).readlines()
	f00 = float(lines[0].split()[1])
	f10 = np.asarray([float(lines[1].split()[1]),float(lines[1].split()[2]),float(lines[1].split()[3])])
	f01 = float(lines[2].split()[1])
	f02 = float(lines[3].split()[1])
	f11 = np.asarray([float(lines[4].split()[1]),float(lines[4].split()[2]),float(lines[4].split()[3])])
	f20 = np.matrix([[float(lines[5].split()[1]),float(lines[5].split()[2]),float(lines[5].split()[3])],
		[float(lines[6].split()[1]),float(lines[6].split()[2]),float(lines[6].split()[3])],
		[float(lines[7].split()[1]),float(lines[7].split()[2]),float(lines[7].split()[3])]])
	return f00, f10, f01, f02, f11, f20



def load_derivative(der, station, component, event):
	from obspy.core import read
	
	tr = read(event+"/derivates/"+der+"/"+station+"*"+component+".sac")[0]
	tr_array = np.asarray(tr.data)
	return tr_array


def fft(tr):
	amp = tr.copy()
	amp.detrend('demean')
	amp.taper(0.05,type='cosine')
	freq = np.fft.rfftfreq(amp.stats.npts, amp.stats.delta)
	amp.data = np.abs(np.fft.rfft(amp.data))#/(0.5*tr.stats.npts)
	return freq, amp

def full_fft(tr):
	sp = np.fft.fft(tr.data)
	freq = np.fft.fftfreq(tr.times().shape[-1])
	return freq, sp


def seismicmoment2par(f00, f10, f01, f02, f11, f20):
	M0 = f00
	qc = f10 / M0
	tauc = f01 / M0
	delta_t = np.sqrt(f02 / M0)
	v0 = f11 / f02

	fxx, fxy, fxz = f20[0,0], f20[0,1], f20[0,2]
	fyx, fyy, fyz = f20[1,0], f20[1,1], f20[1,2]
	fzx, fzy, fzz = f20[2,0], f20[2,1], f20[2,2]

	W =  np.matrix([[fxx/M0, fxy/M0, fxz/M0], [fyx/M0, fyy/M0, fyz/M0],[ fzx/M0, fzy/M0, fzz/M0]])


	r1 = np.array([1,0,0])
	r2 = np.array([0,1,0])
	r3 = np.array([0,0,1])

	r = np.matrix([r1,r2,r3])

	l1 = np.sqrt(np.dot((r1.T * W ),r1) [0,0])/2.
	l2 = np.sqrt(np.dot((r2.T * W ),r2) [0,0])/2.
	l3 = np.sqrt(np.dot((r3.T * W ),r3) [0,0])/2.

	l = np.array([l1,l2,l3])


	return M0, qc, tauc, delta_t, v0, l

	

def origin_time_and_location(cmtfile):
	import datetime

	lines = open(cmtfile).readlines()
	year = int(lines[0].split()[1])
	month = int(lines[0].split()[2])
	day = int(lines[0].split()[3])
	hour = int(lines[0].split()[4])
	minute = int(lines[0].split()[5])
	sec = int(lines[0].split()[6].split(".")[0])
	microsec = 10000*int(lines[0].split()[6].split(".")[1])
	microsec = 0
	
	origin_time = datetime.datetime(year, month, day, hour, minute, sec, microsec)
	
	lat = float(lines[4].split()[1])
	lon = float(lines[5].split()[1])
	depth =  float(lines[6].split()[1])
	
	return origin_time, lat, lon, depth



def calc_f1_dS1_XXX(event_code, station, component,  f1):
	import glob
	derivative_list = ["dSdx", "dSdy", "dSdz"]#, "dS2dx2", "dS2dy2","dS2dz2", "dSdxdy", "dSdxdz", "dSdydz"]


	# load derivatives
	freq_der = {}
	dS = {}
	for der in derivative_list:
		#print "./CMTSOLUTION_" + event_code + "/derivatives/"+der+"/"+station+"."+component+".txt"
		filename = glob.glob("../CMTSOLUTION_" + event_code + "/derivatives/"+der+"/"+station+"."+component+".txt")[0]
		lines = open(filename)
		freq_der[der], dS[der] = np.genfromtxt(filename, delimiter='\t', unpack=True, dtype=np.complex128)
	dS1 = np.matrix([dS["dSdx"],dS["dSdy"],dS["dSdz"]])
	f1_dS = np.dot(f1, dS1)[0]
	
	
	return f1_dS

def calc_f1_dS1(event_code, station, component,  f1):
	import glob
	derivative_list = ["dSdx", "dSdy", "dSdz"]

	# load derivatives
	freq_der = {}
	dS = {}
	for der in derivative_list:
		#print "./CMTSOLUTION_" + event_code + "/derivatives/"+der+"/"+station+"."+component+".txt"
		filename = glob.glob("../CMTSOLUTION_" + event_code + "/derivatives/"+der+"/"+station+"."+component+".txt")[0]
		lines = open(filename)
		freq_der[der], dS[der] = np.genfromtxt(filename, delimiter='\t', unpack=True, dtype=np.complex128)

	f1_dS = []
	for i in range(0,len(freq_der["dSdx"])):
		a = 0.0
		a = dS["dSdx"][i]*f1[0] + dS["dSdy"][i]*f1[1] + dS["dSdz"][i]*f1[2] 
		f1_dS.append(a)


	#dS1 = np.matrix([dS["dSdx"],dS["dSdy"],dS["dSdz"]])
	#f1_dS = np.dot(f1, dS1)[0]
	
	
	return np.asarray(f1_dS)



def calc_f2_dS2(event_code, station, component,  f20):
	import glob
	derivative_list = ["dS2dx2", "dS2dy2","dS2dz2", "dSdxdy", "dSdxdz", "dSdydz"]

	# load derivatives
	freq_der = {}
	dS = {}
	for der in derivative_list:
		filename = glob.glob("./CMTSOLUTION_" + event_code + "/derivatives/"+der+"/"+station+"."+component+".txt")[0]
		lines = open(filename)
		freq_der[der], dS[der] = np.genfromtxt(filename, delimiter='\t', unpack=True, dtype=np.complex128)

	f2_dS2 = []
	for i in range(0,len(freq_der["dS2dx2"])):
		a = 0.0
		a = 	dS["dS2dx2"][i]*f20[0,0] + dS["dSdxdy"][i]*f20[1,0] + dS["dSdxdz"][i]*f20[2,0] +	\
			dS["dSdxdy"][i]*f20[0,1] + dS["dS2dy2"][i]*f20[1,1] + dS["dSdydz"][i]*f20[2,1] +	\
			dS["dSdxdz"][i]*f20[0,2] + dS["dSdydz"][i]*f20[1,2] + dS["dS2dz2"][i]*f20[2,2] 
	
		f2_dS2.append(a)
		
	
#	S = np.matrix( [[dS["dS2dx2"],dS["dSdxdy"],dS["dSdxdz"]], 
#			[dS["dSdxdy"],dS["dS2dy2"],dS["dSdydz"]], 	
#			[dS["dSdxdz"],dS["dSdydz"],dS["dS2dz2"]]])
#	print a
#	print b
#	f2_dS2 = np.trace(np.inner(a,b))
#	print f2_dS2

	return np.asarray(f2_dS2)











	
