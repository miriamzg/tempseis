import numpy as np
import scipy
import math
from numpy import linalg as LA
import matplotlib.pylab as plt
import sys
import os
from math import sin, cos, sqrt, atan2, radians, degrees

def distance(lat1, lon1, lat2, lon2):
	# approximate radius of earth in km
	R = 6371.0

	lat1 = radians(lat1)
	lon1 = radians(lon1)
	lat2 = radians(lat2)
	lon2 = radians(lon2)

	dlon = lon2 - lon1
	dlat = lat2 - lat1

	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
	c = 2 * atan2(sqrt(a), sqrt(1 - a))

	distance_deg =   degrees(c)
	distance_km = R * c

#	print("Result:", distance)
#	print("Should be:", 278.546, "km")

	return distance_deg, distance_km
#==============================================================

def filter_trace(trace, Tmin, Tmax):
	trace.detrend("demean")
	trace.taper(0.05, type="hann")
	trace.filter('bandpass', freqmin = 1/Tmax, freqmax=1/Tmin, corners=4, zerophase=True)
	trace.taper(0.05, type="hann")
	return trace

#--------------------------------------------------------------------

def full_fft(tr):
	sp = np.fft.fft(tr.data)
	omega = 2*np.pi*np.fft.fftfreq(tr.times().shape[-1]) / tr.stats.delta  # ANGULAR frequency
	return omega, sp

#----------------------------------------------------------------------

def get_fm(cmtfile):
	lines = open(cmtfile).readlines()
	Expo = float(lines[6].split("Expo=")[1].split()[0])
	mrr  = float(lines[6].split()[3])*10**Expo
	mtt  = float(lines[6].split()[4])*10**Expo
	mpp  = float(lines[6].split()[5])*10**Expo
	mrt  = float(lines[6].split()[6])*10**Expo
	mrp  = float(lines[6].split()[7])*10**Expo
	mtp  = float(lines[6].split()[8])*10**Expo

	m = [mrr, mtt, mpp, mrt, mrp, mtp]


	return m

#----------------------------------------------------------------------

def full_fft_test(tr):
	sp = np.fft.fft(tr.data)
	freq = np.fft.fftfreq(tr.times().shape[-1])
	return freq, sp
#--------------------------------------------------------------------

def shift_stream(tr, timeshift):
	tr_new = tr.copy()
	tr_new.detrend(type="demean")
	tr_new.taper(type="hann", max_percentage=0.05)
	dt = tr.stats.delta
	npts = tr.stats.npts
	i_timeshift = abs(int(timeshift / dt))
	if timeshift <= 0.0:
			for i in range(0,i_timeshift):
				tr_new.data[i] = 0.0
			j = 0
			for i in range(i_timeshift, npts):
				tr_new.data[i] = tr.data[j]
				j += 1
	else:
			j = i_timeshift
			for i in range(0,npts-i_timeshift):
				tr_new.data[i] = tr.data[j]
				j += 1
			for i in range(npts-i_timeshift, npts):
				tr_new.data[i] = 0.0

	return tr_new

#--------------------------------------------------------------------

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

#--------------------------------------------------------------------

def convolve_trace(S1, duration):
	from scipy import signal

	S1_new = S1.copy()
	dt = S1.stats.delta
	tt = np.arange(-S1_new.times()[-1]-dt, S1_new.times()[-1]+dt, dt)

	T = tt[1] - tt[0]

	t_function = [0.0] * len(tt)
	for i in range(0,len(tt)):
		t = tt[i]
		if t < - duration/2.0 or t > duration / 2.0:
			t_function[i] = 0.0
		else:
			t_function[i] = 1.0

	area = 0.0
	for i in range(0,len(tt)-1):
		area += ((t_function[i] + t_function[i+1]) * T) / 2.0

	for i in range(0,len(tt)):
		t_function[i] = t_function[i] / area


	S1_double = [0.0]*len(tt)

	j = 0
	for i in range(0,len(tt)):
		if tt[i] < 0:
			S1_double[i] = 0.0
		else:
			S1_double[i] = S1_new.data[j]
			j += 1



	convolution = signal.fftconvolve(S1_double, t_function, mode="same") * T
	#convolution = np.convolve(S1_double, t_function, mode="full") * T 

	plt.figure(1, figsize=(11.69, 8.27))
	plt.subplots_adjust(hspace=0.4, bottom=0.1)
	plt.subplot(211)
	plt.title("Source function")
	plt.plot(tt, t_function, color="black")
	plt.xlim(-30,30)
	plt.xlabel("time (s)")

	plt.subplot(212)
	plt.plot(tt, S1_double, color="black", linewidth=2, zorder=0, label="Point source (S1)")
	plt.plot(tt, convolution, color="red", linewidth=1, zorder=1, label="Convolved S1 (S2)")
	plt.xlim(-100,600)
	plt.legend(loc=2, fontsize=12)
	plt.xlabel("time (s)")
	plt.savefig("convolved.eps")
	plt.close()

	j = 0
	for i in range(int(len(convolution)/2.0), len(convolution)):
		S1_new.data[j] = convolution[i]
		j += 1

	return S1_new

#--------------------------------------------------------------------

def calculate_W(Lmax, Lmin, phi, dip, strike):
	# Calculate W matrix from Lmax, Lmin and phi
	phi=np.deg2rad(phi)
	dip=np.deg2rad(dip)
	strike=np.deg2rad(strike)

	#print Lmax, Lmin, phi, dip, strike
	M = [[Lmax, 0, 0],[0, Lmin, 0],[0,0,0]]
	# passo dagli autovettori nel sdr coordinato con l'ellisse a quello xy del piano di faglia
	# rotate I around z by phi angle
	I = [[1,0,0],[0,1,0],[0,0,1]] 		# versors in the reference system of the ellipse
	R_phi = [[np.cos(phi), -np.sin(phi),0],[np.sin(phi), np.cos(phi),0],[0,0,1]]
	R_dip = [[1,0,0],[0,np.cos(-dip), -np.sin(-dip)],[0,np.sin(-dip),np.cos(-dip)]]
	teta = ((np.pi/2.)-strike)
	R_strike = [[np.cos(teta), -np.sin(teta),0],[np.sin(teta),np.cos(teta), 0],[0,0,1]]
	
	u = np.dot(I, R_strike)
	u = np.dot(u, R_dip)
	u = np.dot(u, R_phi)
	S = u

	Sinv = LA.inv(S)
	
	W = np.dot(np.dot(S,M),Sinv)
	W = np.around(W, decimals=2)

	return W


def calculate_R(dip, strike):
	# Calculate R matrix
#	phi=np.deg2rad(phi)
	dip=np.deg2rad(dip)
	strike=np.deg2rad(strike)

	I = [[1,0,0],[0,1,0],[0,0,1]] 		# versors in the reference system of the ellipse
#	R_phi = [[np.cos(phi), -np.sin(phi),0],[np.sin(phi), np.cos(phi),0],[0,0,1]]
	R_dip =[[1	,0				,0],
			[0	,np.cos(-dip)	,-np.sin(-dip)],
			[0	,np.sin(-dip)	,np.cos(-dip)]]

	teta = ((np.pi/2.)-strike)

	R_strike = 	[[np.cos(teta)	, -np.sin(teta)	,0],
				[np.sin(teta)	,np.cos(teta)	, 0],
				[0				,0				,1]]
	
	u = np.dot(I, R_strike)
	u = np.dot(u, R_dip)
	#u = np.dot(u, R_phi)
	R = u
	R = np.around(R, decimals=10)
	return R





#-------------------------------------------------
def halfcos_left(t, t1, smooth_time):
    from numpy import sin, cos, pi
    t = t-t1
    smooth_time = float(smooth_time)
    y = 0.5+(cos((t*pi/smooth_time)+2*pi))/2.
    if t > 0:
        y = 1.0
    if t < -smooth_time:
        y = 0.0
    return y

def halfcos_right(t, t2, smooth_time):
    from numpy import sin, cos, pi
    t = t-t2
    smooth_time = float(smooth_time)
    y = 0.5+cos(t*pi/smooth_time)/2.
    if t < 0:
        y = 1.0
    if t > smooth_time:
        y = 0.0
    return y

def smoothed_box(t, t1, t2, smooth_time):
    y = halfcos_left(t, t1, smooth_time) * halfcos_right(t, t2, smooth_time)
    return y




def create_cutting_trace(tr, t1, t2, smoothing_time):
	f = tr.copy()
	dt = f.stats.delta
	for i in range(0,len(f.data)):
		f.data[i] = smoothed_box(f.times()[i], t1, t2, smoothing_time)

	return f 




def cut_trace(tr, t1, t2, smoothing_time):
	cutting_function = create_cutting_trace(tr, t1, t2, smoothing_time)
	tr_cut = tr.copy()
	tr_cut.data = tr.data*cutting_function.data
	return tr_cut, cutting_function



def trim_trace(tr, t1, t2, smoothing_time, extratime):
	cutting_function = create_cutting_trace(tr, t1, t2, smoothing_time)
	tr_cut = tr.copy()
	tr_cut.data = tr.data*cutting_function.data
	starttime = tr_cut.stats.starttime
	tr_cut.trim(starttime	=	starttime + t1 - extratime,\
				endtime		=	starttime + t2 + extratime, pad=True, fill_value=0.0)
	return tr_cut


#===============================================================

def halfcos_left_abs(t, t1, smooth_time):
    from numpy import sin, cos, pi
    t = t-t1
    smooth_time = float(smooth_time)
    y = 0.5+(cos((t*pi/smooth_time)+2*pi))/2.
    if t > 0:
        y = 1.0
    if t < -smooth_time:
        y = 0.0
    return y


def halfcos_right_abs(t, t2, smooth_time):
    from numpy import sin, cos, pi
    t = t-t2
    smooth_time = float(smooth_time)
    y = 0.5+cos(t*pi/smooth_time)/2.
    if t < 0:
        y = 1.0
    if t > smooth_time:
        y = 0.0
    return y


def smoothed_box_abs(t, origintime, t1, t2, smooth_time):
    y = halfcos_left_abs(origintime + t, t1, smooth_time) * halfcos_right_abs(origintime + t, t2, smooth_time)
    return y


def create_cutting_trace_abs(tr, origintime, t1, t2, smoothing_time):
	f = tr.copy()
	dt = f.stats.delta
	for i in range(0,len(f.data)):
		f.data[i] = smoothed_box_abs(f.times()[i], origintime, t1, t2, smoothing_time)

	return f 



#def cut_trace_abs(tr, starttime, t1, t2, smoothing_time):
#	cutting_function = create_cutting_trace_abs(tr, starttime, t1, t2, smoothing_time)
#	tr_cut = tr.copy()
#	tr_cut.data = tr.data*cutting_function.data
#	return tr_cut, cutting_function

def trim_trace_abs(tr, origintime, t1, t2, smoothing_time, extratime):
	cutting_function = create_cutting_trace_abs(tr, origintime, t1, t2, smoothing_time)
	tr_cut = tr.copy()
#	plt.subplot(311)
#	plt.plot(tr.times(), tr.data, color="black")
#	plt.subplot(312)
#	plt.plot(cutting_function.times(), cutting_function.data, color="red")
	tr_cut.data = tr.data*cutting_function.data
	delta_starttime = origintime - tr_cut.stats.starttime
	#print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>", delta_starttime, starttime, tr_cut.stats.starttime

	tr_cut.trim(starttime	=	t1 - extratime - delta_starttime,\
				endtime		=	t2 + extratime - delta_starttime, pad=True, fill_value=0.0)

#	plt.subplot(313)
#	plt.plot(tr.times(reftime=tr.stats.starttime), tr.data, color="black")
#	plt.plot(tr_cut.times(reftime=tr.stats.starttime), tr_cut.data, color="red")
#	plt.xlim(600,1200)
	dt = tr.stats.delta
#	plt.ylim(min(tr.data[int(600/dt):int(1200/dt)]), max(tr.data[int(600/dt):int(1200/dt)]))
#	plt.savefig("test.png")
#	plt.close()
#	sys.exit()
	return tr_cut



#===============================================================

def create_cutting_trace_zerocross(tr, starttime, t1, t2, smoothing_time):
	f = tr.copy()
	dt = f.stats.delta
	for i in range(0,len(f.data)):
		t = starttime + f.times()[i]
		if t >= t1 and t <= t2:
			f.data[i] = 1.0
		else:
			f.data[i] = 0.0 
	return f 



def trim_trace_zerocross(tr, starttime, t1, t2, smoothing_time, extratime):




	cutting_function = create_cutting_trace_zerocross(tr, starttime, t1, t2, smoothing_time)
	tr_cut = tr.copy()
	tr_cut.data = tr.data*cutting_function.data
	starttime = tr_cut.stats.starttime
	tr_cut.trim(starttime	=	t1 - extratime,\
				endtime		=	t2 + extratime, pad=True, fill_value=0.0)
	return tr_cut













