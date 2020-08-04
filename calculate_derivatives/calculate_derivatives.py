import os
import sys
import numpy as np 
sys.path.insert(0, '../lib/')
from functions import *
from obspy.core import read
import glob
from datetime import datetime
from obspy.core.utcdatetime import UTCDateTime
#from function_derivatives import *


event_code = sys.argv[1]
#channel = "BH"
channel = "MX"

# frequency band for first filtering
Tmin = 17.
Tmax = 300.

sampling_rate = 0.5 # Hz
synthetics_folder = "../database/" + event_code + "/synthetics/"
out_folder = "../database/" + event_code + "/kernels"




def calc_derivatives(Tmin, Tmax, station, comp, sampling_rate, derivatives_folder, out_folder):
	derivatives_folder = derivatives_folder.split("kernels")[0]

	R = 6371.

	lines = open(derivatives_folder.split("kernels")[0]+"/kernels_info.txt", "r").readlines()
	delta_x = float(lines[0].split()[2])
	delta_y = float(lines[1].split()[2])
	delta_z = float(lines[2].split()[2])

	delta_components = ["xp","yp","zp", "xm", "ym", "zm",
			"xpyp","xpym","xmyp","xmym",
			"xpzp","xpzm","xmzp","xmzm",
			"ypzp","ypzm","ymzp","ymzm","point_source"]

	derivative_list = ["dSdx", "dSdy", "dSdz", "dS2dx2", "dS2dy2","dS2dz2", "dSdxdy", "dSdxdz", "dSdydz"]
	tr = {}
	freq = {}
	tr_fft = {}
	dS = {}

	for del_comp in delta_components:
		folder = derivatives_folder  + del_comp + "/"
		tr_tmp = read(folder + "*"+station+".MX"+comp+".sem.sac")[0]
		tr_tmp .detrend("demean")
		tr_tmp.filter('bandpass', freqmin = 1/Tmax, freqmax=1/Tmin, corners=4, zerophase=True) 
		#tr_tmp.interpolate(sampling_rate=sampling_rate, method="linear")
		tr[del_comp, station] = tr_tmp


	for der_comp in derivative_list:
		dS[der_comp] = tr["point_source", station].copy()
		
	# First derivatives
	dS["dSdx"].data = (tr["xp", station].data - tr["xm", station].data) / (2*delta_x)
	dS["dSdy"].data = (tr["yp", station].data - tr["ym", station].data) / (2*delta_y)
	dS["dSdz"].data = (tr["zp", station].data - tr["zm", station].data) / (2*delta_z)

	# Second derivatives
	dS["dS2dx2"].data = (tr["xp", station].data - 2*tr["point_source", station].data + tr["xm", station].data) / ((delta_x)**2)
	dS["dS2dy2"].data = (tr["yp", station].data - 2*tr["point_source", station].data + tr["ym", station].data) / ((delta_y)**2)
	dS["dS2dz2"].data = (tr["zp", station].data - 2*tr["point_source", station].data + tr["zm", station].data) / ((delta_z)**2)

	dS["dSdxdy"].data = (tr["xpyp", station].data - tr["xpym", station].data - tr["xmyp", station].data + tr["xmym", station].data) / (4*delta_x*delta_y)
	dS["dSdxdz"].data = (tr["xpzp", station].data - tr["xpzm", station].data - tr["xmzp", station].data + tr["xmzm", station].data) / (4*delta_x*delta_z)
	dS["dSdydz"].data = (tr["ypzp", station].data - tr["ypzm", station].data - tr["ymzp", station].data + tr["ymzm", station].data) / (4*delta_y*delta_z)

	n = 1
	for der_comp in derivative_list:
		dS[der_comp].interpolate(sampling_rate=sampling_rate, method="linear")
		dS[der_comp].write(out_folder+"/"+station+"_"+comp+"_"+ der_comp + ".sac", format="SAC")
		plt.plot(dS[der_comp].times(), dS[der_comp].data, color="black")
		dS[der_comp].plot(outfile=out_folder+ "/" + station+"_"+comp+"_"+ der_comp + ".png")
		#plt.savefig(derivatives_folder+ "kernels/"+station+"_"+comp+"_"+ wavetype + "_" + der_comp +"_" + str(delta_x)+"km.png")
		plt.close()




	return
#=================================================================================================================================


os.system("mkdir " + out_folder)

station_comp = []
filelist = glob.glob(synthetics_folder+"point_source/*sem.sac")
for fl in filelist:
	station = fl.split("/")[-1].split(".")[1]
	comp = fl.split("/")[-1].split(".")[2][2:3]
	station_comp.append([station, comp])


#=====================================================================
#	calculate derivatives
#=====================================================================
for i in range(0,len(station_comp)):
	station = station_comp[i][0]
	comp = station_comp[i][1]

	if comp in ["Z","R","T"]:
	#	if not os.path.isfile(kernels_folder + "/kernels/" + station + "_" + comp + "_dSdx.sac"):
			print "Calculating derivatives ", station, comp


	#		os.system("python function_derivatives.py " + str(Tmin)+ " " + str(Tmax) + " "+station + " " + comp + " "\
	#														+ str(sampling_rate) + " "+ synthetics_folder + " " + out_folder )

			calc_derivatives(Tmin,Tmax,station,comp,sampling_rate,synthetics_folder,out_folder)

