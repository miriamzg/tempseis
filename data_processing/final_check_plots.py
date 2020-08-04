import os
import sys
import matplotlib.pylab as plt
import glob
import numpy as np
from scipy import signal
import matplotlib
from pylab import *


def onclick(event):
	global xx, yy
	y = event.ydata
	x = event.xdata
	yy.append(y)
	xx.append(x)
	if len(yy) == 1:
		fig.canvas.mpl_disconnect(cid)
		plt.close()



event_code = "CMTSOLUTION_201808190019A_GCMT"

# Frequency band
Tmin_p = 25
Tmax_p = 100

Tmin_s = 25
Tmax_s = 100

Tmin_r = 45  # 125
Tmax_r = 100  #	180


folder = "../database/" + event_code + "/" + "fortran_format_" +\
		 str(Tmin_p) + "_"  + str(Tmax_p) + "_"+ str(Tmin_s) + "_"+ \
		 str(Tmax_s) + "_"+ str(Tmin_r) + "_"+ str(Tmax_r) 



point_source_folder = folder + "/point_source/"
observed_folder = folder + "/observed_data/"
out_file = folder + "/station2use.txt"


use_all = "no" # type "no", if you want to check all the traces one by one


os.system("mkdir " + folder + "/final_check")

out = open(out_file,"w")
out.write("Station\tComp\tType\tUse it?\n")
filelist = glob.glob(point_source_folder+"*")
n = 0



for i in range(0,len(filelist)):
	print "----------------"
	ps_file = filelist[i]
	sta = ps_file.split("/")[-1].split("_")[0]
	comp = ps_file.split("/")[-1].split("_")[1]
	wavetype = ps_file.split("/")[-1].split("_")[2]
	if (comp == "T" and wavetype == "P"):
		pass
	else:
		print ps_file
		obs_file = glob.glob(observed_folder + sta + "_" + comp + "_" + wavetype + "_ff")[0]
		print obs_file
		n+=1
		print sta, comp, wavetype, n , "/", len(filelist)

		lines = open(ps_file).readlines()
		ff, iimag, rreal = [],[],[]
		ps_complex = []
		for i in range(2, len(lines)):
			f = float(lines[i].split()[0])
			real  = float(lines[i].split()[1])
			imm  = float(lines[i].split()[2])
			c = complex(real, imm)
			ff.append(f)
			rreal.append(real)
			iimag.append(imm)
			ps_complex.append(c)

		fig = plt.figure(1, figsize=(11.69, 8.27))
		plt.subplot(311)
		plt.plot(ff[0:int(len(ff)/2.)], rreal[0:int(len(ff)/2.)], color="black")

		plt.subplot(312)
		plt.plot(ff[0:int(len(ff)/2.)], iimag[0:int(len(ff)/2.)], color="black")


		lines = open(obs_file).readlines()
		ff, iimag, rreal = [],[],[]
		obs_complex = []
		for i in range(2, len(lines)):
			f = float(lines[i].split()[0])
			real  = float(lines[i].split()[1])
			imm  = float(lines[i].split()[2])
			c = complex(real, imm)
			ff.append(f)
			rreal.append(real)
			iimag.append(imm)
			obs_complex.append(c)
		
		plt.subplot(311)
		plt.plot(ff[0:int(len(ff)/2.)], rreal[0:int(len(ff)/2.)], color="red")
		#plt.xlim(-0.5,0.5)
		plt.xscale("log")

		plt.subplot(312)
		plt.plot(ff[0:int(len(ff)/2.)], iimag[0:int(len(ff)/2.)], color="red")
		#plt.xlim(-0.5,0.5)
		plt.xscale("log")

		

		npoints = len(ff)
		ps_inv =  np.fft.ifft(ps_complex, n=npoints)
		obs_inv =  np.fft.ifft(obs_complex, n=npoints)

		corr =  round(np.corrcoef(ps_inv.real, obs_inv.real)[0, 1],2)

		xx = np.arange(-10,10,0.1)
		yy1 = 0.000001*np.sin(xx)
		yy2 = 0.0000005*np.sin(xx)

		rratio = []
		ii = []
		for i in range(0,len(obs_inv.real)):
			if abs(obs_inv.real[i]) >= max(abs(obs_inv.real)) / 50.:
				rratio.append((ps_inv.real[i])/(obs_inv.real[i]))
				ii.append(i)

		amp_ratio = round(np.mean(rratio),2)

		print amp_ratio


		plt.subplot(313)
		plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, label="Observed", color="black", zorder=0, linewidth=3)
		plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, label="Point source", color="red", zorder=1)
	#	plt.xlim(0,800)
	#	plt.axhline(max(abs(obs_inv.real)) / 50.)
			
		if corr <= 0.75:
			plt.annotate("correlation coeff.: "+ str(corr),  xy=(20, 150),  xycoords="axes pixels", fontsize=14, color="red")
		else:
			plt.annotate("correlation coeff.: "+ str(corr),  xy=(20, 150),  xycoords="axes pixels", fontsize=14, color="black")
		if amp_ratio <= 0.75:	
			plt.annotate("amplitude coeff.: "+ str(amp_ratio),  xy=(20, 128),  xycoords="axes pixels", fontsize=14, color="red")
		else:
			plt.annotate("amplitude coeff.: "+ str(amp_ratio),  xy=(20, 128),  xycoords="axes pixels", fontsize=14, color="black")

		npts = len(ps_inv.real)
		h = int(npts/2.)

		if use_all != "yes":
			xx, yy = [], []
			cid = fig.canvas.callbacks.connect('button_press_event', onclick)
#mng = plt.get_current_fig_manager()
#mng.resize(*mng.window.maxsize())



		fill([-100,h,h,-100], [-1000,-1000,1000,1000], 'red', alpha=0.15)
		fill([npts,h,h,npts], [-1000,-1000,1000,1000], 'lime', alpha=0.15)

		plt.ylim(1.1*min(min(obs_inv.real),min(ps_inv.real)),1.1*max(max(obs_inv.real),max(ps_inv.real)))
		plt.xlim(0,npts)
		plt.suptitle("Station: "+sta + " Comp: " + comp + " Wavetype: "+wavetype)
		plt.savefig(folder + "/final_check/" + sta+"_"+comp+"_"+wavetype+".png")

		
		

		if use_all != "yes":

			plt.show()
			x = xx[0]
			y = yy[0]

			if x >= h:
				out.write(sta+"\t"+comp+"\t"+wavetype+"\tY\n")
			else:
				out.write(sta+"\t"+comp+"\t"+wavetype+"\tN\n")
		
		#plt.savefig(folder + "/final_check/" + sta+"_"+comp+"_"+wavetype+".png")
		plt.close()
	#	sys.exit()



out.close()





