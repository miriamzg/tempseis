import os
import sys
import numpy as np 
import glob
import matplotlib.pylab as plt
from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-2,2))


def time_windowing(list):
	dt = 1.0
	integral = []
	somma = 0.0
	for i in range(0,len(list.real)):
		#print trRsynt.data[i]
		somma += abs(list[i].real)
		integral.append(somma)
	
	ymin = (max(integral) - min(integral)) * 0.05
	ymax = (max(integral) - min(integral)) * 0.95 	

	istart = 0.0
	iend = 0.0
	for i in range(0,len(integral)):
		if integral[i] >= ymin :
			istart = i
			break
	for i in range(0,len(integral)):
		if integral[i] >= ymax:
			iend = i
			break

	tstart = istart * dt
	tend = iend * dt 

	return tstart, tend

	


#inversion_code = "test7.3"

results_folder = sys.argv[1]


observed_folder = results_folder + "/observed/"
predicted_folder = results_folder + "/best_predicted/"
point_source_folder = "../NA/src/rfi_subs/point_source/"




comp_list = ["Z","R","T"]




station_comp = []
filelist = glob.glob(observed_folder + "*_P_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])

station_comp = sorted(station_comp)

fig = plt.figure(1, figsize=(24, 18))
plt.subplots_adjust(bottom=0.05,
					 left=0.05, right = 0.98,
					 top=0.95, hspace=0.5,
					 wspace=0.25)


ncolumn = 6
n = 1
for i in range(0,len(station_comp)):
	#try:
		station = station_comp[i][0]
		comp = station_comp[i][1]
	
		wavetype = "P"

		observed = observed_folder + station + "_" + comp + "_"+wavetype + "_observed.asc"
		predicted = predicted_folder + station + "_" + comp + "_" + wavetype + "_best_predicted.asc"
		point_source = point_source_folder + station + "_" + comp + "_"+wavetype + "_ps"

		lines =  open(observed).readlines()
		oobsf, oobsr, oobsi, cc = [],[],[],[]
		for i in range(0,len(lines)):
			obsf = float(lines[i].split()[0])
			obsr = float(lines[i].split()[1])
			obsi = float(lines[i].split()[2])
			oobsf.append(obsf/6.28)
			oobsr.append(obsr)
			oobsi.append(obsi)
			c = complex(obsr, obsi)
			cc.append(c)
		obs_complex = np.asarray(cc, dtype=complex)

		lines =  open(point_source).readlines()
		ppsf, ppsr, ppsi, cc = [],[],[],[]
		for i in range(2,len(lines)):
			psf = float(lines[i].split()[0])
			psr = float(lines[i].split()[1])
			psi = float(lines[i].split()[2])
			ppsf.append(psf/6.28)
			ppsr.append(psr)
			ppsi.append(psi)
			c = complex(psr, psi)
			cc.append(c)
		ps_complex = np.asarray(cc, dtype=complex)		


	 	lines = open(predicted).readlines()
		ppredf, ppredr, ppredi, cc = [],[],[],[]
		for i in range(0,len(lines)):
			predf = float(lines[i].split()[0])
			predr = float(lines[i].split()[1])
			predi = float(lines[i].split()[2])
			ppredf.append(predf/6.28)
			ppredr.append(predr)
			ppredi.append(predi)
			c = complex(predr, predi)
			cc.append(c)
		
		pred_complex = np.asarray(cc, dtype=complex)


	 	npoints = len(obs_complex)
	 	obs_inv =  np.fft.ifft(obs_complex, n=npoints)
	 	ps_inv =  np.fft.ifft(ps_complex, n=npoints)
	 	pred_inv =  np.fft.ifft(pred_complex, n=npoints)

	 	mft = 0.0
	 	for j in range(0,len(obs_inv)):
	 		(obs_inv[j].real - pred_inv[j].real)**2
	 		mft += (obs_inv[j].real - pred_inv[j].real)**2
	 	mft = np.sqrt(mft)
	 	mft = format(mft, '1.1e')

	 	tstart, tend = time_windowing(obs_inv)

	 	a = 1+int((len(station_comp)+1)/float(ncolumn))
	 	b = ncolumn
	 	c = n
	 	ax1 = fig.add_subplot(111)
	 	print a,b,c
	 	ax = plt.subplot(a,b,c)
	 	ax.set(xlim=(300, 600))
	 	ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
	 	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, color="black", label="Observed", linewidth=1.5 , zorder=1)	 	
	 	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, color="0.6", label="Point source", linewidth=2.0 , zorder=0)
	 	plt.plot(np.arange(0,len(obs_inv),1), pred_inv.real, color="red", label="Predicted", linewidth=1.0 , zorder=2 )
	 	plt.grid(True)
	 #	plt.ylabel(wavetype + " waves")
	 	ax1.spines['top'].set_color('none')
	 	ax1.spines['bottom'].set_color('none')
	 	ax1.spines['left'].set_color('none')
	 	ax1.spines['right'].set_color('none')
	 	ax1.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	 	 
         #if (c - 1) % b   != 0 :
         #ax.set_yticklabels([])
	 	if c <= ( len(station_comp) - b ) : # this is stupid as it will need changing for each different thing 
  			print c 
 			ax.set_xticklabels([])
		#ax.set_xticklabels([])
		plt.xlim(tstart-30, tend+90)
		if len(station) == 3:
			plt.annotate(station+" "+ comp, xy=(0.75, 0.85), xycoords='axes fraction', fontsize=11)
			plt.annotate( str(mft),  xy=(0.75, 0.15), xycoords='axes fraction', fontsize=11 )
		if len(station) == 4:
			plt.annotate(station+" "+ comp, xy=(0.75, 0.85), xycoords='axes fraction', fontsize=11)
			plt.annotate( str(mft),  xy=(0.75, 0.15), xycoords='axes fraction', fontsize=11 )
	  	n += 1


#	except:
#		pass
#plt.show()
plt.grid(True)
#plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
          #
          #for a in ax[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
#ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
fig.text(0.5, 0.025, 'Time (s)', ha='center', va='center',fontsize=20 )
fig.text(0.025, 0.5, 'Displacement (m) ', ha='center', va='center', rotation='vertical',fontsize=20 )
plt.savefig(results_folder + "/seismograms_PTEST2.pdf")
#plt.show()
plt.close()


#===================================================================
station_comp = []
filelist = glob.glob(observed_folder + "*_S_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])

station_comp = sorted(station_comp)

fig = plt.figure(1, figsize=(26, 20))
plt.subplots_adjust(bottom=0.05,
                    left=0.05, right = 0.98,
                    top=0.95, hspace=0.7,
                    wspace=0.25)



ncolumn = 6
n = 1
for i in range(0,len(station_comp)):
	try:
		station = station_comp[i][0]
		comp = station_comp[i][1]
	
		wavetype = "S"

		observed = observed_folder + station + "_" + comp + "_"+wavetype + "_observed.asc"
		predicted = predicted_folder + station + "_" + comp + "_" + wavetype + "_best_predicted.asc"
		point_source = point_source_folder + station + "_" + comp + "_"+wavetype + "_ps"


		lines =  open(observed).readlines()
		oobsf, oobsr, oobsi, cc = [],[],[],[]
		for i in range(0,len(lines)):
			obsf = float(lines[i].split()[0])
			obsr = float(lines[i].split()[1])
			obsi = float(lines[i].split()[2])
			oobsf.append(obsf/6.28)
			oobsr.append(obsr)
			oobsi.append(obsi)
			c = complex(obsr, obsi)
			cc.append(c)
		obs_complex = np.asarray(cc, dtype=complex)

		
		lines =  open(point_source).readlines()
		ppsf, ppsr, ppsi, cc = [],[],[],[]
		for i in range(2,len(lines)):
			psf = float(lines[i].split()[0])
			psr = float(lines[i].split()[1])
			psi = float(lines[i].split()[2])
			ppsf.append(psf/6.28)
			ppsr.append(psr)
			ppsi.append(psi)
			c = complex(psr, psi)
			cc.append(c)
		ps_complex = np.asarray(cc, dtype=complex)	

	 	lines = open(predicted).readlines()
		ppredf, ppredr, ppredi, cc = [],[],[],[]
		for i in range(0,len(lines)):
			predf = float(lines[i].split()[0])
			predr = float(lines[i].split()[1])
			predi = float(lines[i].split()[2])
			ppredf.append(predf/6.28)
			ppredr.append(predr)
			ppredi.append(predi)
			c = complex(predr, predi)
			cc.append(c)
		
		pred_complex = np.asarray(cc, dtype=complex)


	 	npoints = len(obs_complex)
	 	obs_inv =  np.fft.ifft(obs_complex, n=npoints)
	 	pred_inv =  np.fft.ifft(pred_complex, n=npoints)
	 	ps_inv =  np.fft.ifft(ps_complex, n=npoints)


	 	mft = 0.0
	 	for j in range(0,len(obs_inv)):
	 		(obs_inv[j].real - pred_inv[j].real)**2
	 		mft += (obs_inv[j].real - pred_inv[j].real)**2
	 	mft = np.sqrt(mft)
	 	mft = format(mft, '1.1e')

	 	tstart, tend = time_windowing(obs_inv)

	 	a = 1+int((len(station_comp)+1)/float(ncolumn))
	 	b = ncolumn
	 	c = n
	 	ax1 = fig.add_subplot(111)
	 	print a,b,c
	 	ax = plt.subplot(a,b,c)
	 	ax.set(xlim=(300, 600))
	 	ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
	 	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, color="black", label="Observed", linewidth=1.5 , zorder=1)	 	
	 	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, color="0.6", label="Point source", linewidth=2.0 , zorder=0)
	 	plt.plot(np.arange(0,len(obs_inv),1), pred_inv.real, color="red", label="Predicted", linewidth=1.0 , zorder=2 )
	 	plt.grid(True)
	 #	plt.ylabel(wavetype + " waves")
	 	ax1.spines['top'].set_color('none')
	 	ax1.spines['bottom'].set_color('none')
	 	ax1.spines['left'].set_color('none')
	 	ax1.spines['right'].set_color('none')
	 	ax1.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	 	 
         #if (c - 1) % b   != 0 :
         #ax.set_yticklabels([])
	 	if c <= ( len(station_comp) - b ) : # this is stupid as it will need changing for each different thing 
  			print c 
 			ax.set_xticklabels([])
		#ax.set_xticklabels([])
		plt.xlim(tstart-30, tend+90)
		if len(station) == 3:
			plt.annotate(station+" "+ comp, xy=(0.75, 0.85), xycoords='axes fraction', fontsize=11)
			plt.annotate( str(mft),  xy=(0.75, 0.14), xycoords='axes fraction', fontsize=11 )
		if len(station) == 4:
			plt.annotate(station+" "+ comp, xy=(0.75, 0.85), xycoords='axes fraction', fontsize=11)
			plt.annotate( str(mft),  xy=(0.75, 0.14), xycoords='axes fraction', fontsize=11 )
	  	n += 1
	except:
		pass
fig.text(0.5, 0.02, 'time (s)', ha='center', va='center',fontsize=20 )
fig.text(0.02, 0.5, 'displacement (m) ', ha='center', va='center', rotation='vertical',fontsize=20 )
plt.savefig(results_folder + "/seismograms_STEST2.pdf")
#plt.show()
plt.close()





