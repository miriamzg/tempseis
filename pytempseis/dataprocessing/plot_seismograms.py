import os
import sys
import numpy as np 
import glob
import matplotlib.pylab as plt

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

results_folder = /Users/TheStuffofAlice/Documents/UCL/4th_year/Masters_Project/TEMPSEIS_package_v1.2_WORKING/database/CMTSOLUTION_022601C_GCMT/fortran_format_20_70_20_100_45_100 
#result_file = result_folder + "/rfi_models"


observed_folder = results_folder + "/observed_data/"
#predicted_folder = results_folder + "/best_predicted/"
point_source_folder = results_folder + "/point_source/"




comp_list = ["Z","R","T"]




station_comp = []
filelist = glob.glob(observed_folder + "*_P_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])

station_comp = sorted(station_comp)

fig = plt.figure(1, figsize=(11.69, 8.27))
plt.subplots_adjust(bottom=0.02,
					 left=0.02, right = 0.98,
					 top=0.98, hspace=0.05,
					 wspace=0.05)


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
	 	print a,b,c
	 	ax = plt.subplot(a,b,c)
	 	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, color="black", label="Observed", linewidth=1.5 , zorder=1)	 	
	 	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, color="0.6", label="Point source", linewidth=2.0 , zorder=0)
	 	plt.plot(np.arange(0,len(obs_inv),1), pred_inv.real, color="red", label="Predicted", linewidth=1.0 , zorder=2 )
	 #	plt.ylabel(wavetype + " waves")

	 	ax.set_yticklabels([])
		#ax.set_xticklabels([])
		plt.xlim(tstart-30, tend+90)
		if len(station) == 3:
			plt.annotate(station+" "+ comp, xy=(tend+15, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
		if len(station) == 4:
			plt.annotate(station+" "+ comp, xy=(tend+25, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
	  	n += 1
#	except:
#		pass
plt.savefig(results_folder + "/seismograms_P.png")
plt.close()




#===================================================================
station_comp = []
filelist = glob.glob(observed_folder + "*_S_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])

station_comp = sorted(station_comp)

fig = plt.figure(1, figsize=(11.69, 8.27))
plt.subplots_adjust(bottom=0.02,
					 left=0.02, right = 0.98,
					 top=0.98, hspace=0.05,
					 wspace=0.05)


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
	 	print a,b,c
	 	ax = plt.subplot(a,b,c)
	 	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, color="black", label="Observed", linewidth=1.5 , zorder=1)	 	
	 	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, color="0.6", label="Point source", linewidth=2.0 , zorder=0)
	 	plt.plot(np.arange(0,len(obs_inv),1), pred_inv.real, color="red", label="Predicted", linewidth=1.0 , zorder=2 )

	 #	plt.ylabel(wavetype + " waves")

	 	ax.set_yticklabels([])
		ax.set_xticklabels([])
		plt.xlim(tstart-30, tend+90)
		if len(station) == 3:
			plt.annotate(station+" "+ comp, xy=(tend+15, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
		if len(station) == 4:
			plt.annotate(station+" "+ comp, xy=(tend+25, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
	  	n += 1
	except:
		pass
plt.savefig(results_folder + "/seismograms_S.png")
plt.close()


station_comp = []
filelist = glob.glob(observed_folder + "*_W_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])

station_comp = sorted(station_comp)

fig = plt.figure(1, figsize=(11.69, 8.27))
plt.subplots_adjust(bottom=0.02,
					 left=0.02, right = 0.98,
					 top=0.98, hspace=0.05,
					 wspace=0.05)


ncolumn = 6
n = 1
for i in range(0,len(station_comp)):
	try:
		station = station_comp[i][0]
		comp = station_comp[i][1]
	
		wavetype = "W"

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
	 	print a,b,c
	 	ax = plt.subplot(a,b,c)
	 	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, color="black", label="Observed", linewidth=2.0 , zorder=0)	 	
	 	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, color="gray", label="Point source", linewidth=1.5 , zorder=1)
	 	plt.plot(np.arange(0,len(obs_inv),1), pred_inv.real, color="red", label="Predicted", linewidth=1.0 , zorder=2 )
	 #	plt.ylabel(wavetype + " waves")

	 	ax.set_yticklabels([])
		ax.set_xticklabels([])
		plt.xlim(tstart-200, tend+200)
		if len(station) == 3:
			plt.annotate(station+" "+ comp, xy=(tend+30, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
		if len(station) == 4:
			plt.annotate(station+" "+ comp, xy=(tend+20, 0.2*max(obs_inv.real)), fontsize=8)
			plt.annotate("mft " + str(mft),  xy=(tend, 0.8*min(obs_inv.real)), fontsize=8 )
	  	n += 1
	except:
		pass
	 	
plt.savefig(results_folder + "/seismograms_W.png")
plt.close()




	 	#sys.exit()










