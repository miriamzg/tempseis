import os
import sys
import numpy as np 
sys.path.insert(0, '../lib/')
from functions import *
from obspy.core import read
import glob
from datetime import datetime
from obspy.core.utcdatetime import UTCDateTime

event_code = sys.argv[1]
real = "yes"

channel = "MX" # for synthetic test




if real == "yes":
	channel = "BH"

# frequency band for first filtering
Tmin = 17.
Tmax = 300.


# Frequency bands
# p waves
Tmin_p = 20
Tmax_p = 70

# s waves
Tmin_s = 20
Tmax_s = 100

# surface waves
Tmin_r = 45
Tmax_r = 100

# time before and after the time window
extratime = 800. #800

sampling_rate = 0.5 # Hz

comp_list = ["Z","R","T"]
wavetype_list = ["P","S","W"]

data_folder   	= "../database/" + event_code + "/processed_data/noisy/" 
if real == "yes":
	data_folder   	= "../database/" + event_code + "/data_ready2use/" 
arrivals_file 	= "../database/" + event_code + "/picking_times_" + str(Tmin_p)+"_"+str(Tmax_p)+"_"+str(Tmin_s)+"_"+str(Tmax_s)+"_"+ str(Tmin_r)+"_"+str(Tmax_r)+".txt"
kernels_folder 	= "../database/" + event_code + "/kernels/"
ps_folder 		= "../database/" + event_code + "/synthetics/point_source/"  
out_folder 		= "../database/" + event_code + "/fortran_format_" + str(Tmin_p) + "_"  + str(Tmax_p) + "_"+ \
			str(Tmin_s) + "_"+ str(Tmax_s) + "_"+ str(Tmin_r) + "_"+ str(Tmax_r) 


#os.system("rm -r " + out_folder)
os.system("mkdir " + out_folder)

# read cut times from picking time file
station_comp = []
cut_file = arrivals_file
c_lines = open(cut_file).readlines()
cut_times = {}
for i in range(7,len(c_lines)):
	station = c_lines[i].split()[0] 
	comp = c_lines[i].split()[1].split(channel)[1]
	starttime 	= UTCDateTime(c_lines[i].split()[2])
	begin 		= float(c_lines[i].split()[3])
	origintime 	= UTCDateTime(c_lines[i].split()[4])
	p_start 	= UTCDateTime(c_lines[i].split()[5])
	p_end 		= UTCDateTime(c_lines[i].split()[6])
	s_start 	= UTCDateTime(c_lines[i].split()[7])
	s_end 		= UTCDateTime(c_lines[i].split()[8])
	r_start 	= UTCDateTime(c_lines[i].split()[9])
	r_end 		= UTCDateTime(c_lines[i].split()[10])
	cut_times[station,comp] = [starttime, begin, origintime, p_start, p_end, s_start, s_end, r_start, r_end]
	station_comp.append([station, comp])



#=====================================================================
#	Cut kernels and filtering
#=====================================================================
os.system("mkdir " + kernels_folder + "cut")
derivative_list = ["dSdx", "dSdy", "dSdz", "dS2dx2", "dS2dy2","dS2dz2", "dSdxdy", "dSdxdz", "dSdydz"]
for i in range(0,len(station_comp)):
	station = station_comp[i][0]
	comp = station_comp[i][1]
	origintime = cut_times[station, comp][2]
	for der in derivative_list:
		for wavetype in wavetype_list:
			#if not os.path.isfile(kernels_folder+"/kernels/"+station+"_"+comp+"_"+wavetype+"_"+der+".sac"):

				print "Cutting and filtering traces ", station, comp, der, wavetype
				starttime = cut_times[station, comp][0]
				if (comp == "T" and wavetype == "P"):
					pass
				else:
					if wavetype == "P":
						Tmin = Tmin_p
						Tmax = Tmax_p
						t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
					if wavetype == "S":
						Tmin = Tmin_s
						Tmax = Tmax_s
						t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]
					#if wavetype == "W":
						#Tmin = Tmin_r
						#Tmax = Tmax_r
						#t1, t2 = cut_times[station, comp][7], cut_times[station, comp][8]

					tr_tmp = read(kernels_folder+station+"_"+comp+"_"+der+".sac")[0]
					tr_cut_f = tr_tmp.copy()
					tr_cut_f = filter_trace(tr_tmp, float(Tmin), float(Tmax))
					smoothing_time = 0.1 #Tmax
					tr_cut_f = trim_trace_abs(tr_cut_f, origintime, t1, t2, smoothing_time, extratime)

					tr_cut_f.write(kernels_folder+"cut/"+station+"_"+comp+"_"+wavetype+"_"+der+".sac", format="SAC")


# =====================================================================
# convert derivatives into ascii (in frequency domain)
# =====================================================================
os.system("mkdir " + out_folder + "/derivatives/")

derivative_list = ["dSdx", "dSdy", "dSdz", "dS2dx2", "dS2dy2","dS2dz2", "dSdxdy", "dSdxdz", "dSdydz"]
for i in range(0,len(station_comp)):
		station = station_comp[i][0]
		comp = station_comp[i][1]
		print "Converting derivatives ", station, comp
		for der in derivative_list:
			for wavetype in wavetype_list:
				if (comp == "T" and wavetype == "P"):
					pass
				else:
				
					filename = glob.glob(kernels_folder + "cut/" + station+"_"+comp+"_"+wavetype+"_"+der+".sac")[0]

					trace = read(filename)[0]
					trace.interpolate(sampling_rate=sampling_rate, method="linear")

					omega, sp = full_fft(trace)
					out = open(out_folder + "/derivatives/" + station + "_"+comp+"_"+wavetype+"_" + der, "w")
					out.write(str(len(sp))+"\n")
					out.write("omega\t\treal\t\t\timaginary\n")
					for i in range(0,len(sp)):
						f = omega[i]
						a = np.real(sp[i])
						b = np.imag(sp[i])
						out.write(str('%.8f' % f)+"\t"+str('%.8e' % a)+"\t\t"+str('%.8e' % b)+"\n")
					out.close()

#=====================================================================
#Convert observed seismograms into asci (in frequency domain)
#=====================================================================

os.system("mkdir " + data_folder + "cut")
os.system("rm " + out_folder + "/observed_data/*")
os.system("mkdir " + out_folder + "/observed_data")

for i in range(0,len(station_comp)):
		station = station_comp[i][0]
		comp = station_comp[i][1]
		origintime = cut_times[station, comp][2]
		filelist = glob.glob(data_folder + "*." + station + "*" + channel + comp + "*.int.noisy")
		if real == "yes":
			filelist = glob.glob(data_folder + "*." + station + "*" + channel + comp )
		starttime = cut_times[station, comp][0]
		for wavetype in wavetype_list:
			if (comp == "T" and wavetype == "P"):
				pass
			else:
				
				if wavetype == "P":
					Tmin = Tmin_p
					Tmax = Tmax_p
					t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
				if wavetype == "S":
					Tmin = Tmin_s
					Tmax = Tmax_s
					t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]
				#if wavetype == "W":
					#Tmin = Tmin_r
					#Tmax = Tmax_r
					#t1, t2 = cut_times[station, comp][7], cut_times[station, comp][8]

				print "Preparing real data ", station, comp, wavetype, Tmin, Tmax

				if len(filelist) == 1:
					filename = filelist[0]	

					trace_tmp = read(filename)[0]
					start1 = trace_tmp.stats.starttime
					

					#trace = shift_stream(trace, timeshift)
					trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
			#		
					
				#	f = create_cutting_trace(trace_filt_tmp, t1, t2, 100.)
				#	trace_filt = trace_filt_tmp.copy()
				#	trace_filt.data = trace_filt_tmp.data * f.data
				#	print starttime, trace_filt_tmp.stats.starttime
				#	print "observed", starttime, t1, t2
					trace_filt = trim_trace_abs(trace_filt_tmp, starttime, t1, t2, float(Tmax), extratime)
				#	trace_filt.interpolate(sampling_rate=sampling_rate, method="linear")
					trace_filt.write(data_folder + "cut/" + station + "_"+comp+"_"+wavetype + ".sac", format="SAC")
					#-------test-------------------------
				#	tr_obs = read("/home/andrea/Dropbox/Work/NERC_Project_UCL/Real_data_inversion/Data/CMTSOLUTION_201310251710A_SYNT_50.0/data_ready2use_3D/II.ABKT.S3.MX" + comp)[0]
				#	plt.plot(tr_obs.times(),tr_obs.data, color="black", linewidth=3)
				#	plt.plot(trace_filt.times(), trace_filt.data, color="black", linewidth=3)
				#	plt.xlim(200,1000)
				#	plt.savefig("test.png")
				#	sys.exit()
					#------------------------------------

		
					omega, sp = full_fft(trace_filt)

					out = open(out_folder + "/observed_data/" + station + "_"+comp+"_"+wavetype+"_ff", "w")
					out.write(str(len(sp))+"\n")
					out.write("omega\t\treal\t\t\timaginary\n")
					for i in range(0,len(sp)):
						f = omega[i]
						a = np.real(sp[i])
						b = np.imag(sp[i])
						out.write(str('%.8f' % f)+"\t"+str('%.8e' % a)+"\t\t"+str('%.8e' % b)+"\n")
					out.close()

#=====================================================================
# 	# Convert point source seismogram to asci (in frequency domain)
#=====================================================================
os.system("mkdir " + ps_folder + "cut")
os.system("rm " + out_folder + "/point_source/*")
os.system("mkdir " + out_folder + "/point_source")
#ps_folder = "../../Synthetic_tests/source_dimensions/" + plane_code + "/point_source/"   
for i in range(0,len(station_comp)):
		station = station_comp[i][0]
		comp = station_comp[i][1]
		starttime = cut_times[station, comp][0]
		begin = cut_times[station, comp][1]
		origintime = cut_times[station, comp][2]
		for wavetype in wavetype_list:
 			if (comp == "T" and wavetype == "P"):
 				pass
 			else:
				
				if wavetype == "P":
					Tmin = Tmin_p
					Tmax = Tmax_p
					t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
				if wavetype == "S":
					Tmin = Tmin_s
					Tmax = Tmax_s
					t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]
				#if wavetype == "W":
					#Tmin = Tmin_r
					#Tmax = Tmax_r
					#t1, t2 = cut_times[station, comp][7], cut_times[station, comp][8]

				print "Preparing point source ", station, comp, wavetype, Tmin, Tmax


				filename = ps_folder + "*" + station + ".MX"+comp+".sem.sac"
				trace_tmp = read(filename)[0]
				trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
				trace_filt_tmp.interpolate(sampling_rate=sampling_rate, method="linear")
				trace_filt = trim_trace_abs(trace_filt_tmp, origintime, t1, t2, float(Tmax), extratime)
				trace_filt.write(ps_folder + "cut/" + station + "_"+comp+"_"+wavetype + ".sac", format="SAC")

				omega, sp = full_fft(trace_filt)
				out = open(out_folder + "/point_source/" + station + "_"+comp+"_"+wavetype+"_ps", "w")
				out.write(str(len(sp))+"\n")
				out.write("omega\t\treal\t\t\timaginary\n")
				for i in range(0,len(sp)):
					f = omega[i]
					a = np.real(sp[i])
					b = np.imag(sp[i])
					out.write(str('%.8f' % f)+"\t"+str('%.8e' % a)+"\t\t"+str('%.8e' % b)+"\n")
				out.close()




# os.system("rm ../NA/data/rfi_files/OBS/*" )
# os.system("rm ../NA/src/rfi_subs/point_source/*" )
# os.system("rm ../NA/src/rfi_subs/kernels/*" )

# os.system("cp " + out_folder + "/observed_data/* ../NA/data/rfi_files/OBS/" )
# os.system("cp " + out_folder + "/point_source/*  ../NA/src/rfi_subs/point_source/" )
# os.system("cp " + out_folder + "/derivatives/*   ../NA/src/rfi_subs/kernels/" )


# #=====================================================================
# #	write station file
# #=====================================================================
# print "Writing station file..."
# filelist = glob.glob("../NA/data/rfi_files/OBS/*")
# out = open("rfi.in","w")
# out.write("#\n# Input file for receiver function inversion specific information\n#\n")
# out.write("rfi_param                              /* input model   */\n")
# out.write("rfi_models                             /* output models */\n")
# out.write(str(len(filelist))+"\t\t\t/* nwave */\n")

# for fl in filelist:
# 	wavetype = fl.split("/")[-1].split("_")[2]
# 	data_name = fl.split("/")[-1].split("_ff")[0]
# 	if wavetype == "P":
# 		weight = "1.0"
# 	if wavetype == "S":
# 		weight = "1.0"
# 	if wavetype == "W":
# 		weight = "1.0"
# 	out.write(data_name + "\n")
# 	out.write(weight+ "\n")
# out.write("1                                       /* iwrite_models */")
# out.close()





#=====================================================================
# Plot and check final traces in the NA folder
#=====================================================================
print "Plotting and checking final traces..."
obs_folder 			= out_folder + "/observed_data/"
point_source_folder = out_folder + "/point_source/"
kernels_folder 		= out_folder + "/kernels/"

plot_folder = out_folder + "/comparision_plots"
os.system("mkdir " + plot_folder)

filelist = glob.glob(obs_folder+"*_ff")
for fl in filelist:
	print fl
	plt.figure(1, figsize=(11.69, 8.27))
	filename = fl.split("/")[-1]

	ps_file = point_source_folder + filename.replace("ff","ps")
	kernel_files = glob.glob(kernels_folder + filename.split("ff")[0]+"*")


	lines = open(fl).readlines()
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
	plt.plot(ff, rreal, color="black")

	plt.subplot(312)
	plt.plot(ff, iimag, color="black")

	lines = open(ps_file).readlines()
	ff, iimag, rreal = [],[],[]
	ps_complex  = []
	for i in range(2, len(lines)):
		f = float(lines[i].split()[0])
		real  = float(lines[i].split()[1])
		imm  = float(lines[i].split()[2])
		c = complex(real, imm)
		ff.append(f)
		rreal.append(real)
		iimag.append(imm)
		ps_complex.append(c)


	npoints = len(ff)
	obs_inv =  np.fft.ifft(obs_complex, n=npoints)
	ps_inv =  np.fft.ifft(ps_complex, n=npoints)

	plt.subplot(311)
	plt.plot(ff, rreal, color="red")
	plt.xlim(-0.5,0.5)

	plt.subplot(312)
	plt.plot(ff, iimag, color="red")
	plt.xlim(-0.5,0.5)

	plt.subplot(313)
	plt.plot(np.arange(0,len(obs_inv),1), obs_inv.real, label="Observed", color="black")
	plt.plot(np.arange(0,len(ps_inv),1), ps_inv.real, label="Point source", color="red")
	#plt.xlim(2500,2600)
	plt.legend(loc=2)

	# for k in kernel_files:
	# 	lines = open(k).readlines()
	# 	ff, iimag, rreal = [],[],[]
	# 	for i in range(2, len(lines)):
	# 		f = float(lines[i].split()[0])
	# 		real  = float(lines[i].split()[1])
	# 		imm  = float(lines[i].split()[2])

	# 		ff.append(f)
	# 		rreal.append(real)
	# 		iimag.append(imm)

	# 	plt.subplot(211)
	# 	plt.plot(ff, rreal)
	# 	plt.xlim(-0.5,0.5)

	# 	plt.subplot(212)
	# 	plt.plot(ff, iimag)
	# 	plt.xlim(-0.5,0.5)



	plt.savefig(plot_folder+"/"+filename+".png")
	plt.close()

# #=============================================================







