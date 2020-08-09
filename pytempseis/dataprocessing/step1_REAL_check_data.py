import os
import sys
from obspy.core import read
import glob
import matplotlib.pylab as plt
from pylab import *
sys.path.insert(0, '../lib/')


def onclick(event):
	global xx, yy
	y = event.ydata
	x = event.xdata
	yy.append(y)
	xx.append(x)
	if len(yy) == 1:
		fig.canvas.mpl_disconnect(cid)
		plt.close()

def filter_trace(trace, Tmin, Tmax):
	trace.detrend("demean")
	trace.taper(0.05, type="hann")
	trace.filter('bandpass', freqmin = 1/Tmax, freqmax=1/Tmin, corners=4, zerophase=True)
	trace.taper(0.05, type="hann")
	return trace


event_code = sys.argv[1]

Data_folder = "../database/" + event_code + "/processed_data/"
Synt_folder = "../database/" + event_code + "/synthetics/point_source"


plot_folder = Data_folder + "/first_check_plots/"
os.system("mkdir " + plot_folder)

Tmin = 35.
Tmax = 150.

xlim_min = 0.
xlim_max = 6000.

output_file ="../database/" + event_code + "/first_check.txt"
out = open(output_file,"w")
out.write("Sta\tComp\tQuality (Y: ok, N: bad, I: inverted)\n")
#-------------------------------------------------------
sta_lines = open("../STATIONS").readlines()
station_list = []
for i in range(0,len(sta_lines)):
#	if sta_lines[i].split()[0] == "ANTO":
		station_list.append(sta_lines[i].split()[0])

components_list = ["Z","R","T"]

for i in range(0,len(station_list)):
	station = station_list[i]
	for c in components_list:
		print "Station: ", station, " Component: ", c, "   ", i+1, "/", len(station_list)
		if len(glob.glob(Data_folder+"/*" + station + "*BH" + c)) == 1:
		#if len(glob.glob(Data_folder+"/*" + station + "*MX" + c + ".sem.sac_post_processed_shift.noisy")) == 1:
		#if len(glob.glob(Data_folder+"/*" + station + "*MX" + c + ".sem.sac_post_processed_shift.noisy")) == 1:
			filename_real = glob.glob(Data_folder+"/*" + station + "*BH" + c)[0]
			#filename_real = glob.glob(Data_folder+"/*" + station + "*MX" + c + ".sem.sac_post_processed_shift.noisy")[0]
			filename_synt = glob.glob(Synt_folder+"/*" + station + "*MX" + c + ".sem.sac")[0]


			tr_r = read(filename_real)[0]
			tr_r = filter_trace(tr_r, Tmin, Tmax)

			tr_s = read(filename_synt)[0]
			tr_s = filter_trace(tr_s, Tmin, Tmax)




			xlim_mid = (xlim_max + xlim_min) / 2.



			xx, yy = [], []
			fig, ax = plt.subplots()
			fig = plt.figure(1, figsize=(11.69, 8.27))
			plt.subplot(211)
			#plt.plot(tr_r.times(reftime=tr_r.stats.starttime), tr_r.data, color="black")
			#plt.plot(tr_s.times(reftime=tr_r.stats.starttime), tr_s.data, color="red")
			plt.plot(tr_r.times(reftime=tr_r.stats.starttime), tr_r.data, color="black", linewidth=3, zorder=0)
			plt.plot(tr_s.times(reftime=tr_r.stats.starttime), tr_s.data, color="red",zorder=10)
			plt.axvline(xlim_mid, linestyle=":")
			fill([0,xlim_mid,xlim_mid,0], [0,0,1000,1000], 'lime', alpha=0.2, edgecolor='r')
			fill([0,xlim_mid,xlim_mid,0], [0,0,-1000,-1000], 'red', alpha=0.2, edgecolor='r')
			fill([xlim_mid,80000,80000,xlim_mid], [-1000,-1000,1000,1000], 'yellow', alpha=0.2, edgecolor='r')
			plt.ylim(1.5*min(tr_r.data),1.5*max(tr_r.data))

			plt.annotate('SAVE IT', xy=(0.1,0.75) ,xycoords='axes fraction', fontsize=20)
			plt.annotate('BIN IT', xy=(0.1,0.25) ,xycoords='axes fraction', fontsize=20)
			plt.annotate('FLIP IT', xy=(0.75,0.75) ,xycoords='axes fraction', fontsize=20)
			
			plt.xlim(xlim_min, xlim_max)

			plt.subplot(212)
		#	plt.plot(tr_r.times(reftime=tr_r.stats.starttime), -tr_r.data, color="black")
		#	plt.plot(tr_s.times(reftime=tr_r.stats.starttime), tr_s.data, color="red")
			plt.plot(tr_r.times(reftime=tr_r.stats.starttime), -tr_r.data, color="black")
			plt.plot(tr_s.times(reftime=tr_r.stats.starttime), tr_s.data, color="red")
			plt.axvline(xlim_mid, linestyle=":")
			plt.xlim(xlim_min, xlim_max)


		
			cid = fig.canvas.callbacks.connect('button_press_event', onclick)
			mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())

			plt.suptitle(tr_r.id)
			plt.show()

			x = xx[0]
			y = yy[0]

			if y > 0 and x < xlim_mid:
				line = station + "\t" + c + "\tY\n"
				out.write(line)
				status = "SAVED"
			if y < 0 and x < xlim_mid:
				line = station + "\t" + c + "\tN\n"
				out.write(line)
				status = "CANCELLED"
			if x > xlim_mid:
				line = station + "\t" + c + "\tI\n"
				out.write(line)
				status = "FLIPPED"
			
			print status
			# plot and save
			fig = plt.figure(1, figsize=(11.69, 8.27))
			plt.subplot(211)
			plt.plot(tr_r.times(), tr_r.data, color="black", linewidth=3, zorder=0)
			plt.plot(tr_s.times(), tr_s.data, color="red",zorder=10)
			plt.ylim(1.5*min(tr_r.data),1.5*max(tr_r.data))
			plt.xlim(0, xlim_max)

			plt.subplot(212)
			plt.plot(tr_r.times(), -tr_r.data, color="black")
			plt.plot(tr_s.times(), tr_s.data, color="red")
			plt.xlim(0, xlim_max)
						
			plt.suptitle(tr_r.id+"\n"+status)
			plt.savefig(plot_folder+tr_r.id+".png")
			plt.close()

		else:
			print "file not found"


















