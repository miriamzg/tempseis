import os
import sys
import matplotlib.pylab as plt
from obspy.core import read
from obspy.taup import TauPyModel
import glob
import obspy
#import geopy.distance
import numpy as np
from obspy.clients.iris import Client
sys.path.insert(0, '../lib/')
from functions import *
from obspy.core.util import AttribDict


real = "yes"

event_code = sys.argv[1]
channel = "MX"

folder = "../database/" + event_code 
data_folder = folder + "/processed_data/"

cmt_file = folder + "/" + event_code


if real == "yes":
    channel = "BH"
    data_folder = folder + "/data_ready2use/"

Tmin_p = 20 #25
Tmax_p = 70 #60

Tmin_s = 20
Tmax_s = 100

Tmin_r = 45  # 125
Tmax_r = 100  #	180

pretime_p 	= 100
posttime_p 	= 100

pretime_s 	= 100
posttime_s 	= 100

pretime_r 	= 2*Tmax_r   	# 500
posttime_r 	= 4*Tmax_r	# 500


plot_folder = folder + "/plots_picking_"+ str(Tmin_p) + "_" + str(Tmax_p) + "_" + str(Tmin_s) + "_" + str(Tmax_s) + "_" + str(Tmin_r) + "_" + str(Tmax_r) + "/"
os.system("mkdir " + plot_folder)


lines = open(cmt_file).readlines()
ev_lat = float(lines[4].split()[1])
ev_lon = float(lines[5].split()[1])
ev_dep = float(lines[6].split()[1])


taup_model = TauPyModel(model="iasp91")

sta_lines = open("../STATIONS").readlines()

file_list = sorted(glob.glob(data_folder +  "*"+channel+"Z*.corr.int"))
file_list += sorted(glob.glob(data_folder + "*"+channel+"R*.corr.int"))
file_list += sorted(glob.glob(data_folder + "*"+channel+"T*.corr.int"))

if real=="yes":
    file_list = []
#file_list = glob.glob(data_folder +  "*ANTO*"+channel+"*")
    file_list = sorted(glob.glob(data_folder +  "*"+channel+"Z"))
    file_list += sorted(glob.glob(data_folder + "*"+channel+"R"))
    file_list += sorted(glob.glob(data_folder + "*"+channel+"T"))



out = open(folder + "/picking_times_"+str(Tmin_p)+"_"+str(Tmax_p)+"_"+str(Tmin_s)+"_"+str(Tmax_s)+"_"+ str(Tmin_r)+"_"+str(Tmax_r)+".txt","w")
out.write("Period bands p: " + str(Tmin_p) + "\t" + str(Tmax_p)+"\n")
out.write("Period bands s: " + str(Tmin_s) + "\t" + str(Tmax_s)+"\n")
out.write("Period bands r: " + str(Tmin_r) + "\t" + str(Tmax_r)+"\n")
out.write("Ev lat: " + str(ev_lat)+"\n")
out.write("Ev lon: " + str(ev_lon)+"\n")
out.write("Ev dep: " + str(ev_dep)+"\n")
out.write("Station\tChannel\tstarttime\tbegin\torigintime\tp start\tp end\ts start\ts end\tr start\tr end\n")
for file in file_list:
	#try:
        tr = read(file)[0]

        tr_p = tr.copy()
        tr_p = filter_trace(tr_p, float(Tmin_p), float(Tmax_p))

        tr_s = tr.copy()
        tr_s = filter_trace(tr_s, float(Tmin_s), float(Tmax_s))

        tr_r = tr.copy()
        tr_r = filter_trace(tr_r, float(Tmin_r), float(Tmax_r))




        env_p = obspy.signal.filter.envelope(tr_p.data)
        env_s = obspy.signal.filter.envelope(tr_s.data)
        env_surface = obspy.signal.filter.envelope(tr_r.data)
        tr.stats['snr'] = []
        tr.stats['cutpoints_in_s'] = []
        tr.stats['cutpoints_noise_in_s'] = []
        tr.stats['phase'] = []
        starttime = tr.stats.starttime
        begin = float(tr.stats.sac.get("b"))
        origintime = starttime - begin
        endtime = tr.stats.endtime
        st_lat =  tr.stats["sac"]["stla"]
        st_lon =  tr.stats["sac"]["stlo"]

        station = tr.stats.station
        for i in range(0,len(sta_lines)):
            if station == sta_lines[i].split()[0]:
                print sta_lines[i]
                st_lat = float(sta_lines[i].split()[2])
                st_lon = float(sta_lines[i].split()[3])
		#print st_lat, st_lon

        delta = tr.stats.delta
        channel = tr.stats.channel
		#dist_km = geopy.distance.vincenty((ev_lat, ev_lon), (st_lat, st_lon)).km
		
#		result = client.distaz(stalat=st_lat, stalon=st_lon, evtlat=ev_lat,evtlon=ev_lon)
#		baz = result['backazimuth']
#		az = result['azimuth']
#		dist_km = result['distancemeters']/1000.
#		dist_deg = result['distance']
#		print "Distance (deg): ", dist_deg
        dist_deg, dist_km = distance(st_lat, st_lon, ev_lat, ev_lon)
        print "Distance (deg): ", dist_deg
        if dist_deg < 140. and dist_deg > 10.:

			#print station, channel, st_lat, st_lon, dist_km, dist_deg, baz, az


			#---
			# First rough p and s time picking using travel time in model iasp91
			#---
            arrivals_p = taup_model.get_travel_times(source_depth_in_km=ev_dep, phase_list=['P','Pdiff','pP','PcP','sP','PKP','PKS','PKKP','PKP','PS'],distance_in_degree=obspy.geodetics.kilometer2degrees(dist_km))
            p_arrival = tr.stats['p_arrival'] = int(arrivals_p[0].time)
            extra_p_arrival= []
            extra_p_arrival_name=[]
            for i in range(1,len(arrivals_p)):
                #print i
                extra_p_arrival.append(int(arrivals_p[i].time))
                #print extra_p_arrival
                extra_p_arrival_name.append ((arrivals_p[i].name))
            #print extra_p_arrival_name
           
                
            arrivals_s = taup_model.get_travel_times(source_depth_in_km=ev_dep, phase_list= ['S','Sdiff','pS','SP','sS','PS','SKS','SP','SKP'],distance_in_degree=obspy.geodetics.kilometer2degrees(dist_km))
            s_arrival = tr.stats['s_arrival'] = int(arrivals_s[0].time)
            extra_s_arrival= []
            extra_s_arrival_name=[]
            for i in range(1,len(arrivals_s)):
                #print i
                extra_s_arrival.append(int(arrivals_s[i].time))
                #print extra_s_arrival
                extra_s_arrival_name.append ((arrivals_s[i].name))
#print extra_s_arrival_name

#
#s_arrival1 = int(arrivals_s[1].time)
#s_arrival1_name =(arrivals_s[1].name)
#print s_arrival1_name
#s_arrival2 = int(arrivals_s[2].time)
#s_arrival2_name =(arrivals_s[2].name)
#print s_arrival2_name
#s_arrival3 = int(arrivals_s[3].time)
#s_arrival3_name =(arrivals_s[3].name)
#print s_arrival3_name
#s_arrival4 = int(arrivals_s[4].time)
#s_arrival4_name =(arrivals_s[4].name)
#print s_arrival4_name
#s_arrival5 = int(arrivals_s[5].time)
#s_arrival5_name =(arrivals_s[5].name)
#print s_arrival5_name

#s_arrival6 = int(arrivals_s[6].time)
#s_arrival6_name =(arrivals_s[6].name)
#print s_arrival6_name
                    
            
            
            
            
            p_start = abs(begin) + p_arrival - pretime_p
            p_end = abs(begin) + p_arrival + posttime_p
            i_p_start = int(p_start / delta)
            i_p_end = int(p_end / delta)

            s_start = abs(begin) + s_arrival - pretime_s
            s_end = abs(begin) + s_arrival + posttime_s
            i_s_start = int(s_start / delta)
            i_s_end = int(s_end / delta)




			#---
			# Refine p and s picking using envelope inside the time window 
			#---
            if i_p_start < 0:
                i_p_start = 0
            if i_s_start < 0:
                i_s_start = 0

            p_time = np.argmax(abs(tr_p.data)[i_p_start:i_p_end])*delta + p_start
            p_max = max(abs(tr_p.data)[i_p_start:i_p_end])
            p_start_final = p_time - pretime_p
            p_end_final = p_time + posttime_p

            s_time = np.argmax(abs(tr_s.data)[i_s_start:i_s_end])*delta + s_start
            s_max = max(abs(tr_s.data)[i_s_start:i_s_end])
            s_start_final = s_time - pretime_s
            s_end_final = s_time + posttime_s



			#---
			# Surface picking using envelope
			#---
			# rough determination of arrival time with min and max vs velocity
            r_start1 = abs(begin) + dist_km / 4.9
            r_start2 = abs(begin) + dist_km / 3.0

            ir_start1 = int(r_start1 / delta)
            ir_start2 = int(r_start2 / delta)

            time_rayleigh = r_start1+np.argmax(env_surface[ir_start1:ir_start2])*delta

            r_start_final = s_end_final


            for i in reversed(range(0,int(time_rayleigh/delta)  )):
                max_env = max(env_surface)
                if env_surface[i] <= 0.1 * max_env:
                    r_start_final = i*delta
                    break

            for i in range(int(time_rayleigh/delta) , len(env_surface) ):
                max_env = max(env_surface)
                if env_surface[i] <= 0.02 * max_env:
                    r_end_final = i*delta
                    break

			# if the starting time has not been found, then put the time window equals to 0
            if r_start_final == s_end_final:
                r_start_final = r_end_final


			#---
			# convert relative times to absolute
			#---
            p_start_final_abs =  starttime + p_start_final
            p_end_final_abs =  starttime + p_end_final
            s_start_final_abs =  starttime + s_start_final
            s_end_final_abs =  starttime + s_end_final
            r_start_final_abs =  starttime + r_start_final
            r_end_final_abs =  starttime + r_end_final

			#---
			# Write the picking time into a file
			#---
			#print file, tr.stats.starttime
			#sys.exit()

            line = station + "\t" + channel + "\t" + str(starttime) + "\t" + str(begin) + "\t" + str(origintime) + "\t"+ str(p_start_final_abs) + 					"\t" + str(p_end_final_abs) +  "\t" + str(s_start_final_abs) + 					"\t" + str(s_end_final_abs) + "\t" + str(r_start_final_abs) + "\t" + str(r_end_final_abs) + "\n"
            out.write(line)


			#======================================================
			# cut traces
            extratime = 800.
            tr_cut_p = trim_trace_abs(tr_p, tr_p.stats.starttime, p_start_final_abs, p_end_final_abs, Tmax_p, extratime)
            tr_cut_s = trim_trace_abs(tr_s, tr_s.stats.starttime, s_start_final_abs, s_end_final_abs, Tmax_s, extratime)
            tr_cut_r = trim_trace_abs(tr_r, tr_r.stats.starttime, r_start_final_abs, r_end_final_abs, Tmax_r, extratime)

		#	tr_cut_p = trim_trace_abs(tr_p, tr_p.stats.starttime, p_start_final_abs, p_end_final_abs, Tmax_p, extratime)
		#	tr_cut_s = trim_trace_abs(tr_s, tr_s.stats.starttime, s_start_final_abs, s_end_final_abs, Tmax_s, extratime)
		#	tr_cut_r = trim_trace_abs(tr_r, tr_r.stats.starttime, r_start_final_abs, r_end_final_abs, Tmax_r, extratime)

			
            plt.figure(1, figsize=(11.69, 8.27))
            plt.subplot(311)
            plt.plot(tr_p.times(), tr_p.data, color="black")
            plt.plot(tr_cut_p.times()+p_start_final-extratime, tr_cut_p.data, color="red", linewidth=2)
            plt.axvline(origintime + p_arrival-50, color="gray",linewidth=1)
            plt.axvline(origintime + p_arrival+100,color="gray",linewidth=1)
            	#	plt.axvline(pP_arrival,color = "blue",linestyle=":",linewidth=1)
            	#	plt.axvline(PcP_arrival,color = "cyan",linestyle=":",linewidth=1)
            		#plt.axvline(sP_arrival,color = "magenta",linestyle=":",linewidth=1 )
            plt.ylim(-p_max*1.1, p_max*1.1)
            plt.xlim(p_arrival-400,  r_end_final+2000)
            plt.xlim(0,6000)
            plt.axvline(p_time, color="red",linestyle=":")
            plt.axvline(p_start_final, color="red",linewidth=2)
            plt.axvline(p_end_final, color="red",linewidth=2)
            plt.ylabel("P waves")
            for i in range(1,len(extra_p_arrival)):
                plt.axvline(extra_p_arrival[i],color = "blue",linestyle=":",linewidth=2.5)
                plt.text(extra_p_arrival[i],p_max*1.5, str(extra_p_arrival_name[i]), rotation=90)
                print extra_p_arrival_name[i]
			

			#phases = [ "
            #for arrivals_p_sub in arrivals_p:
#plt.axvline(int(arrivals_p_sub.time), linewidth = 1, color = 'g')

		#	[plt.axvline(arrivals_p_sub,linewidth = 1, color = 'g') for arrivals_p_sub in arrivals_p]

            plt.subplot(312)
            plt.plot(tr_s.times(), tr_s.data, color="black",zorder=0)
            plt.plot(tr_cut_s.times()+s_start_final-extratime, tr_cut_s.data, color="red", linewidth=2)
		#	plt.plot(tr_cut_s.times(), tr_cut_s.data, color="red",zorder=2)
            plt.axvline(s_arrival-50, color="gray",linewidth=1)
            plt.axvline(s_arrival+150, color="gray",linewidth=1)
            plt.xlim(p_arrival-400,  r_end_final+2000)
            plt.ylim(-s_max*1.1, s_max*1.1)
            plt.axvline(s_time, color="red",linestyle=":")
            plt.axvline(s_start_final, color="red",linewidth=2)
            plt.axvline(s_end_final, color="red",linewidth=2)
            plt.ylabel("S waves")
            for i in range(1,len(extra_s_arrival)):
                plt.axvline(extra_s_arrival[i],color = "blue",linestyle=":",linewidth=2.5)
                plt.text(extra_s_arrival[i],s_max*1.2, str(extra_s_arrival_name[i]), rotation=90)
                print extra_s_arrival_name[i]
            	#	[plt.axvline(arrivals_s_sub,linewidth = 1, color = 'g') for arrivals_s_sub in arrivals_s]

#for arrivals_s_sub in arrivals_s:                                               plt.axvline(int(arrivals_s_sub.time), linewidth = 1, color = 'g')
            

            plt.subplot(313)
            plt.plot(tr_r.times(), tr_r.data, color="black")
            plt.plot(tr_cut_r.times()+r_start_final-extratime, tr_cut_r.data, color="red", linewidth=2)
            plt.axvline(r_start_final, color="red",linewidth=2)
            plt.axvline(r_end_final, color="red",linewidth=2)
            plt.plot(tr_r.times(), env_surface, color="black", linestyle=":")
            plt.xlim(p_arrival-400,  r_end_final+2000)
            plt.axvline(r_start1, color="black", linestyle="--")
            plt.axvline(r_start2, color="black", linestyle="--")
            plt.ylabel("surf. waves")



            plt.suptitle("Station: " + station + " Channel: " + channel)
            plt.savefig(plot_folder + "/" +station + "_" + channel + "_picking.png")
            plt.close()
			#sys.exit()
	#except:
	#	pass
out.close()
#	sys.exit()






