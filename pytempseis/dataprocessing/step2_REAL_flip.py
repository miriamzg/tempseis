import os
import sys
from obspy.core import read
import glob
import matplotlib.pylab as plt
sys.path.insert(0, '../lib/')

event_code = sys.argv[1]
data_folder = "../database/"+event_code+"/processed_data/"
output_folder = "../database/"+event_code+"/data_ready2use/"
channel = "BH"




os.system("mkdir " + output_folder)
os.system("rm " + output_folder + "*.sac")



sta_lines = open("../database/"+event_code + "/first_check.txt").readlines()
station_list = []
for i in range(1,len(sta_lines)):
	station = sta_lines[i].split()[0]
	comp = sta_lines[i].split()[1]
	status = sta_lines[i].split()[2]
	station_list.append([station, comp, status])


for i in range(0,len(station_list)):
	try: 
		station = station_list[i][0]
		comp = station_list[i][1]
		status = station_list[i][2]
		
		print station, comp, status
		
		filename = glob.glob(data_folder + "*" + station + "*.00."+channel+comp)[0]
	#	filename = glob.glob(data_folder+"/*" + station + "*MX" + comp + ".sem.sac_post_processed_shift.noisy")[0]


		tr = read(filename)[0]
		tri = tr.copy()
		

		


		if status == "I":
			for j in range(0,len(tr.data)):
				tri.data[j] = -tr.data[j]
			tri.write(output_folder + tr.id, format="SAC")

		if status == "Y":
			tr.write(output_folder + tr.id, format="SAC")
		
		if status == "X":
			pass

	except IndexError:
		pass


	



	
	
	

	
	











