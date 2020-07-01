import os
import sys
import glob
sys.path.insert(0, '../../lib/')
from functions import *


#folder = "../../database/CMTSOLUTION_201310251710A_GCMT/"
folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/database/CMTSOLUTION_201310251710A_GCMT/"
#cmt_file =  folder + "CMTSOLUTION_201310251710A_GCMT"
cmt_file =  folder + "CMTSOLUTION_201310251710A_GCMT"
use_selection_file = "yes"


selection_file = folder + "fortran_format_25_60_25_100_45_100/station2use.txt"






station_list = []
lines = open("../../STATIONS").readlines()
for i in range(0,len(lines)):
	station = lines[i].split()[0]
	st_lat = float(lines[i].split()[2])
	st_lon = float(lines[i].split()[3])
	station_list.append([station, st_lat, st_lon])


lines = open(cmt_file).readlines()
ev_lat = float(lines[4].split()[1])
ev_lon = float(lines[5].split()[1])


stations2skip = [] #["GAR", "INCN"]



if use_selection_file == "yes":
	data2use = []
	lines = open(selection_file).readlines()
	for i in range(1,len(lines)):
		sta = lines[i].split()[0]
		comp = lines[i].split()[1]
		wavetype = lines[i].split()[2]
		use = lines[i].split()[3]
		data2use.append([sta, comp, wavetype, use])



#=====================================================================
#	write station file
#=====================================================================
print "Writing station file..."
filelist = glob.glob("../NA/data/rfi_files/OBS/*")
out = open("rfi.in","w")
out.write("#\n# Input file for receiver function inversion specific information\n#\n")
out.write("rfi_param                              /* input model   */\n")
out.write("rfi_models                             /* output models */\n")
out.write(str(len(filelist))+"\t\t\t/* nwave */\n")
n = 0
use_it = "Y"
for fl in filelist:
	station = fl.split("/")[-1].split("_")[0]
	wavetype = fl.split("/")[-1].split("_")[2]
	data_name = fl.split("/")[-1].split("_ff")[0]
	comp = fl.split("/")[-1].split("_")[1]
	print data_name, comp, wavetype


	for i in range(0,len(data2use)):
		if station == data2use[i][0] and comp == data2use[i][1] and wavetype == data2use[i][2]:
			use_it = data2use[i][3]
			break
	
	if use_it == "Y":

		for i in range(0,len(station_list)):
			if station == station_list[i][0]:
				st_lat = station_list[i][1]
				st_lon = station_list[i][2]
				break
		
		dist = distance(st_lat, st_lon, ev_lat, ev_lon)[0]
		print st_lat, st_lon, dist
		if station not in stations2skip:
		#		if (wavetype == "P" or wavetype == "S") and dist > 10. and dist < 50.:
                        if (wavetype == "P" or wavetype == "S"):
				print "yes"
				weight = "1.0"
				out.write(data_name + "\n")
				out.write(weight+ "\n")
				n += 1

		#	if wavetype == "S" and comp == "Z":
		#		weight = "1.0"
		#		out.write(data_name + "\n")
		#		out.write(weight+ "\n")
		#		n += 1
		#	if wavetype == "W" and dist > 20. and dist < 60.:  # 40 120
                        if wavetype == "W":
				weight = "1.0"
				print "yes"
				out.write(data_name + "\n")
				out.write(weight+ "\n")
				n += 1
			else: 
				print "no"

out.write("1                                       /* iwrite_models */")
out.close()


out = open("rfi.in_new","w")
lines = open("rfi.in","r").readlines()
for i in range(0,len(lines)):
	
	if i == 5:
		out.write(str(n)+"\t\t\t/* nwave */\n")
	else:
		out.write(lines[i])
out.close()
os.system("mv rfi.in_new rfi.in")



