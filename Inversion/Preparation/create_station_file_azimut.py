# -*- coding: utf-8 -*-
import os
import sys
import glob
sys.path.insert(0, '../../lib/')
from functions import *
import numpy as np
from matplotlib import cm
from adjustText import adjust_text


def get_station_coords(station):
	station_list = []
	lines = open("../../STATIONS").readlines()
	for i in range(0,len(lines)):
		sta = lines[i].split()[0]
		st_lat = float(lines[i].split()[2])
		st_lon = float(lines[i].split()[3])
		station_list.append([sta, st_lat, st_lon])


	for j in range(0,len(station_list)):
		if station == station_list[j][0]:
			st_lat = station_list[j][1]
			st_lon = station_list[j][2]
			break

	return st_lat, st_lon

def use_it_or_not(station, comp, wavetype, data2use):
	use_it = "N"
	for i in range(0,len(data2use)):
			if station == data2use[i][0] and comp == data2use[i][1] and wavetype == data2use[i][2]:
				use_it = data2use[i][3]
				break
	return use_it

def calculate_initial_compass_bearing(pointA, pointB):

    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")
    
    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
                                       * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360
    
    return compass_bearing


folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/database/CMTSOLUTION_201810131110A_GCMT/"
cmt_file =  folder + "CMTSOLUTION_201810131110A_GCMT"
use_selection_file = "yes"


dist_min_P = 20.
dist_max_P = 90.
dist_min_S = 20.
dist_max_S = 90.
dist_min_W = 40.
dist_max_W = 120.

selection_file = folder + "fortran_format_18_40_18_60_45_100/station2use.txt"

use_selection_file = "yes"







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

use_it = "N"

for wt in ["P","S"]:
	print "Preparing " + wt + " waves"
	bbz = []
	if wt == "P":
		dist_min = dist_min_P
		dist_max = dist_max_P
	if wt == "S":
		dist_min = dist_min_S
		dist_max = dist_max_S
        #if wt == "W":
        #dist_min = dist_min_W
#dist_max = dist_max_W

	bbaz = []
	ddist = []
	ss_list = []
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	ax.set_theta_zero_location('W',offset=-90)
	ax.set_theta_direction(-1)
	for fl in sorted(filelist):
		station = fl.split("/")[-1].split("_")[0]
		wavetype = fl.split("/")[-1].split("_")[2]
		data_name = fl.split("/")[-1].split("_ff")[0]
		comp = fl.split("/")[-1].split("_")[1]

		st_lat, st_lon = get_station_coords(station)
		dist = distance(st_lat, st_lon, ev_lat, ev_lon)[0]
		st = (st_lat, st_lon)
		ev = (ev_lat, ev_lon)
		baz = calculate_initial_compass_bearing(ev, st)
		use_it = use_it_or_not(station, comp, wavetype, data2use)
		#use_it = "Y"

		if wavetype == wt and use_it == "Y" and dist >= dist_min and dist <= dist_max:
			#if station not in ss_list:
				ss_list.append(station)
				bbaz.append(baz)
				ddist.append(dist)
				ax.scatter(math.radians(baz), dist, marker="^", c="g",s=70,zorder = 2)
				fig.suptitle( wt + " waves")
				ax.set_ylim(0,90)
				ax.set_yticks(np.arange(10,90,10))
				#ax.annotate(station, xy=(math.radians(baz), dist))
				n+=1


	baz_bins = np.arange(0,370,10.)
	counts, bins = np.histogram(bbaz, bins=baz_bins)

	plt.savefig("az_coverageTEST_"+wt+".pdf")
	plt.close()


	plt.hist(ddist, bins=100)
	plt.savefig("dist_hist"+wt+".png")
	plt.close()

	plt.hist(bbaz, bins=baz_bins)
	for i in range(0,len(counts)):
		b = (bins[i] + bins[i+1])/2.
		plt.scatter(b, counts[i])
	plt.savefig("baz_hist_"+wt+".png")
	plt.close()


	fig = plt.figure()
	ax2 = fig.add_subplot(111, projection='polar')
	ax2.set_theta_zero_location('W',offset=-90)
	ax2.set_theta_direction(-1)
	ax2.set_xticks([0,3.14/2., 3.14, 3*3.14/2.])
	for fl in sorted(filelist):
		data_name = fl.split("/")[-1]
		station = fl.split("/")[-1].split("_")[0]
		wavetype = fl.split("/")[-1].split("_")[2]
		data_name = fl.split("/")[-1].split("_ff")[0]
		comp = fl.split("/")[-1].split("_")[1]

		st_lat, st_lon = get_station_coords(station)
		dist = distance(st_lat, st_lon, ev_lat, ev_lon)[0]
		if wt == "P":
			dist_min = dist_min_P
			dist_max = dist_max_P
		if wt == "S":
			dist_min = dist_min_S
			dist_max = dist_max_S
		if wt == "W":
			dist_min = dist_min_W
			dist_max = dist_max_W

		st = (st_lat, st_lon)
		ev = (ev_lat, ev_lon)
		baz = calculate_initial_compass_bearing(ev, st)
		use_it = use_it_or_not(station, comp, wavetype, data2use)
	#	use_it = "Y"
		
		
		if wavetype == wt and use_it == "Y" and dist >= dist_min and dist <= dist_max:
			for i in range(0,len(bins)-1):
				baz_min = bins[i]
				baz_max = bins[i+1]
				

				if baz_min not in bbz:
					bbz.append(baz_min)
					plt.axvline(math.radians(baz_min), color="0.8", zorder=0)

				if baz >= baz_min and baz <= baz_max:
					density_sta = 1/float(counts[i])
					print counts[i], density_sta
					break

			#print station, comp, wavetype, dist, baz, density_sta
			c = plt.scatter(math.radians(baz), dist, c=density_sta, vmin=0.0, vmax=1.0, cmap=cm.jet, zorder=10, marker="^", s=40)
			ax2.annotate(station, xy=(math.radians(baz), dist))
			out.write(data_name + "\n")
			weight = str(round(density_sta,3))
			out.write(weight+ "\n")
	plt.ylim(0,dist_max)

	plt.colorbar(c)	

	plt.savefig("tmp_"+wt+".png")
	plt.close()

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


os.system("cp rfi.in ../NA/data/rfi_files")

