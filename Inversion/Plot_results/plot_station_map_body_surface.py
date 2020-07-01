import os
import matplotlib.pylab as plt
import numpy as np
import sys
from matplotlib import cm
import glob
from mpl_toolkits.basemap import Basemap
from operator import itemgetter
import math
from obspy.imaging.beachball import beach

def coord_from_dist_angle(lon1, lat1, d, brng):
	event_location = [lon1, lat1]
	#p1 = [lat1, lon1]
	#dist =  vincenty(p1,p2, miles=False)

	lat1 = np.deg2rad(event_location[1])
	lon1 = np.deg2rad(event_location[0])

	R = 6371.0 #Radius of the Earth
	brng = np.deg2rad(brng) #Bearing is 90 degrees converted to radians.
	#d = 15 #Distance in km

	#lat2  52.20444 - the lat result I'm hoping for
	#lon2  0.36056 - the long result I'm hoping for.

	#lat1 = math.radians(52.20472) #Current lat point converted to radians
	#lon1 = math.radians(0.14056) #Current long point converted to radians

	lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
	     math.cos(lat1)*math.sin(d/R)*math.cos(brng))

	lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
	             math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

	lat2 = round(np.rad2deg(lat2),3)
	lon2 = round(np.rad2deg(lon2),3)

	return lon2, lat2


R = 6371.0

results_folder = sys.argv[1]
observed_folder = results_folder + "/observed/"

station_comp = []
filelist = glob.glob(observed_folder + "*_observed.asc")
for i in range(0,len(filelist)):
	station = filelist[i].split("/")[-1].split("_")[0]
	comp = filelist[i].split("/")[-1].split("_")[1]
	station_comp.append([station, comp])


ev_lon = float(sys.argv[3])
ev_lat = float(sys.argv[2])
#ev_lon = 153.43
#ev_lat = 52.71

station_file = "../../STATIONS"
station_list_avail = []
for i in range(0,len(station_comp)):
	station = station_comp[i][0]
	lines = open(station_file).readlines()
	#print station
	for i in range(0,len(lines)):

		if station == lines[i].split()[0]:
			lat = float(lines[i].split()[2])
			lon = float(lines[i].split()[3])
			station_list_avail.append([station, lon, lat])	

fig,ax = plt.subplots()

#plt.figure(figsize=[10,10])
map = Basemap(projection='aeqd',lat_0=ev_lat,lon_0=ev_lon,resolution='l')
#map = Basemap(projection='robin',lon_0=ev_lon,resolution='c')
map.drawcoastlines(linewidth=0.5)
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))
map.fillcontinents(color='0.8')



map.scatter(ev_lon, ev_lat, latlon=True, s=400, color="red", edgecolor="black",linewidth=1, marker="*", zorder=10)

for i in range(0,len(station_list_avail)):
	station = station_list_avail[i][0]
	st_lat = float(station_list_avail[i][2])
	st_lon = float(station_list_avail[i][1])
	map.scatter(st_lon, st_lat, marker="^", latlon=True, facecolor="lime", edgecolor="black", linewidth=1, s=200, zorder=10)
	#plt.annotate(code, xy=(x-250000, y-300000), color="red", size=7.5, fontweight="bold")






plt.savefig(results_folder + "/Station_map.png")
plt.savefig(results_folder + "/Station_map.eps")
plt.close()


#======================================================

fig = plt.figure(1,figsize=[15,15])
ax = fig.add_subplot(111)
map = Basemap(projection='aeqd',lat_0=ev_lat,lon_0=ev_lon,resolution='c')
map.drawcoastlines(linewidth=1)
map.fillcontinents(color='wheat',lake_color='aqua')
map.drawmapboundary(fill_color='skyblue')

llon2, llat2 = [],[]
llon3, llat3 = [],[]
for d in range(0,360):
	dist_km = R*math.radians(40.)
	lon2, lat2 = coord_from_dist_angle(ev_lon, ev_lat, dist_km, d)

	dist_km = R*math.radians(120.)
	lon3, lat3 = coord_from_dist_angle(ev_lon, ev_lat, dist_km, d)

	llon2.append(lon2)
	llat2.append(lat2)
	llon3.append(lon3)
	llat3.append(lat3)
map.plot(llon2, llat2, latlon=True,  zorder=20, color="red", linestyle="--", linewidth=4)
map.plot(llon3, llat3, latlon=True,  zorder=20, color="red", linestyle="--", linewidth=4)



for i in range(0,len(station_list_avail)):
	station = station_list_avail[i][0]
	st_lat = float(station_list_avail[i][2])
	st_lon = float(station_list_avail[i][1])
	map.scatter(st_lon, st_lat, marker="^", latlon=True, facecolor="lime", edgecolor="black", linewidth=2, s=400, zorder=10)
	#plt.annotate(code, xy=(x-250000, y-300000), color="red", size=7.5, fontweight="bold")

fm = [-6.29e26, 5.37e24, 6.24e26, -1.38e26, -6.37e25, 2.29e25]
x, y = map(-120,5)
beach1 = beach(fm, xy=(x,y), width=5000000, facecolor="C4")
beach1.set_zorder(100)
ax = plt.gca()
ax.add_collection(beach1) 

map.plot([-111,5],[ev_lon,ev_lat])
map.scatter(ev_lon, ev_lat, latlon=True, s=800, color="white", edgecolor="red",linewidth=3, marker="*", zorder=10)
plt.savefig(results_folder + "/Azimuthal_coverage.png")
plt.close()

#=======================================================================
#++++++++++++++++++++++++++++++++++++++++++

results_file = results_folder + "/station_mft.asc"
lines = open(results_file).readlines()
mmft, iindex = [],[]
for i in range(0,len(lines)):
	if lines[i].split()[0] == "Tot":
		mft = float(lines[i].split()[2])
		mmft.append(mft)
		iindex.append(i)

i_min = min(enumerate(mmft), key=itemgetter(1))[0]
line_best = iindex[i_min]



for c in ["Z"]:
	plt.figure(1,figsize=[10,10])
	map = Basemap(projection='aeqd',lat_0=ev_lat,lon_0=ev_lon,resolution='l')
	#map = Basemap(projection='robin',lon_0=ev_lon,resolution='c')
	map.drawcoastlines(linewidth=0.2)
	map.scatter(ev_lon, ev_lat, latlon=True, s=400, color="red", edgecolor="black",linewidth=1, marker="*", zorder=10)

	ssmft, ssta, llon, llat = [],[],[],[]
	for n in range(2,1000):
		if lines[line_best - n].split()[0][0]!="*":
			sta = lines[line_best - n].split()[1].split("_")[0]
			comp = lines[line_best - n].split()[1].split("_")[1].split("ff")[0]
			s_mft = float(lines[line_best - n].split()[3])
			if comp == c:
				for j in range(0,len(station_list_avail)):
					if sta == station_list_avail[j][0]:
						st_lon = float(station_list_avail[j][1])
						st_lat = float(station_list_avail[j][2])
						llon.append(st_lon)
						llat.append(st_lat)
						ssmft.append(float(s_mft))
						ssta.append(sta)
						x, y = map(st_lon, st_lat)
						plt.annotate(sta, xy=(x-500000,y-1000000), color="black",
				 			size=5, fontweight="bold")

						break
		
			

				

		else:
			break

	print "CICCIO" , min(ssmft)
#	sys.exit()
	#print min(ssmft)
	#print "------------"
	#print min(ssmft)*10.
	s = map.scatter(llon, llat, marker="^", latlon=True,
							c=ssmft, linewidth=1,
							s=100, zorder=10,
							vmin = min(ssmft), vmax=min(ssmft)*20,
							cmap=cm.jet)
	


	plt.colorbar(s, label="misfit")
	plt.suptitle("Station misfit " + c + " component")
	plt.savefig(results_folder + "/"+c+"_misfit_map.png", dpi=200)
	plt.close()





