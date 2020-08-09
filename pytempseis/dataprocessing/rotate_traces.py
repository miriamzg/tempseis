import os
import sys
import glob
import pdb
from obspy.core import read
import obspy
from obspy.clients.iris import Client

client = Client()



pdb.set_trace()
data_path = "../Users/TheStuffofAlice/Documents/UCL/4th_year/Masters_Project/TEMPSEIS_package_v1.2_WORKING/database/CMTSOLUTION_201505301123A_GCMT/synthetics/point_source/"


ev_lon = 140.5600
ev_lat = 27.9400
ev_dep = 680.7000


filelist = glob.glob("/Users/TheStuffofAlice/Documents/UCL/4th_year/Masters_Project/TEMPSEIS_package_v1.2_WORKING/database/CMTSOLUTION_201505301123A_GCMT/synthetics/point_source/"+"/*MXZ*.sac")


for fl in filelist:
	Zfile = fl
	print Zfile
 	Nfile = Zfile.replace("MXZ","MXN")
 	Efile = Zfile.replace("MXZ","MXE")
 	Rfile = Zfile.replace("MXZ","MXR")
 	Tfile = Zfile.replace("MXZ","MXT")


	if os.path.isfile(Nfile) and os.path.isfile(Efile):
		st = read(Zfile)
		st += read(Nfile)
		st += read(Efile)

		#print Zfile, Rfile, Tfile
		Ztr = st[0]
		st_lat =  Ztr.stats["sac"]["stla"]
		st_lon =  Ztr.stats["sac"]["stlo"]
		station = Ztr.stats.station
		delta = Ztr.stats.delta
		channel = Ztr.stats.channel
		#dist_km = geopy.distance.vincenty((ev_lat, ev_lon), (st_lat, st_lon)).km
		
		result = client.distaz(stalat=st_lat, stalon=st_lon, evtlat=ev_lat,evtlon=ev_lon)
		baz = result['backazimuth']
		az = result['azimuth']
		dist_km = result['distancemeters']/1000.
		dist_deg = result['distance']
		print station, channel, st_lat, st_lon, dist_km, dist_deg, baz, az


		line = "sac > macro rotate.macro " + Nfile + " " + Efile + " " + str(baz) + " " + Rfile + " " + Tfile
#		print line
		os.system(line)

		Rtr = read(Rfile)[0]
		Rtr.stats["channel"] = "MXR"
		Rtr.write(Rfile, format="SAC")

		Ttr = read(Tfile)[0]
		Ttr.stats["channel"] = "MXT"
		Ttr.write(Tfile, format="SAC")

		#sys.exit()
