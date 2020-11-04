import os
import sys
import glob
from obspy.core import read
from obspy.clients.iris import Client

event_code = sys.argv[1]
database = sys.argv[2]


data_path = f"{database}/{event_code}/synthetics/point_source/"
client = Client()


ev_lon = 140.5600
ev_lat = 27.9400
ev_dep = 680.7000


filelist = glob.glob(f"{data_path}/*MXZ*.sac")


for fl in filelist:
    Zfile = fl
    print(Zfile)
    Nfile = Zfile.replace("MXZ", "MXN")
    Efile = Zfile.replace("MXZ", "MXE")
    Rfile = Zfile.replace("MXZ", "MXR")
    Tfile = Zfile.replace("MXZ", "MXT")

    if os.path.isfile(Nfile) and os.path.isfile(Efile):
        st = read(Zfile)
        st += read(Nfile)
        st += read(Efile)

        Ztr = st[0]
        st_lat = Ztr.stats["sac"]["stla"]
        st_lon = Ztr.stats["sac"]["stlo"]
        station = Ztr.stats.station
        delta = Ztr.stats.delta
        channel = Ztr.stats.channel

        result = client.distaz(
            stalat=st_lat, stalon=st_lon, evtlat=ev_lat, evtlon=ev_lon
        )
        baz = result["backazimuth"]
        az = result["azimuth"]
        dist_km = result["distancemeters"] / 1000.0
        dist_deg = result["distance"]
        print(station, channel, st_lat, st_lon, dist_km, dist_deg, baz, az)

        line = (
            "sac > macro rotate.macro "
            + Nfile
            + " "
            + Efile
            + " "
            + str(baz)
            + " "
            + Rfile
            + " "
            + Tfile
        )
        os.system(line)

        Rtr = read(Rfile)[0]
        Rtr.stats["channel"] = "MXR"
        Rtr.write(Rfile, format="SAC")

        Ttr = read(Tfile)[0]
        Ttr.stats["channel"] = "MXT"
        Ttr.write(Tfile, format="SAC")
