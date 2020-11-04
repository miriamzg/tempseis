import os
import sys
import glob
from obspy.core import read
from obspy.clients.iris import Client
from pytempseis import _lib


def rotate_traces(event_code, database):
    data_path = f"{database}/{event_code}/synthetics/point_source/"
    filelist = glob.glob(f"{data_path}/*MXZ*.sac")

    cmt_finite_fault = f"{database}/{event_code}/{event_code}"
    cmt_lines = open(cmt_finite_fault).readlines()
    for line in cmt_lines:
        if line.split()[0] == "latitude":
            ev_lat = float(line.split()[1])
        elif line.split()[0] == "longitude":
            ev_lon = float(line.split()[1])
        else:
            continue

    client = Client()
    for Zfile in filelist:
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
            channel = Ztr.stats.channel

            result = client.distaz(
                stalat=st_lat, stalon=st_lon, evtlat=ev_lat, evtlon=ev_lon
            )
            baz = result["backazimuth"]
            az = result["azimuth"]
            dist_km = result["distancemeters"] / 1000.0
            dist_deg = result["distance"]
            print(station, channel, st_lat, st_lon, dist_km, dist_deg, baz, az)

            os.system(
                f"sac > macro {_lib}/rotate.macro {Nfile} {Efile} {baz} {Rfile} {Tfile}"
            )

            Rtr = read(Rfile)[0]
            Rtr.stats["channel"] = "MXR"
            Rtr.write(Rfile, format="SAC")

            Ttr = read(Tfile)[0]
            Ttr.stats["channel"] = "MXT"
            Ttr.write(Tfile, format="SAC")


if __name__ == "__main__":
    event_code = sys.argv[1]
    database = sys.argv[2]
    rotate_traces(event_code, database)
