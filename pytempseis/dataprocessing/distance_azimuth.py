from obspy.core import read
from obspy.taup import TauPyModel
import glob
import obspy
from pytempseis.functions import distance, filter_trace


real = True

event_code = "CMTSOLUTION_201809061549A_GCMT"

folder = f"../database/{event_code}"
cmt_file = f"{folder}/{event_code}"

channel = "BH" if real else "MX"
data_folder = f"{folder}/data_ready2use/" if real else f"{folder}/processed_data/"

Tmin_p = 20  # 25
Tmax_p = 70  # 60

Tmin_s = 20
Tmax_s = 100

Tmin_r = 45  # 125
Tmax_r = 100  # 180

pretime_p = 100
posttime_p = 100

pretime_s = 100
posttime_s = 100

pretime_r = 2 * Tmax_r  # 500
posttime_r = 4 * Tmax_r  # 500


lines = open(cmt_file).readlines()
ev_lat = float(lines[4].split()[1])
ev_lon = float(lines[5].split()[1])
ev_dep = float(lines[6].split()[1])


taup_model = TauPyModel(model="iasp91")

sta_lines = open("../STATIONS").readlines()

if real:
    file_list = []
    selected_stations = [
        "AIS",
        "ERM",
        "INCN",
        "JOHN",
        "KDAK",
        "MA2",
        "MBWA",
        "PMG",
        "TAOE",
        "TAU",
        "WAKE",
        "YAK",
        "YSS",
    ]
    for i in range(0, len(selected_stations)):
        station = str(selected_stations[i])
        print(station)
        file_list += glob.glob(f"{data_folder}*{station}*{channel}Z")


for file in file_list:
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
    tr.stats["snr"] = []
    tr.stats["cutpoints_in_s"] = []
    tr.stats["cutpoints_noise_in_s"] = []
    tr.stats["phase"] = []
    starttime = tr.stats.starttime
    begin = float(tr.stats.sac.get("b"))
    origintime = starttime - begin
    endtime = tr.stats.endtime
    st_lat = tr.stats["sac"]["stla"]
    st_lon = tr.stats["sac"]["stlo"]

    station = tr.stats.station
    for i in range(0, len(sta_lines)):
        if station == sta_lines[i].split()[0]:
            print(sta_lines[i])
            st_lat = float(sta_lines[i].split()[2])
            st_lon = float(sta_lines[i].split()[3])

    delta = tr.stats.delta
    channel = tr.stats.channel

    dist_deg, dist_km = distance(st_lat, st_lon, ev_lat, ev_lon)
    distance_in_degree = obspy.geodetics.kilometer2degrees(dist_km)
    print(file)
    print(f"\nDistance (deg): {dist_deg}")
    print(f"\nDistance (km): {dist_km}")
