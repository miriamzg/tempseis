import os
import glob
from pytempseis.functions import distance
import numpy as np
from argparse import ArgumentParser


def get_station_coords(station, stations_file):
    station_list = []
    lines = open(stations_file).readlines()
    for i in range(0, len(lines)):
        sta = lines[i].split()[0]
        st_lat = float(lines[i].split()[2])
        st_lon = float(lines[i].split()[3])
        station_list.append([sta, st_lat, st_lon])

    for j in range(0, len(station_list)):
        if station == station_list[j][0]:
            st_lat = station_list[j][1]
            st_lon = station_list[j][2]
            break

    return st_lat, st_lon


def use_it_or_not(station, comp, wavetype, data2use):
    use_it = False
    for i in range(0, len(data2use)):
        if (
            station == data2use[i][0]
            and comp == data2use[i][1]
            and wavetype == data2use[i][2]
        ):
            use_it = data2use[i][3]
            break
    return False if use_it == "N" else True


def calculate_initial_compass_bearing(pointA, pointB):

    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = np.radians(pointA[0])
    lat2 = np.radians(pointB[0])
    diffLong = np.radians(pointB[1] - pointA[1])

    x = np.sin(diffLong) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (
        np.sin(lat1) * np.cos(lat2) * np.cos(diffLong)
    )

    initial_bearing = np.arctan2(x, y)
    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


parser = ArgumentParser(description="Builds station plots and homti.in")
parser.add_argument("database", type=str, help="Path to database")
parser.add_argument("event_code", type=str, help="GCMT event code")
parser.add_argument("filtering", type=str)
args = parser.parse_args()

folder = os.path.join(args.database, args.event_code)
cmt_file = os.path.join(folder, args.event_code)
stations_file = os.path.join(args.database, "STATIONS")

dist_min_P = 20.0
dist_max_P = 90.0
dist_min_S = 20.0
dist_max_S = 90.0
dist_min_W = 40.0
dist_max_W = 120.0

lines = open(cmt_file).readlines()
ev_lat = float(lines[4].split()[1])
ev_lon = float(lines[5].split()[1])

selection_file = os.path.join(folder, args.filtering, "station2use.txt")
data2use = []
lines = open(selection_file).readlines()
for i in range(1, len(lines)):
    sta = lines[i].split()[0]
    comp = lines[i].split()[1]
    wavetype = lines[i].split()[2]
    use = lines[i].split()[3]
    data2use.append([sta, comp, wavetype, use])


# =====================================================================
# 	write station file
# =====================================================================
print("Writing station file...")
filelist = glob.glob(os.path.join(folder, args.filtering, "observed_data", "*"))
outlines = []
outlines.append("#\n# Input file for higher order MT inversion specific information\n#\n")
outlines.append("homti_param                              /* input model   */\n")
outlines.append("homti_models                             /* output models */\n")
nwave = 0

for wt in ["P", "S"]:
    print(f"Preparing {wt} waves")
    bbz = []
    if wt == "P":
        dist_min = dist_min_P
        dist_max = dist_max_P
    else:
        dist_min = dist_min_S
        dist_max = dist_max_S

    bbaz = []
    for fl in sorted(filelist):
        station = fl.split("/")[-1].split("_")[0]
        wavetype = fl.split("/")[-1].split("_")[2]
        data_name = fl.split("/")[-1].split("_ff")[0]
        comp = fl.split("/")[-1].split("_")[1]

        st_lat, st_lon = get_station_coords(station, stations_file)
        dist = distance(st_lat, st_lon, ev_lat, ev_lon)[0]
        st = (st_lat, st_lon)
        ev = (ev_lat, ev_lon)
        baz = calculate_initial_compass_bearing(ev, st)
        use_it = use_it_or_not(station, comp, wavetype, data2use)

        if wavetype == wt and use_it and dist >= dist_min and dist <= dist_max:
            bbaz.append(baz)
            nwave += 1

    baz_bins = np.arange(0, 370, 10.0)
    counts, bins = np.histogram(bbaz, bins=baz_bins)

    for fl in sorted(filelist):
        station = fl.split("/")[-1].split("_")[0]
        wavetype = fl.split("/")[-1].split("_")[2]
        data_name = fl.split("/")[-1].split("_ff")[0]
        comp = fl.split("/")[-1].split("_")[1]

        st_lat, st_lon = get_station_coords(station, stations_file)
        dist = distance(st_lat, st_lon, ev_lat, ev_lon)[0]
        st = (st_lat, st_lon)
        ev = (ev_lat, ev_lon)
        baz = calculate_initial_compass_bearing(ev, st)
        use_it = use_it_or_not(station, comp, wavetype, data2use)

        if wavetype == wt and use_it and dist >= dist_min and dist <= dist_max:
            for i in range(len(bins) - 1):
                baz_min = bins[i]
                baz_max = bins[i + 1]

                if baz_min not in bbz:
                    bbz.append(baz_min)

                if baz >= baz_min and baz <= baz_max:
                    density_sta = 1 / float(counts[i])
                    print(counts[i], density_sta)
                    break
            outlines.append(f"{data_name}\n")
            weight = str(round(density_sta, 3))
            outlines.append(f"{weight}\n")

outlines.insert(3, f"{nwave}\t\t\t\t\t/* nwave */\n")
outlines.append("1                                       /* iwrite_models */")

with open("homti.in", "w") as outfile:
    outfile.writelines(outlines)
