import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import cm
import glob
import cartopy.crs as ccrs
from operator import itemgetter
from obspy.imaging.beachball import beach


def coord_from_dist_angle(lon1, lat1, d, brng):
    event_location = [lon1, lat1]

    lat1 = np.deg2rad(event_location[1])
    lon1 = np.deg2rad(event_location[0])

    R = 6371.0  # Radius of the Earth
    brng = np.deg2rad(brng)  # Bearing is 90 degrees converted to radians.

    lat2 = np.arcsin(
        np.sin(lat1) * np.cos(d / R) + np.cos(lat1) * np.sin(d / R) * np.cos(brng)
    )

    lon2 = lon1 + np.arctan2(
        np.sin(brng) * np.sin(d / R) * np.cos(lat1),
        np.cos(d / R) - np.sin(lat1) * np.sin(lat2),
    )

    lat2 = round(np.rad2deg(lat2), 3)
    lon2 = round(np.rad2deg(lon2), 3)

    return lon2, lat2


R = 6371.0

results_folder = sys.argv[1]
observed_folder = f"{results_folder}/observed/"

station_comp = []
filelist = glob.glob(f"{observed_folder}/*_observed.asc")
for file in filelist:
    station = file.split("/")[-1].split("_")[0]
    comp = file.split("/")[-1].split("_")[1]
    station_comp.append([station, comp])


ev_lon = float(sys.argv[3])
ev_lat = float(sys.argv[2])

station_file = "../../database/STATIONS"
station_list_avail = []
for station, comp in station_comp:
    with open(station_file) as file:
        for line in file:
            if station == line.split()[0]:
                lat = float(line.split()[2])
                lon = float(line.split()[3])
                station_list_avail.append([station, lon, lat])
                break

fig, ax = plt.subplots()

map = plt.axes(
    projection=ccrs.AzimuthalEquidistant(
        central_latitude=ev_lat, central_longitude=ev_lon
    )
)
map.coastlines()
map.set_global()
map.scatter(
    ev_lon,
    ev_lat,
    s=400,
    color="red",
    edgecolor="black",
    linewidth=1,
    marker="*",
    zorder=10,
)

for i in range(0, len(station_list_avail)):
    station = station_list_avail[i][0]
    st_lat = float(station_list_avail[i][2])
    st_lon = float(station_list_avail[i][1])
    map.scatter(
        st_lon,
        st_lat,
        marker="^",
        facecolor="lime",
        edgecolor="black",
        linewidth=1,
        s=200,
        zorder=10,
    )


plt.savefig(f"{results_folder}/Station_map.png")
plt.savefig(f"{results_folder}/Station_map.eps")
plt.close()


# ======================================================

fig = plt.figure(1, figsize=[15, 15])
ax = fig.add_subplot(111)
map = plt.axes(
    projection=ccrs.AzimuthalEquidistant(
        central_latitude=ev_lat, central_longitude=ev_lon
    )
)
map.coastlines()
map.stock_img()
map.set_global

llon2, llat2 = [], []
llon3, llat3 = [], []
for d in range(0, 360):
    dist_km = R * np.radians(40.0)
    lon2, lat2 = coord_from_dist_angle(ev_lon, ev_lat, dist_km, d)

    dist_km = R * np.radians(120.0)
    lon3, lat3 = coord_from_dist_angle(ev_lon, ev_lat, dist_km, d)

    llon2.append(lon2)
    llat2.append(lat2)
    llon3.append(lon3)
    llat3.append(lat3)
map.plot(llon2, llat2, zorder=20, color="red", linestyle="--", linewidth=4)
map.plot(llon3, llat3, zorder=20, color="red", linestyle="--", linewidth=4)


for i in range(0, len(station_list_avail)):
    station = station_list_avail[i][0]
    st_lat = float(station_list_avail[i][2])
    st_lon = float(station_list_avail[i][1])
    map.scatter(
        st_lon,
        st_lat,
        marker="^",
        facecolor="lime",
        edgecolor="black",
        linewidth=2,
        s=400,
        zorder=10,
    )

fm = [-6.29e26, 5.37e24, 6.24e26, -1.38e26, -6.37e25, 2.29e25]
beach1 = beach(fm, xy=(-120, 5), width=5000000, facecolor="C4")
beach1.set_zorder(100)
ax = plt.gca()
ax.add_collection(beach1)

map.plot([-111, 5], [ev_lon, ev_lat])
map.scatter(
    ev_lon,
    ev_lat,
    s=800,
    color="white",
    edgecolor="red",
    linewidth=3,
    marker="*",
    zorder=10,
)
plt.savefig(f"{results_folder}/Azimuthal_coverage.png")
plt.close()

# =======================================================================
# ++++++++++++++++++++++++++++++++++++++++++

results_file = f"{results_folder}/station_mft.asc"
lines = open(results_file).readlines()
mmft, iindex = [], []
for i in range(0, len(lines)):
    if lines[i].split()[0] == "Tot":
        mft = float(lines[i].split()[2])
        mmft.append(mft)
        iindex.append(i)

i_min = min(enumerate(mmft), key=itemgetter(1))[0]
line_best = iindex[i_min]


for c in ["Z"]:
    plt.figure(1, figsize=[10, 10])
    map = plt.axes(
        projection=ccrs.AzimuthalEquidistant(
            central_latitude=ev_lat, central_longitude=ev_lon
        )
    )
    map.coastlines()
    map.scatter(
        ev_lon,
        ev_lat,
        s=400,
        color="red",
        edgecolor="black",
        linewidth=1,
        marker="*",
        zorder=10,
    )

    ssmft, ssta, llon, llat = [], [], [], []
    for n in range(2, 1000):
        if lines[line_best - n].split()[0][0] != "*":
            sta = lines[line_best - n].split()[1].split("_")[0]
            comp = lines[line_best - n].split()[1].split("_")[1].split("ff")[0]
            s_mft = float(lines[line_best - n].split()[3])
            if comp == c:
                for j in range(0, len(station_list_avail)):
                    if sta == station_list_avail[j][0]:
                        st_lon = float(station_list_avail[j][1])
                        st_lat = float(station_list_avail[j][2])
                        llon.append(st_lon)
                        llat.append(st_lat)
                        ssmft.append(float(s_mft))
                        ssta.append(sta)
                        x, y = map(st_lon, st_lat)
                        plt.annotate(
                            sta,
                            xy=(x - 500000, y - 1000000),
                            color="black",
                            size=5,
                            fontweight="bold",
                        )

                        break

        else:
            break

    print("CICCIO", min(ssmft))
    s = map.scatter(
        llon,
        llat,
        marker="^",
        c=ssmft,
        linewidth=1,
        s=100,
        zorder=10,
        vmin=min(ssmft),
        vmax=min(ssmft) * 20,
        cmap=cm.turbo,
    )

    plt.colorbar(s, label="misfit")
    plt.suptitle(f"Station misfit {c} component")
    plt.savefig(f"{results_folder}/{c}_misfit_map.png", dpi=200)
    plt.close()
