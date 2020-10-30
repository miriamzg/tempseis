import os
import sys
from obspy.core import read
import glob


def flip_seismograms(event_code, database):
    data_folder = f"{database}/{event_code}/processed_data/"
    output_folder = f"{database}/{event_code}/data_ready2use/"
    channel = "BH"

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    with open(f"{database}/{event_code}/first_check.txt") as file:
        next(file)
        station_list = []
        for line in file:
            station = line.split()[0]
            comp = line.split()[1]
            status = line.split()[2]
            station_list.append([station, comp, status])

    for line in station_list:
        station, comp, status = line
        print(station, comp, status)

        filename = glob.glob(f"{data_folder}*{station}*.00.{channel}{comp}")[0]

        tr = read(filename)[0]
        tri = tr.copy()

        if status == "I":
            for j in range(0, len(tr.data)):
                tri.data[j] = -tr.data[j]
            tri.write(output_folder + tr.id, format="SAC")

        if status == "Y":
            tr.write(output_folder + tr.id, format="SAC")

        if status == "X":
            pass


if __name__ == "__main__":
    event_code = sys.argv[1]
    database = sys.argv[2]
    flip_seismograms(event_code, database)
