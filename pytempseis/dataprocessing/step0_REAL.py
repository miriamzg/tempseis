import os
import sys
import glob
from obspy.core import read


Event_code = sys.argv[1]
database = sys.argv[2]

if not os.path.exists(database):
    raise FileNotFoundError(database)

Data_folder = f"{database}/{Event_code}/raw_data/"
if not os.path.exists(Data_folder):
    raise FileNotFoundError(Data_folder)

out_folder = f"{database}/{Event_code}/processed_data/"
if not os.path.exists(out_folder):
    print(f"Making output directory {out_folder}")
    os.mkdir(out_folder)

perl_dir = f"{os.path.dirname(os.path.abspath(__file__))}/../../lib/seis_process/bin"

sampling_rate = 0.5

# -------------------------------------------------------
sta_lines = open(f"{database}/STATIONS").readlines()
station_list = []
for i in range(0, len(sta_lines)):
    station_list.append(sta_lines[i].split()[0])

# =================================================
# Post processing real data
# =================================================
# extract from seed to sac
for i in range(0, len(station_list)):
    try:
        station_name = str(station_list[i])
        print(station_name)
        rdseed_file1 = glob.glob(f"{Data_folder}/{station_name}.*.mseed")[0]
        print(rdseed_file1)
        rdseed_file2 = glob.glob(f"{Data_folder}/*-{station_name}.*.dataless")[0]
        print(rdseed_file2)
        command = f"rdseed -d -o 1 -p -f {rdseed_file1} -g {rdseed_file2} -q {Data_folder}> /dev/null"
        os.system(command)
    except IndexError:
        pass

# apply instrument correction
CMT_folder = f"{Data_folder}/../"
CMT_file = CMT_folder + Event_code
command = f"perl {perl_dir}/process_data.pl -x corr -i -m {CMT_file} -s 2 -t 5/500 {Data_folder}*.SAC"
os.system(command)

# rotate 12 to RT
command = f"perl {perl_dir}/rotate.pl {Data_folder}*BH1*.corr"
os.system(command)

# rotate NE to RT
command = f"perl {perl_dir}/rotate.pl {Data_folder}*BHE*.corr"
os.system(command)

# Rename files
channel = "BH"
n = 1
filelist = []
for station in station_list:
    Zfile_list = glob.glob(f"{Data_folder}*{station}.00.{channel}*Z*.corr")

    if len(Zfile_list) != 0:
        Zfile = Zfile_list[0]
        Nfile = Zfile.replace(channel + "Z", channel + "1")
        Efile = Zfile.replace(channel + "Z", channel + "2")
        Rfile = Zfile.replace(channel + "Z", channel + "R")
        Tfile = Zfile.replace(channel + "Z", channel + "T")
        # check if it exists
        if not os.path.isfile(Nfile):
            Nfile = Zfile.replace(channel + "Z", channel + "N")
        if not os.path.isfile(Efile):
            Efile = Zfile.replace(channel + "Z", channel + "E")

        if (
            os.path.isfile(Nfile)
            and os.path.isfile(Efile)
            and os.path.isfile(Zfile)
            and os.path.isfile(Rfile)
            and os.path.isfile(Tfile)
        ):
            Ntr = read(Nfile)[0]
            Etr = read(Efile)[0]
            Ztr = read(Zfile)[0]
            Rtr = read(Rfile)[0]
            Ttr = read(Tfile)[0]

            Ntr.interpolate(sampling_rate=sampling_rate, method="linear")
            Etr.interpolate(sampling_rate=sampling_rate, method="linear")
            Ztr.interpolate(sampling_rate=sampling_rate, method="linear")
            Rtr.interpolate(sampling_rate=sampling_rate, method="linear")
            Ttr.interpolate(sampling_rate=sampling_rate, method="linear")

            Ztr.write(f"{out_folder}/{Ztr.id}", format="SAC")
            filelist.append(Ztr.id)
            Ntr.write(
                f"{out_folder}/{Ntr.id.replace(channel + '1', channel + 'N')}",
                format="SAC",
            )
            Etr.write(
                f"{out_folder}/{Etr.id.replace(channel + '2', channel + 'E')}",
                format="SAC",
            )
            Rtr.write(f"{out_folder}/{Rtr.id}", format="SAC")
            Ttr.write(f"{out_folder}/{Ttr.id}", format="SAC")
            n += 1
