import glob
import os
import sys
import numpy as np
from pytempseis.functions import shift_stream, filter_trace
from obspy.core import read
import matplotlib.pyplot as plt

event_code = sys.argv[1]
database = sys.argv[2]

raw_data_folder = f"{database}/{event_code}/raw_data/"
processed_data_folder = f"{database}/{event_code}/processed_data/"
cmt_finite_fault = f"{database}/{event_code}/{event_code}"
sampling_rate = 0.5  # Hz

if not os.path.exists(processed_data_folder):
    os.system(f"mkdir {processed_data_folder}")

# ------- compute time centroid and count point sources ----------
tt = []
cmt_lines = open(cmt_finite_fault).readlines()
n_points = 0
for i in range(0, len(cmt_lines)):
    if cmt_lines[i].split()[0] == "time":
        time_shift = float(cmt_lines[i].split()[2])
        tt.append(time_shift)
        n_points += 1


time_centroid = np.median(tt)
print(f"time centroid: {time_centroid}")
filelist = glob.glob(f"{raw_data_folder}*.sac")
for i in range(0, len(filelist)):
    raw_filename = filelist[i]
    filename = filelist[i].split("/")[-1]

    tr = read(raw_filename)[0]
    tr_new = tr.copy()

    for i in range(0, len(tr_new.data)):
        tr_new.data[i] = tr_new.data[i] / float(n_points)
    tr_new = shift_stream(tr_new, time_centroid)

    tr_new = filter_trace(tr_new, 17.0, 300.0)
    tr_new.interpolate(sampling_rate=sampling_rate, method="cubic")

    tr_new.write(f"{processed_data_folder}{filename}.corr.int", format="SAC")
    print(tr)


plt.hist(tt, bins=100)
plt.savefig("hist.png")
