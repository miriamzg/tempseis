import os
import sys
import matplotlib.pylab as plt
import glob
from obspy.core import read
sys.path.insert(0, '../lib/')
from functions import *


file1 = "../database/CMTSOLUTION_201310251710A_SYNT_50/processed_data/II.AAK.MXZ.sem.sac.corr.int"
file2 = "../database/CMTSOLUTION_201310251710A_SYNT_50/synthetics/point_source/II.AAK.MXZ.sem.sac"


#file2 = "/home/andrea/Dropbox/Work/NERC_Project_UCL/Real_data_inversion/Data/\
#CMTSOLUTION_201310251710A_SYNT_50.0/II.AAK.MXZ.sem.sac_post_processed_shift"

Tmin = 25.
Tmax = 100.

tr1 = read(file1)[0]
tr2 = read(file2)[0]

tr1 = filter_trace(tr1, float(Tmin), float(Tmax))
tr2 = filter_trace(tr2, float(Tmin), float(Tmax))

#tr1 = shift_stream(tr1, 3.75)

plt.plot(tr1.times(), tr1.data, color="black",	label="file1")
plt.plot(tr2.times(), tr2.data, color="red",	label="file2")
plt.xlim(560,570)
plt.xticks(np.arange(560,570, 1.))
plt.axvline(554, color="C1")
plt.axvline(549, color="C2")
plt.axvline(544, color="C3")
plt.axvline(539, color="C4")
plt.axvline(534, color="C5")

plt.ylim(-0.000035, 0.00005)
plt.legend()
plt.savefig("test.png")
plt.close()






