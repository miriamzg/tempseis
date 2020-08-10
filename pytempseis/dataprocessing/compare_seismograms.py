import matplotlib.pylab as plt
import numpy as np
from obspy.core import read
from pytempseis.functions import filter_trace


file1 = "../database/CMTSOLUTION_201310251710A_SYNT_50/processed_data/II.AAK.MXZ.sem.sac.corr.int"
file2 = "../database/CMTSOLUTION_201310251710A_SYNT_50/synthetics/point_source/II.AAK.MXZ.sem.sac"

Tmin = 25.0
Tmax = 100.0

tr1 = read(file1)[0]
tr2 = read(file2)[0]

tr1 = filter_trace(tr1, float(Tmin), float(Tmax))
tr2 = filter_trace(tr2, float(Tmin), float(Tmax))

plt.plot(tr1.times(), tr1.data, color="black", label="file1")
plt.plot(tr2.times(), tr2.data, color="red", label="file2")
plt.xlim(560, 570)
plt.xticks(np.arange(560, 570, 1.0))
plt.axvline(554, color="C1")
plt.axvline(549, color="C2")
plt.axvline(544, color="C3")
plt.axvline(539, color="C4")
plt.axvline(534, color="C5")

plt.ylim(-0.000035, 0.00005)
plt.legend()
plt.savefig("test.png")
plt.close()
