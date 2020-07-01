import matplotlib.pylab as plt
import os
import sys
import numpy as np



station = sys.argv[1].split("/")[-1].split("_Z_ff")[0]


synt_file = "rfi_files/SYNT/"+station+"_Z_ff"
obs_file = "rfi_files/SYNT/"+station+"_Z_obs"


freq_synt, real_synt, imag_synt, mft = [],[],[],[]
lines = open(synt_file).readlines()
for i in range(2,len(lines)):
	f_synt = float(lines[i].split()[0])
	r_synt = float(lines[i].split()[1])
	i_synt = float(lines[i].split()[2])
	if f_synt >= 0.0:
		freq_synt.append(f_synt)
		real_synt.append(r_synt)
		imag_synt.append(i_synt)

freq_obs, real_obs, imag_obs = [],[],[]
lines = open(obs_file).readlines()
for i in range(2,len(lines)):
	f_obs = float(lines[i].split()[0])
	r_obs = float(lines[i].split()[1])
	i_obs = float(lines[i].split()[2])
	if f_obs >= 0.0:
		freq_obs.append(f_obs)
		real_obs.append(r_obs)
		imag_obs.append(i_obs)

mft = 0.0
for i in range(0,len(freq_obs)):
	m = (real_synt[i] - real_obs[i])**2
	mft += m



plt.subplot(211)
plt.plot(freq_obs, real_obs, color="red", label="obs", linewidth=2, zorder=0)
plt.plot(freq_synt, real_synt, color="blue", label="synt", zorder=1)
plt.xscale("log")
plt.xlim(0.06,1.5)


plt.subplot(212)
plt.plot(freq_obs, imag_obs, color="red", label="obs", linewidth=2, zorder=0)
plt.plot(freq_synt, imag_synt, color="blue", label="synt", zorder=1)
plt.xscale("log")
plt.xlim(0.06,1.5)

plt.legend(loc=3)
plt.savefig(station+"_tmp.png")
plt.close()














