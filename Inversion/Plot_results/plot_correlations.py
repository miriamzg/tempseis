import os
import matplotlib.pylab as plt
import numpy as np
import sys
from matplotlib.pylab import cm
from operator import itemgetter



result_folder = sys.argv[1]
#result_folder = "/Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/CMTSOLUTION_201505301123A_GCMT/aspect_subhor_0.5" 
result_file = result_folder + "/rfi_models"
#fault_file = "../results/" + code + "/201604151625A_TEST0_info.txt"

#Amax_t = 20.
#Amin_t = 5
#Duration_t = 7.5
#phi_t = 0.0
#v_abs_t = 4.0
#v_ang_t = 0.0

perc_threshold = 20.  # %

# lines = open(fault_file, "r").readlines()
# phi_t = float(lines[11].split()[1])
# Amax_t = float(lines[12].split()[1])
# Amin_t = float(lines[13].split()[1])
# Duration_t = float(lines[14].split()[1])
# dip_t =  float(lines[7].split()[1])
# strike_t =  float(lines[8].split()[1])
# E_t = Amin_t / Amax_t
	
mmft = []
lines = open(result_file).readlines()
short_lines = len(lines) - 2000
for i in range(3,short_lines):
	if len(lines[i].split()) != 0:
		if lines[i].split()[0] == "model:":
			misfit = float(lines[i].split("value:")[1])
			mmft.append(misfit)

min_mft = min(mmft)
mft_threshold = min_mft *(1+perc_threshold/100.)


mmisfit, AAmax, AAmin, pphi, sstrike, iindex, EE, ddip, dduration = [],[],[],[],[],[],[],[],[]
AAmax_a, AAmin_a, pphi_a, dduration_a = [],[],[],[]
vv_abs, vv_ang = [],[]
vv_abs_a, vv_ang_a = [],[]
rrc = []
ttc = []
mmft = []

lines = open(result_file).readlines()
mod_initial = int(lines[0].split()[0])
mod_per_iter = int(lines[1].split()[0])
n_iter = int(lines[2].split()[0])
n = 0
for i in range(3,short_lines):
	if len(lines[i].split()) != 0:
		if lines[i].split()[0] == "model:":
			n += 1
			iindex.append(n)
			misfit = float(lines[i].split("value:")[1])

			rc_x, rc_y, rc_z  = float(lines[i+1].split()[2]),float(lines[i+1].split()[3]),float(lines[i+1].split()[4])
			tc = float(lines[i+2].split()[2])
			duration = float(lines[i+3].split()[1]) * np.sqrt(3.)
			v_abs = float(lines[i+4].split()[2])
			v_ang = float(lines[i+4].split()[5])
			Amax = (float(lines[i+5].split()[1])) * np.sqrt(3.)
			Amin = float(lines[i+5].split()[3]) * np.sqrt(3.)
			phi =  float(lines[i+5].split()[5])

			strike =  float(lines[i+6].split()[1])
			dip = float(lines[i+6].split()[3])
			


			if misfit <= mft_threshold:
				mmisfit.append([misfit,n])
				mmft.append(misfit)
				#EE.append(E)
				ttc.append(tc)
				AAmax.append(Amax)
				AAmin.append(Amin)
				pphi.append(phi)
				sstrike.append(strike)
				ddip.append(dip)
				dduration.append(duration)
				vv_abs.append(v_abs)
				vv_ang.append(v_ang)


			AAmax_a.append(Amax)
			AAmin_a.append(Amin)
			pphi_a.append(phi)
			dduration_a.append(duration)
			vv_abs_a.append(v_abs)
			vv_ang_a.append(v_ang)

i_best = min(enumerate(mmft), key=itemgetter(1))[0] 
Amax_best =  AAmax[i_best]
Amin_best =  AAmin[i_best]
Phi_best =   pphi[i_best]
Duration_best =  dduration[i_best]
v_abs_best = vv_abs[i_best]
v_ang_best = vv_ang[i_best]


vmin = min(mmft)
vmax = min(mmft)*(1+perc_threshold/100.)
cmap = cm.jet

print vmin, vmax

font = {'size'   : 15}

plt.rc('font', **font)

fc_true = "black"
ec_true = "lime"
s_true = 500
lw_true = 3

fig = plt.figure(figsize=(18,15))
plt.subplots_adjust(left=0.1, right = 0.9, top=0.9, bottom=0.1, hspace=0.25, wspace=0.25)
plt.subplot(6,6,2)
cp = plt.scatter(AAmax, AAmin, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap,zorder=1)
plt.scatter(AAmax_a, AAmin_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.ylabel("A min", fontsize=20)
plt.xlim(min(AAmax_a), max(AAmax_a))
plt.ylim(min(AAmin_a),max(AAmax_a))
plt.scatter(Amax_best, Amin_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amax_t, Amin_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cbar = fig.colorbar(cp, cax=cbar_ax)
cbar.set_label(label='Misfit',size=30)
cbar.ax.tick_params(labelsize=20, size=4, width=4) 

plt.subplot(6,6,8)
cp = plt.scatter(AAmax, pphi, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap,zorder=1)
plt.scatter(AAmax_a, pphi_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.ylabel("Phi (deg)", fontsize=20)
plt.xlim(min(AAmax_a), max(AAmax_a))
plt.ylim(min(pphi_a), max(pphi_a))
#plt.xlim(0,30)
plt.scatter(Amax_best, Phi_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red",zorder=2)
#plt.scatter(Amax_t, phi_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)

plt.subplot(6,6,9)
cp = plt.scatter(AAmin, pphi, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmin_a, pphi_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(AAmin_a),max(AAmin_a))
plt.ylim(min(pphi_a),max(pphi_a))
plt.scatter(Amin_best, Phi_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amin_t, phi_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
#plt.xlim(0,30)

plt.subplot(6,6,14)
cp = plt.scatter(AAmax, dduration, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmax_a, dduration_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.ylabel("Duration (s)", fontsize=20)
plt.xlim(min(AAmax_a),max(AAmax_a))
plt.ylim(min(dduration_a),max(dduration_a))
plt.scatter(Amax_best, Duration_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amax_t, Duration_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
#plt.xlim(0,30)

plt.subplot(6,6,15)
cp = plt.scatter(AAmin, dduration, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmin_a, dduration_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(AAmin_a),max(AAmin_a))
plt.ylim(min(dduration_a),max(dduration_a))
plt.scatter(Amin_best, Duration_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amin_t, Duration_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
#plt.xlim(0,30)

plt.subplot(6,6,16)
cp = plt.scatter(pphi, dduration, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(pphi_a, dduration_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(pphi_a),max(pphi_a))
plt.ylim(min(dduration_a),max(dduration_a))
plt.scatter(Phi_best, Duration_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(phi_t, Duration_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)

plt.subplot(6,6,20)
plt.ylabel("|v| (km/s)", fontsize=20)
plt.scatter(AAmax, vv_abs, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmax_a, vv_abs_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(Amax_best, v_abs_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amax_t, v_abs_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(AAmax_a), max(AAmax_a))
plt.ylim(min(vv_abs_a), max(vv_abs_a))

plt.subplot(6,6,21)
cp = plt.scatter(AAmin, vv_abs, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmin_a, vv_abs_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(AAmin_a),max(AAmin_a))
plt.ylim(min(vv_abs_a),max(vv_abs_a))
plt.scatter(Amin_best, v_abs_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amin_t, v_abs_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)

plt.subplot(6,6,22)
cp = plt.scatter(pphi, vv_abs, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(pphi_a, vv_abs_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(pphi_a),max(pphi_a))
plt.ylim(min(vv_abs_a),max(vv_abs_a))
plt.scatter(Phi_best, v_abs_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(phi_t, v_abs_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)

plt.subplot(6,6,23)
cp = plt.scatter(dduration, vv_abs, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(dduration_a, vv_abs_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.xlim(min(dduration_a),max(dduration_a))
plt.ylim(min(vv_abs_a),max(vv_abs_a))
plt.scatter(Duration_best, v_abs_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Duration_t, v_abs_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)

plt.subplot(6,6,26)
plt.scatter(AAmax, vv_ang, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmax_a, vv_ang_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(Amax_best, v_ang_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amax_t, v_ang_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(AAmax_a), max(AAmax_a))
plt.ylim(min(vv_ang_a), max(vv_ang_a))
plt.xlabel("Amax (km)", fontsize=20)
plt.ylabel("v angle (deg)", fontsize=20)

plt.subplot(6,6,27)
plt.scatter(AAmin, vv_ang, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(AAmin_a, vv_ang_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(Amin_best, v_ang_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Amin_t, v_ang_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(AAmin_a), max(AAmin_a))
plt.ylim(min(vv_ang_a), max(vv_ang_a))
plt.xlabel("Amin (km)", fontsize=20)

plt.subplot(6,6,28)
plt.scatter(pphi, vv_ang, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(pphi_a, vv_ang_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(Phi_best, v_ang_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(phi_t, v_ang_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(pphi_a), max(pphi_a))
plt.ylim(min(vv_ang_a), max(vv_ang_a))
#plt.xlim(-20,20)
#plt.ylim(-20,20)
plt.xlabel("Phi (deg)", fontsize=20)

plt.subplot(6,6,29)
plt.scatter(dduration, vv_ang, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(dduration_a, vv_ang_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(Duration_best, v_ang_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(Duration_t, v_ang_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(dduration_a), max(dduration_a))
plt.ylim(min(vv_ang_a), max(vv_ang_a))
plt.xlabel("Duration (s)", fontsize=20)

plt.subplot(6,6,30)
plt.scatter(vv_abs, vv_ang, c=mmft,s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
plt.scatter(vv_abs_a, vv_ang_a, color="0.8",s=5, linewidth=0,zorder=0)
plt.scatter(v_abs_best, v_ang_best, linewidth=4, marker="o", s=200, facecolors="none", edgecolors="red")
#plt.scatter(v_abs_t, v_ang_t, marker="*", facecolors=fc_true, edgecolors=ec_true, linewidth=lw_true,  s=s_true)
plt.xlim(min(vv_abs_a), max(vv_abs_a))
plt.ylim(min(vv_ang_a), max(vv_ang_a))
plt.xlabel("|v| (km/s)", fontsize=20)





plt.subplot(6,6,1)
plt.ylabel("A min", fontsize=20)
plt.hist(AAmin_a, bins=30, orientation='horizontal', color="black")
plt.xscale("log")
plt.ylim(min(AAmin_a), max(AAmin_a))

plt.subplot(6,6,7)
plt.ylabel("Phi (deg)", fontsize=20)
plt.hist(pphi_a, bins=30, orientation='horizontal', color="black")
plt.xscale("log")
plt.ylim(min(pphi_a), max(pphi_a))

plt.subplot(6,6,13)
plt.ylabel("Duration (s)", fontsize=20)
plt.hist(dduration_a, bins=30, orientation='horizontal', color="black")
plt.xscale("log")
plt.ylim(min(dduration_a), max(dduration_a))

plt.subplot(6,6,19)
plt.ylabel("|v| (km/s)", fontsize=20)
plt.hist(vv_abs_a, bins=30, orientation='horizontal', color="black")
plt.xscale("log")
plt.ylim(min(vv_abs_a), max(vv_abs_a))

plt.subplot(6,6,25)
plt.ylabel("v angle (deg)", fontsize=20)
plt.hist(vv_ang_a, bins=30, orientation='horizontal', color="black")
plt.xscale("log")
plt.ylim(min(vv_ang_a), max(vv_ang_a))

plt.subplot(6,6,32)
plt.xlabel("Amax (km)", fontsize=20)
plt.hist(AAmax_a, bins=30, color="black")
plt.xlim(min(AAmax_a), max(AAmax_a))
plt.yscale("log")

plt.subplot(6,6,33)
plt.xlabel("Amin (km)", fontsize=20)
plt.hist(AAmin_a, bins=30, color="black")
plt.xlim(min(AAmin_a), max(AAmin_a))
plt.yscale("log")


plt.subplot(6,6,34)
plt.xlabel("Phi (deg)", fontsize=20)
plt.hist(pphi_a, bins=30, color="black")
plt.xlim(min(pphi_a), max(pphi_a))
plt.yscale("log")

plt.subplot(6,6,35)
plt.xlabel("Duration (s)", fontsize=20)
plt.hist(dduration_a, bins=30, color="black")
plt.xlim(min(dduration_a), max(dduration_a))
plt.yscale("log")

plt.subplot(6,6,36)
plt.xlabel("|v| (km/s)", fontsize=20)
plt.hist(vv_abs_a, bins=30, color="black")
plt.xlim(min(vv_abs_a), max(vv_abs_a))
plt.yscale("log")





plt.savefig(result_folder+"/Correlations.png")
plt.close()
