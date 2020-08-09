import matplotlib.pylab as plt
import numpy as np
import sys
from matplotlib import cm
from operator import itemgetter
import operator

result_folder = sys.argv[1]
result_file = f"{result_folder}/rfi_models"
scalar_moment_file = f"{result_folder}/scalar_moment.txt"

perc_treshold = 5.0

synt = False

mmisfit = []
AAmax = []
AAmin = []
pphi = []
sstrike = []
iindex = []
EE = []
ddip = []
dduration = []
AArea = []
AAmax_real = []
AAmin_real = []
vv_abs = []
vv_ang = []
rrc_mod = []
rrc_ang = []
dduration_real = []
rrc = []
ttc = []
mmft = []
nmodel = []
Sstress_drop = []
Sstress_drop_Mcguire = []
Sstress_drop_Amax = []
c = 1

# defining Mo from where it is saved in a text file
Mo = open(scalar_moment_file).readlines()
Mo = int(Mo[0])

lines = open(result_file).readlines()
mod_initial = int(lines[0].split()[0])
mod_per_iter = int(lines[1].split()[0])
n_iter = int(lines[2].split()[0])
n = 0
for i in range(3, len(lines)):
    if len(lines[i].split()) != 0:
        if lines[i].split()[0] == "model:":

            iindex.append(n)
            misfit = float(lines[i].split("value:")[1])

            rc_mod, rc_ang = (
                float(lines[i + 1].split()[2]),
                float(lines[i + 1].split()[3]),
            )
            tc = float(lines[i + 2].split()[2])
            duration = float(lines[i + 3].split()[1])
            v_abs = float(lines[i + 4].split()[2])
            v_ang = float(lines[i + 4].split()[5])
            Amax = float(lines[i + 5].split()[1])
            Amin = float(lines[i + 5].split()[3])
            phi = float(lines[i + 5].split()[5])

            strike = float(lines[i + 6].split()[1])
            dip = float(lines[i + 6].split()[3])

            mmisfit.append([misfit, n])
            mmft.append(misfit)
            # EE.append(E)
            rrc_mod.append(rc_mod)
            rrc_ang.append(rc_ang)
            ttc.append(tc)
            AAmax.append(Amax)
            AAmax_real.append(Amax * np.sqrt(3.0))
            AAmin.append(Amin)
            AAmin_real.append(Amin * np.sqrt(3.0))

            AArea.append(Amax * np.sqrt(3.0) * Amin * np.sqrt(3.0))
            Sstress_drop.append(
                ((c * Mo) / ((np.sqrt(Amax * Amin * 1000 * 1000 * 3)) ** (3))) / 1000000
            )
            pphi.append(phi)
            sstrike.append(strike)
            ddip.append(dip)
            dduration.append(duration)
            dduration_real.append(duration * np.sqrt(3.0))
            vv_abs.append(v_abs)
            vv_ang.append(v_ang)
            nmodel.append(n)
            n += 1

i_best = min(enumerate(mmft), key=itemgetter(1))[0]
min_mft = min(mmft)
mft_treshold = min_mft * (1 + perc_treshold / 100.0)

dduration_bests = []
vv_abs_bests = []
vv_angle_bests = []
AAmax_bests = []
AAmin_bests = []
pphi_bests = []
ttc_bests = []
stress_drop_bests = []
stress_drop_Mcguire_bests = []
stress_drop_Amax_bests = []
iindex_bests = []
j = 0
for i in range(3, (len(lines) - 1000)):
    if len(lines[i].split()) != 0:
        if lines[i].split()[0] == "model:":

            misfit = float(lines[i].split("value:")[1])
            if misfit <= mft_treshold:

                rc_mod, rc_ang = (
                    float(lines[i + 1].split()[2]),
                    float(lines[i + 1].split()[3]),
                )
                tc = float(lines[i + 2].split()[2])
                duration = float(lines[i + 3].split()[1]) * np.sqrt(3.0)
                v_abs = float(lines[i + 4].split()[2])
                v_ang = float(lines[i + 4].split()[5])
                Amax = (float(lines[i + 5].split()[1])) * np.sqrt(3.0)
                Amin = float(lines[i + 5].split()[3]) * np.sqrt(3.0)
                phi = float(lines[i + 5].split()[5])
                stress_drop = (
                    (c * Mo) / ((np.sqrt((Amax * Amin * 1000 * 1000))) ** (3)) / 1000000
                )
                stress_drop_Mcguire = ((c * Mo) / (np.pi * Amax * Amin * Amin)) / 1e15
                stress_drop_Amax = (
                    (c * Mo) / ((np.sqrt((Amax * Amax * 1000 * 1000))) ** (3)) / 1000000
                )
                ttc_bests.append([tc, misfit])
                dduration_bests.append([duration, misfit])
                vv_abs_bests.append([v_abs, misfit])
                vv_angle_bests.append([v_ang, misfit])
                AAmax_bests.append([Amax, misfit])
                AAmin_bests.append([Amin, misfit])
                pphi_bests.append([phi, misfit])
                stress_drop_bests.append([stress_drop, misfit])
                stress_drop_Mcguire_bests.append([stress_drop_Mcguire, misfit])
                stress_drop_Amax_bests.append([stress_drop_Amax, misfit])

                iindex_bests.append(j)
                j += 1


ttc_bests_sorted = sorted(
    ttc_bests, key=operator.itemgetter(1), reverse=False
)  # ordino i misfit dal peggiore al migliore
dduration_bests_sorted = sorted(
    dduration_bests, key=operator.itemgetter(1), reverse=False
)  # ordino i misfit dal peggiore al migliore
vv_abs_bests_sorted = sorted(
    vv_abs_bests, key=operator.itemgetter(1), reverse=False
)  # ordino i misfit dal peggiore al migliore
vv_angle_bests_sorted = sorted(
    vv_angle_bests, key=operator.itemgetter(1), reverse=False
)
AAmax_bests_sorted = sorted(AAmax_bests, key=operator.itemgetter(1), reverse=False)
AAmin_bests_sorted = sorted(AAmin_bests, key=operator.itemgetter(1), reverse=False)
pphi_bests_sorted = sorted(pphi_bests, key=operator.itemgetter(1), reverse=False)
stress_drop_bests_sorted = sorted(
    stress_drop_bests, key=operator.itemgetter(1), reverse=False
)
stress_drop_bests_Mcguire_sorted = sorted(
    stress_drop_Mcguire_bests, key=operator.itemgetter(1), reverse=False
)
stress_drop_bests_Amax_sorted = sorted(
    stress_drop_Amax_bests, key=operator.itemgetter(1), reverse=False
)

best_duration = round(dduration_bests_sorted[0][0], 2)
best_v_abs = round(vv_abs_bests_sorted[0][0], 2)
best_v_ang = round(vv_angle_bests_sorted[0][0], 2)
best_Amax = round(AAmax_bests_sorted[0][0], 2)
best_Amin = round(AAmin_bests_sorted[0][0], 2)
best_phi = round(pphi_bests_sorted[0][0], 2)
best_tc = round(ttc_bests_sorted[0][0], 2)
best_stress_drop = round(stress_drop_bests_sorted[0][0], 2)
best_stress_drop_Mcguire = round(stress_drop_bests_Mcguire_sorted[0][0], 2)
best_stress_drop_Amax = round(stress_drop_bests_Amax_sorted[0][0], 2)

min_duration = round(min(dduration_bests, key=operator.itemgetter(0))[0], 2)
max_duration = round(max(dduration_bests, key=operator.itemgetter(0))[0], 2)
min_v_abs = round(min(vv_abs_bests, key=operator.itemgetter(0))[0], 2)
max_v_abs = round(max(vv_abs_bests, key=operator.itemgetter(0))[0], 2)
min_v_ang = round(min(vv_angle_bests, key=operator.itemgetter(0))[0], 2)
max_v_ang = round(max(vv_angle_bests, key=operator.itemgetter(0))[0], 2)
min_Amax = round(min(AAmax_bests, key=operator.itemgetter(0))[0], 2)
max_Amax = round(max(AAmax_bests, key=operator.itemgetter(0))[0], 2)
min_Amin = round(min(AAmin_bests, key=operator.itemgetter(0))[0], 2)
max_Amin = round(max(AAmin_bests, key=operator.itemgetter(0))[0], 2)
min_phi = round(min(pphi_bests, key=operator.itemgetter(0))[0], 2)
max_phi = round(max(pphi_bests, key=operator.itemgetter(0))[0], 2)
min_tc = round(min(ttc_bests, key=operator.itemgetter(0))[0], 2)
max_tc = round(max(ttc_bests, key=operator.itemgetter(0))[0], 2)
min_stress_drop = round(min(stress_drop_bests, key=operator.itemgetter(0))[0], 2)
max_stress_drop = round(max(stress_drop_bests, key=operator.itemgetter(0))[0], 2)
min_stress_drop_Mcguire = round(
    min(stress_drop_Mcguire_bests, key=operator.itemgetter(0))[0], 2
)
max_stress_drop_Mcguire = round(
    max(stress_drop_Mcguire_bests, key=operator.itemgetter(0))[0], 2
)
min_stress_drop_Amax = round(
    min(stress_drop_Amax_bests, key=operator.itemgetter(0))[0], 2
)
max_stress_drop_Amax = round(
    max(stress_drop_Amax_bests, key=operator.itemgetter(0))[0], 2
)


print("***")
print(f"Duration:  {best_duration} (  {min_duration} -  {max_duration} )")
print(f"v abs:  {best_v_abs} (  {min_v_abs} -  {max_v_abs} )")
print(f"v ang:  {best_v_ang} (  {min_v_ang} -  {max_v_ang} )")
print(f"A max:  {best_Amax} (  {min_Amax} -  {max_Amax} )")
print(f"A min:  {best_Amin} (  {min_Amin} -  {max_Amin} )")
print(f"Phi:  {best_phi} (  {min_phi} -  {max_phi} )")
print(f"Stress Drop: {best_stress_drop} ({min_stress_drop} - {max_stress_drop} )")
print(
    f"Stress Drop Mcguire : {best_stress_drop_Mcguire} ( {min_stress_drop_Mcguire} - {max_stress_drop_Mcguire} )"
)
print(
    f"Stress Drop A_max : {best_stress_drop_Amax} ( {min_stress_drop_Amax} - {max_stress_drop_Amax} )"
)
print("***")


# Some statistics
rc_mean_module = round(np.median(rrc_mod), 1)
rc_mean_angle = round(np.median(rrc_ang), 1)

out = open(f"{result_folder}/best_parameters.txt", "w")

out.write(f"*** threshold: {str(perc_treshold)}% ***\n")
out.write(f"Time centroid: {str(best_tc)}     \n")
out.write(f"spatial centroid (delta):\n")
out.write(f"rc module: \t{str(rc_mean_module)} \n")
out.write(f"rc angle: \t{str(rc_mean_angle)} \n")
out.write(
    f"Duration: \t{str(best_duration)}\t( {str(min_duration)} - {str(max_duration)} )\n"
)
out.write(f"Amax: \t\t{str(best_Amax)}\t( {str(min_Amax)} - {str(max_Amax)} )\n")
out.write(f"Amin: \t\t{str(best_Amin)}\t( {str(min_Amin)} - {str(max_Amin)} )\n")
out.write(f"Phi: \t\t{str(best_phi)}\t( {str(min_phi)} - {str(max_phi)} )\n")
out.write(f"V abs: \t\t{str(best_v_abs)}\t( {str(min_v_abs)} - {str(max_v_abs)} )\n")
out.write(f"V ang: \t\t{str(best_v_ang)}\t( {str(min_v_ang)} - {str(max_v_ang)} )\n")
out.write(
    f"Stress Drop: \t\t{str(best_stress_drop)}\t( {str(min_stress_drop)} - {str(max_stress_drop)} )\n"
)
out.write(
    f"Stress Drop Mcguire: \t\t{str(best_stress_drop_Mcguire)}\t( {str(min_stress_drop_Mcguire)} - {str(max_stress_drop_Amax)} )\n"
)
out.write(
    f"Stress Drop A max: \t\t{str(best_stress_drop_Amax)}\t( {str(min_stress_drop_Amax)} - {str(max_stress_drop_Amax)} )\n"
)
out.write("***")
out.close()


# ==================================================================
plt.figure(1, figsize=(8.27, 11.69))
plt.subplot(311)
plt.scatter(iindex, rrc_mod, color="black", s=2, linewidth=0)
plt.xlim(min(iindex), max(iindex))
plt.ylabel("Rc module (km)")
plt.xlabel("# model")

plt.subplot(312)
plt.scatter(iindex, rrc_ang, color="black", s=2, linewidth=0)
plt.xlim(min(iindex), max(iindex))
plt.ylabel("rc angle (deg)")
plt.xlabel("# model")

plt.subplot(313)
plt.scatter(iindex, ttc, color="black", s=2, linewidth=0)
plt.xlim(min(iindex), max(iindex))
plt.ylabel("tc (s)")
plt.xlabel("# model")

plt.savefig(result_folder + "/centroid_evolution.png")
plt.close()

# ==================================================================
plt.figure(1, figsize=(8.27, 11.69))
plt.subplot(711)
plt.scatter(iindex, zip(*mmisfit)[0], color="black", s=2, linewidth=0)
plt.ylabel("Misfit")
plt.xlabel("# model")
plt.xlim(min(iindex), iindex[-1000])
plt.ylim(0.1 * min(zip(*mmisfit)[0]), 10 * max(zip(*mmisfit)[0]))

plt.subplot(712)
plt.scatter(iindex, AAmax_real, color="black", s=2, linewidth=0)
plt.scatter(iindex_bests, zip(*AAmax_bests)[0], color="red", s=2, linewidth=0, zorder=2)
plt.axhline(best_Amax, color="0.5", zorder=9)
plt.axhline(min_Amax, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_Amax, color="0.5", linestyle="--", zorder=9)
plt.xlim(min(iindex), iindex[-1000])
plt.ylabel("$A_{max}$ (km)")
plt.xlabel("# model")

plt.subplot(713)
plt.scatter(iindex, AAmin_real, color="black", s=2, linewidth=0)
plt.scatter(
    iindex_bests, zip(*AAmin_bests)[0], color="red", s=2, linewidth=0, zorder=10
)
plt.ylabel("$A_{min}$ (km)")
plt.xlabel("# model")
plt.xlim(min(iindex), iindex[-1000])
plt.axhline(best_Amin, color="0.5", zorder=9)
plt.axhline(min_Amin, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_Amin, color="0.5", linestyle="--", zorder=9)

plt.subplot(714)
plt.scatter(iindex, pphi, color="black", s=2, linewidth=0)
plt.scatter(iindex_bests, zip(*pphi_bests)[0], color="red", s=2, linewidth=0, zorder=10)
plt.ylabel("$\Phi$ (deg)")
plt.xlabel("# model")
plt.axhline(best_phi, color="0.5", zorder=9)
plt.axhline(min_phi, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_phi, color="0.5", linestyle="--", zorder=9)
plt.xlim(min(iindex), iindex[-1000])

plt.subplot(715)
plt.scatter(iindex, dduration_real, color="black", s=2, linewidth=0)
plt.scatter(
    iindex_bests, zip(*dduration_bests)[0], color="red", s=2, linewidth=0, zorder=10
)
plt.ylabel("$\Delta$ t (s)")
plt.xlabel("# model")
plt.xlim(min(iindex), iindex[-1000])
plt.axhline(best_duration, color="0.5", zorder=9)
plt.axhline(min_duration, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_duration, color="0.5", linestyle="--", zorder=9)

plt.subplot(716)
plt.scatter(iindex, vv_abs, color="black", s=2, linewidth=0)
plt.scatter(
    iindex_bests, zip(*vv_abs_bests)[0], color="red", s=2, linewidth=0, zorder=10
)
plt.ylabel("|v| (km/s)")
plt.xlabel("# model")
plt.xlim(min(iindex), iindex[-1000])
plt.axhline(best_v_abs, color="0.5", zorder=9)
plt.axhline(min_v_abs, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_v_abs, linestyle="--", zorder=9)

plt.subplot(717)
plt.scatter(iindex, vv_ang, color="black", s=2, linewidth=0)
plt.scatter(
    iindex_bests, zip(*vv_angle_bests)[0], color="red", s=2, linewidth=0, zorder=10
)
plt.ylabel("v angle (deg)")
plt.xlabel("# model")
plt.xlim(min(iindex), iindex[-1000])
plt.axhline(best_v_ang, color="0.5", zorder=9)
plt.axhline(min_v_ang, color="0.5", linestyle="--", zorder=9)
plt.axhline(max_v_ang, color="0.5", linestyle="--", zorder=9)

plt.savefig(result_folder + "/parameters_evolution.png")
plt.close()

# ==========================================================
plt.figure(1, figsize=(8.27, 11.69))
font = {"size": 15}
plt.rc("font", **font)

plt.subplot(611)
plt.subplots_adjust(
    left=None, bottom=0.08, right=None, top=0.98, wspace=0.5, hspace=0.5
)
plt.hist(AAmax_real, bins=100, color="C0")
plt.xlabel("$Amax (km$)", size=18)
plt.ylabel("# models")

plt.subplot(612)
plt.hist(AAmin_real, bins=100, color="C1")
plt.xlabel("Amin (km)", size=18)
plt.ylabel("# models")

plt.subplot(613)
plt.hist(pphi, bins=100, color="C2")
plt.ylabel("# models")
plt.xlabel("Phi angle (degrees)", size=18)

plt.subplot(614)
plt.hist(dduration_real, bins=100, color="C3")
plt.ylabel("# models")
plt.xlabel("Duration (seconds)", size=18)

plt.subplot(615)
plt.hist(vv_abs, bins=100, color="C4")
plt.ylabel("# models")
plt.xlabel("|v| (km/s)", size=18)

plt.subplot(616)
plt.hist(vv_ang, bins=100, color="C5")
plt.ylabel("# models")
plt.xlabel("v angle (degrees)", size=18)

plt.savefig(result_folder + "/hist1.png")
plt.close()

# =======================================================
plt.figure(1, figsize=(8.27, 11.69))
plt.subplot(611)
plt.subplots_adjust(
    left=None, bottom=0.08, right=None, top=0.98, wspace=0.5, hspace=0.5
)
plt.hist(AAmax_real, bins=100, color="C0")
plt.xlabel("Amax (km)")
plt.ylabel("# models")
plt.yscale("log")

plt.subplot(612)
plt.hist(AAmin_real, bins=100, color="C1")
plt.xlabel("Amin (km)")
plt.ylabel("# models")
plt.yscale("log")

plt.subplot(613)
plt.hist(pphi, bins=100, color="C2")
plt.ylabel("# models")
plt.xlabel("Phi angle (degrees)")
plt.yscale("log")

plt.subplot(614)
plt.hist(dduration_real, bins=100, color="C3")
plt.ylabel("# models")
plt.xlabel("Duration (seconds)")
plt.yscale("log")

plt.subplot(615)
plt.hist(vv_abs, bins=100, color="C4")
plt.ylabel("# models")
plt.xlabel("|v| (km/s)")
plt.yscale("log")

plt.subplot(616)
plt.hist(vv_ang, bins=100, color="C5")
plt.ylabel("# models")
plt.xlabel("v angle (degrees)")
plt.yscale("log")

plt.savefig(result_folder + "/hist2.png")
plt.close()

# =======================================
fig, axs = plt.subplots(4, 2, figsize=(15, 20))
ax1 = plt.subplot(421)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    AAmax_real,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$A_{max}$ (km)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "a)", fontweight="bold", transform=ax1.transAxes, size=20)

###################

ax2 = plt.subplot(422)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    AAmin_real,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$A_{min}$ (km)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "b)", fontweight="bold", transform=ax2.transAxes, size=20)


ax3 = plt.subplot(423)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    pphi,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$\Phi$ ($^\circ$)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "c)", fontweight="bold", transform=ax3.transAxes, size=20)

ax4 = plt.subplot(424)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    dduration_real,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$\Delta$ t (s)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "d)", fontweight="bold", transform=ax4.transAxes, size=20)


ax5 = plt.subplot(425)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    vv_abs,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$v_0$ (km s $^{-1}$)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "e)", fontweight="bold", transform=ax5.transAxes, size=20)

ax6 = plt.subplot(426)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    vv_ang,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$\Theta$ ($^\circ$)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "f)", fontweight="bold", transform=ax6.transAxes, size=20)

ax6 = plt.subplot(427)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3
)
it = plt.scatter(
    Sstress_drop,
    zip(*mmisfit)[0],
    c=nmodel,
    s=30,
    vmin=min(nmodel),
    vmax=max(nmodel),
    cmap=cm.rainbow,
    linewidth=0.5,
)
plt.xlabel("$\Delta \sigma$ (MPa)")
plt.ylabel("Misfit")
plt.yscale("log")
plt.ylim(min(zip(*mmisfit)[0]), max(zip(*mmisfit)[0]))
plt.text(-0.05, 1.05, "g)", fontweight="bold", transform=ax6.transAxes, size=20)

ax7 = plt.subplot(428)
ax7.set_visible(False)

cbaxes = fig.add_axes([0.2, 0.03, 0.6, 0.02])
cb = plt.colorbar(it, cax=cbaxes, label="iteration number", orientation="horizontal")
plt.savefig(f"{result_folder}/combined.png")
plt.close()

# ==================================================================
