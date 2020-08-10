import os
import sys
import numpy as np
from pytempseis.functions import full_fft, filter_trace, trim_trace_abs
from obspy.core import read
import glob
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt

event_code = sys.argv[1]
real = True
channel = "BH" if real else "MX"

# frequency band for first filtering
Tmin = 17.0
Tmax = 300.0

# Frequency bands
# p waves
Tmin_p = 20
Tmax_p = 70

# s waves
Tmin_s = 20
Tmax_s = 100

# surface waves
Tmin_r = 45
Tmax_r = 100

# time before and after the time window
extratime = 800.0  # 800
sampling_rate = 0.5  # Hz

comp_list = ["Z", "R", "T"]
wavetype_list = ["P", "S", "W"]

if real:
    data_folder = f"../database/{event_code}/data_ready2use/"
else:
    data_folder = f"../database/{event_code}/processed_data/noisy/"

arrivals_file = (
    f"../database/{event_code}/picking_times_{Tmin_p}_{Tmax_p}_{Tmin_s}_{Tmax_s}_{Tmin_r}_{Tmax_r}.txt"
)
kernels_folder = f"../database/{event_code}/kernels/"
ps_folder = f"../database/{event_code}/synthetics/point_source/"
out_folder = (
    f"../database/{event_code}/fortran_format_{Tmin_p}_{Tmax_p}_{Tmin_s}_{Tmax_s}_{Tmin_r}_{Tmax_r}"
)

os.system(f"mkdir {out_folder}")

# read cut times from picking time file
station_comp = []
cut_file = arrivals_file
c_lines = open(cut_file).readlines()
cut_times = {}
for i in range(7, len(c_lines)):
    station = c_lines[i].split()[0]
    comp = c_lines[i].split()[1].split(channel)[1]
    starttime = UTCDateTime(c_lines[i].split()[2])
    begin = float(c_lines[i].split()[3])
    origintime = UTCDateTime(c_lines[i].split()[4])
    p_start = UTCDateTime(c_lines[i].split()[5])
    p_end = UTCDateTime(c_lines[i].split()[6])
    s_start = UTCDateTime(c_lines[i].split()[7])
    s_end = UTCDateTime(c_lines[i].split()[8])
    r_start = UTCDateTime(c_lines[i].split()[9])
    r_end = UTCDateTime(c_lines[i].split()[10])
    cut_times[station, comp] = [
        starttime,
        begin,
        origintime,
        p_start,
        p_end,
        s_start,
        s_end,
        r_start,
        r_end,
    ]
    station_comp.append([station, comp])


# =====================================================================
# 	Cut kernels and filtering
# =====================================================================
os.system(f"mkdir {kernels_folder}cut")
derivative_list = [
    "dSdx",
    "dSdy",
    "dSdz",
    "dS2dx2",
    "dS2dy2",
    "dS2dz2",
    "dSdxdy",
    "dSdxdz",
    "dSdydz",
]
for i in range(0, len(station_comp)):
    station = station_comp[i][0]
    comp = station_comp[i][1]
    origintime = cut_times[station, comp][2]
    for der in derivative_list:
        for wavetype in wavetype_list:
            print(f"Cutting and filtering traces, {station}, {comp}, {der}, {wavetype}")
            starttime = cut_times[station, comp][0]
            if comp == "T" and wavetype == "P":
                pass
            else:
                if wavetype == "P":
                    Tmin = Tmin_p
                    Tmax = Tmax_p
                    t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
                if wavetype == "S":
                    Tmin = Tmin_s
                    Tmax = Tmax_s
                    t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]

                tr_tmp = read(
                    kernels_folder + station + "_" + comp + "_" + der + ".sac"
                )[0]
                tr_cut_f = tr_tmp.copy()
                tr_cut_f = filter_trace(tr_tmp, float(Tmin), float(Tmax))
                smoothing_time = 0.1  # Tmax
                tr_cut_f = trim_trace_abs(
                    tr_cut_f, origintime, t1, t2, smoothing_time, extratime
                )

                tr_cut_f.write(
                    f"{kernels_folder}cut/{station}_{comp}_{wavetype}_{der}.sac",
                    format="SAC",
                )


# =====================================================================
# convert derivatives into ascii (in frequency domain)
# =====================================================================
os.system(f"mkdir {out_folder}/derivatives/")

derivative_list = [
    "dSdx",
    "dSdy",
    "dSdz",
    "dS2dx2",
    "dS2dy2",
    "dS2dz2",
    "dSdxdy",
    "dSdxdz",
    "dSdydz",
]
for i in range(0, len(station_comp)):
    station = station_comp[i][0]
    comp = station_comp[i][1]
    print(f"Converting derivatives {station}, {comp}")
    for der in derivative_list:
        for wavetype in wavetype_list:
            if comp == "T" and wavetype == "P":
                pass
            else:

                filename = glob.glob(
                    f"{kernels_folder}cut/{station}_{comp}_{wavetype}_{der}.sac"
                )[0]

                trace = read(filename)[0]
                trace.interpolate(sampling_rate=sampling_rate, method="linear")

                omega, sp = full_fft(trace)
                out = open(
                    f"{out_folder}/derivatives/{station}_{comp}_{wavetype}_{der}",
                    "w",
                )
                out.write(f"len(sp)\n")
                out.write(f"omega\t\treal\t\t\timaginary\n")
                for i in range(0, len(sp)):
                    f = omega[i]
                    a = np.real(sp[i])
                    b = np.imag(sp[i])
                    out.write(f"{f:.8f}\t{a:.8e}\t\t{b:.8e}\n")
                out.close()

# =====================================================================
# Convert observed seismograms into asci (in frequency domain)
# =====================================================================

os.system(f"mkdir {data_folder}cut")
os.system(f"rm {out_folder}/observed_data/*")
os.system(f"mkdir {out_folder}/observed_data")

for i in range(0, len(station_comp)):
    station = station_comp[i][0]
    comp = station_comp[i][1]
    origintime = cut_times[station, comp][2]
    if real:
        filelist = glob.glob(f"{data_folder}*.{station}*{channel}{comp}")
    else:
        filelist = glob.glob(f"{data_folder}*.{station}*{channel}{comp}*.int.noisy")
    starttime = cut_times[station, comp][0]
    for wavetype in wavetype_list:
        if comp == "T" and wavetype == "P":
            pass
        else:

            if wavetype == "P":
                Tmin = Tmin_p
                Tmax = Tmax_p
                t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
            if wavetype == "S":
                Tmin = Tmin_s
                Tmax = Tmax_s
                t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]

            print(f"Preparing real data , {station}, {comp}, {wavetype}, {Tmin}, {Tmax}")

            if len(filelist) == 1:
                filename = filelist[0]

                trace_tmp = read(filename)[0]
                start1 = trace_tmp.stats.starttime

                trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
                trace_filt = trim_trace_abs(
                    trace_filt_tmp, starttime, t1, t2, float(Tmax), extratime
                )
                trace_filt.write(
                    f"{data_folder}cut/{station}_{comp}_{wavetype}.sac",
                    format="SAC",
                )
                omega, sp = full_fft(trace_filt)
                out = open(
                    f"{out_folder}/observed_data/{station}_{comp}_{wavetype}_ff",
                    "w",
                )
                out.write(f"{len(sp)}\n")
                out.write("omega\t\treal\t\t\timaginary\n")
                for i in range(0, len(sp)):
                    f = omega[i]
                    a = np.real(sp[i])
                    b = np.imag(sp[i])
                    out.write(f"{f:.8f}\t{a:.8e}\t\t{b:.8e}\n")
                out.close()

# =====================================================================
# 	# Convert point source seismogram to asci (in frequency domain)
# =====================================================================
os.system(f"mkdir {ps_folder}cut")
os.system(f"rm {out_folder}/point_source/*")
os.system(f"mkdir {out_folder}/point_source")
for i in range(0, len(station_comp)):
    station = station_comp[i][0]
    comp = station_comp[i][1]
    starttime = cut_times[station, comp][0]
    begin = cut_times[station, comp][1]
    origintime = cut_times[station, comp][2]
    for wavetype in wavetype_list:
        if comp == "T" and wavetype == "P":
            pass
        else:

            if wavetype == "P":
                Tmin = Tmin_p
                Tmax = Tmax_p
                t1, t2 = cut_times[station, comp][3], cut_times[station, comp][4]
            if wavetype == "S":
                Tmin = Tmin_s
                Tmax = Tmax_s
                t1, t2 = cut_times[station, comp][5], cut_times[station, comp][6]

            print(f"Preparing point source , {station}, {comp}, {wavetype}, {Tmin}, {Tmax}")

            filename = f"{ps_folder}*{station}.MX{comp}.sem.sac"
            trace_tmp = read(filename)[0]
            trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
            trace_filt_tmp.interpolate(sampling_rate=sampling_rate, method="linear")
            trace_filt = trim_trace_abs(
                trace_filt_tmp, origintime, t1, t2, float(Tmax), extratime
            )
            trace_filt.write(
                f"{ps_folder}cut/{station}_{comp}_{wavetype}.sac",
                format="SAC",
            )

            omega, sp = full_fft(trace_filt)
            out = open(
                f"{out_folder}/point_source/{station}_{comp}_{wavetype}_ps",
                "w",
            )
            out.write(str(len(sp)) + "\n")
            out.write("omega\t\treal\t\t\timaginary\n")
            for i in range(0, len(sp)):
                f = omega[i]
                a = np.real(sp[i])
                b = np.imag(sp[i])
            out.write(f"{f:.8f}\t{a:.8e}\t\t{b:.8e}\n")
            out.close()

# =====================================================================
# Plot and check final traces in the NA folder
# =====================================================================
print(f"Plotting and checking final traces...")
obs_folder = f"{out_folder}/observed_data/"
point_source_folder = f"{out_folder}/point_source/"
kernels_folder = f"{out_folder}/kernels/"

plot_folder = f"{out_folder}/comparision_plots"
os.system(f"mkdir {plot_folder}")

filelist = glob.glob(f"{obs_folder}*_ff")
for fl in filelist:
    print(fl)
    plt.figure(1, figsize=(11.69, 8.27))
    filename = fl.split("/")[-1]

    ps_file = point_source_folder + filename.replace("ff", "ps")
    kernel_files = glob.glob(kernels_folder + filename.split("ff")[0] + "*")

    lines = open(fl).readlines()
    ff, iimag, rreal = [], [], []
    obs_complex = []
    for i in range(2, len(lines)):
        f = float(lines[i].split()[0])
        real = float(lines[i].split()[1])
        imm = float(lines[i].split()[2])
        c = complex(real, imm)
        ff.append(f)
        rreal.append(real)
        iimag.append(imm)
        obs_complex.append(c)

    plt.subplot(311)
    plt.plot(ff, rreal, color="black")

    plt.subplot(312)
    plt.plot(ff, iimag, color="black")

    lines = open(ps_file).readlines()
    ff, iimag, rreal = [], [], []
    ps_complex = []
    for i in range(2, len(lines)):
        f = float(lines[i].split()[0])
        real = float(lines[i].split()[1])
        imm = float(lines[i].split()[2])
        c = complex(real, imm)
        ff.append(f)
        rreal.append(real)
        iimag.append(imm)
        ps_complex.append(c)

    npoints = len(ff)
    obs_inv = np.fft.ifft(obs_complex, n=npoints)
    ps_inv = np.fft.ifft(ps_complex, n=npoints)

    plt.subplot(311)
    plt.plot(ff, rreal, color="red")
    plt.xlim(-0.5, 0.5)

    plt.subplot(312)
    plt.plot(ff, iimag, color="red")
    plt.xlim(-0.5, 0.5)

    plt.subplot(313)
    plt.plot(
        np.arange(0, len(obs_inv), 1), obs_inv.real, label="Observed", color="black"
    )
    plt.plot(
        np.arange(0, len(ps_inv), 1), ps_inv.real, label="Point source", color="red"
    )
    plt.legend(loc=2)
    plt.savefig(plot_folder + "/" + filename + ".png")
    plt.close()

# #=============================================================
