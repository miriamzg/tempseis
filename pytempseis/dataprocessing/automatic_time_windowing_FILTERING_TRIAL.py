import os
import sys
import matplotlib.pylab as plt
from obspy.core import read
from obspy.taup import TauPyModel
import glob
import obspy
import numpy as np
from pytempseis.functions import filter_trace, distance, trim_trace_abs


class WaveArrivals:
    def __init__(self, Tmin, Tmax, pretime, posttime, phaselist=[]):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.pretime = pretime
        self.posttime = posttime
        self.phaselist = phaselist

    def body_arrivals(self, taup, ev_dep, dist_km):
        arrivals = taup_model.get_travel_times(
            source_depth_in_km=ev_dep,
            phase_list=self.phaselist,
            distance_in_degree=obspy.geodetics.kilometer2degrees(dist_km),
        )
        self.arrival = tr.stats["p_arrival"] = int(arrivals[0].time)
        self.extra_p_arrival = [int(arrivals[i].time) for i in range(1, len(arrivals))]
        self.extra_p_arrival_name = [
            arrivals[i].name for i in range(1, len(arrivals))
        ]


p_waves = WaveArrivals(
    20,
    70,
    100,
    100,
    [
        "P",
        "Pdiff",
        "pP",
        "PcP",
        "sP",
        "PKP",
        "PKS",
        "PKKP",
        "PKP",
        "PS",
    ],
)
s_waves = WaveArrivals(
    20, 100, 100, 100, ["S", "Sdiff", "pS", "SP", "sS", "PS", "SKS", "SP", "SKP"]
)
r_waves = WaveArrivals(45, 100, 200, 400)
waves = {"p": p_waves, "s": s_waves, "r": r_waves}

id_string = "_".join([f"{waves[w].Tmin}_{waves[w].Tmax}" for w in waves])

real = True

event_code = sys.argv[1]
database = sys.argv[2]

folder = f"{database}/{event_code}"
cmt_file = f"{folder}/{event_code}"

if real:
    channel = "BH"
    data_folder = f"{folder}/data_ready2use/"
else:
    channel = "MX"
    data_folder = f"{folder}/processed_data/"

plot_folder = f"{folder}/plots_picking_{id_string}/"
if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)

lines = open(cmt_file).readlines()
ev_lat = float(lines[4].split()[1])
ev_lon = float(lines[5].split()[1])
ev_dep = float(lines[6].split()[1])

taup_model = TauPyModel(model="iasp91")

sta_lines = open(f"{database}/STATIONS").readlines()

if real:
    file_list = sorted(glob.glob(f"{data_folder}*{channel}Z"))
    file_list += sorted(glob.glob(f"{data_folder}*{channel}R"))
    file_list += sorted(glob.glob(f"{data_folder}*{channel}T"))
else:
    file_list = sorted(glob.glob(f"{data_folder}*{channel}Z*.corr.int"))
    file_list += sorted(glob.glob(f"{data_folder}*{channel}R*.corr.int"))
    file_list += sorted(glob.glob(f"{data_folder}*{channel}T*.corr.int"))

out = open(
    f"{folder}/picking_times_{id_string}.txt",
    "w",
)
for w, wave in waves.items():
    out.write(f"Period bands {w}: {wave.Tmin}\t{wave.Tmax}\n")
out.write(f"Ev lat: {ev_lat}\n")
out.write(f"Ev lon: {ev_lon}\n")
out.write(f"Ev dep: {ev_dep}\n")
out.write(
    "Station\tChannel\tstarttime\tbegin\torigintime\tp start\tp end\ts start\ts end\tr start\tr end\n"
)
for file in file_list:
    tr = read(file)[0]

    for wave in waves.values():
        trace = tr.copy()
        wave.tr = filter_trace(trace, wave.Tmin, wave.Tmax)
        wave.env = obspy.signal.filter.envelope(wave.tr.data)

    tr.stats["snr"] = []
    tr.stats["cutpoints_in_s"] = []
    tr.stats["cutpoints_noise_in_s"] = []
    tr.stats["phase"] = []
    starttime = tr.stats.starttime
    begin = float(tr.stats.sac.get("b"))
    origintime = starttime - begin
    endtime = tr.stats.endtime
    st_lat = tr.stats["sac"]["stla"]
    st_lon = tr.stats["sac"]["stlo"]

    station = tr.stats.station
    for line in sta_lines:
        if station == line.split()[0]:
            print(line)
            st_lat = float(line.split()[2])
            st_lon = float(line.split()[3])

    delta = tr.stats.delta
    channel = tr.stats.channel
    dist_deg, dist_km = distance(st_lat, st_lon, ev_lat, ev_lon)
    print(f"Distance (deg): {dist_deg}")
    if dist_deg < 140.0 and dist_deg > 10.0:
        for w, wave in waves:
            if w == "r":
                continue

            wave.body_arrivals(taup_model, ev_dep, dist_km)
            wave.start = abs(begin) + wave.arrival - wave.pretime
            wave.end = abs(begin) + wave.arrival + wave.posttime
            wave.i_start = max(0, int(wave.start / delta))
            wave.i_end = int(wave.end / delta)

            # ---
            # Refine p and s picking using envelope inside the time window
            # ---

            wave.time = np.argmax(abs(wave.tr.data)[wave.i_start: wave.i_end]) * delta + wave.start
            wave.max = max(abs(wave.tr.data)[wave.i_start: wave.i_end])
            wave.start_final = wave.time - wave.pretime
            wave.end_final = wave.time + wave.posttime

        # ---
        # Surface picking using envelope
        # ---
        # rough determination of arrival time with min and max vs velocity
        r_start1 = abs(begin) + dist_km / 4.9
        r_start2 = abs(begin) + dist_km / 3.0

        ir_start1 = int(r_start1 / delta)
        ir_start2 = int(r_start2 / delta)

        time_rayleigh = r_start1 + np.argmax(env_surface[ir_start1:ir_start2]) * delta

        r_start_final = s_end_final

        for i in reversed(range(0, int(time_rayleigh / delta))):
            max_env = max(env_surface)
            if env_surface[i] <= 0.1 * max_env:
                r_start_final = i * delta
                break

        for i in range(int(time_rayleigh / delta), len(env_surface)):
            max_env = max(env_surface)
            if env_surface[i] <= 0.02 * max_env:
                r_end_final = i * delta
                break

        # if the starting time has not been found, then put the time window equals to 0
        if r_start_final == s_end_final:
            r_start_final = r_end_final

        # ---
        # convert relative times to absolute
        # ---
        p_start_final_abs = starttime + p_start_final
        p_end_final_abs = starttime + p_end_final
        s_start_final_abs = starttime + s_start_final
        s_end_final_abs = starttime + s_end_final
        r_start_final_abs = starttime + r_start_final
        r_end_final_abs = starttime + r_end_final

        # ---
        # Write the picking time into a file
        # ---
        line = f"{station}\t{channel}l\t{starttime}\t{begin}\t{origintime}\t{p_start_final_abs}\t{p_end_final_abs}\t{s_start_final_abs}\t{s_end_final_abs}\t{r_start_final_abs}\t{r_end_final_abs}\n"
        out.write(line)

        # ======================================================
        # cut traces
        extratime = 800.0
        tr_cut_p = trim_trace_abs(
            tr_p,
            tr_p.stats.starttime,
            p_start_final_abs,
            p_end_final_abs,
            Tmax_p,
            extratime,
        )
        tr_cut_s = trim_trace_abs(
            tr_s,
            tr_s.stats.starttime,
            s_start_final_abs,
            s_end_final_abs,
            Tmax_s,
            extratime,
        )
        tr_cut_r = trim_trace_abs(
            tr_r,
            tr_r.stats.starttime,
            r_start_final_abs,
            r_end_final_abs,
            Tmax_r,
            extratime,
        )

        plt.figure(1, figsize=(11.69, 8.27))
        plt.subplot(311)
        plt.plot(tr_p.times(), tr_p.data, color="black")
        plt.plot(
            tr_cut_p.times() + p_start_final - extratime,
            tr_cut_p.data,
            color="red",
            linewidth=2,
        )
        plt.axvline(origintime + p_arrival - 50, color="gray", linewidth=1)
        plt.axvline(origintime + p_arrival + 100, color="gray", linewidth=1)
        plt.ylim(-p_max * 1.1, p_max * 1.1)
        plt.xlim(p_arrival - 400, r_end_final + 2000)
        plt.xlim(0, 6000)
        plt.axvline(p_time, color="red", linestyle=":")
        plt.axvline(p_start_final, color="red", linewidth=2)
        plt.axvline(p_end_final, color="red", linewidth=2)
        plt.ylabel("P waves")
        for i in range(1, len(extra_p_arrival)):
            plt.axvline(extra_p_arrival[i], color="blue", linestyle=":", linewidth=2.5)
            plt.text(
                extra_p_arrival[i],
                p_max * 1.5,
                str(extra_p_arrival_name[i]),
                rotation=90,
            )
            print(extra_p_arrival_name[i])

        plt.subplot(312)
        plt.plot(tr_s.times(), tr_s.data, color="black", zorder=0)
        plt.plot(
            tr_cut_s.times() + s_start_final - extratime,
            tr_cut_s.data,
            color="red",
            linewidth=2,
        )
        plt.axvline(s_arrival - 50, color="gray", linewidth=1)
        plt.axvline(s_arrival + 150, color="gray", linewidth=1)
        plt.xlim(p_arrival - 400, r_end_final + 2000)
        plt.ylim(-s_max * 1.1, s_max * 1.1)
        plt.axvline(s_time, color="red", linestyle=":")
        plt.axvline(s_start_final, color="red", linewidth=2)
        plt.axvline(s_end_final, color="red", linewidth=2)
        plt.ylabel("S waves")
        for i in range(1, len(extra_s_arrival)):
            plt.axvline(extra_s_arrival[i], color="blue", linestyle=":", linewidth=2.5)
            plt.text(
                extra_s_arrival[i],
                s_max * 1.2,
                str(extra_s_arrival_name[i]),
                rotation=90,
            )
            print(extra_s_arrival_name[i])

        plt.subplot(313)
        plt.plot(tr_r.times(), tr_r.data, color="black")
        plt.plot(
            tr_cut_r.times() + r_start_final - extratime,
            tr_cut_r.data,
            color="red",
            linewidth=2,
        )
        plt.axvline(r_start_final, color="red", linewidth=2)
        plt.axvline(r_end_final, color="red", linewidth=2)
        plt.plot(tr_r.times(), env_surface, color="black", linestyle=":")
        plt.xlim(p_arrival - 400, r_end_final + 2000)
        plt.axvline(r_start1, color="black", linestyle="--")
        plt.axvline(r_start2, color="black", linestyle="--")
        plt.ylabel("surf. waves")

        plt.suptitle(f"Station: {station} Channel: {channel}")
        plt.savefig(f"{plot_folder}/{station}_{channel}_picking.png")
        plt.close()
out.close()
