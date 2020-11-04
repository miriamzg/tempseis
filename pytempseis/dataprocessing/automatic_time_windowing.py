import os
import sys
import matplotlib.pylab as plt
from obspy.core import read
from obspy.taup import TauPyModel
import glob
import obspy
import numpy as np
from copy import copy
from pytempseis.functions import filter_trace, distance, trim_trace_abs


class WaveArrivals:
    def __init__(self, Tmin, Tmax, pretime, posttime, phaselist=[], wavetype="body"):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.pretime = pretime
        self.posttime = posttime
        self.phaselist = phaselist
        if wavetype not in ["body", "p", "s", "surface"]:
            raise ValueError("wavetype must be either 'body', 'p', 's', or 'surface'")
        self.wavetype = wavetype

    def body_arrivals(self, taup, ev_dep, dist_km):
        arrivals = taup.get_travel_times(
            source_depth_in_km=ev_dep,
            phase_list=self.phaselist,
            distance_in_degree=obspy.geodetics.kilometer2degrees(dist_km),
        )
        self.arrival = int(arrivals[0].time)
        self.extra_arrival = [int(arrivals[i].time) for i in range(1, len(arrivals))]
        self.extra_arrival_name = [arrivals[i].name for i in range(1, len(arrivals))]

    def body_windows(self, begin, delta):
        self.start = abs(begin) + self.arrival - self.pretime
        self.end = abs(begin) + self.arrival + self.posttime
        self.i_start = max(0, int(self.start / delta))
        self.i_end = int(self.end / delta)

    def body_refine_windows(self, delta):
        """
        Refine p and s picking using envelope inside the time window
        """
        self.time = (
            np.argmax(abs(self.tr.data)[self.i_start: self.i_end]) * delta + self.start
        )
        self.max = max(abs(self.tr.data)[self.i_start: self.i_end])
        self.start_final = self.time - self.pretime
        self.end_final = self.time + self.posttime

    def rayleigh_windows(self, delta, swave_end, begin, dist_km):
        """
        Surface picking using envelope
        rough determination of arrival time with min and max vs velocity
        """
        self.r_start1 = abs(begin) + dist_km / 4.9
        self.r_start2 = abs(begin) + dist_km / 3.0

        ir_start1 = int(self.r_start1 / delta)
        ir_start2 = int(self.r_start2 / delta)

        self.time = self.r_start1 + np.argmax(self.env[ir_start1:ir_start2]) * delta
        self.start_final = copy(swave_end)

        for i in reversed(range(int(self.time / delta))):
            self.max = max(self.env)
            if self.env[i] <= 0.1 * self.max:
                self.start_final = i * delta
                break

        for i in range(int(self.time / delta), len(self.env)):
            self.max = max(self.env)
            if self.env[i] <= 0.02 * self.max:
                self.end_final = i * delta
                break

        # if the starting time has not been found, then put the time window equals to 0
        if self.start_final == swave_end:
            self.start_final = self.end_final

    def absolute_times(self):
        self.start_final_abs = self.tr.stats.starttime + self.start_final
        self.end_final_abs = self.tr.stats.starttime + self.end_final

    def cut_trace(self):
        self.tr_cut = trim_trace_abs(
            self.tr,
            self.tr.stats.starttime,
            self.start_final_abs,
            self.end_final_abs,
            self.Tmax,
            self.extratime,
        )

    def plot_traces(self, xlim=[0, 6000], ylabel=""):
        plt.plot(self.tr.times(), self.tr.data, color="black")
        plt.plot(
            self.tr_cut.times() + self.start_final - self.extratime,
            self.tr_cut.data,
            color="red",
            linewidth=2,
        )
        if self.wavetype != "surface":
            self.plot_arrivals()
        plt.ylim(-self.max * 1.1, self.max * 1.1)
        plt.xlim(xlim)
        plt.axvline(self.time, color="red", linestyle=":")
        plt.axvline(self.start_final, color="red", linewidth=2)
        plt.axvline(self.end_final, color="red", linewidth=2)
        plt.ylabel(ylabel)

    def plot_arrivals(self):
        plt.axvline(self.origintime + self.arrival - 50, color="gray", linewidth=1)
        plt.axvline(self.origintime + self.arrival + 100, color="gray", linewidth=1)
        if hasattr(self, "extra_arrival") and hasattr(self, "extra_arrival_name"):
            for i in range(len(self.extra_arrival)):
                plt.axvline(
                    self.extra_arrival[i], color="blue", linestyle=":", linewidth=2.5
                )
                plt.text(
                    self.extra_arrival[i],
                    self.max * 1.2,
                    str(self.extra_arrival_name[i]),
                    rotation=90,
                )


def build_figure(p_waves, s_waves, r_waves, station, channel, outfile):
    """
    Builds the final figure. *_waves are WaveArrival instances
    """
    plt.figure(1, figsize=(11.69, 8.27))
    plt.subplot(311)
    p_waves.plot_traces(xlim=[0, 6000], ylabel="P waves")

    plt.subplot(312)
    s_waves.plot_traces(
        xlim=[p_waves.arrival - 400, r_waves.end_final + 2000], ylabel="S waves"
    )

    plt.subplot(313)
    r_waves.plot_traces(
        xlim=[p_waves.arrival - 400, r_waves.end_final + 2000], ylabel="Surface waves"
    )

    plt.suptitle(f"Station: {station} Channel: {channel}")
    plt.savefig(outfile)
    plt.close()


def get_datafiles(folder, real=True):
    if real:
        channel = "BH"
        data_folder = f"{folder}/data_ready2use/"
        file_list = sorted(glob.glob(f"{data_folder}*{channel}Z"))
        file_list += sorted(glob.glob(f"{data_folder}*{channel}R"))
        file_list += sorted(glob.glob(f"{data_folder}*{channel}T"))
    else:
        channel = "MX"
        data_folder = f"{folder}/processed_data/"
        file_list = sorted(glob.glob(f"{data_folder}*{channel}Z*.corr.int"))
        file_list += sorted(glob.glob(f"{data_folder}*{channel}R*.corr.int"))
        file_list += sorted(glob.glob(f"{data_folder}*{channel}T*.corr.int"))
    return file_list


def automatic_time_windowing(event_code, database, waves, real=True):
    id_string = "_".join([f"{int(wave.Tmin)}_{int(wave.Tmax)}" for wave in waves])
    folder = f"{database}/{event_code}"
    cmt_file = f"{folder}/{event_code}"

    plot_folder = f"{folder}/plots_picking_{id_string}/"
    if not os.path.exists(plot_folder):
        os.mkdir(plot_folder)

    lines = open(cmt_file).readlines()
    ev_lat = float(lines[4].split()[1])
    ev_lon = float(lines[5].split()[1])
    ev_dep = float(lines[6].split()[1])

    taup_model = TauPyModel(model="iasp91")

    sta_lines = open(f"{database}/STATIONS").readlines()
    file_list = get_datafiles(folder, real)

    p_waves, s_waves, r_waves = waves
    with open(f"{folder}/picking_times_{id_string}.txt", "w") as out:
        for wave in waves:
            out.write(f"Period bands {wave.wavetype}: {wave.Tmin}\t{wave.Tmax}\n")
        out.write(f"Ev lat: {ev_lat}\n")
        out.write(f"Ev lon: {ev_lon}\n")
        out.write(f"Ev dep: {ev_dep}\n")
        out.write(
            "Station\tChannel\tstarttime\tbegin\torigintime\tp start\tp end\ts start\ts end\tr start\tr end\n"
        )
        for file in file_list:
            tr = read(file)[0]
            tr.stats["snr"] = []
            tr.stats["cutpoints_in_s"] = []
            tr.stats["cutpoints_noise_in_s"] = []
            tr.stats["phase"] = []
            starttime = tr.stats.starttime
            begin = float(tr.stats.sac.get("b"))
            origintime = starttime - begin
            endtime = tr.stats.endtime
            extratime = 8000.0

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
                for wave in waves:
                    trace = tr.copy()
                    wave.tr = filter_trace(trace, wave.Tmin, wave.Tmax)
                    wave.env = obspy.signal.filter.envelope(wave.tr.data)

                    wave.begin = begin
                    wave.origintime = origintime
                    wave.extratime = extratime
                    if wave.wavetype == "surface":
                        wave.rayleigh_windows(delta, s_waves.end_final, begin, dist_km)
                    else:
                        wave.body_arrivals(taup_model, ev_dep, dist_km)
                        wave.body_windows(begin, delta)
                        wave.body_refine_windows(delta)

                    wave.absolute_times()
                    wave.cut_trace()

                line = "\t".join(
                    [f"{wave.start_final_abs}\t{wave.end_final_abs}" for wave in waves]
                )
                out.write(
                    f"{station}\t{channel}\t{starttime}\t{begin}\t{origintime}\t{line}\n"
                )

                build_figure(
                    p_waves,
                    s_waves,
                    r_waves,
                    station,
                    channel,
                    f"{plot_folder}/{station}_{channel}_picking.png",
                )


if __name__ == "__main__":
    event_code = sys.argv[1]
    database = sys.argv[2]

    p_waves = WaveArrivals(
        25,
        60,
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
        wavetype="p",
    )
    s_waves = WaveArrivals(
        25,
        100,
        100,
        100,
        ["S", "Sdiff", "pS", "SP", "sS", "PS", "SKS", "SP", "SKP"],
        wavetype="s",
    )
    r_waves = WaveArrivals(45, 100, 200, 400, wavetype="surface")
    waves = [p_waves, s_waves, r_waves]
    real = True
    automatic_time_windowing(event_code, database, waves, real)
