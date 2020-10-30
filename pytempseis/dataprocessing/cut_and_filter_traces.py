import os
import sys
import shutil
import numpy as np
from pytempseis.functions import full_fft, filter_trace, trim_trace_abs
from obspy.core import read
import glob
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt


def process(eventcode, database, periods, real=True):
    id_string = "_".join(
        [
            f"{periods[T]}"
            for T in ["Tmin_p", "Tmax_p", "Tmin_s", "Tmax_s", "Tmin_r", "Tmax_r"]
        ]
    )

    (
        data_folder,
        kernels_folder,
        ps_folder,
        out_folder,
        arrivals_file,
    ) = setup_directories(eventcode, database, id_string, real)

    wavetype_list = ["P", "S", "W"]
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

    traces = Traces(data_folder, kernels_folder, ps_folder, arrivals_file, out_folder, periods)

    for station, comp in traces.cut_times:
        for wavetype in wavetype_list:
            if comp == "T" and wavetype == "P":
                pass
            else:
                for der in derivative_list:
                    print(f"Cutting and filtering traces, {station}, {comp}, {wavetype}")
                    traces.cut_filter_kernels(station, comp, der, wavetype)

                    print(f"Converting derivatives {station}, {comp}, {wavetype}")
                    traces.derivatives_to_ascii(station, comp, der, wavetype)

                print(f"Preparing real data , {station}, {comp}, {wavetype}")
                traces.observed_to_ascii(station, comp, wavetype, real)

                print(f"Preparing point source , {station}, {comp}, {wavetype}")
                traces.pointsource_to_ascii(station, comp, wavetype)

    plot_final(out_folder)


class Traces:
    def __init__(
        self, datafolder, kernelfolder, psfolder, arrivalsfile, outfolder, periods
    ):
        self.datafolder = datafolder
        self.kernelfolder = kernelfolder
        self.psfolder = psfolder
        self.arrivalsfile = arrivalsfile
        self.outfolder = outfolder
        self.periods = periods

        self._get_cut_times()
        self.sampling_rate = 0.5
        self.extratime = 800

    def cut_filter_kernels(self, station, comp, der, wavetype):
        # =====================================================================
        # Cut kernels and filtering
        # =====================================================================
        origintime = self.cut_times[station, comp][2]
        Tmin, Tmax, t1, t2 = self._select_period_and_cut(wavetype)

        tr_tmp = read(f"{self.kernelfolder}{station}_{comp}_{der}.sac")[0]
        tr_cut_f = tr_tmp.copy()
        tr_cut_f = filter_trace(tr_tmp, float(Tmin), float(Tmax))
        smoothing_time = 0.1  # Tmax
        tr_cut_f = trim_trace_abs(
            tr_cut_f, origintime, t1, t2, smoothing_time, self.extratime
        )

        tr_cut_f.write(
            f"{self.self.kernelfolder}cut/{station}_{comp}_{wavetype}_{der}.sac",
            format="SAC",
        )

    def derivatives_to_ascii(self, station, comp, der, wavetype):
        # =====================================================================
        # convert derivatives into ascii (in frequency domain)
        # =====================================================================
        filename = glob.glob(
            f"{self.kernelfolder}cut/{station}_{comp}_{wavetype}_{der}.sac"
        )[0]

        trace = read(filename)[0]
        trace.interpolate(sampling_rate=self.sampling_rate, method="linear")

        omega, sp = full_fft(trace)
        self._write_fft_to_file(
            omega,
            sp,
            f"{self.outfolder}/derivatives/{station}_{comp}_{wavetype}_{der}",
        )

    def observed_to_ascii(self, station, comp, wavetype, real=True):
        # =====================================================================
        # Convert observed seismograms into asci (in frequency domain)
        # =====================================================================
        if real:
            filelist = glob.glob(f"{self.datafolder}*.{station}*{channel}{comp}")
        else:
            filelist = glob.glob(f"{self.datafolder}*.{station}*{channel}{comp}*.int.noisy")

        Tmin, Tmax, t1, t2 = self._select_period_and_cut(wavetype)
        starttime = self.cut_times[station, comp][0]
        if len(filelist) == 1:
            filename = filelist[0]

            trace_tmp = read(filename)[0]

            trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
            trace_filt = trim_trace_abs(
                trace_filt_tmp, starttime, t1, t2, float(Tmax), self.extratime
            )
            trace_filt.write(
                f"{self.datafolder}cut/{station}_{comp}_{wavetype}.sac",
                format="SAC",
            )
            omega, sp = full_fft(trace_filt)
            self._write_fft_to_file(
                omega,
                sp,
                f"{self.outfolder}/observed_data/{station}_{comp}_{wavetype}_ff",
            )

    def pointsource_to_ascii(self, station, comp, wavetype):
        # =====================================================================
        # Convert point source seismogram to asci (in frequency domain)
        # =====================================================================
        Tmin, Tmax, t1, t2 = self._select_period_and_cut(wavetype)
        origintime = self.cut_times[station, comp][2]

        filename = f"{self.psfolder}*{station}.MX{comp}.sem.sac"
        trace_tmp = read(filename)[0]
        trace_filt_tmp = filter_trace(trace_tmp, float(Tmin), float(Tmax))
        trace_filt_tmp.interpolate(
            sampling_rate=self.sampling_rate, method="linear"
        )
        trace_filt = trim_trace_abs(
            trace_filt_tmp, origintime, t1, t2, float(Tmax), self.extratime
        )
        trace_filt.write(
            f"{self.psfolder}cut/{station}_{comp}_{wavetype}.sac",
            format="SAC",
        )

        omega, sp = full_fft(trace_filt)
        self._write_fft_to_file(
            omega,
            sp,
            f"{self.outfolder}/point_source/{station}_{comp}_{wavetype}_ps",
        )

    def _get_cut_times(self):
        # read cut times from picking time file
        self.cut_times = {}
        with open(self.arrivalsfile) as f:
            for _ in range(7):
                next(f)
            for line in f:
                station = line.split()[0]
                comp = line.split()[1].split(channel)[1]
                starttime = UTCDateTime(line.split()[2])
                begin = float(line.split()[3])
                origintime = UTCDateTime(line.split()[4])
                p_start = UTCDateTime(line.split()[5])
                p_end = UTCDateTime(line.split()[6])
                s_start = UTCDateTime(line.split()[7])
                s_end = UTCDateTime(line.split()[8])
                r_start = UTCDateTime(line.split()[9])
                r_end = UTCDateTime(line.split()[10])
                self.cut_times[station, comp] = [
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

    def _select_period_and_cut(self, wavetype):
        if wavetype == "P":
            Tmin = self.periods["Tmin_p"]
            Tmax = self.periods["Tmax_p"]
            t1, t2 = self.cut_times[station, comp][3], self.cut_times[station, comp][4]
        if wavetype == "S":
            Tmin = self.periods["Tmin_s"]
            Tmax = self.periods["Tmax_s"]
            t1, t2 = self.cut_times[station, comp][5], self.cut_times[station, comp][6]
        return Tmin, Tmax, t1, t2

    def _write_fft_to_file(self, omega, sp, file):
        with open(file) as out:
            out.write(f"{len(sp)}\n")
            out.write("omega\t\treal\t\t\timaginary\n")
            for i in range(0, len(sp)):
                f = omega[i]
                a = np.real(sp[i])
                b = np.imag(sp[i])
                out.write(f"{f:.8f}\t{a:.8e}\t\t{b:.8e}\n")


event_code = sys.argv[1]
database = sys.argv[2]

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


def setup_directories(event_code, database, id_string, real=True):
    def _clean_directory(folder):
        if os.path.exists(folder):
            os.mkdir(folder)
        else:
            shutil.rmtree(folder)
            os.mkdir(folder)

    _root = f"{database}/{event_code}"

    if real:
        data_folder = f"{_root}/data_ready2use/"
    else:
        data_folder = f"{_root}/processed_data/noisy/"
    out_folder = f"{_root}/fortran_format_{id_string}"
    kernels_folder = f"{_root}/kernels/"
    ps_folder = f"{_root}/synthetics/point_source/"

    folders_to_make = [
        out_folder,
        f"{data_folder}cut",
        f"{kernels_folder}cut",
        f"{ps_folder}cut",
        f"{out_folder}/derivatives/",
        f"{out_folder}/observed_data",
        f"{out_folder}/point_source",
    ]

    for folder in folders_to_make:
        _clean_directory(folder)

    arrivals_file = f"{_root}/picking_times_{id_string}.txt"
    if not os.path.exists(arrivals_file):
        raise FileNotFoundError(arrivals_file)

    return data_folder, kernels_folder, ps_folder, out_folder, arrivals_file


def read_fft_file(file):
    with open(file) as file:
        ff, iimag, rreal = [], [], []
        cmplx = []
        for _ in range(2):
            next(file)
        for line in file:
            f = float(line.split()[0])
            real = float(line.split()[1])
            imm = float(line.split()[2])
            c = complex(real, imm)
            ff.append(f)
            rreal.append(real)
            iimag.append(imm)
            cmplx.append(c)
    return ff, rreal, iimag, cmplx


def plot_final(out_folder):
    # =====================================================================
    # Plot and check final traces in the NA folder
    # =====================================================================
    print("Plotting and checking final traces...")
    obs_folder = f"{out_folder}/observed_data/"
    point_source_folder = f"{out_folder}/point_source/"
    plot_folder = f"{out_folder}/comparision_plots"
    os.mkdir(f"{plot_folder}")

    filelist = glob.glob(f"{obs_folder}*_ff")
    for fl in filelist:
        print(fl)
        plt.figure(1, figsize=(11.69, 8.27))
        filename = fl.split("/")[-1]

        ps_file = point_source_folder + filename.replace("ff", "ps")

        obs_ff, obs_rreal, obs_iimag, obs_complex = read_fft_file(fl)
        ps_ff, ps_rreal, ps_iimag, ps_complex = read_fft_file(ps_file)

        npoints = len(ps_ff)
        obs_inv = np.fft.ifft(obs_complex, n=npoints)
        ps_inv = np.fft.ifft(ps_complex, n=npoints)

        plt.subplot(311)
        plt.plot(obs_ff, obs_rreal, color="black")
        plt.plot(ps_ff, ps_rreal, color="red")
        plt.xlim(-0.5, 0.5)

        plt.subplot(312)
        plt.plot(obs_ff, obs_iimag, color="black")
        plt.plot(ps_ff, ps_iimag, color="red")
        plt.xlim(-0.5, 0.5)

        plt.subplot(313)
        plt.plot(
            np.arange(0, len(obs_inv), 1), obs_inv.real, label="Observed", color="black"
        )
        plt.plot(
            np.arange(0, len(ps_inv), 1), ps_inv.real, label="Point source", color="red"
        )
        plt.legend(loc=2)
        plt.savefig(f"{plot_folder}/{filename}.png")
        plt.close()
