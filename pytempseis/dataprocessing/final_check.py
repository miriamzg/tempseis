import os
import sys
import matplotlib.pylab as plt
import glob
import numpy as np
from pylab import fill
from pytempseis.functions import PointPicker, read_fft_file


def final_check(event_code, database, id_string, use_all=False):
    folder = f"{database}/{event_code}/fortran_format_{id_string}"

    point_source_folder = f"{folder}/point_source/"
    observed_folder = f"{folder}/observed_data/"
    out_file = f"{folder}/station2use.txt"
    os.mkdir(f"{folder}/final_check")

    out = open(out_file, "w")
    out.write("Station\tComp\tType\tUse it?\n")
    filelist = glob.glob(f"{point_source_folder}*")

    for n, ps_file in enumerate(filelist, 1):
        print("----------------")
        sta = ps_file.split("/")[-1].split("_")[0]
        comp = ps_file.split("/")[-1].split("_")[1]
        wavetype = ps_file.split("/")[-1].split("_")[2]
        if comp == "T" and wavetype == "P":
            pass
        if wavetype == "W":
            pass
        else:
            print(ps_file)
            obs_file = glob.glob(f"{observed_folder}{sta}_{comp}_{wavetype}_ff")[0]
            print(obs_file)
            print(sta, comp, wavetype, n, "/", len(filelist))

            obs_ff, obs_rreal, obs_iimag, obs_complex = read_fft_file(obs_file)
            ps_ff, ps_rreal, ps_iimag, ps_complex = read_fft_file(ps_file)

            trim = len(ps_ff) // 2

            fig = plt.figure(1, figsize=(11.69, 8.27))
            plt.subplot(311)
            plt.plot(ps_ff[0:trim], ps_rreal[0:trim], color="black")
            plt.plot(obs_ff[0:trim], obs_rreal[0:trim], color="red")
            plt.xscale("log")

            plt.subplot(312)
            plt.plot(ps_ff[0:trim], ps_iimag[0:trim], color="black")
            plt.plot(obs_ff[0:trim], obs_iimag[0:trim], color="red")
            plt.xscale("log")

            npoints = len(obs_ff)
            ps_inv = np.fft.ifft(ps_complex, n=npoints)
            obs_inv = np.fft.ifft(obs_complex, n=npoints)

            corr = round(np.corrcoef(ps_inv.real, obs_inv.real)[0, 1], 2)

            rratio = []
            for i in range(0, len(obs_inv.real)):
                if abs(obs_inv.real[i]) >= max(abs(obs_inv.real)) / 50.0:
                    rratio.append((ps_inv.real[i]) / (obs_inv.real[i]))

            amp_ratio = round(np.mean(rratio), 2)

            plt.subplot(313)
            plt.plot(
                np.arange(0, len(obs_inv), 1),
                obs_inv.real,
                label="Observed",
                color="black",
                zorder=0,
                linewidth=3,
            )
            plt.plot(
                np.arange(0, len(ps_inv), 1),
                ps_inv.real,
                label="Point source",
                color="red",
                zorder=1,
            )

            if corr <= 0.67 or corr >= 1.29:
                corr_txt_colour = "red"
            else:
                corr_txt_colour = "black"
            plt.annotate(
                f"correlation coeff.: {corr}",
                xy=(20, 150),
                xycoords="axes pixels",
                fontsize=14,
                color=corr_txt_colour,
            )

            if amp_ratio <= 0.67 or amp_ratio >= 1.29:
                amp_txt_colour = "red"
            else:
                amp_txt_colour = "black"
            plt.annotate(
                f"amplitude coeff.: {amp_ratio}",
                xy=(20, 128),
                xycoords="axes pixels",
                fontsize=14,
                color=amp_txt_colour,
            )

            npts = len(ps_inv.real)
            h = npts // 2
            fill([-100, h, h, -100], [-1000, -1000, 1000, 1000], "red", alpha=0.15)
            fill([npts, h, h, npts], [-1000, -1000, 1000, 1000], "lime", alpha=0.15)

            plt.ylim(
                1.1 * min(min(obs_inv.real), min(ps_inv.real)),
                1.1 * max(max(obs_inv.real), max(ps_inv.real)),
            )
            plt.xlim(0, npts)
            plt.suptitle(f"Station: {sta} Comp: {comp} Wavetype: {wavetype}")
            plt.savefig(f"{folder}/final_check/{sta}_{comp}_{wavetype}.png")

            if not use_all:
                pointpick = PointPicker(fig)
                plt.show()
                x = pointpick.xx[0]
                y = pointpick.yy[0]

                if x >= h:
                    out.write(f"{sta}\t{comp}\t{wavetype}\tY\n")
                else:
                    out.write(f"{sta}\t{comp}\t{wavetype}\tN\n")

            plt.close()

    out.close()


if __name__ == "__main__":
    event_code = sys.argv[1]
    database = sys.argv[2]

    # Frequency band
    Tmin_p = 25
    Tmax_p = 60

    Tmin_s = 25
    Tmax_s = 100

    Tmin_r = 45
    Tmax_r = 100

    id_string = f"{Tmin_p}_{Tmax_p}_{Tmin_s}_{Tmax_s}_{Tmin_r}_{Tmax_r}"
    use_all = False

    final_check(event_code, database, id_string, use_all)
