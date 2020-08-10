import matplotlib.pylab as plt
from obspy.core import read
import os
import glob


synt_folder = "../database/CMTSOLUTION_201310251710A_SYNT_50/processed_data/"
noise_folder = "../ambient_noise/"
out_folder = "../database/CMTSOLUTION_201310251710A_SYNT_50/processed_data/noisy_data/"

sampling_rate = 0.5


os.system(f"mkdir {out_folder}")
os.system(f"mkdir {out_folder}/noise_plots")

station_list = []
fl = glob.glob(f"{synt_folder}*.int")
for f in fl:
    station = f.split("/")[-1].split(".")[1]
    net = f.split("/")[-1].split(".")[0]
    station_list.append([station, net])
# ====================================================


n = 0
for c in ["Z", "R", "T"]:
    for i in range(0, len(station_list)):
        station = station_list[i][0]
        net = station_list[i][1]

        synt_file = glob.glob(
            f"{synt_folder}{net}.{station}.MX{c}*.corr.int"
        )

        if c == "Z":
            c1 = "Z"
            c2 = "Z"

        if c == "R":
            c1 = "1"
            c2 = "N"
        if c == "T":
            c1 = "2"
            c2 = "E"

        noise_file = glob.glob(f"{noise_folder}*.{station}*BH{c1}*")
        if len(noise_file) == 0:
            noise_file = glob.glob(f"{noise_folder}*.{station}*BH{c2}*")
        if len(noise_file) == 0:
            noise_file = glob.glob(f"{noise_folder}*.AAK*BH{c1}*")
        if len(noise_file) == 0:
            noise_file = glob.glob(f"{noise_folder}*.AAK*BH{c2}*")
        if len(noise_file) > 1:
            noise_file = noise_file[0]

        print(synt_file)
        print(noise_file)
        print("**********")

        if len(noise_file) == 1:
            n += 1

            synt_tr = read(synt_file[0])[0]
            noise_tr = read(noise_file[0])[0]

            plt.subplot(311)
            plt.plot(synt_tr.times(), synt_tr.data, color="black", linewidth=2)
            plt.plot(noise_tr.times(), noise_tr.data, color="red", linewidth=2)
            plt.xlim(0, 7000)

            synt_tr.interpolate(sampling_rate=sampling_rate, method="linear")
            noise_tr.interpolate(sampling_rate=sampling_rate, method="linear")

            synt_noisy = synt_tr.copy()

            for i in range(0, len(synt_tr.data)):
                synt_noisy.data[i] = synt_tr.data[i] + noise_tr.data[i]

            synt_noisy.write(
                out_folder + synt_file[0].split("/")[-1] + ".noisy", format="SAC"
            )

            plt.subplot(312)
            plt.plot(synt_noisy.times(), synt_noisy.data, color="black")
            plt.xlim(0, 7000)
            plt.savefig(f"{out_folder}/noise_plots/{station}_c.png")
            plt.close()

print(f"n traces found: {n}")
