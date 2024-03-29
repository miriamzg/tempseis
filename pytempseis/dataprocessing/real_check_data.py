import os
from obspy.core import read
import glob
import matplotlib.pylab as plt
from pytempseis.functions import filter_trace, PointPicker


def check_and_pick(tr_r, tr_s, xlim_max=6000.0, xlim_min=0.0):
    xlim_mid = (xlim_max + xlim_min) / 2.0

    fig = plt.figure(figsize=(11.69, 8.27))
    plt.subplot(211)
    plt.plot(
        tr_r.times(reftime=tr_r.stats.starttime),
        tr_r.data,
        color="black",
        linewidth=3,
        zorder=0,
    )
    plt.plot(
        tr_s.times(reftime=tr_r.stats.starttime),
        tr_s.data,
        color="red",
        zorder=10,
    )
    plt.axvline(xlim_mid, linestyle=":")
    plt.fill(
        [0, xlim_mid, xlim_mid, 0],
        [0, 0, 1000, 1000],
        "lime",
        alpha=0.2,
        edgecolor="r",
    )
    plt.fill(
        [0, xlim_mid, xlim_mid, 0],
        [0, 0, -1000, -1000],
        "red",
        alpha=0.2,
        edgecolor="r",
    )
    plt.fill(
        [xlim_mid, 80000, 80000, xlim_mid],
        [-1000, -1000, 1000, 1000],
        "yellow",
        alpha=0.2,
        edgecolor="r",
    )
    plt.ylim(1.5 * min(tr_r.data), 1.5 * max(tr_r.data))

    plt.annotate(
        "SAVE IT", xy=(0.1, 0.75), xycoords="axes fraction", fontsize=20
    )
    plt.annotate(
        "BIN IT", xy=(0.1, 0.25), xycoords="axes fraction", fontsize=20
    )
    plt.annotate(
        "FLIP IT", xy=(0.75, 0.75), xycoords="axes fraction", fontsize=20
    )

    plt.xlim(xlim_min, xlim_max)

    plt.subplot(212)
    plt.plot(
        tr_r.times(reftime=tr_r.stats.starttime), -tr_r.data, color="black"
    )
    plt.plot(tr_s.times(reftime=tr_r.stats.starttime), tr_s.data, color="red")
    plt.axvline(xlim_mid, linestyle=":")
    plt.xlim(xlim_min, xlim_max)

    pointpick = PointPicker(fig)

    plt.suptitle(tr_r.id)
    plt.show()

    return pointpick.xx[0], pointpick.yy[0]


def plot_traces(tr_r, tr_s, status="", plot_folder="./", xlim_max=6000):
    plt.figure(figsize=(11.69, 8.27))
    plt.subplot(211)
    plt.plot(tr_r.times(), tr_r.data, color="black", linewidth=3, zorder=0)
    plt.plot(tr_s.times(), tr_s.data, color="red", zorder=10)
    plt.ylim(1.5 * min(tr_r.data), 1.5 * max(tr_r.data))
    plt.xlim(0, xlim_max)

    plt.subplot(212)
    plt.plot(tr_r.times(), -tr_r.data, color="black")
    plt.plot(tr_s.times(), tr_s.data, color="red")
    plt.xlim(0, xlim_max)

    plt.suptitle(tr_r.id + "\n" + status)
    plt.savefig(os.path.join(plot_folder, f"{tr_r.id}.png"))
    plt.close()


def check_data(event_code, database, Tmin=35.0, Tmax=150.0):
    Data_folder = f"{database}/{event_code}/processed_data/"
    Synt_folder = f"{database}/{event_code}/synthetics/point_source"
    plot_folder = f"{Data_folder}/first_check_plots/"
    if not os.path.exists(plot_folder):
        os.mkdir(plot_folder)

    xlim_min = 0.0
    xlim_max = 6000.0
    xlim_mid = (xlim_max + xlim_min) / 2.0

    output_file = f"{database}/{event_code}/first_check.txt"
    out = open(output_file, "w")
    out.write("Sta\tComp\tQuality (Y: ok, N: bad, I: inverted)\n")
    # -------------------------------------------------------
    with open(f"{database}/STATIONS") as f:
        station_list = []
        for line in f:
            station_list.append(line.split()[0])

    components_list = ["Z", "R", "T"]

    for i, station in enumerate(station_list):
        for c in components_list:
            print(f"Station: {station} Component: {c}    {i + 1} / {len(station_list)}")
            if len(glob.glob(f"{Data_folder}/*{station}*BH{c}")) == 1:
                filename_real = glob.glob(f"{Data_folder}/*{station}*BH{c}")[0]
                filename_synt = glob.glob(f"{Synt_folder}/*{station}*MX{c}.sem.sac")[0]

                tr_r = read(filename_real)[0]
                tr_r = filter_trace(tr_r, Tmin, Tmax)

                tr_s = read(filename_synt)[0]
                tr_s = filter_trace(tr_s, Tmin, Tmax)

                x, y = check_and_pick(tr_r, tr_s, xlim_max, xlim_min)

                if y > 0 and x < xlim_mid:
                    status = "SAVED"
                if y < 0 and x < xlim_mid:
                    status = "BINNED"
                if x > xlim_mid:
                    status = "FLIPPED"

                print(status)
                line = f"{station}\t{c}\t{status}\n"
                out.write(line)
                plot_traces(tr_r, tr_s, status, plot_folder, xlim_max)

            else:
                print("file not found")


if __name__ == "__main__":
    import yaml

    with open("parameters.yml", "r") as file:
        params = yaml.full_load(file)
        eventcode = params["eventcode"]
        database = params["database"]
        Tmin = params["periods"]["Tmin"]
        Tmax = params["periods"]["Tmax"]

    check_data(eventcode, database, Tmin, Tmax)
