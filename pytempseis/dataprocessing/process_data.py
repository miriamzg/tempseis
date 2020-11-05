from pytempseis.dataprocessing import (
    preprocessing,
    check_data,
    flip_seismograms,
    WaveArrivals,
    automatic_time_windowing,
    calculate_derivatives,
    cut_and_filter,
    final_check,
    compute_centroid_time,
    rotate_traces,
)
from pytempseis.functions import build_id_string
import yaml


def real_preprocessing(eventcode, database):
    response = input("Step 0: Extract SAC files from seed? (y/n)\t")
    while response not in ["y", "n"]:
        response = input("Please provide valid response (y/n)\t")
    if response == "y":
        preprocessing(eventcode, database)

    print("Step 1: First check")
    check_data(eventcode, database)

    print("Step2: Flip some seismograms")
    flip_seismograms(eventcode, database)


def synth_preprocessing(eventcode, database):
    print("Step 1: Rotate traces")
    rotate_traces(eventcode, database)

    print("Step 2: Shift to centroid time")
    compute_centroid_time(eventcode, database)


if __name__ == "__main__":
    with open("parameters.yml", "r") as file:
        params = yaml.full_load(file)
        eventcode = params['eventcode']
        database = params['database']
        real = params["real"]
        periods = params["periods"]
        p_phases = params["p_phases"]
        s_phases = params["s_phases"]

    id_string = build_id_string(periods)

    if real:
        real_preprocessing(eventcode, database)
    else:
        synth_preprocessing(eventcode, database)

    print("Step3: Automatic windowing")
    p_waves = WaveArrivals(
        periods["Tmin_p"],
        periods["Tmax_p"],
        p_phases,
        wavetype="p",
    )
    s_waves = WaveArrivals(
        periods["Tmin_s"],
        periods["Tmax_s"],
        s_phases,
        wavetype="s",
    )
    r_waves = WaveArrivals(periods["Tmin_r"], periods["Tmax_r"], wavetype="surface")
    waves = [p_waves, s_waves, r_waves]
    automatic_time_windowing(eventcode, database, waves, id_string, real)

    print("Step4: Calculate derivatives")
    calculate_derivatives(eventcode, database, periods["Tmin"], periods["Tmax"])

    print("Step5: Cut and filter traces")
    cut_and_filter(eventcode, database, periods, real)

    response = input("Step6: Do a final check? (y/n)\t")
    while response not in ["y", "n"]:
        response = input("Please provide valid response (y/n)\t")
    if response == "y":

        response = input("Step6: Do you want to check all traces manually? (y/n)\t")
        while response not in ["y", "n"]:
            response = input("Please provide valid response (y/n)\t")
        if response == "y":
            use_all = False
        else:
            use_all = True
        final_check(eventcode, database, id_string, use_all)
