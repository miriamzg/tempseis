from pytempseis.dataprocessing import (
    preprocessing,
    check_data,
    flip_seismograms,
    WaveArrivals,
    automatic_time_windowing,
    calculate_derivatives,
    cut_and_filter,
    final_check,
)

from argparse import ArgumentParser

parser = ArgumentParser(
    description="Prepares real and sythentic data and kernels for NA"
)
parser.add_argument("eventcode", type=str, help="GCMT code for event to be processed")
parser.add_argument(
    "database",
    type=str,
    help="Path to database.  This folder should contain a folder called <eventcode>",
)
parser.add_argument("--synth", action="store_true", help="Synthetic experiment")

args = parser.parse_args()

response = input("Step 0: Extract SAC files from seed? (y/n)")
while response not in ["y", "n"]:
    response = input("Please provide valid response (y/n)")
if response == "y":
    preprocessing(args.eventcode, args.database)

print("Step 1: First check")
check_data(args.eventcode, args.database)

print("Step2: Flip some seismograms")
flip_seismograms(args.eventcode, args.database)

print("Step3: Automatic windowing")
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
automatic_time_windowing(args.eventcode, args.database, waves, not args.synth)

print("Step4: Calculate derivatives")
calculate_derivatives(args.eventcode, args.database, 17, 300)

print("Step5: Cut and filter traces")
periods = {
    "Tmin": 17.0,  # first filtering
    "Tmax": 300.0,
    "Tmin_p": 25,  # p waves
    "Tmax_p": 60,
    "Tmin_s": 25,  # s waves
    "Tmax_s": 100,
    "Tmin_r": 45,  # surface waves
    "Tmax_r": 100,
}
cut_and_filter(args.eventcode, args.database, periods, not args.synth)

response = input("Step6: Do a final check? (y/n)")
while response not in ["y", "n"]:
    response = input("Please provide valid response (y/n)")
if response == "y":
    id_string = "_".join(
        [
            f"{periods[T]}"
            for T in ["Tmin_p", "Tmax_p", "Tmin_s", "Tmax_s", "Tmin_r", "Tmax_r"]
        ]
    )
    response = input("Step6: Do you want to check all traces manually? (y/n)")
    while response not in ["y", "n"]:
        response = input("Please provide valid response (y/n)")
    if response == "y":
        use_all = True
    else:
        use_all = False
    final_check(args.eventcode, args.database, id_string, use_all)
