import numpy as np
from numpy import linalg as LA
import os
from obspy.core import read


def load_derivatives(station):
    import glob

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

    # load derivatives
    freq_der = {}
    dS = {}
    for der in derivative_list:
        filename = glob.glob(f"derivatives_mineos/{station}_Z_{der}.txt")[0]
        freq_der[der], dS[der] = np.genfromtxt(
            filename, delimiter="\t", unpack=True, dtype=np.complex128
        )

    return freq_der, dS


def filter(trace, Tmin, Tmax):
    trace.detrend("demean")
    trace.taper(0.05, type="cosine")
    trace.filter(
        "bandpass", freqmin=1 / Tmax, freqmax=1 / Tmin, corners=4, zerophase=True
    )

    return trace


def full_fft(tr):
    sp = np.fft.fft(tr.data)
    freq = np.fft.fftfreq(tr.times().shape[-1])
    return freq, sp


def calculate_W(Lmax, Lmin, phi, dip, strike):
    # Calculate W matrix from Lmax, Lmin and phi
    print(Lmax, Lmin, phi, dip, strike)
    phi = np.deg2rad(phi)
    dip = np.deg2rad(dip)
    strike = np.deg2rad(strike)

    M = [[Lmax, 0, 0], [0, Lmin, 0], [0, 0, 0]]
    # passo dagli autovettori nel sdr coordinato con l'ellisse a quello xy del piano di faglia
    # rotate I around z by phi angle
    I = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]  # versors in the reference system of the ellipse
    R_phi = [[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]]

    print(f"dip: {dip}")
    R_dip = [
        [1, 0, 0],
        [0, np.cos(-dip), -np.sin(-dip)],
        [0, np.sin(-dip), np.cos(-dip)],
    ]
    teta = (np.pi / 2.0) - strike
    print(f"teta: {teta}\t{strike}")
    R_strike = [
        [np.cos(teta), -np.sin(teta), 0],
        [np.sin(teta), np.cos(teta), 0],
        [0, 0, 1],
    ]

    u = np.dot(I, R_strike)
    u = np.dot(u, R_dip)
    u = np.dot(u, R_phi)

    S = u
    Sinv = LA.inv(S)

    W = np.dot(np.dot(S, M), Sinv)
    print(W)
    return W


source_code = "test0.1"
fault_file = "../../../prepare_observed_synt_data/" + source_code + "_data.txt"
fl = open(fault_file)
lines = fl.readlines()
station_list = []
for i in range(0, len(lines)):

    if lines[i].split()[0] == "Event_code:":
        event_code = lines[i].split()[1]
    if lines[i].split()[0] == "CMT_file:":
        cmtfile = lines[i].split()[1]
    if lines[i].split()[0] == "rc_lon:":
        rc_lon = float(lines[i].split()[1])
    if lines[i].split()[0] == "rc_lat:":
        rc_lat = float(lines[i].split()[1])
    if lines[i].split()[0] == "rc_depth:":
        rc_dep = float(lines[i].split()[1])
    if lines[i].split()[0] == "dip_angle:":
        dip = float(lines[i].split()[1])
    if lines[i].split()[0] == "strike_angle:":
        strike = float(lines[i].split()[1])
    if lines[i].split()[0] == "Max_axis:":
        Amax = float(lines[i].split()[1])
    if lines[i].split()[0] == "Min_axis:":
        Amin = float(lines[i].split()[1])
    if lines[i].split()[0] == "Stations:":
        for j in range(1, 100):
            if lines[i + j].split()[0][0] != "*":
                station_list.append(lines[i + j].split()[0])
            else:
                break
fl.close()

station_list = ["COCO"]

rc = [rc_lon, rc_lat, rc_dep]

Lmax = (Amax ** 2) / 4.0
Lmin = (Amin ** 2) / 4.0

Lmax = 100.0
Lmin = 10.0
phi = 10.0
strike = 262
dip = 26

dip_rad = np.deg2rad(dip)
strike_rad = np.deg2rad(strike)
phi_rad = np.deg2rad(phi)

model = "prem_noocean.txt"
c = "Z"

for station in station_list:
    freq_der, dS = load_derivatives(station)
    dS2 = [
        [dS["dS2dx2"], dS["dSdxdy"], dS["dSdxdz"]],
        [dS["dSdxdy"], dS["dS2dy2"], dS["dSdydz"]],
        [dS["dSdxdz"], dS["dSdydz"], dS["dS2dz2"]],
    ]
    S1_trace = read(f"point_source/NA.{station}..LHZ.sac.disp_S1")[0]
    S1_freq, S1_fft = full_fft(S1_trace)
    npoints = len(S1_trace.data)

    inv_S1_calc = np.fft.ifft(S1_fft, n=npoints)

    W = calculate_W(Lmax, Lmin, phi, dip, strike)
    E = 0.5 * np.tensordot(W, dS2) / 2.
    S2_predicted = S1_fft + E
    inv_S2_pred = np.fft.ifft(S2_predicted, n=npoints)

    predicted_trace = S1_trace.copy()
    for i in range(0, len(inv_S2_pred)):
        predicted_trace.data[i] = np.real(inv_S2_pred[i])
    predicted_trace.write(
        "../../data/rfi_files/NA_SRF/" + predicted_trace.id + ".sac.disp_S1",
        format="SAC",
    )

os.chdir("../../data/")
