import matplotlib.pylab as plt


station = "CO_0"
observed_file = f"../../data/rfi_files/OBS/{station}_Z_fft"
predicted_real_file = f"./{station}_predicted_real.txt"
predicted_imag_file = f"./{station}_predicted_imag.txt"
point_source_file = f"point_source/{station}_Z_fft"


lines = open(observed_file).readlines()
obs_freq, obs_real, obs_imag = [], [], []
for i in range(2, len(lines)):
    obs_freq.append(float(lines[i].split()[0]))
    obs_real.append(float(lines[i].split()[1]))
    obs_imag.append(float(lines[i].split()[2]))

lines = open(point_source_file).readlines()
psource_freq, psource_real, psource_imag = [], [], []
for i in range(2, len(lines)):
    psource_freq.append(float(lines[i].split()[0]))
    psource_real.append(float(lines[i].split()[1]))
    psource_imag.append(float(lines[i].split()[2]))


lines = open(predicted_real_file).readlines()
pred_freq_real, pred_real = [], []
for i in range(1, len(lines)):
    pred_freq_real.append(float(lines[i].split()[0]))
    pred_real.append(float(lines[i].split()[1]))

lines = open(predicted_imag_file).readlines()
pred_freq_imag, pred_imag = [], []
for i in range(1, len(lines)):
    pred_freq_imag.append(float(lines[i].split()[0]))
    pred_imag.append(float(lines[i].split()[1]))


fig = plt.figure(1, figsize=(8.27, 11.69))
plt.subplot(511)
plt.plot(
    obs_freq[0: len(obs_freq) // 2],
    obs_real[0: len(obs_freq) // 2],
    color="red",
)
plt.plot(
    psource_freq[0: len(obs_freq) // 2],
    psource_real[0: len(obs_freq) // 2],
    color="black",
)
plt.xscale("log")
plt.xlim(0.01, 0.2)

plt.subplot(512)
plt.plot(
    obs_freq[0: len(obs_freq) // 2 - 1],
    obs_real[0: len(obs_freq) // 2 - 1],
    color="red",
)
plt.plot(
    pred_freq_real[0: len(obs_freq) // 2 - 1],
    pred_real[0: len(obs_freq) // 2 - 1],
    color="black",
)
plt.xscale("log")
plt.xlim(0.01, 0.2)

plt.savefig("tmp.png")
