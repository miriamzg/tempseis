import matplotlib.pylab as plt



station = "COCO"
lines = open("../point_source/"+station+"_Z_fft").readlines()
ps_freq, ps_real, ps_imag = [],[],[]
for i in range(2,len(lines)):
	if float(lines[i].split()[0]) >= 0.0:
		ps_freq.append(float(lines[i].split()[0]))
		ps_real.append(float(lines[i].split()[1]))
		ps_imag.append(float(lines[i].split()[2]))


lines = open("../../../data/rfi_files/OBS/"+station+"_Z_fft").readlines()
obs_freq, obs_real, obs_imag = [],[],[]
for i in range(2,len(lines)):
	if float(lines[i].split()[0]) >= 0.0:
		obs_freq.append(float(lines[i].split()[0]))
		obs_real.append(float(lines[i].split()[1]))
		obs_imag.append(float(lines[i].split()[2]))

lines = open("./predicted").readlines()
pred_freq, pred_real, pred_imag = [],[],[]
for i in range(2,len(lines)):
	if float(lines[i].split()[0]) >= 0.0:
		pred_freq.append(float(lines[i].split()[0]))
		pred_real.append(float(lines[i].split()[1]))
		pred_imag.append(float(lines[i].split()[2]))


plt.subplot(211)
plt.plot(ps_freq, ps_real, color="black")
plt.plot(obs_freq, obs_real, color="red")
plt.plot(pred_freq, pred_real, color="blue")
plt.xscale("log")
plt.xlim(0.01, 0.2)

plt.subplot(212)
plt.plot(ps_freq, ps_imag, color="black")
plt.plot(obs_freq, obs_imag, color="red")
plt.plot(pred_freq, pred_imag, color="blue")
plt.xscale("log")
plt.xlim(0.01, 0.2)

plt.savefig("trace.png")
plt.close()


