import matplotlib.pylab as plt



filename = "testAmax"
lines = open(filename).readlines()
AAmax, mmft = [], []
for i in range(0,len(lines)):
	AAmax.append(float(lines[i].split()[0]))
	mmft.append(float(lines[i].split()[1]))

plt.plot(AAmax, mmft)
#plt.xlim(0,20)
#plt.ylim(0,1E-34)
plt.savefig("testAmax.png")
plt.close()
