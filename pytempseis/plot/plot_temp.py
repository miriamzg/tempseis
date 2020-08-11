import matplotlib.pylab as plt
import sys


station2plot = sys.argv[1]
comp = sys.argv[2]
wavetype = sys.argv[3]

obs = f"rfi_files/TEMP/{station2plot}_{comp}_{wavetype}_observed.asc"
pred = f"rfi_files/TEMP/{station2plot}_{comp}_{wavetype}_predicted.asc"

for n in range(0, 1000):
    obsf, obsr, obsi = [], [], []
    predf, predr, predi = [], [], []

    lines = open(obs).readlines()
    for i in range(0, len(lines)):
        if float(lines[i].split()[0]) >= 0:
            obsf.append(float(lines[i].split()[0]))
            obsr.append(float(lines[i].split()[1]))
            obsi.append(float(lines[i].split()[2]))

    lines = open(pred).readlines()
    for i in range(0, len(lines)):
        if float(lines[i].split()[0]) >= 0:
            predf.append(float(lines[i].split()[0]))
            predr.append(float(lines[i].split()[1]))
            predi.append(float(lines[i].split()[2]))

    xmin, xmax = 0.1, 0.3

    plt.figure(1, figsize=(8.27, 11.69))
    plt.subplots_adjust(bottom=None, top=0.92, hspace=0.3, wspace=0.2)
    plt.subplot(211)
    plt.plot(obsf, obsr, color="black", zorder=0, linewidth=2)
    plt.plot(predf, predr, color="red")
    plt.ylim(min(obsr) * 1.1, max(obsr) * 1.1)
    plt.xlim(xmin, xmax)
    plt.title(f"{comp} component - real part", y=1.05)

    plt.subplot(212)
    plt.plot(obsf, obsi, color="black", zorder=0, linewidth=2)
    plt.plot(obsf, predi, color="red")
    plt.xlim(xmin, xmax)
    plt.ylim(min(obsi) * 1.1, max(obsi) * 1.1)
    plt.title(f"{comp} component - imag part", y=1.05)

    plt.suptitle(f"Station: {station2plot}")
    plt.savefig("tmp.png")
    plt.close()
