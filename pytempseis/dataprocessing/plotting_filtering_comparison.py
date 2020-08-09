import os
import sys
import matplotlib.pylab as plt
import glob
import numpy as np
from scipy import signal
import matplotlib
from pylab import *

event_code = "CMTSOLUTION_201505301123A_GCMT"
x = [25, 25,20 ]
y = [ 150,150,70]

x1 = [25, 25,20, ]
y1 = [ 200,150,100]
ps_inv_list = []

obs_inv_list = []

ps_inv = []

obs_inv = []

station = []
station_list = []

component = []
component_list = []

wavetype2 = []
wavetype_list = []

for i in range(0,len(x)):
# Frequency band
    Tmin_p = x[i]
    Tmax_p = y[i]

    Tmin_s = x1[i]
    Tmax_s = y1[i]

    Tmin_r = 45  # 125
    Tmax_r= 100  #    180


    folder = "../database/" + event_code + "/" + "fortran_format_" +\
        str(Tmin_p) + "_"  + str(Tmax_p) + "_"+ str(Tmin_s) + "_"+ \
        str(Tmax_s) + "_"+ str(Tmin_r) + "_"+ str(Tmax_r)



    point_source_folder = folder + "/point_source/"
    observed_folder = folder + "/observed_data/"
    out_file = folder + "/station2use.txt"


    use_all = "no" # type "no", if you want to check all the traces one by one


    os.system("mkdir " + "../database/" + event_code  + "/filtering_testing")

    out = open(out_file,"w")
    out.write("Station\tComp\tType\tUse it?\n")
    filelist = glob.glob(point_source_folder+"*")
    n = 0
    

    for i in range(0,len(filelist)):
        print "----------------"
        ps_file = filelist[i]
        sta=(ps_file.split("/")[-1].split("_")[0])
        comp=(ps_file.split("/")[-1].split("_")[1])
        wavetype=(ps_file.split("/")[-1].split("_")[2])
        if (comp == "T" and wavetype == "P"):
            pass
        else:
            print ps_file
            obs_file = glob.glob(observed_folder + sta + "_" + comp + "_" + wavetype + "_ff")[0]
            print obs_file
            n+=1
            print sta, comp, wavetype, n , "/", len(filelist)
            
            lines = open(ps_file).readlines()
            ff, iimag, rreal = [],[],[]
            ps_complex = []
            for i in range(2, len(lines)):
                f = float(lines[i].split()[0])
                real  = float(lines[i].split()[1])
                imm  = float(lines[i].split()[2])
                c = complex(real, imm)
                ff.append(f)
                rreal.append(real)
                iimag.append(imm)
                ps_complex.append(c)


            fig = plt.figure(1, figsize=(11.69, 8.27))

            lines = open(obs_file).readlines()
            ff, iimag, rreal = [],[],[]
            obs_complex = []
            for i in range(2, len(lines)):
                f = float(lines[i].split()[0])
                real  = float(lines[i].split()[1])
                imm  = float(lines[i].split()[2])
                c = complex(real, imm)
                ff.append(f)
                rreal.append(real)
                iimag.append(imm)
                obs_complex.append(c)

            npoints = len(ff)
            ps_inv.append(np.fft.ifft(ps_complex, n=npoints))
            obs_inv.append( np.fft.ifft(obs_complex, n=npoints))
            station.append(sta)
            component.append(comp)
            wavetype2.append(wavetype)


    ps_inv_list.append(ps_inv)
    obs_inv_list.append(obs_inv)

#print obs_inv_list[3]



for i in range(0,len(filelist)):
        uplim = []
        lowlim = []
        if wavetype2[i] == "S":
            
            #print wavetype2[i]
            uplim.append(str(y1[0]))
            uplim.append(str(y1[1]))
            uplim.append(str(y1[2]))
            #uplim.append(str(y1[3]))
                
            lowlim.append(str(x1[0]))
            lowlim.append(str(x1[1]))
            lowlim.append(str(x1[2]))
            #lowlim.append(str(x1[3]))
        elif wavetype2[i] == "P":
            uplim.append(str(y[0]))
            uplim.append(str(y[1]))
            uplim.append(str(y[2]))
            #uplim.append(str(y[3]))
                
            lowlim.append(str(x[0]))
            lowlim.append(str(x[1]))
            lowlim.append(str(x[2]))
#lowlim.append(str(x[3]))

        else:
                
            uplim.append("100")
            lowlim.append("45")
            uplim.append("100")
            lowlim.append("45")
            uplim.append("100")
            lowlim.append("45")
            uplim.append("100")
            lowlim.append("45")



#print uplim
        #print lowlim
        fig = plt.figure(1, figsize=(11.69, 8.27))
        plt.suptitle("Station: "+station[i] + " Comp: " + component[i] + " Wavetype: "+wavetype2[i])
        plt.subplot(411)
        plt.plot(np.arange(0,len(obs_inv_list[0][i]),1), obs_inv_list[0][i].real, label="Observed", color="black",zorder=0, linewidth=3)
        plt.plot(np.arange(0,len(ps_inv_list[0][i]),1), ps_inv_list[0][i].real, label="Point source", color="red", zorder=1)
        plt.ylabel("Filter:"+uplim[0] +"s" + lowlim[0] +"s")
        plt.subplot(412)
        plt.plot(np.arange(0,len(obs_inv_list[1][i]),1), obs_inv_list[1][i].real, label="Observed", color="black", zorder=0, linewidth=3)
        plt.plot(np.arange(0,len(ps_inv_list[1][i]),1), ps_inv_list[1][i].real, label="Point source", color="red", zorder=1)
        plt.ylabel("Filter:"+uplim[1] +"s" + lowlim[1] +"s")
        plt.subplot(413)
        plt.plot(np.arange(0,len(obs_inv_list[2][i]),1), obs_inv_list[2][i].real, label="Observed", color="black", zorder=0, linewidth=3)
        plt.plot(np.arange(0,len(ps_inv_list[2][i]),1), ps_inv_list[2][i].real, label="Point source", color="red", zorder=1)
        plt.ylabel("Filter:"+uplim[2] +"s" + lowlim[2] +"s")
#plt.subplot(414)
#plt.plot(np.arange(0,len(obs_inv_list[3][i]),1), obs_inv_list[3][i].real, label="Observed", color="black", zorder=0, linewidth=3)
#plt.plot(np.arange(0,len(ps_inv_list[3][i]),1), ps_inv_list[3][i].real, label="Point source", color="red", zorder=1)
#plt.ylabel("Filter:"+uplim[3] +"s" + lowlim[3] +"s")
        #plt.xlabel()
        plt.savefig("../database/" + event_code + "/filtering_testing/" + station[i]+"_"+component[i]+"_"+wavetype2[i]+"2.png")
        plt.close()
    #plt.show()
    #   plt.close()

    
    
out.close()
    
# plt.show()




#for i in range(0,len(filelist)):
#  plt.subplot(411)
#plt.plot(np.arange(0,len(obs_inv_list[0][i]),1), obs_inv_list[0][i].real, label="Observed", color="black", #zorder=0, linewidth=3)
#  plt.plot(np.arange(0,len(ps_inv_list[0][i]),1), ps_inv_list[0][i].real, label="Point source", color="red", zorder=1)
#  plt.subplot(412)
#  plt.plot(np.arange(0,len(obs_inv_list[1][i]),1), obs_inv_list[1][i].real, label="Observed", color="black", zorder=0, linewidth=3)
#  plt.plot(np.arange(0,len(ps_inv_list[1][i]),1), ps_inv_list[1][i].real, label="Point source", color="red", zorder=1)
#  plt.subplot(413)
#  plt.plot(np.arange(0,len(obs_inv_list[2][i]),1), obs_inv_list[2][i].real, label="Observed", color="black", zorder=0, linewidth=3)
#  plt.plot(np.arange(0,len(ps_inv_list[2][i]),1), ps_inv_list[2][i].real, label="Point source", color="red", zorder=1)
#  plt.subplot(414)
#  plt.plot(np.arange(0,len(obs_inv_list[3][i]),1), obs_inv_list[3][i].real, label="Observed", color="black", zorder=0, linewidth=3)
#   plt.plot(np.arange(0,len(ps_inv_list[3][i]),1), ps_inv_list[3][i].real, label="Point source", color="red", zorder=1)
#   plt.show()
#   plt.close()

