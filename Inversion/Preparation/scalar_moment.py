import os
import matplotlib.pylab as plt
import numpy as np
import sys
from matplotlib import cm
from operator import itemgetter
import operator
import glob
import re
from math import sqrt



Event_code = "CMTSOLUTION_201810131110A_GCMT"
Data_folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/database/" + Event_code
#inversion_code = "unbounded"
#Results_folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Results/"     + Event_code + "/" + inversion_code
Result_folder = sys.argv[1]


os.chdir(Data_folder)
filename ="CMTSOLUTION_201810131110A_GCMT"

# all of this is just pulling out the moments infomration from the CMTsolution file ( I THINK THIS IS A STUPID WAY TO DO THIS
f = open(filename,"r")
lines = f.readlines()
l=[]
for i in range(7, 13):
    s = lines[i]
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    final_list = [float(x) for x in re.findall(match_number, s)]
    l.append(final_list)


Mrr = l[0]
Mtt = l[1]
Mpp = l[2]
Mrt = l[3]
Mrp = l[4]
Mtp = l[5]

#this is just putting the exponetials back together becuase it did not like the way the exponetials were written becuase they had a plus in front so I removed it and then put the bits back together
Mrr = Mrr[0] * (10**Mrr[1]) * 0.00000010 #converted to Nm
Mtt = Mtt[0] * (10**Mtt[1]) * 0.00000010
Mpp = Mpp[0] * (10**Mpp[1]) * 0.00000010
Mrt = Mrt[0] * (10**Mrt[1]) * 0.00000010
Mrp = Mrp[0] * (10**Mrp[1]) * 0.00000010
Mtp = Mtp[0] * (10**Mtp[1]) * 0.00000010

#now I am going to attempt to work out the moment magnitude using the formula from the shearer textbook

Mo = (1/sqrt(2)) * (sqrt((Mrr**2)+(Mtt**2)+(Mpp**2)+(Mrt**2)+(Mrp**2)+(Mtp**2)))

os.chdir(Result_folder)
f = open('scalar_moment.txt', 'wb')
f.write("%d" % Mo)
f.close()

#Amin = 3
#Amax = 6
#c =1
#stress_drop = (( c* Mo)/  ((sqrt((Amax*Amin*1000*1000))) **(3)))/ 1000000



#os.chdir(Results_folder)
#f = open("best_parameters.txt","r")
#lines = f.readlines()
#Amin_all =  lines[7].split()
#Amin = Amin_all[1]
#Amin = float(Amin)
#Amax_all =  lines[6].split()
#Amax = Amax_all[1]
#Amax = float(Amax)
#print Amin
#print Amax


# calculating the stress drop and writing it to a file
#c =1
#area =Amin * Amax
#root_area  = sqrt (area)
#stress_drop = (c* scalar_moment)/ root_area
#print stress_drop
#f = open('stress_drop.txt', 'wb')
#f.write("The stress drop is %d" % stress_drop)
#f.close()

