import os
import sys
import glob
import re
from math import sqrt


Event_code = "CMTSOLUTION_201505301123A_GCMT"
Data_folder = "/Users/TheStuffofAlice/Documents/UCL/4th_year/Masters_Project/TEMPSEIS_package_v1.2_WORKING/database/" + Event_code


os.chdir(Data_folder)
filename ="CMTSOLUTION_201505301123A_GCMT"

# all of this is just pulling out the moments infomration from the CMTsolution file ( I THINK THIS IS A STUPID WAY TO DO THIS
f = open(filename,"r")
lines = f.readlines()
l=[]
for i in range(7, 13):
    s = lines[i]
    print s
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
Mrr = Mrr[0] * (10**Mrr[1])
Mtt = Mtt[0] * (10**Mtt[1])
Mpp = Mpp[0] * (10**Mpp[1])
Mrt = Mrt[0] * (10**Mrt[1])
Mrp = Mrp[0] * (10**Mrp[1])
Mtp = Mtp[0] * (10**Mtp[1])

#now I am going to attempt to work out the moment magnitude using the formula from the shearer textbook

Mo = (1/sqrt(2)) * (sqrt((Mrr**2)+(Mtt**2)+(Mpp**2)+(Mrt**2)+(Mrp**2)+(Mtp**2)))
print Mo
f = open('scalar_moment.txt', 'wb')
f.write("%d" % Mo)
f.close()

