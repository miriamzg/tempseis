#!/usr/bin/env python
# coding: utf-8

# In[47]:


import os
import sys
import re

def contains_word(s, w):
    return (' ' + w + ' ') in (' ' + s + ' ')


event_code = "CMTSOLUTION_201505301123A_GCMT" #change this when running the script
fortran_format = '/fortran_format_20_70_20_100_45_100/' #this also needs to be changed 
#results_folder = "/Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/"      + event_code + "/" + inversion_code
Data_folder = "/Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/database/" + event_code
from collections import OrderedDict
lines = open(Data_folder + fortran_format + '/station2use.txt').readlines()
station_list = []

for i in range(1,len(lines)-1):
    sta = lines[i].split()[0]
    station_list.append(sta)
station_list = list(set(station_list))
station_list = sorted(station_list)
#station_list = station_list.sort()
with open(Data_folder + fortran_format + 'station_remove_list.txt', 'w') as f:
    for item in station_list:
        f.write("%s\n" % item)

for i in range(0,len(station_list)):
    bad_words = station_list[i]
    with open(Data_folder + fortran_format + 'station2use.txt') as oldfile, open(Data_folder + fortran_format + '%s.text' % bad_words, "w") as newfile:
            for line in oldfile:
                 if contains_word(line.split()[0], bad_words) == False: 
                    newfile.write(line)


