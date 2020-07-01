#!/bin/bash

#  looping_newdat_paramiters.sh

file="CMTSOLUTION_122200D_GCMT"

echo " Duration , Error, Amax , Error, Amin, Error, Phi, Error, Vabs, Error, Vang, Error , Stress Drop, Error, " > combined.csv
echo "File," > side.csv
echo "File" > side2.csv

echo "Misfit" > misfit.csv

for f in $(find /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat -name 'best_parameters.txt')
do
SUBSTRING=$(echo $f| cut -d'/' -f 10)
echo $SUBSTRING >> side.csv
awk -F: 'NR >= 4  { print $2 }' $f | awk -F"[()]" '{print $1 ,",", $2,","}' | awk '{printf("%s ", $0)}' | awk -F"[,,]"  '{print ",", $0}' >> combined.csv


done

for f in $(find /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat -name 'best_model.asc')
do
SUBSTRING=$(echo $f| cut -d'/' -f 10)
echo $SUBSTRING  >> side2.csv
tail -n 2 $f | awk -F: '{if ($2 !="" ) print $2}' >> misfit.csv
done
paste side.csv combined.csv | sort  > /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/paramiters1.csv

paste -d, side2.csv misfit.csv | sort | awk '{print $2}'  > /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/paramiters2.csv

cd /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/

paste -d, paramiters1.csv paramiters2.csv  > /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/newdat.csv

awk -F, '{print $1,",",$2,",",$4,",",$6,",",$8,",",$10,",",$12,",",$14,",",$23 "\n","",",", "("$3")",",","("$5")",",","("$7")",",","("$9")",",","("$11")",",","("$13")",",","("$15")"}' /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/newdat.csv > /Users/TheStuffofAlice/Dropbox/TEMPSEIS_package_v1.2_WORKING/Inversion/Results/$file/newdat/newdat2.csv

#
#  Created by Alice on 02/08/2019.
#
