import os
import sys


event_code = "CMTSOLUTION_201810131110A_GCMT"

inversion_code = "test8_31.18"
lat = " 52.71"
lon = " 153.43"



os.system("mkdir /Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Results/" + event_code)
results_folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Results/"     + event_code + "/" + inversion_code
Data_folder = "/Users/miriamgauntlett/TEMPSEIS_PACKAGE/database/" + event_code

os.system("mkdir " + results_folder)
os.system("mkdir " + results_folder + "/observed/")
os.system("mkdir " + results_folder + "/best_predicted/")

os.chdir ("/Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Preparation")
os.system("cp /Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Preparation/az_coverageTEST_P.pdf /Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Results/"     + event_code + "/" + inversion_code)
os.system("cp /Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Preparation/az_coverageTEST_S.pdf /Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Results/"     + event_code + "/" + inversion_code)
os.system("python scalar_moment.py " + results_folder)

#sys.exit()

os.system( "")
os.chdir("/Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/NA/src")
os.system("make clean")
os.system("make all")
os.chdir("../data")
os.system("../bin/rfi_na")

# #sys.exit()

os.system("cp rfi_files/NA_MDL/rfi_models " + results_folder)
os.system("cp rfi_files/NA_MDL/rfi_param " + results_folder)
os.system("cp rfi_files/rfi.in " + results_folder)
os.system("cp na.in  " + results_folder)
os.system("mv rfi_files/FINAL/*observed.asc " + results_folder + "/observed/")
os.system("mv rfi_files/FINAL/*best_predicted.asc " + results_folder + "/best_predicted/")
os.system("cp station_mft.asc " + results_folder)
os.system("cp best_model.asc " + results_folder)


os.chdir ("/Users/miriamgauntlett/TEMPSEIS_PACKAGE/Inversion/Plot_results")
os.system("python plot_results_subplot.py " + results_folder)
os.system("python plot_fault_results.py " + results_folder)            # uncomment for real data test
# os.system("python plot_fault_results_SYNT.py " + results_folder)  # uncomment for synthetic test
os.system("python plot_correlations.py " + results_folder)
os.system("python plot_station_map_body_surface.py " + results_folder + lat + lon)
#os.system("python plot_observed_vs_best_predicted_body_surface.py" + results_folder)
os.system("python plot_seismograms2.py "+ results_folder)

