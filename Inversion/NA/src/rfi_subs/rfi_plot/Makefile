#
#	Make discription file for RFI plot programs
#

#
# set top level directory
DEMO = ../../..

#
include $(DEMO)/src/makena.macros   # Edit this file to change location or
                                    # type of compilers (e.g. to gcc, g77)
GLIB = graphic/lib/glib.a
LDIR    = -L./graphic/xplot/lib
RFI = ../theo.o ../qlayer.o ../read_data.o \
      ../forward_modelling.o ../calcmisfit.o ../display_model.o \
      ../display_final.o ../lofw.o ../four1.o ../ran3.o ../write_srf.o

# make all programs 

all  : $(X)

yes  : glib synrf synrf_noise plot_model plot_rf-x plot_rf-p $(MACH)

no   : glib synrf synrf_noise plot_model plot_rf-p $(MACH)


.f.o  :
	$(F77) $(FFLAGS) -c $*.f

synrf : synrf.o 
	$(F77) -o $(DEMO)/bin/synrf synrf.o ../theo.o \
                  ../qlayer.o ../ran3.o ../four1.o ../lofw.o ../write_srf.o 

synrf_noise : synrf_noise.o 
	$(F77) -o $(DEMO)/bin/synrf_noise synrf_noise.o ../theo_noise.o \
		  ../qlayer.o ../ran3.o ../four1.o ../lofw.o ../write_srf.o
	cd $(DEMO)/data/rfi_files/rfi_plot; ./run_synns; \
	mv NA_SRF/rec_func.syn ORF/rec_func ;\
	cp ORF/rec_func ORF/rec_func_demo;

plot_model : plot_model.o
	$(F77) -o $(DEMO)/bin/plot_model plot_model.o \
		../indexx.o ../lofw.o $(GLIB)
#
plot_rf-x : plot_rf.o
	$(F77) $(LDIR) -O plot_rf.o  \
	../read_vmodelNA.o $(RFI) \
	-lxpakw $(XFLAGS) -o $(DEMO)/bin/plot_rf-x
#
plot_rf-p : plot_rf.o
	$(F77) $(LDIR) -O plot_rf.o  \
	../read_vmodelNA.o $(RFI) \
	-lpspak -o $(DEMO)/bin/plot_rf-p

glib : 	
	cd graphic; make -f Makefile
	cd graphic/xplot; make -f Makefile

bigendian:
	cd $(DEMO)/data/rfi_files/NA_SRF ;\
	\rm -f rec_func_demo ;\
	ln -s ../../demo_results/rec_func_500_20_2_sun rec_func_demo ;\
        cd ../NA_ORF ;\
	\rm -f rec_func_demo ;\
	ln -s rec_func_demo_sun rec_func_demo


litendian: 
	cd $(DEMO)/data/rfi_files/NA_SRF ;\
	\rm -f rec_func_demo ;\
	ln -s ../../demo_results/rec_func_500_20_2_dec rec_func_demo ;\
        cd ../NA_ORF ;\
	\rm -f rec_func_demo ;\
	ln -s rec_func_demo_dec rec_func_demo

sgi: bigendian

lin: litendian

sun: bigendian

dec: litendian

