#
#	Make discription file for Receiver function example subroutines
#

#FFLAGS = -O4
#F77 = /opt/SUNWspro/bin/f77

DEMO = ../..

include $(DEMO)/src/makena.macros   # Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

DIR = .
RFI = $(DIR)/theo.o $(DIR)/qlayer.o $(DIR)/read_data.o $(DIR)/read_seed.o \
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o $(DIR)/display_model.o \
      $(DIR)/display_final.o $(DIR)/lofw.o $(DIR)/four1.o \
      $(DIR)/output_summary.o $(DIR)/read_vmodelNA.o $(DIR)/write_data.o \
      $(DIR)/ran3.o $(DIR)/theo_noise.o $(DIR)/indexx.o $(DIR)/write_srf.o\
      $(DIR)/read_data_fortran.o

# make rfi_subs 

rfi_subs: $(RFI) 
	#cd ./rfi_plot; make -f Makefile

.f.o:
	$(F77) $(FFLAGS) -c $*.f
.c.o:
	$(CC) $(CFLAGS) -c $*.c
