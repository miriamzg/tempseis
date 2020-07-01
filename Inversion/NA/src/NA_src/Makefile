#
#	Make file for program NA.f
#

# location of fortran compiler

F77 = /opt/local/SUNWspro/bin/f77

# compiler options

FFLAGS = -O

#---- start user specific make --------------------------------------------

# edit the following to link in user routines

# user supplied subroutines called by user interface routines
# must be edited for user problem

DIR = rfi_subs
RFI = $(DIR)/theo.o $(DIR)/qlayer.o $(DIR)/read_data.o $(DIR)/read_seed.o \
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o $(DIR)/display_model.o \
      $(DIR)/display_final.o $(DIR)/lofw.o $(DIR)/cft.o \
      $(DIR)/output_summary.o $(DIR)/read_vmodelNA.o $(DIR)/write_dataNA.o \
      $(DIR)/sac.SunOS5.a

USERFILES = lib/rfi.o $(RFI)

#---- end user specific make ----------------------------------------------

all: user na_subs na 

na:
	$(F77) $(FFLAGS) -o ../bin/na na.f \
	lib/na_subs.o $(USERFILES) 
	cp na.f Backups

na_subs:
	$(F77) $(FFLAGS) -c na_subs.f
	cp na_subs.f Backups
	mv na_subs.o lib

#---- start user specific make --------------------------------------------

# edit the following to compile user routines

user:
	$(F77) $(FFLAGS) -c rfi.f
	mv rfi.o lib
	cp rfi.f Backups

#---- end user specific make ----------------------------------------------
