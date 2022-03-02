NA_DIR = /usr/local/opt/Demo-NA# path to NA installation
NA_OBJ = \
		 $(NA_DIR)/src/lib/na.o \
		 $(NA_DIR)/src/lib/utils.o
include $(NA_DIR)/src/makena.macros

MAIN = homti_na

SUBS_DIR = homti_subs
SUBS_OBJ = \
         $(SUBS_DIR)/forward_modelling.o \
         $(SUBS_DIR)/calcmisfit.o \
         $(SUBS_DIR)/display_model.o \
         $(SUBS_DIR)/display_final.o \
         $(SUBS_DIR)/read_vmodelNA.o \
         $(SUBS_DIR)/lofw.o \
         $(SUBS_DIR)/output_summary.o


$(SUBS_OBJ): $(SUBS_DIR)/%.o: $(SUBS_DIR)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

subs: $(SUBS_OBJ)

main: $(SUBS_OBJ)
ifeq ($(strip $(TYPE)), serial)
	$(F77) $(FFLAGS) -o ../bin/$(MAIN) $(MAIN).f \
	$(NA_OBJ) $(SUBS_OBJ)
else ifeq ($(strip $(TYPE)), mpi)
	$(MF77) $(MPI) $(FFLAGS) -o ../bin/$(MAIN) $(MAIN).f \
	$(NA_OBJ) $(SUBS_OBJ)
else
	echo "Wrong TYPE ($(TYPE)), must be either 'serial' or 'mpi'."
endif

all: na subs main

na:
ifeq ($(strip $(TYPE)), serial)
	cd $(NA_DIR)/src; make -f Makefile na
else ifeq ($(strip $(TYPE)), mpi)
	cd $(NA_DIR)/src; make -f Makefile na_mpi
endif

clean:
	rm -f ../bin/*
	rm -f $(SUBS_OBJ)
	cd $(NA_DIR)/src; make -f Makefile clean