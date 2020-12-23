CC      = gcc                   # C compiler
F77     = gfortran              # Fortran compiler
MF77    = mpif77				# MPI fortran compiler
FFLAGS  = -O -ffixed-line-length-300 -fallow-argument-mismatch   # Fortran compiler options
CFLAGS  = -O -Wno-return-type   # C compiler options
MPI     = -DNA_MPIP=1           # MPI option for continuous parallel

MAIN = homti_na
TYPE = serial

SUBS_DIR = homti_subs
SUBS_OBJ = \
         $(SUBS_DIR)/forward_modelling.o \
         $(SUBS_DIR)/calcmisfit.o \
         $(SUBS_DIR)/display_model.o \
         $(SUBS_DIR)/display_final.o \
         $(SUBS_DIR)/read_vmodelNA.o

NA_DIR = /usr/local/opt/NA# path to NA installation
NA_OBJ = \
		 $(NA_DIR)/src/lib/na.o \
		 $(NA_DIR)/src/lib/utils.o \
		 $(NA_DIR)/src/rfi_subs/lofw.o \
		 $(NA_DIR)/src/rfi_subs/output_summary.o

$(SUBS_OBJ): $(SUBS_DIR)/%.o: $(SUBS_DIR)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

subs: $(SUBS_OBJ)

main: $(SUBS_OBJ)
ifeq ($(TYPE), serial)
	$(F77) $(FFLAGS) -o ../bin/$(MAIN) $(MAIN).f \
	$(NA_OBJ) $(SUBS_OBJ)
else ifeq ($(TYPE), mpi)
	$(MF77) $(MPI) $(FFLAGS) -o ../bin/$(MAIN) $(MAIN).f \
	$(NA_OBJ) $(SUBS_OBJ)
else
	echo "Wrong TYPE, must be either 'serial' or 'mpi'."
endif

all: subs main

clean:
	rm -f ../bin/*
	rm -f $(SUBS_OBJ)