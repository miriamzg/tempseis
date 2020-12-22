CC      = gcc                   # C compiler
F77     = gfortran                 # Fortran compiler
FFLAGS  = -O   -ffixed-line-length-300  -fallow-argument-mismatch        # Fortran compiler options
CFLAGS  = -O   -Wno-return-type	        # C compiler options



SUBS_DIR = homnti_subs
SUBS_OBJ = \
         $(SUBS_DIR)/forward_modelling.o \
         $(SUBS_DIR)/calcmisfit.o \
         $(SUBS_DIR)/display_model.o \
         $(SUBS_DIR)/display_final.o \
         $(SUBS_DIR)/read_vmodelNA.o

NA_DIR = /usr/local/opt/NA	# directory to NA installation
NA_OBJ = \
		 $(NA_DIR)/src/lib/na.o \
		 $(NA_DIR)/src/lib/utils.o


main:
	$(F77) $(FFLAGS) -o ../bin/$(USER_MAIN) $(USER_MAIN).F \
	lib/na.o lib/utils.o $(USER_OBJ)
