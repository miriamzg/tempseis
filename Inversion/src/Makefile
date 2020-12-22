CC      = gcc                   # C compiler
F77     = gfortran                 # Fortran compiler
FFLAGS  = -O   -ffixed-line-length-300  -fallow-argument-mismatch        # Fortran compiler options
CFLAGS  = -O   -Wno-return-type	        # C compiler options


MAIN = homti_na

SUBS_DIR = homti_subs
SUBS_OBJ = \
         $(SUBS_DIR)/forward_modelling.o \
         $(SUBS_DIR)/calcmisfit.o \
         $(SUBS_DIR)/display_model.o \
         $(SUBS_DIR)/display_final.o \
         $(SUBS_DIR)/read_vmodelNA.o

NA_DIR = /usr/local/opt/NA# directory to NA installation
NA_OBJ = \
		 $(NA_DIR)/src/lib/na.o \
		 $(NA_DIR)/src/lib/utils.o \
		 $(NA_DIR)/src/rfi_subs/lofw.o \
		 $(NA_DIR)/src/rfi_subs/output_summary.o

$(SUBS_OBJ): $(SUBS_DIR)/%.o: $(SUBS_DIR)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

subs: $(SUBS_OBJ)

main: $(SUBS_OBJ)
	$(F77) $(FFLAGS) -o ../bin/$(MAIN) $(MAIN).f \
	$(NA_OBJ) $(SUBS_OBJ)

all: subs main

clean:
	rm -f ../bin/*
	rm -f $(SUBS_OBJ)