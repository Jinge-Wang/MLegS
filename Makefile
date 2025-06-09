FC = mpiifort 
# Fortran compiler (FC) (ifort: Intel compiler - recommended)

SRC_DIR = src
# Source code directory
OBJ_DIR = obj
# Objective file directory
MOD_DIR = mod
# Module file directory
BIN_DIR = bin
# Bin file directory
OUTPUT_DIR = output

#FFLG = -g -traceback -check all -fp-stack-check -qopenmp -warn nonused
#FFLG = -g -traceback -qopenmp -warn nonused
#FFLG = -r8 -O3 -qopenmp -warn nounused
FFLG = -r8 -O3 -qopenmp -mkl=parallel -warn nounused #-fltconsistency
# Fortran compiler flag
FINC = -I${MKLROOT}/include/fftw 
#FINC = -I/opt/intel/oneapi/mkl/latest/include/fftw
# Fortran include files
FLIB = -mkl
#FLIB =  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 \
#        -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread \
#        -lm -ldl
# Fortran libraries

SRC_F = $(wildcard $(SRC_DIR)/*.f90)
# List of source files.

# OBJ_C = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_C))
OBJ_F = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRC_F)) $(wildcard $(OBJ_DIR)/*.o)
# List of object files

ifneq ($(OBJ_DIR),)
  $(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))
endif
# If the obj folder not exists, make this directory

ifneq ($(MOD_DIR),)
  $(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
  FFLG+= -module $(MOD_DIR)
endif
# If the obj folder not exists, make this directory

ifneq ($(BIN_DIR),)
  $(shell test -d $(BIN_DIR) || mkdir -p $(BIN_DIR))
endif
# If the obj folder not exists, make this directory

ifneq ($(OUTPUT_DIR),)
  $(shell test -d $(OUTPUT_DIR) || mkdir -p $(OUTPUT_DIR))
endif
# If the output folder not exists, make this directory

EXE_F = addperturb_non \
	addperturb_split \
	evp_parametric \
	evp_print \
	init \
	vort \
	vort9 \
	postproc_mpi \
# Program scripts in the f90 folder (ADD/REMOVE THE PROGRAM LISTS HERE)
# Using 'make [program_name]' will create the executable program file in the bin folder
# 'make new' will wipe out all obj, mod and exec files and re-compile the first exe_f file

.PHONY: clean new

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(FC) $(FFLG) $(FINC) $(FLIB) -c $< -o $@
# Create an object file from c source files
# Module files go to the MOD_DIR directory

all: $(OBJ_F)

$(EXE_F): $(OBJ_F)
	$(FC) $(FFLG) $(FINC) $(FLIB) ./f90/$@.f90 $^ -o $(BIN_DIR)/$@_exec
# Produce an executable file labeled by _exec

clean:
	rm -f ./*.o ./*.mod ./*.dat ./*.DAT ./*.info
	find ./$(OBJ_DIR)/ ! -name 'fm*' -exec rm -f {} \;
	find ./$(MOD_DIR)/ ! -name 'fm*' -exec rm -f {} \;
	rm -rf ./$(BIN_DIR)/*
	rm -rf ./$(OUTPUT_DIR)/*.dat
	rm -rf ./$(OUTPUT_DIR)/*.info
# Remove all objects, modules and executable files

swipe:
	rm -f ./*.dat ./*.DAT ./*.info
	rm -rf ./$(OUTPUT_DIR)/*.dat
	rm -rf ./$(OUTPUT_DIR)/*.DAT 
	rm -rf ./$(OUTPUT_DIR)/*.info
	rm -rf ./$(OUTPUT_DIR)/especData_MK/*.output

new:
	make clean
	make
# Clean and produce an executable file from the first program script in EXE_F

#include Makefile.dependencies
include Makefile.depend
# Reference to module dependencies
