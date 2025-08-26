# Auto-detect OS
UNAME_S := $(shell uname -s)

# Default settings for Linux with Intel compilers
ifeq ($(UNAME_S),Linux)
	FC = mpiifort
	FFLG = -r8 -O3 -mkl=parallel -warn nounused #-fltconsistency
	FINC = -I${MKLROOT}/include/fftw
	# FINC = -I/opt/intel/oneapi/mkl/latest/include/fftw
	FLIB = -mkl
	# FLIB =  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 \
	#        -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread \
	#        -lm -ldl
	MODULE_FLAG = -module
	# Set OpenMP flags
	FFLG += -qopenmp
	# Add debug flags and sanitization options
	# FFLG += -g -traceback -check all -fp-stck-check

# macOS settings with homebrew packages
else ifeq ($(UNAME_S),Darwin)
	FC = mpifort
	FFLG = -fdefault-real-8 -fdefault-double-8 -O3 -ffree-line-length-none
	FINC = -I$(shell brew --prefix fftw)/include -I$(shell brew --prefix openblas)/include
	FLIB = -L$(shell brew --prefix fftw)/lib -lfftw3 -L$(shell brew --prefix openblas)/lib -lopenblas -lgcc
	MODULE_FLAG = -J
	# Set OpenMP flags
	FFLG += -fopenmp
	# Add debug flags and sanitization options
	FFLG += -g
	# Set mpifort to use gfortran
	export OMPI_FC = gfortran
endif

# Directory structure definitions
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod
BIN_DIR = bin
OUTPUT_DIR = output

# List of source files
SRC_F = $(wildcard $(SRC_DIR)/*.f90)

# List of object files
OBJ_F = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRC_F)) $(wildcard $(OBJ_DIR)/*.o)

# Create directories
$(shell mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR) $(OUTPUT_DIR))

# Add module directory to compiler flags
FFLG += $(MODULE_FLAG) $(MOD_DIR)

# Scan the f90 folder for all .f90 scripts
F90_DIR = f90
F90_SCRIPTS = $(basename $(notdir $(wildcard $(F90_DIR)/*.f90)))
EXE_F = $(F90_SCRIPTS)

# # Program scripts in the f90 folder
# EXE_F = addperturb_non \
# 	addperturb_split \
# 	evp_parametric \
# 	evp_print \
# 	init \
# 	vort \
# 	vort9 \
# 	postproc_mpi \
# 	bsnsq_ivp \ 
# 	test_prodct \

# Phony targets
.PHONY: all clean new swipe $(EXE_F)

# Default target to compile all executables
all: $(addprefix $(BIN_DIR)/,$(addsuffix _exec,$(EXE_F)))

# Create an object file from fortran source files
# Explicit dependency on the directories ensures they exist
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLG) $(FINC) $(FLIB) -c $< -o $@

# Produce executable files
$(BIN_DIR)/%_exec: f90/%.f90 $(OBJ_F)
	$(FC) $(FFLG) $(FINC) $(FLIB) $^ -o $@

# Individual targets for each executable
$(EXE_F): %: $(BIN_DIR)/%_exec

# Remove all objects, modules and executable files
# Preserve lib_fm files in obj and mod directories
clean:
	rm -f ./*.o ./*.mod ./*.dat ./*.DAT ./*.info
	find ./$(OBJ_DIR)/ -type f ! -name 'fm*' -exec rm -f {} \;
	find ./$(MOD_DIR)/ -type f ! -name 'fm*' -exec rm -f {} \;
	rm -rf ./$(BIN_DIR)/*
	rm -rf ./$(OUTPUT_DIR)/*.dat
	rm -rf ./$(OUTPUT_DIR)/*.info

# Remove all data files but keep compiled code
swipe:
	rm -f ./*.dat ./*.DAT ./*.info
	rm -rf ./$(OUTPUT_DIR)/*.dat
	rm -rf ./$(OUTPUT_DIR)/*.DAT 
	rm -rf ./$(OUTPUT_DIR)/*.info
	rm -rf ./$(OUTPUT_DIR)/especData_MK/*.output

# Clean and produce all executable files
new: clean all

# Include module dependencies
include Makefile.depend
