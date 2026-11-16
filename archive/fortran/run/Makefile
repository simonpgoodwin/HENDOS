# HENDOS - the ARCHIVED Fortran version of the N-body code
#
# Run make in this directory (HENDOS/archive/fortran/run):
#   make            build hermite_fortran here
#   make clean      remove it and the build directory
#
# This version is no longer the default (that is the C code in
# HENDOS/src, built in HENDOS/run). It is kept because it is the only
# version with the Barnes-Hut tree. It reads the same parameter and
# initial conditions files and writes the same output files.
# NOTE: it has never been compiled - see ../README.md.

FC      = gfortran
FFLAGS  = -O3 -march=native -fopenmp
# flag that tells the compiler where to put .mod files
# (gfortran: -J   ifort/ifx: -module )
MODFLAG = -J

SRC     = ../src
BUILD   = build

# modules in the order they must be compiled (each uses ones before it)
MODS    = nbvars nbsetup nbtree nbac nbforce nbstep nbfew nbbin nbout \
          nbrestart nbana hermite
OBJS    = $(addprefix $(BUILD)/,$(addsuffix .o,$(MODS)))

all: hermite_fortran

hermite_fortran: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# which modules each file uses
$(BUILD)/nbsetup.o: $(BUILD)/nbvars.o
$(BUILD)/nbtree.o:  $(BUILD)/nbvars.o
$(BUILD)/nbac.o:    $(BUILD)/nbvars.o
$(BUILD)/nbforce.o: $(BUILD)/nbvars.o $(BUILD)/nbtree.o $(BUILD)/nbac.o
$(BUILD)/nbstep.o:  $(BUILD)/nbvars.o $(BUILD)/nbtree.o $(BUILD)/nbac.o $(BUILD)/nbforce.o
$(BUILD)/nbbin.o:   $(BUILD)/nbvars.o $(BUILD)/nbtree.o $(BUILD)/nbac.o $(BUILD)/nbforce.o \
                    $(BUILD)/nbstep.o $(BUILD)/nbfew.o
$(BUILD)/nbout.o:   $(BUILD)/nbvars.o $(BUILD)/nbforce.o $(BUILD)/nbbin.o
$(BUILD)/nbrestart.o: $(BUILD)/nbvars.o $(BUILD)/nbsetup.o $(BUILD)/nbac.o $(BUILD)/nbout.o
$(BUILD)/nbana.o:  $(BUILD)/nbvars.o $(BUILD)/nbout.o
$(BUILD)/hermite.o: $(BUILD)/nbvars.o $(BUILD)/nbsetup.o $(BUILD)/nbtree.o \
                    $(BUILD)/nbac.o $(BUILD)/nbforce.o $(BUILD)/nbstep.o \
                    $(BUILD)/nbbin.o $(BUILD)/nbout.o $(BUILD)/nbrestart.o \
                    $(BUILD)/nbana.o

$(BUILD)/%.o: $(SRC)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG)$(BUILD) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD) hermite_fortran

.PHONY: all clean
