# HENDOS - HErmite N-body Dynamics Of Systems
#
# Run make in this directory (HENDOS/run):
#   make            build the executables hermite and makeics here
#   make hermite_c  build the C version of the main code (needs only a C
#                   compiler; see ../csrc/CSRC_README.md)
#   make clean      remove the executables and the build directory
#
# The source code is in ../src (initial conditions generators in
# ../src/ics). Object and module files go in ./build, out of the way.
# Output from runs goes in ./data (or wherever datadir in params.nml says).

FC      = gfortran
FFLAGS  = -O3 -march=native -fopenmp
# flag that tells the compiler where to put .mod files
# (gfortran: -J   ifort/ifx: -module )
MODFLAG = -J

SRC     = ../src
ICS     = $(SRC)/ics
BUILD   = build

# modules in the order they must be compiled (each uses ones before it)
MODS    = nbvars nbsetup nbtree nbac nbforce nbstep nbfew nbbin nbout \
          nbrestart nbana hermite
OBJS    = $(addprefix $(BUILD)/,$(addsuffix .o,$(MODS)))

# the C version of the main code (../csrc)
CC      = gcc
CFLAGS  = -O3 -march=native -fopenmp -std=gnu11 -Wall
CSRC    = ../csrc
CFILES  = $(CSRC)/main.c $(CSRC)/params.c $(CSRC)/io.c $(CSRC)/force.c \
          $(CSRC)/ac.c $(CSRC)/few.c $(CSRC)/groups.c $(CSRC)/restart.c

all: hermite makeics

hermite_c: $(CFILES) $(CSRC)/hendos.h
	$(CC) $(CFLAGS) -o $@ $(CFILES) -lm

hermite: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# initial conditions generator - the program makeics.f90 plus the files
# of routines it calls (they use no modules)
ICSRC   = $(ICS)/icrandom.f90 $(ICS)/plummersubs.f90 $(ICS)/fractalsubs.f90 \
          $(ICS)/masssubs.f90 $(ICS)/usersubs.f90 $(ICS)/makeics.f90

makeics: $(ICSRC)
	$(FC) $(FFLAGS) -o $@ $(ICSRC)

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
	rm -rf $(BUILD) hermite hermite_c makeics

.PHONY: all clean
