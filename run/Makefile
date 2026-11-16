# HENDOS - HErmite N-body Dynamics Of Systems
#
# Run make in this directory (HENDOS/run):
#   make            build the executables hermite and makeics here
#   make hermite    the N-body code only (needs only a C compiler)
#   make makeics    the initial conditions generator only (needs Fortran)
#   make clean      remove the executables
#
# The N-body code is in ../src (C); the initial conditions generator is
# in ../src/ics (Fortran). Output from runs goes in ./data (or wherever
# datadir in params.nml says).
# The old Fortran version of the N-body code (the only one with the
# Barnes-Hut tree) is kept in ../archive/fortran, with its own Makefile.

CC      = gcc
CFLAGS  = -O3 -march=native -fopenmp -std=gnu11 -Wall
FC      = gfortran
FFLAGS  = -O3 -march=native -fopenmp

SRC     = ../src
ICS     = $(SRC)/ics

# the N-body code: one header and these files (see ../src/README.md)
CFILES  = $(SRC)/main.c $(SRC)/params.c $(SRC)/io.c $(SRC)/force.c \
          $(SRC)/ac.c $(SRC)/few.c $(SRC)/groups.c $(SRC)/restart.c

all: hermite makeics

hermite: $(CFILES) $(SRC)/hendos.h
	$(CC) $(CFLAGS) -o $@ $(CFILES) -lm

# initial conditions generator - the program makeics.f90 plus the files
# of routines it calls (they use no modules)
ICSRC   = $(ICS)/icrandom.f90 $(ICS)/plummersubs.f90 $(ICS)/fractalsubs.f90 \
          $(ICS)/masssubs.f90 $(ICS)/usersubs.f90 $(ICS)/makeics.f90

makeics: $(ICSRC)
	$(FC) $(FFLAGS) -o $@ $(ICSRC)

clean:
	rm -f hermite makeics

.PHONY: all clean
