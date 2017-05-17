INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= prim2matlab

TRANSFORM	= fftw3
MODSOBJ		= io.o meshs.o mpi.o parameters.o \
		  timestep.o transform.o variables.o velocity.o

#COMPILER	= g95 -C  
#COMPFLAGS	= -cpp -c -O3
#COMPILER	= ifort -i-dynamic #-C #-static 
#COMPFLAGS	= -cpp -c -O3 -heap-arrays 1024 -mcmodel=medium
#COMPILER	= pgf90 #-C
#COMPFLAGS	= -Mpreprocess -c -fast #-mcmodel=medium
#COMPILER	= pathf90 #-pg -C
#COMPFLAGS	= -cpp -c -O3 -OPT:Ofast -march=opteron -fno-second-underscore


COMPILER	= gfortran
COMPFLAGS	= -ffree-line-length-none -x f95-cpp-input -c -O3 \
		  -I/usr/include \
                  #-C #-pg
LIBS		= \
		  -L/home/ash/lib \
		  cheby.o -lfftw3 -llapack -lnetcdff -lnetcdf \
		  # -lblas -lcurl

#------------------------------------------------------------------------
all : 	$(MODSOBJ) $(PROGDIR)main.f90
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)main.f90
	$(COMPILER) -o ./main.out main.o $(MODSOBJ) $(LIBS)

install : main.out
	if test ! -d $(INSTDIR); then mkdir -p $(INSTDIR); fi
	mv ./main.out $(INSTDIR)
	date > $(INSTDIR)/main.info
	echo $(HOSTNAME) >> $(INSTDIR)/main.info
	pwd >> $(INSTDIR)/main.info
	echo $(COMPILER) $(COMPFLAGS) >> $(INSTDIR)/main.info
	grep "define _N" parallel.h >> $(INSTDIR)/main.info
	cut -d! -f1 $(PROGDIR)parameters.f90 | grep = | \
	   cut -d: -f3  >> $(INSTDIR)/main.info

util : 	$(MODSOBJ) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) $(COMPFLAGS) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) -o ./$(UTIL).out $(UTIL).o $(MODSOBJ) $(LIBS)

#------------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.il core *.out

#------------------------------------------------------------------------
io.o : $(PROGDIR)io.f90 velocity.o 
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)io.f90

meshs.o : $(PROGDIR)meshs.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)cheby.f
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)meshs.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o : $(PROGDIR)parameters.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

variables.o : $(PROGDIR)variables.f90 meshs.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 timestep.o transform.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

