FC = gfortran
FFLAGS = -O3 -std=f95 -fall-intrinsics
#-Wall 


OBJECTS = rnglib.o ranlib_poisson.o LabFuncs.o utils.o expt.o  modulation.o DDrate.o like.o

all: test GenerateEvents

test: $(OBJECTS)
	$(FC) $(FFLAGS) -o test test.f90 $(OBJECTS)

GenerateEvents: $(OBJECTS)
	$(FC) $(FFLAGS) -o GenerateEvents GenerateEvents.f90 $(OBJECTS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm -f *.o
	rm -f *.mod
	rm test
	rm GenerateEvents