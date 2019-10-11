FC = gfortran
FFLAGS = -O3 -std=f95 -fall-intrinsics
#-Wall 


OBJECTS = rnglib.o ranlib_poisson.o LabFuncs.o utils.o expt.o  modulation.o DDrate.o like.o

all: test GenerateEvents GenerateEvents_Asimov GridLike testSpectra

test: $(OBJECTS)
	$(FC) $(FFLAGS) -o test test.f90 $(OBJECTS)

testSpectra: $(OBJECTS)
	$(FC) $(FFLAGS) -o testSpectra testSpectra.f90 $(OBJECTS)


GridLike: $(OBJECTS)
	$(FC) $(FFLAGS) -o GridLike GridLike.f90 $(OBJECTS)

GenerateEvents: $(OBJECTS)
	$(FC) $(FFLAGS) -o GenerateEvents GenerateEvents.f90 $(OBJECTS)

GenerateEvents_Asimov: $(OBJECTS)
	$(FC) $(FFLAGS) -o GenerateEvents_Asimov GenerateEvents_Asimov.f90 $(OBJECTS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm -f *.o
	rm -f *.mod
	rm -f test
	rm -f GenerateEvents
	rm -f GenerateEvents_Asimov
	rm -f GridLike
	rm -f GridLike2
	rm -f testSpectra 