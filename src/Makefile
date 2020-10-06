FC = gfortran
FFLAGS = -O3 -std=f95 -fall-intrinsics
#-Wall 


OBJECTS = rnglib.o ranlib_poisson.o LabFuncs.o utils.o expt.o  modulation.o DDrate.o like.o cdf.o d1mach.o dzero.o stat.o

testLikes: $(OBJECTS)
	$(FC) $(FFLAGS) -o ../testLikes testLikes.f90 $(OBJECTS)

calcContour: $(OBJECTS)
	$(FC) $(FFLAGS) -o ../calcContour calcContour.f90 $(OBJECTS)

testSpectra: $(OBJECTS)
	$(FC) $(FFLAGS) -o ../testSpectra testSpectra.f90 $(OBJECTS)

GenerateEvents: $(OBJECTS)
	$(FC) $(FFLAGS) -o ../GenerateEvents GenerateEvents.f90 $(OBJECTS)

GenerateEvents_Asimov: $(OBJECTS)
	$(FC) $(FFLAGS) -o ../GenerateEvents_Asimov GenerateEvents_Asimov.f90 $(OBJECTS)

all: GenerateEvents GenerateEvents_Asimov testSpectra calcContour testLikes

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^

clean:
	rm -f *.o
	rm -f *.mod
	rm -f ../GenerateEvents
	rm -f ../GenerateEvents_Asimov
	rm -f ../testSpectra 
	rm -f ../calcContour
	rm -f ../testLikes
