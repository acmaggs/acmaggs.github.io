#CXX=icc
OBJ= net.o main.o param.o mc.o recordings.o  energy.o wiring.o output.o 
CXXFLAGS=-g -Wall -pedantic   -Wno-unused  -Wno-shadow  -O
#MATLIB=-L/misc/stubbsopt/matlab12/bin/glnx86 -L/misc/stubbsopt/matlab12/extern/lib/glnx86  -lmx -lut -lmat -lmi -lm

#CPPFLAGS=-I/misc/stubbsopt/matlab12/extern/include 

all: simul


simul: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(MATLIB)


clean:
	rm -f *.o

depend:
	makedepend *.cc *.h
backup:tar
	scp elec2b.tar.gz stubbs:
	scp elec2b.tar.gz hogarth:
	scp elec2b.tar.gz hockney:
tar:
	tar cvfz elec2b.tar.gz *.cc *.h Makefile License
# DO NOT DELETE
