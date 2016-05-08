#Makefile
CC := h5c++

OMPFLAG := -openmp

#To get compile and linker flags with g++
MPIraw := $(shell mpic++ -show ./test.false)
HDF5raw := $(shell h5c++ -show ./test.false)

MPInoprog = $(subst g++ ,,${MPIraw})
MPIlinker = $(shell echo ${MPInoprog} | awk 'BEGIN {FS="./test.false"} {print $$1}')
MPIcompile = $(shell echo ${MPInoprog} | awk 'BEGIN {FS="./test.false"} {print $$2}')

HDF5noprog = $(subst g++ ,,${HDF5raw})
HDF5linker = $(shell echo ${HDF5noprog} | awk 'BEGIN {FS="./test.false"} {print $$1}')
HDF5compile = $(shell echo ${HDF5noprog} | awk 'BEGIN {FS="./test.false"} {print $$2}')

CompileFlags = ${MPIcompile} -std=c++11 -ggdb ${HDF5compile}
LinkerFlags = ${MPIlinker} ${HDF5linker}

EXEC := fluid

info:
	@echo "Usage:\tmake all\n\tmake clean"

all: fluid

fluid: main.o Simulation.o
	$(CC) $(CompileFlags) main.o Simulation.o -o fluid ${LinkerFlags}
#	$(RM) main.o Simulation.o

main.o: main.cpp
	${CC} $(CompileFlags) -c ./main.cpp ${LinkerFlags}

Simulation.o: Simulation.cpp Simulation.h
	$(CC) $(CompileFlags) -c ./Simulation.cpp ${LinkerFlags}  

clean:
	$(RM) $(EXEC) main.o Simulation.o
