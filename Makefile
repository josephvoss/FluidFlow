#Makefile
CC := mpic++
OMPFLAG := -openmp
FLAGS := -std=c++11 -ggdb

EXEC := fluid

info:
	@echo "Usage:\tmake all\n\tmake clean"

all: fluid

fluid: main.o Simulation.o
	$(CC) $(FLAGS) main.o Simulation.o -o fluid
#:	$(RM) main.o Simulation.o

main.o: main.cpp
	${CC} $(FLAGS) -c ./main.cpp

Simulation.o:
	$(CC) $(FLAGS) -c ./Simulation.cpp

clean:
	$(RM) $(EXEC) main.o Simulation.o
