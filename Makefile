#Makefile
CC := mpic++
OMPFLAG := -openmp

EXEC := fluid

info:
	@echo "Usage:\tmake all\n\tmake clean"

all: fluid

fluid: main.o Simulation.o
	$(CC) main.o Simulation.o -o fluid
	$(RM) main.o Simulation.o

main.o: main.c
	${CC} -c ./main.cpp

Simulation.o:
	$(CC) -c ./Simulation.cpp

clean:
	$(RM) $(EXEC) main.o Simulation.o
