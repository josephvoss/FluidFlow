#Makefile
ERR = $(shell which icpc >/dev/null; echo $$?)
ifeq "$(ERR)" "0"
    CC := icc
    OMPFLAG := -openmp
else	
    CC := gcc
    OMPFLAG := -fopenmp
endif

EXECS := main
OBJ := $(EXECS:=.o)
SRC := $(OBJ:.o=.c)

main: main.c
	${CC} $(OMPFLAG) -O0 -o graphTransversal $@.c -lm

clean:
	$(RM) $(EXECS)
