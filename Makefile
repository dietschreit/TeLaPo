CCOMP = gcc
#OPT = -m64 -O3 -Wall -lm # linux (Ubuntu) version
OPT = -m64 -O3 -Wall # MacOS X version


INCLUDE = include
SRC = src
BIN = bin

vpath %.h $(INCLUDE)
vpath %.o $(BIN)
vpath %.c $(SRC)
.SUFFIXES:
PROGRAMS = tela
VLEARN_OBJ = main.o getArgs.o output.o read_dcd.o time.o observables.o pair_correlation_function.o statistics.o entropy.o bricks.o checks.o random_knuth.o tetra_rw.o tetra_fww.o tetra_saw1.o tetra_saw2.o tetra_saw3.o tetra_fsaw3.o tetra_bsaw3.o tetra_fbsaw3.o intra_potential.o SASA.o


DEPENDFILE = .depend
SOURCES = $(wildcard $(SRC)/*.c)

all: depend $(PROGRAMS)

depend: $(SOURCES)
	$(CCOMP) -I$(INCLUDE) -MM $(SOURCES) > $(DEPENDFILE)

main: $(VLEARN_OBJ)
	$(CCOMP) -o main $(patsubst %.o,$(BIN)/%.o,$(VLEARN_OBJ))


%.o: %.c
	$(CCOMP) -I$(INCLUDE) -o $(BIN)/$@ -c $< $(OPT)

clean :
	rm -f $(BIN)/*.o
	for i in . $(SRC) $(INCLUDE); do rm -f $$i/*~; done
	for i in $(PROGRAMS); do rm -f $$i; done
