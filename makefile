CC:= g++

CCFLAGS:=  -pthread -Wall -fPIC -pedantic -Wno-long-long -g -O3 -fopenmp

CCOFLAGS:= -fexpensive-optimizations

CCFLAGS+=$(CCOFLAGS)

GLOBALDEPEND:= $(wildcard include/*.hpp)


TEST:= test/RunTests

EXTRA_LIBS:= extra_libs

RLLIB:= RLlib/libRLlib.a

INCLUDE:= -Iinclude -IRLlib/include -Iarpack

LIBS:=  -lm -LRLlib -lRLlib -lconfig++
LIBS+= -liomp5

SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o) 


MAINS := Flux
MAINSOBJ:= $(MAINS:%=bin/%.o)
MAINSNAME := Flux



.PHONY: clean doc all

all: bin $(RLLIB) $(OBJ) $(MAINS) $(TEST)


bin: 
	mkdir -p bin

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) bin/%.o
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $^ $(LIBS)  -o $(MAINSNAME)

$(OBJ): bin/%.o: src/%.cc include/%.hh $(GLOBALDEPEND)
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@

$(TEST): $(MAINS)
	@$(MAKE) -C test makerun

clean: backup_clean
	@echo Cleaning up...
	@rm -rf bin/*
	@rm -rf $(MAINSOBJ)
	@$(MAKE) -C RLlib clean
	@$(MAKE) -C test clean

backup_clean: 
	@echo Removing tilde and hashtag files...
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;

$(RLLIB):
	$(MAKE) -C RLlib
