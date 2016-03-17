
include ../GNUmakefile.common.mk


EXE := dock
OBJ += dock.o
OBJ += load.o data.o rmsd.o util.o util_print.o

default: $(EXE)

$(EXE): $(OBJ)
	$(CXX_HOST) $(LIBPATH) $(OBJ) -o $@ $(LINKFLAGS)

dock.o: dock.C
	$(CXX_HOST) $(CXXFLAGS_HOST) -c $< -o $@

include ../GNUmakefile.common.rules

test: dock
	./dock


