IDIR = header
CXX = g++
LDFLAGS=-lm

CXXFLAGS = -std=c++14 -O2 -Wall -I$(IDIR)
DEBUGFLAGS = -std=c++14 -g -Wall -I$(IDIR)

ODIR = obj
SDIR = src
OUTDIR = .
DEBUGODIR = obj/Debug
DEBUGDIR = bin/Debug

GENMOD = circle_circle_intersection auxiliary_functions

_DEPS =
DEPS = $(patsubst %, $(IDIR)/%.h, $(_DEPS) $(GENMOD))

_OBJ = main
OBJ = $(patsubst %, $(ODIR)/%.o, $(_OBJ) $(GENMOD))

DEBUGOBJ = $(patsubst %, $(DEBUGODIR)/%.o, $(_OBJ) $(GENMOD))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(LDFLAGS)

$(OUTDIR)/ayudantia.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

$(DEBUGODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(DEBUGFLAGS)

$(DEBUGDIR)/AyudantiaDebug: $(DEBUGOBJ)
	$(CXX) -o $@ $^ $(DEBUGFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

