IDIR = header
CXX = g++
LDFLAGS=-lm

CXXFLAGS = -std=c++11 -O2 -Wall -I$(IDIR)
DEBUGFLAGS = -std=c++11 -g -Wall -I$(IDIR)

ODIR = obj
SDIR = src
OUTDIR = .
DEBUGODIR = obj/Debug
DEBUGDIR = .

GENMOD = circle_circle_intersection auxiliary_functions cech_scale

_DEPS =
DEPS = $(patsubst %, $(IDIR)/%.h, $(_DEPS) $(GENMOD))

_OBJ = main
OBJ = $(patsubst %, $(ODIR)/%.o, $(_OBJ) $(GENMOD))

DEBUGOBJ = $(patsubst %, $(DEBUGODIR)/%.o, $(_OBJ) $(GENMOD))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(LDFLAGS)

$(OUTDIR)/Cech-scale: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

$(DEBUGODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(DEBUGFLAGS)

$(DEBUGDIR)/debug: $(DEBUGOBJ)
	$(CXX) -o $@ $^ $(DEBUGFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(DEBUGODIR)/*.o   *~ core $(INCDIR)/*~

