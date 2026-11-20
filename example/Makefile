# PTAPI - examples
#   make              quickstart, exampleTSP and the three examples of version 1
#   make quickstart   only the quick start (docs/getting-started.md)
#   make clean
# One example for each group of functions: see functions/ (make run).
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -pthread
HDR      = $(wildcard ../include/*.h)
PROGRAMS = quickstart exampleTSP mainSSP mainIPMTC mainRCPMS

all: $(PROGRAMS)

quickstart: quickstart.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $< -o $@
exampleTSP: exampleTSP.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $< -o $@
main%: main%.cpp %.cpp %.h $(HDR)
	$(CXX) $(CXXFLAGS) $< $*.cpp -o $@

clean:
	rm -f $(PROGRAMS)

.PHONY: all clean
