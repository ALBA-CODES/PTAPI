# PTAPI v2 - tuning the PT with irace (example: SSP)
#   make         builds the program that irace runs
#   make check   one call of the target runner, as irace does it
#   make tune    runs irace (needs R and the irace package)
#   make clean
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -pthread
IRACE    ?= $(shell Rscript -e 'cat(system.file("bin", "irace", package="irace"))' 2>/dev/null)

ptSSP: ptSSP.cpp ../example/SSP.cpp ../example/SSP.h $(wildcard ../include/*.h)
	$(CXX) $(CXXFLAGS) ptSSP.cpp ../example/SSP.cpp -o $@

check: ptSSP
	./target-runner 1 1 1234 ../instances/SSP/F1001.txt --TEMP_INIT 0.2 --TEMP_FIM 1 --N_REPLICAS 16 --MCL 400

tune: ptSSP
	@test -n "$(IRACE)" || { echo "irace not found: install it in R with install.packages(\"irace\")"; exit 1; }
	$(IRACE) --scenario scenario.txt

clean:
	rm -f ptSSP irace.Rdata

.PHONY: check tune clean
