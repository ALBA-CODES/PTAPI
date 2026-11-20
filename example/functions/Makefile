# PTAPI v2 - one example for each group of functions of the API
#   make          build all the examples
#   make run      build and run them all (each one returns 0 if it behaved as described)
#   make clean
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -pthread
HDR      = tsp.h $(wildcard ../../include/*.h)
EXAMPLES = exampleSeed exampleStop exampleEvaluate exampleTemperatures exampleInitialSolutions \
           exampleTempUpdate exampleCalibrateTemperatures exampleCalibrateLadder exampleMoves \
           exampleCheckpoint exampleSwapScheme exampleStatistics exampleV1Compatibility \
           exampleAutoUpdate exampleTiming exampleExperiment exampleComparison exampleOptions

all: $(EXAMPLES)

%: %.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $< -o $@

run: all
	@for e in $(EXAMPLES); do \
		echo "=== $$e ==="; \
		if [ $$e = exampleCheckpoint ]; then ./$$e demo; else ./$$e; fi || { echo "FAILED: $$e"; exit 1; }; \
		echo; \
	done; echo "all the examples ran"

clean:
	rm -f $(EXAMPLES) exampleCheckpoint.ckpt
	rm -rf exampleExperiment.out exampleComparison.out

.PHONY: all run clean
