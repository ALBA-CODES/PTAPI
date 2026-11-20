# PTAPI - performance comparison between the version 1 (../v1/include)
# and v2 (../include). Always run it before accepting a change in the core.
#
#   make bench      synthetic problem (cheap evaluation: shows the cost of the API)
#   make bench-sd   real problem of ../v1/SD (expensive evaluation)
#
# Each line shows: wall time, CPU time, evaluations, evaluations/s.
# Compare the versions on an idle machine, or at least run them interleaved.
CXX  ?= g++
OLD  ?= ../v1
B    = build
INST ?= $(OLD)/instances/J50C25_1.txt

$(B)/old/bench $(B)/new/bench: bench.cpp $(wildcard ../include/*.h)
	rm -rf $(B)/old $(B)/new && mkdir -p $(B)/old $(B)/new
	cp -r $(OLD)/SD $(OLD)/include $(B)/old/ && cp bench.cpp $(B)/old/
	cp -r $(OLD)/SD ../include     $(B)/new/ && cp bench.cpp $(B)/new/
	cd $(B)/old && $(CXX) bench.cpp SD/SD.cpp -std=c++20 -O3 -pthread -o bench
	cd $(B)/new && $(CXX) bench.cpp SD/SD.cpp -std=c++20 -O3 -pthread -o bench

# replicas MCL PTL threads [time] [updateType] [updateInterval] [tempMin] [tempMax]
bench: $(B)/old/bench $(B)/new/bench
	@for cfg in "8 500 2000 8" "8 500 2000 24" "30 500 500 24" "16 20 20000 16"; do \
		for r in 1 2 3; do for v in old new; do printf "%-4s tsp:100 [%s] " $$v "$$cfg"; $(B)/$$v/bench tsp:100 $$cfg 0 0 50 1 60 | cut -d' ' -f1-4; done; done; done

bench-sd: $(B)/old/bench $(B)/new/bench
	@for cfg in "8 300 60 8" "8 300 60 24" "30 300 16 24"; do \
		for r in 1 2 3; do for v in old new; do printf "%-4s SD [%s] " $$v "$$cfg"; $(B)/$$v/bench $(INST) $$cfg | cut -d' ' -f1-5; done; done; done

clean:
	rm -rf $(B)

.PHONY: bench bench-sd clean
