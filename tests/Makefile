# PTAPI - tests
#   make test    build and run the unit tests
#   make tsan    the same under ThreadSanitizer (data races)
#   make asan    the same under AddressSanitizer + UBSan (memory errors)
#   make all     test + tsan + asan
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -Wshadow -pthread
SRC      = test_ptapi.cpp test_second_unit.cpp
HDR      = $(wildcard ../include/*.h)

test: test_ptapi
	./test_ptapi

tsan: test_ptapi_tsan
	./test_ptapi_tsan

asan: test_ptapi_asan
	./test_ptapi_asan

all: test tsan asan

test_ptapi: $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -O2 $(SRC) -o $@

test_ptapi_tsan: $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -O1 -g -fsanitize=thread $(SRC) -o $@

test_ptapi_asan: $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined $(SRC) -o $@

clean:
	rm -f test_ptapi test_ptapi_tsan test_ptapi_asan

.PHONY: test tsan asan all clean
