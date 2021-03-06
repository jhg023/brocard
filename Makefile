CXXFLAGS=-std=c++17 -march=native -mtune=native -fomit-frame-pointer -fno-stack-protector -fno-rtti -flto

all: brocard bench test

brocard:
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -Ofast -lflint -fopenmp Brocard.cpp -o bin/brocard

bench:
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -O0 -lcelero -lflint bench/factorial_bench.cpp -o bin/factorial_bench
	$(CXX) $(CXXFLAGS) -Ofast -lcelero bench/count_trailing_zeros_bench.cpp -o bin/count_trailing_zeros_bench
	$(CXX) $(CXXFLAGS) -Ofast -lcelero -lflint bench/mulmod_bench.cpp -o bin/mulmod_bench
	$(CXX) $(CXXFLAGS) -Ofast -lcelero -lflint bench/jacobi_bench.cpp -o bin/jacobi_bench

clean:
	rm -f bin/*

test:
	$(CXX) $(CXXFLAGS) -Ofast -lgtest -lflint test/mulmod_test.cpp -o bin/mulmod_test
	$(CXX) $(CXXFLAGS) -Ofast -lgtest -lflint test/jacobi_test.cpp -o bin/jacobi_test


.PHONY: all test clean bench brocard
