CXXFLAGS=-std=c++17 -march=native -mtune=native -fomit-frame-pointer -fno-stack-protector -fno-rtti -flto

all: brocard bench

brocard:
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -Ofast -lflint -fopenmp Brocard.cpp -o bin/brocard

bench:
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -O0 -lcelero -lflint bench/factorial_bench.cpp -o bin/factorial_bench
	$(CXX) $(CXXFLAGS) -Ofast -lcelero bench/count_trailing_zeros_bench.cpp -o bin/count_trailing_zeros_bench
	$(CXX) $(CXXFLAGS) -Ofast -lcelero -lflint bench/mulmod_bench.cpp -o bin/mulmod_bench

clean:
	rm -f bin/*


.PHONY: all test clean bench brocard
