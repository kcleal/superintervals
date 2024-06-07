CXXFLAGS=	-g -Wall -O3 -std=c++17 -march=native
LIBS=		-lz
EXE=		run-tests run-cpp-libs
INCLUDE= -I../matrylist -I./3rd-party
# GCC use -ftree-vectorize -fopt-info-vec-optimized
ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

all:$(EXE)

run-tests: tests.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o run-tests

run-cpp-libs: bench.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o $@

clean:
		rm -fr *.o a.out *.dSYM $(EXE)
