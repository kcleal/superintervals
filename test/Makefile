CXXFLAGS=	-g -Wall -O3 -std=c++17 #-march=native
#CXXFLAGS=	-g -Wall -O1 -std=c++17 -march=native
#CFLAGS=		-g -Wall -O2 -Wc++-compat
LIBS=		-lz
EXE=		run-cpp-libs run_tests bedcov
INCLUDE= -I../sortedintersect -I./3rd-party
ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

all:$(EXE)

bedcov: bedcov-siList.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o siList

tests: test_sotedintersect.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o run_tests

run-cpp-libs: bench.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o $@

clean:
		rm -fr *.o a.out *.dSYM $(EXE)
