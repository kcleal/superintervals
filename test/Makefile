CXXFLAGS=	-g -Wall -O3 -std=c++17 -march=native
#CXXFLAGS=	-g -Wall -O1 -std=c++17 -march=native
#CFLAGS=		-g -Wall -O2 -Wc++-compat
LIBS=		-lz
EXE=		run-tests run-cpp-libs  #bedcov
INCLUDE= -I../matrylist -I./3rd-party
ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

all:$(EXE)

#bedcov: bedcov-siList.cpp
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o siList

run-tests: tests.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o run-tests

run-cpp-libs: bench.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o $@

clean:
		rm -fr *.o a.out *.dSYM $(EXE)
