CXXFLAGS+=	-O3 -std=c++17 -march=native -mtune=native -fomit-frame-pointer -fno-exceptions -fno-rtti #-DSI_NOSIMD
LIBS=
EXE=		run-tests run-cpp-libs
INCLUDE= -I../src -I./3rd-party
# GCC use -ftree-vectorize -fopt-info-vec-optimized

# to build coitrees use:
# cd 3rd-party/coitrees; RUSTFLAGS="-Ctarget-cpu=native" cargo run --release --example bed-intersect


ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

all:$(EXE)

build_clibs:
	$(CC) $(CFLAGS) -c ./3rd-party/intervaldb.c -o intervaldb.o
	$(CC) $(CFLAGS) -c ./3rd-party/cgranges.c -o cgranges.o

run-tests: tests.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) -o run-tests

run-cpp-libs: bench.cpp build_clibs
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< cgranges.o intervaldb.o $(LIBS) -o $@

clean:
	rm -fr *.o a.out *.dSYM $(EXE)
