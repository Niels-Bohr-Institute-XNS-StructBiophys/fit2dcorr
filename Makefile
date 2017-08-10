# GNU Make v3.81 supports else ifeq, v3.80 not !
#
# Optimization

# gcc c(++) compilers
ifeq ($(OLEVEL),DEBUG)
	OPT:=
else
	OPT:=-O$(OLEVEL)
endif

ifeq ($(CC),gcc)
	OPENMP=-fopenmp
	ifeq ($(OOPTIONS),AVX)
		OPT:=$(OPT) -march=corei7-avx -mtune=corei7-avx
	endif
endif
ifeq ($(CC),g++)
	OPENMP=-fopenmp
	ifeq ($(OOPTIONS),AVX)
		OPT:=$(OPT) -march=corei7-avx -mtune=corei7-avx
	endif
endif

# intel c(++) compilers
ifeq ($(CC),icc)
	OPENMP=-openmp
	ifeq ($(OOPTIONS),-xP)
		OPT:=$(OPT) -xP
	endif
endif
ifeq ($(CC),icpc)
	OPENMP=-openmp
	ifeq ($(OOPTIONS),-xP)
		OPT:=$(OPT) -xP
	endif
endif

# compiler flags
CFLAGS:=-Wall $(OPENMP)
ifeq ($(OLEVEL),DEBUG)
	CFLAGS:=-g3 $(CFLAGS)
endif

# libraries
LFLAGS:=-lm -ltiff
LFLAGS:=$(LFLAGS) -lc

all: fit2dcorr
fit2dcorr: fit2dcorr.o
	@echo $(CFLAGS)
	@echo $(LFLAGS)
	$(CC) $(OPT) $(CFLAGS) -o fit2dcorr fit2dcorr.o $(LFLAGS)
fit2dcorr.o: fit2dcorr.cpp
	$(CC) $(OPT) $(CFLAGS) -c fit2dcorr.cpp
clean:
	rm -f *.o fit2dcorr
