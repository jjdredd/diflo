WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/../include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm -lgsl -lgslcblas
# for the below libs we should make an .a lib with ar ?
INPUT_F = ../input/data.o
ANALYS_F = ../analysis/fft.o ../analysis/dist.o ../analysis/momentum.o \
../analysis/handedness.o

all: rp_test.elf

rp_test.elf: rp_test.cpp inputlib analys
	g++ $(CXXFLAGS) $(LDFLAGS) rp_test.cpp \
	$(INPUT_F) $(ANALYS_F) \
	-o rp_test.elf

inputlib:
	$(MAKE) -C ../input

analys:
	$(MAKE) -C ../analysis

clean:
	$(MAKE) -C ../input clean
	$(MAKE) -C ../analysis clean
	rm rp_test.elf
