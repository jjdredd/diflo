WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm
# for the below libs we should make an .a lib with ar ?
INPUT_F = input/data.o
ANALYS_F = analysis/fft.o analysis/dist.o analysis/momentum.o \
analysis/handedness.o

all: diflo.elf sdnopt.elf multip.elf

sdnopt.elf: CFLAGS += -I/home/ltph/usubov/usr/include
sdnopt.elf: CXXFLAGS += -I/home/ltph/usubov/usr/include
sdnopt.elf: LDFLAGS += -L/home/ltph/usubov/usr/lib -lgsl -lgslcblas
sdnopt.elf: sdnopt.cpp inputlib analys
	g++ $(CXXFLAGS) $(LDFLAGS) sdnopt.cpp \
	$(INPUT_F) $(ANALYS_F) \
	-o sdnopt.elf

multip.elf: multip.cpp $(INPUT_F)
	g++ $(CXXFLAGS) $(LDFLAGS) multip.cpp \
	$(INPUT_F) -o multip.elf

diflo.elf: CXXFLAGS += -I/home/ltph/usubov/usr/include
diflo.elf: LDFLAGS += -L/home/ltph/usubov/usr/lib -lgsl -lgslcblas
diflo.elf: LDFLAGS += -L/home/ltph/usubov/usr/lib -lgsl -lgslcblas
diflo.elf: diflo.cpp inputlib analys
	g++ $(CXXFLAGS) $(LDFLAGS) diflo.cpp \
	$(INPUT_F) $(ANALYS_F) \
	-o diflo.elf

inputlib:
	$(MAKE) -C input

analys:
	$(MAKE) -C analysis

clean:
	$(MAKE) -C input clean
	$(MAKE) -C analysis clean
	rm diflo.elf sdnopt.elf multip.elf
