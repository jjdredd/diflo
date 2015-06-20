WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS)
export CXXFLAGS = $(CFLAGS)
export LDFLAGS = -lm
INPUT_F = input/data.o
ANALYS_F = analysis/fft.o analysis/dist.o analysis/momentum.o

all: diflo.elf diflo.server.elf64

diflo.server.elf64: CFLAGS += -msse4 -I/home/ltph/usubov/usr/include
diflo.server.elf64: CXXFLAGS += -I/home/ltph/usubov/usr/include
diflo.server.elf64: LDFLAGS += -L/home/ltph/usubov/usr/lib -lgsl -lgslcblas
diflo.server.elf64: diflo.cpp inputlib analys
	g++ $(CXXFLAGS) $(LDFLAGS) diflo.cpp \
	$(INPUT_F) $(ANALYS_F) \
	-o diflo.server.elf64

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
	rm diflo.elf diflo.server.elf64
