WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ./include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)
LDFLAGS=-lgsl -lm
INPUT_F=input/data.o
ANALYS_F=analysis/fft.o analysis/dist.o analysis/momentum.o

all: diflo.elf diflo.semi-static.elf

diflo.semi-static.elf: diflo.cpp inputlib analys
	g++ $(CXXFLAGS) diflo.cpp \
	$(INPUT_F) $(ANALYS_F) \
	-o diflo.semi-static.elf /usr/lib/libgsl.a

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
	rm diflo.elf diflo.semi-static.elf
