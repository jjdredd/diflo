WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ./include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)
INPUTLIB=input/data.o input/event.o input/particle.o
ANALYS=analysis/dist.o analysis/fft.o

diflo.elf: diflo.cpp $(INPUTLIB) $(ANALYS)
	g++ $(CXXFLAGS) diflo.cpp $(INPUTLIB) $(ANALYS) -o diflo.elf

$(INPUTLIB):
	$(MAKE) -C input

$(ANALYS):
	$(MAKE) -C analysis

clean:
	$(MAKE) -C input clean
	$(MAKE) -C analysis clean
	rm diflo.elf
