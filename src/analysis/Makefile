WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ../include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)

# just organize everything into an .a lib to avoid mess

all: dist.o fft.o momentum.o

dist.o: dist.cpp
	g++ $(CXXFLAGS) -c -o dist.o dist.cpp

fft.o: fft.cpp
	g++ $(CXXFLAGS) -c -o fft.o fft.cpp

momentum.o: momentum.cpp
	g++ $(CXXFLAGS) -c -o momentum.o momentum.cpp

clean:
	rm dist.o fft.o momentum.o
