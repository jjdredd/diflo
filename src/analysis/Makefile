WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ../include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)

dist.o: dist.cpp
	g++ $(CXXFLAGS) -c -o dist.o dist.cpp

fft.o: fft.cpp
	g++ $(CXXFLAGS) -c -o fft.o fft.cpp

clean:
	rm dist.o
