WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ../include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)

inputlib: data.o event.o particle.o

data.o: data.cpp
	g++ $(CXXFLAGS) -c -o data.o data.cpp

event.o: event.cpp
	g++ $(CXXFLAGS) -c -o event.o event.cpp

particle.o: particle.cpp
	g++ $(CXXFLAGS) -c -o particle.o particle.cpp

clean:
	rm event.o particle.o data.o
