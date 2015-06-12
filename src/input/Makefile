WFLAGS=-W -Wall -Wextra
IFLAGS=-iquote ../include
OPTFLAGS=-march=native -ffast-math -O3 -msse4 -funroll-loops
CFLAGS=$(IFLAGS) $(WFLAGS) $(OPTFLAGS)
CXXFLAGS=$(CFLAGS)

data.o: data.cpp
	g++ $(CXXFLAGS) -c -o data.o data.cpp

clean:
	rm data.o
