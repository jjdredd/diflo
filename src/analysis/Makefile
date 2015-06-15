
# just organize everything into an .a lib to avoid mess ???

all: dist.o fft.o momentum.o

dist.o: dist.cpp
	g++ $(CXXFLAGS) -c -o dist.o dist.cpp

fft.o: fft.cpp
	g++ $(CXXFLAGS) -c -o fft.o fft.cpp

momentum.o: momentum.cpp
	g++ $(CXXFLAGS) -c -o momentum.o momentum.cpp

clean:
	rm dist.o fft.o momentum.o
