
# just organize everything into an .a lib to avoid mess ???

dist.o: dist.cpp
	g++ $(CXXFLAGS) -c -o dist.o dist.cpp

fft.o: fft.cpp
	g++ $(CXXFLAGS) -c -o fft.o fft.cpp

momentum.o: momentum.cpp
	g++ $(CXXFLAGS) -c -o momentum.o momentum.cpp

handedness.o: handedness.cpp
	g++ $(CXXFLAGS) -c -o handedness.o handedness.cpp

clean:
	rm dist.o fft.o momentum.o handedness.o
