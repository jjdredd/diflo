data_all: data.o root_data.o

data.o: data.cpp
	g++ $(CXXFLAGS) -c -o data.o data.cpp

root_data.o: root_data.cpp
	g++ $(CXXFLAGS) -c -o root_data.o root_data.cpp

clean:
	rm data.o root_data.o
