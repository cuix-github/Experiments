CXXFLAGS=-g -Wall -DNDEBUG -fopenmp -ffast-math
CXX=g++-4.6

all: smoke
clean:
	rm -f *.o smoke

main.o: main.cpp smoke.h
smoke.o: smoke.cpp smoke.h

smoke: main.o smoke.o
	$(CXX) -fopenmp main.o smoke.o -o smoke -lglut -lpng -lGLU -lGL
