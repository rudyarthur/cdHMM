CXX=g++

CXXFLAGS = -O3  -std=c++11
CFLAGS = -O3  -std=c++11

all: cline

default: CXXFLAGS = -O3  -std=c++11
default: CFLAGS = -O3  -std=c++11
default: all

debug: CXXFLAGS = -g -std=c++11
debug: CFLAGS =  -g -std=c++11
debug: all

profile: CXXFLAGS = -pg -O3 -std=c++11
profile: CFLAGS =  -pg -O3 -std=c++11
profile: all

cline: hmm.cpp
	$(CXX) $(CXXFLAGS) -o hmm hmm.cpp

test: test.cpp 
	$(CXX) $(CXXFLAGS) -o test test.cpp
clean:
	rm hmm
	rm test 


