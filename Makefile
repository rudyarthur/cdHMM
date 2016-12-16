CXX=g++

CXXFLAGS = -O3  -std=c++11
CFLAGS = -O3  -std=c++11

all: test

IFLAGS = -I./include

default: CXXFLAGS = -O3  -std=c++11
default: CFLAGS = -O3  -std=c++11
default: all

debug: CXXFLAGS = -g -std=c++11
debug: CFLAGS =  -g -std=c++11
debug: all

profile: CXXFLAGS = -pg -O3 -std=c++11
profile: CFLAGS =  -pg -O3 -std=c++11
profile: all

test: test.cpp 
	#$(CXX) $(CXXFLAGS) -o test test.cpp $(IFLAGS)
	$(CXX) $(CXXFLAGS) -o text_analysis text_analysis.cpp $(IFLAGS)
clean:
	rm test 
	rm text_analysis


