CC := g++
CXXFLAGS += -std=c++26 -fmodules-ts

run: all
	./test_rt

all: test_rt

test_rt: test_rt.o ts.o rt.o

test_rt.o: test_rt.cpp ts.o rt.o

ts.o: ts.cpp

rt.o: rt.cpp


