CXX = g++
VERSION = 0.01
NAME = gspan
CXXFLAGS = -std=c++11 -O2 -Wall -Wno-deprecated
LDFLAGS =  -lm -ldl #remove "-static" for Mac,in linux need to compile?
OBJ = gspan.o cachetree.o lpboost.o

all:  base

base:	gspan.o cachetree.o lpboost.o

clean:
	rm -f ${OBJ} lpboost.o ${TARGETS} *~

