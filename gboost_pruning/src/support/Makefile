CXX = g++
VERSION = 0.01
NAME = gspan
CXXFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
#CXXFLAGS = -std=c++11 -O2 -Wall -Wno-deprecated
LDFLAGS =  -lm -ldl #remove "-static" for Mac,in linux need to compile?
TARGETS = lpboost_p_s
OBJ = ../gspan.o ../cachetree.o ../lpboost.o support.o

all:  lpboost

lpboost: ${OBJ}
	${CXX} ${CXXFLAGS} -o ${TARGETS} ${OBJ} ${LDFLAGS} -I ../../glpk-4.8/include ../../glpk-4.8/src/libglpk.a

clean:
	rm -f ${OBJ} lpboost.o ${TARGETS} *~

