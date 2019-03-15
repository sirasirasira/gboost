CXX = g++
VERSION = 0.01
NAME = gspan
CXXFLAGS = -std=c++11 -pg -O2 -Wall -Wno-deprecated
LDFLAGS = -lm -ldl #remove "-static" for Mac,in linux need to compile?
TARGETS = lpboost_ban
OBJ = ../lpboost.o ../gspan.o ../cachetree.o bandit.o

all:  lpboost 

lpboost:  ${OBJ}
	${CXX} ${CXXFLAGS} -o ${TARGETS} ${OBJ} ${LDFLAGS} -I ../../glpk-4.8/include/ ../../glpk-4.8/src/libglpk.a

clean:
	rm -f ${TARGETS} *~ model model_reg
	rm -f -r itr

