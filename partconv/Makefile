CXX=g++

FLAGS= -Wall -g -std=c++11 -fopenmp
INCLUDE= -I./patoh 
LIB= -L./patoh

OBJ= graph.o

all: part evalord

part: $(OBJ) part.cpp  graph.hpp
	$(CXX) $(FLAGS) $(INCLUDE) $(LIB) -o $@ $(OBJ) $@.cpp -lpatoh -lm 
evalord: $(OBJ) evalord.cpp graph.hpp 
	$(CXX) $(FLAGS) $(INCLUDE) $(LIB) -o $@ $(OBJ) $@.cpp -lpatoh -lm

%.o: %.cpp 
	$(CXX) $(FLAGS) $(INCLUDE) -c $< 

.PHONY: clean

clean:
	@rm -rf *.o part evalord *~ 
