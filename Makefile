.SUFFIXES : .cpp .o
CC=gcc
CXX=g++
INC = ./inc
CFLAGS= -g -std=c++14 -O2
EXCUTABLE=bin/graph

#SRC = $(wildcard *.cpp)
OBJ = main.o graph.o
SRC = $(OBJ:.o=.cpp)
OBJ_DIR = $(addprefix obj/,$(OBJ))
SRC_DIR = $(addprefix src/,$(SRC))
BIN_EXCUTE =bin/graph
TARGET=graph
MOVE_OBJ=mv $@ obj/


all: $(TARGET)
	$(BIN_EXCUTE)
	
$(TARGET): $(OBJ)
	echo $(OBJ_DIR)
	$(CXX) -I$(INC) $(CFLAGS) $(OBJ_DIR) -o $(BIN_EXCUTE)

main.o: src/main.cpp inc/graph.h
	$(CXX) -I$(INC) $(CFLAGS) -c src/main.cpp src/graph.cpp 
	$(MOVE_OBJ)

graph.o: src/graph.cpp inc/graph.h
	$(CXX) -I$(INC) $(CFLAGS) -c src/graph.cpp
	$(MOVE_OBJ)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
