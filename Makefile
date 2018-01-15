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
exist=1

all: $(TARGET)
#	$(BIN_EXCUTE)

$(TARGET): $(OBJ)

ifneq (exist, $(shell [ ! -d "obj" ]))
	echo "Directory 'obj' doesn't exist"
	$(shell mkdir obj)
endif

	mv *.o obj/
	echo $(OBJ_DIR)
	$(CXX) -I$(INC) $(CFLAGS) $(OBJ_DIR) -o $(BIN_EXCUTE)

main.o: 
	$(CXX) -I$(INC) $(CFLAGS) -c src/main.cpp src/graph.cpp 

graph.o:
	$(CXX) -I$(INC) $(CFLAGS) -c src/graph.cpp

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
	rm -rf *.o
