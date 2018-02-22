CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3 # -pg
LIBS=-lm -lfftw3
INC=-I /home/username/bin/eigen-3.2.8
# INC=-I/usr/local/lib/ is not explicitly added

MAIN_SRC=main.cc
MAIN_OBJ=$(MAIN_SRC:.cc=.o)
SRC=$(wildcard */*.cc)
OBJ=$(SRC:.cc=.o)
TARGET=../bin/plum

all: $(TARGET)

$(TARGET): $(MAIN_OBJ) $(OBJ)
	$(CXX) -o $@ $(MAIN_OBJ) $(OBJ) $(CXXFLAGS) $(LIBS)
$(MAIN_OBJ): $(MAIN_SRC)
	$(CXX) $(INC) -c $< -o $@ $(CXXFLAGS)
%.o: %.cc
	$(CXX) $(INC) -c $< -o $@ $(CXXFLAGS)
clean: 
	-rm $(MAIN_OBJ) $(OBJ)
run:
	bin/oops < in

