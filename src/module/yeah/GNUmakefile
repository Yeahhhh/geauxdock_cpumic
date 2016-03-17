
SRC_C = $(wildcard ./c/*.c)
SRC_CPP = $(wildcard ./cpp/*.cpp)
SRC_OPENCL += ./opencl/$(wildcard *.c)

OBJ_C := $(SRC_C:.c=.o)
OBJ_CPP := $(SRC_CPP:.cpp=.o)
OBJ_OPENCL := $(SRC_OPENCL:.c=.o)




#HEADPATH += -Ix
#LIBPATH += -Lxx
#LINKFLAGS += -lxx





#CC := icc #-mmic
#CXX := icpc #-mmic
CC := gcc
CXX := g++
#CC := clang
#CXX := clang++

FLAGS := -O3 -fPIC





#all: libc libcpp
all: libc


%.o: %.c
	$(CC) $(FLAGS) $(HEADPATH) -c $< -o $@

%.o: %.cpp
	$(CXX) $(FLAGS) $(HEADPATH) -c $< -o $@

libc: $(OBJ_C)
	ar rcs ../lib/libyeahc.a $(OBJ_C)

libcpp: $(OBJ_C) $(OBJ_CPP)
	ar rcs ../lib/libyeahcpp.a $(OBJ_C) $(OBJ_CPP)

libopencl: $(OBJ_OPENCL)
	ar rcs ../lib/libyeahopencl.a $(OBJ_OPENCL)


echo:
	@echo $(OBJ_C)
	@echo $(OBJ_CPP)

clean:
	$(RM) $(OBJ_C) $(OBJ_CPP)



