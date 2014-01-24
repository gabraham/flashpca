
.PHONY: all

EIGEN=/usr/include/eigen3

all: flashpca2

OBJ = \
   randompca.o \
   flashpca.o \
   data.o \
   util.o

CXXFLAGS = -Iboost -IEigen -I/usr/include/eigen3

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -msse2
else
   CXXFLAGS += -march=native
endif

BOOST = -lboost_system-mt \
   -lboost_iostreams-mt \
   -lboost_filesystem-mt \
   -lboost_program_options
 
debug: LDFLAGS = $(BOOST)
debug: CXXFLAGS += -O0 -ggdb3
debug: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca2 $^ $(LDFLAGS)

flashpca2: LDFLAGS = $(BOOST)
flashpca2: CXXFLAGS += -g -O3 -DNDEBUG \
   -funroll-loops -ftree-vectorize -ffast-math -fopenmp
flashpca2: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca2 $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) flashpca

