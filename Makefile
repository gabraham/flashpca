
.PHONY: all

VERSION=2.0

EIGEN_INC=/usr/local/include/eigen
BOOST_INC=/usr/local/include/boost
BOOST_LIB=/usr/local/lib
SPECTRA_INC=spectra/include

all: flashpca
static: flashpca_x86-64

OBJ = \
   randompca.o \
   flashpca.o \
   svdwide.o \
   svdtall.o \
   data.o \
   util.o

CXXFLAGS += -I${SPECTRA_INC} -I${BOOST_INC} -I${EIGEN_INC}

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -msse2 -DEIGEN_DONT_PARALLELIZE -std=c++11
   BOOST = ${BOOST_LIB}/libboost_program_options.a
else
   CXXFLAGS += -march=native -fopenmp -std=c++0x
   BOOST = -L${BOOST_LIB} -lboost_program_options
endif


debug: LDFLAGS = $(BOOST)
debug: CXXFLAGS += -O0 -ggdb3 -DVERSION=\"$(VERSION)\"
debug: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

flashpca: LDFLAGS = $(BOOST)
flashpca: CXXFLAGS += -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math
flashpca: flashpca.o randompca.o data.o util.o svdwide.o svdtall.o
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

flashpca_x86-64: LDFLAGS = $(BOOST) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
flashpca_x86-64: CXXFLAGS += -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math -static
flashpca_x86-64: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca_x86-64 $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) flashpca flashpca_x86-64
