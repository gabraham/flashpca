
.PHONY: all

VERSION=1.3

EIGEN_INC=/usr/local/include/eigen
BOOST_INC=/usr/local/include/boost
BOOST_LIB=/usr/local/lib
SPECTRA_INC=spectra

all: flashpca predict
static: flashpca_x86-64

OBJ = \
   randompca.o \
   flashpca.o \
   data.o \
   util.o

OBJ2 = \
   predict.o \
   data.o \
   util.o

CXXFLAGS = -I${SPECTRA_INC} -I${BOOST_INC} -I${EIGEN_INC}

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -msse2 -DEIGEN_DONT_PARALLELIZE
else
   CXXFLAGS += -march=native -fopenmp -std=c++0x
endif

BOOST = -L${BOOST_LIB} \
   -lboost_system \
   -lboost_iostreams \
   -lboost_filesystem \
   -lboost_program_options

debug: LDFLAGS = $(BOOST)
debug: CXXFLAGS += -O0 -ggdb3 -DVERSION=\"$(VERSION)\"
debug: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

flashpca: LDFLAGS = $(BOOST)
flashpca: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math
flashpca: flashpca.o randompca.o data.o util.o
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

predict: LDFLAGS = $(BOOST)
predict: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math
predict: $(OBJ2)
	$(CXX) $(CXXFLAGS) -o predict $^ $(LDFLAGS)

flashpca_x86-64: LDFLAGS = $(BOOST) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
flashpca_x86-64: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math -static
flashpca_x86-64: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca_x86-64 $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ2): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OBJ2) flashpca flashpca_x86-64


