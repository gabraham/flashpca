
.PHONY: all

VERSION=1.2

EIGEN_INC=/usr/local/include/eigen
BOOST_INC=/usr/local/include/boost
BOOST_LIB=/usr/local/lib

all: flashpca predict
static: flashpca_x86-64

OBJ = \
   predict.o \
   randompca.o \
   flashpca.o \
   data.o \
   util.o

CXXFLAGS = -I${BOOST_INC} -I${EIGEN_INC}

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -msse2
else
   CXXFLAGS += -march=native
endif

BOOST = -L${BOOST_LIB} \
   -lboost_system \
   -lboost_iostreams \
   -lboost_filesystem \
   -lboost_program_options

debug: LDFLAGS = $(BOOST)
debug: CXXFLAGS += -O0 -ggdb3 -DVERSION=\"$(VERSION)\" -fopenmp
debug: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

flashpca: LDFLAGS = $(BOOST)
flashpca: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math -fopenmp
flashpca: flashpca.o randompca.o data.o util.o
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

predict: LDFLAGS = $(BOOST)
predict: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math -fopenmp
predict: predict.o data.o util.o
	$(CXX) $(CXXFLAGS) -o predict $^ $(LDFLAGS)

flashpca_x86-64: LDFLAGS = $(BOOST) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
flashpca_x86-64: CXXFLAGS += -g -O3 -DNDEBUG -DVERSION=\"$(VERSION)\" \
   -funroll-loops -ftree-vectorize -ffast-math -fopenmp -static
flashpca_x86-64: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca_x86-64 $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) flashpca flashpca_x86-64


