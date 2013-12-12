
.PHONY: all

all: flashpca

OBJ = \
   randompca.o \
   flashpca.o \
   data.o \
   util.o \
   cache.o

CXXFLAGS = -Iboost -IEigen

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -msse2
else
   CXXFLAGS += -march=native
endif

BOOST = -lboost_system-mt \
   -lboost_iostreams-mt \
   -lboost_filesystem-mt
 
debug: LDFLAGS = $(BOOST)
debug: CXXFLAGS += -O0 -ggdb3
debug: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

flashpca: LDFLAGS = $(BOOST)
flashpca: CXXFLAGS += -O3 -DNDEBUG -funroll-loops -ftree-vectorize
flashpca: $(OBJ)
	$(CXX) $(CXXFLAGS) -o flashpca $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) flashpca

