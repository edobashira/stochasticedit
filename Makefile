CXX=g++
CXXFLAGS= -Wno-sign-compare -g -O2

LDFLAGS= -lm -ldl -lfst

all: edit 


edit: edit.o
	$(CXX) $^ -o $@  $(LDFLAGS)

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o *~ *.so edit
  
  
