# define compiler
CC = g++-5
# define compile time flags
CFLAGS = -g3 -Wall -g -DDEBUG
CFLAGS += -fopenmp
CINCLUDES = -I/usr/local/include -I . 

# Libraries 
LFLAGS = -lm 
LFFTFLAGS = -lfftw3 -lfftw3_omp 
LYAMLFLAGS = -lyaml-cpp
CXX = $(CC)
CXXFLAGS = $(CFLAGS)

OBJS = HardRodMain.o HrAniBody.o HrIsoBody.o fftw++.o spGrid.o tGrid.o Propagator.o Ml3dPrint2Screen.o 
OBJS +=	ParseParams.o MayerFncHardRod.o EqDist.o Rho3DMaker.o NlHardRodDr.o CheckBrokenDen.o
OBJS += OPclass.o HRwriter.o

MAIN = ../diffuse.exe

.PHONY: all

all: $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CINCLUDES) $(CFLAGS) $(OBJS) -o $@ $(LFLAGS) $(LFFTFLAGS) $(LYAMLFLAGS) 
	 

.PHONY: clean
clean:
	rm *.o

# $@: target $^: all dependenies $<: first dependency
