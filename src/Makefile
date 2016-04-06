# define compiler
CC = g++-5
# define compile time flags
CFLAGS = -g3 -Wall -g -DDEBUG
#CFLAGS = -fopenmp -O3 -std=c++11
#CFLAGS = -fopenmp -O3 
# define any directories over than /usr/include/ containing header files
CINCLUDES = -I/usr/local/include -I . -I/usr/local/include/yaml-cpp
# define libraries paths over than /usr/include
LFLAGS = -lm 
FFTFLAGS = -lfftw3 -lfftw3_omp 
YAMLFLAGS = -lyaml-cpp
OMPFLAGS = -fopenmp
# define cpp source files (currently not using)
# SRCS = diffuseFT1D.cpp fftw++.cc (currently not using)
# make objects with macro (currently not using)
# OBJS = $(SRCS:.cpp = .o)
# OBJs += $(SRCS:.cc = .o)
CXX = $(CC)
CXXFLAGS = $(CFLAGS)
#CXXFLAGS += $(CINCLUDES)

OBJS = HardRodMain.o HrAniBody.o HrIsoBody.o fftw++.o spGrid.o tGrid.o Propagator.o Ml3dPrint2Screen.o 
OBJS +=	ParseParams.o MayerFncHardRod.o EqDist.o Rho3DMaker.o NlHardRodDr.o CheckBrokenDen.o
OBJS += OPclass.o HRwriter.o
#OBJS += OPclass.o

MAIN = ../diffuse.exe

.PHONY: all

all: $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CINCLUDES) $(CFLAGS) $(OMPFLAGS) $(OBJS) -o $@ $(LFLAGS) $(FFTFLAGS) $(YAMLFLAGS) 
	 
#$(MAIN): $(OBJS) 
#	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(cincludes)

.PHONY: clean
clean:
	rm *.o

# $@: target $^: all dependenies $<: first dependency
