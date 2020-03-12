MAKEFLAGS +=" -j1 -l0.001"

#ifdef TCE_CUDA # not really needed, but keep the guard for now.
 OBJ_OPTIMIZE = fvimatchl32.o fvimatchg32_blocking.o \
	fvimatchg32.o      fvinomatchgeneral.o fvigeneralolap.o \
	BlockingCase.o  main.o  Parameters.o  ParameterTuner.o  TransposeSpec.o model.o

 #LIB_DEFINES = -I./includes $(CUDA_INCLUDE)

LIB_INCLUDES = -I./includes

ifdef TCE_CUDA
LIB_INCLUDES += $(CUDA_INCLUDE)
endif

 #CUDA_FLAGS = -I./includes $(CUDA_INCLUDE)

LIBRARY = libtce.a

#ARCH =  -arch=sm_35 
#CUDA_FLAGS += -O3 -Xcompiler -fPIC -std=c++11 -DNOHTIME -Xptxas --warn-on-spills $(ARCH)

USES_BLAS = 

HEADERS =

include ../../config/makefile.h
include ../../config/makelib.h

#endif
