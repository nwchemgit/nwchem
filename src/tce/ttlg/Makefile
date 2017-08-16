IFLAGS =  -Iincludes/
ARCH =  -arch=sm_35 
#ARCH = -arch=sm_35
#NVFLAGS = -O3 -std=c++11  -DNOHTIME -arch=sm_35
#NVFLAGS = -O3 -std=c++11 $(IFLAGS)  $(CFLAGS) -DNOHTIME -Xptxas --warn-on-spills -Xptxas -dlcm=ca  $(ARCH)
NVFLAGS = -O3 -Xcompiler -fPIC -std=c++11 $(IFLAGS)  $(CFLAGS) -DNOHTIME -Xptxas --warn-on-spills   $(ARCH)
#NVFLAGS = -O3 -std=c++11  -DNOHTIME  -Xptxas --warn-on-spills -Xptxas -dlcm=cg
#NVFLAGS = -O3 -std=c++11  -DNOHTIME  -Xptxas --warn-on-spills -Xptxas -dlcm=ca -arch=sm_35
#NVFLAGS = -O3 -std=c++11  -DNOHTIME  -Xptxas --warn-on-spills -Xptxas -flcm=ca -arch=sm_35

HEADERS = kernels/fvigeneralolap.h
SCR1 =  fvimatchl32.cu fvimatchg32_blocking.cu \
	fvimatchg32.cu      fvinomatchgeneral.cu fvigeneralolap.cu 
#SCR1 =  matrixtranspose.cu fvimatchl32.cu fvimatchg32l120.cu        fvinomatchalldiff.cu  fvinomatchonerepeat.cu  fvinomatchsisame.cu \
	fvimatchg32.cu  fvimatchl32_factor.cu  fvinomatchg32.cu      fvinomatchrepdiff.cu fvigeneral.cu fvigeneralolap.cu 

SCR2 = test.cpp BlockingCase.cpp  main.cpp  Parameters.cpp  ParameterTuner.cpp  TransposeSpec.cpp model.c
OBJS = $(SCR1:.cu=.o) $(SCR2:.cpp=.o) 

MAIN = ttlg
MAIN_L = libttlg.so

LIBS = -lgomp -lm
all: $(MAIN) 
lib: $(MAIN_L)
exe: $(MAIN)
$(MAIN_L) : $(OBJS)
	nvcc -shared $(NVFLAGS) $^ -o $(MAIN_L) $(LIBS)
$(MAIN) : $(OBJS)
	nvcc $(NVFLAGS) $^ -o $(MAIN) $(LIBS)
    
%.o : %.cpp
	nvcc $(NVFLAGS) -c $< -o $@

%.o : %.cu
	nvcc $(NVFLAGS) -c $< -o $@
fvigeneralolap.cu : kernels/fvigeneralolap.h


clean:
	rm -f *.o $(MAIN) $(MAIN_L)

