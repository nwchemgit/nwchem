#TODO: Double check the GPU part
INCLUDE(FindCUDA)

if(CUDA_TOOLKIT_INCLUDE) 
    list(APPEND TCE_INCLUDES ${CUDA_TOOLKIT_INCLUDE})
else()
    message(WARNING "CUDA_TOOLKIT_INCLUDE not set.")
endif()

set(TTLG_CUDA_SOURCE_FILES 
fvimatchl32.cu 
fvimatchg32.cu 
fvimatchg32_blocking.cu 
fvinomatchgeneral.cu
fvigeneralolap.cu)


set(TTLG_CPP_SOURCE_FILES 
 model.c
 main.cpp
 Parameters.cpp  
 BlockingCase.cpp 
 ParameterTuner.cpp 
 TransposeSpec.cpp)

set_source_files_properties(${TTLG_CPP_SOURCE_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

cuda_add_library(ttlg STATIC ${TTLG_CPP_SOURCE_FILES} ${TTLG_CUDA_SOURCE_FILES})
target_include_directories(ttlg PUBLIC ${PROJECT_SOURCE_DIR}/src/tce/ttlg/includes)
target_compile_options(ttlg PRIVATE -O3 -Xcompiler -fPIC -std=c++11 -DNOHTIME -Xptxas --warn-on-spills -arch=sm_35)
