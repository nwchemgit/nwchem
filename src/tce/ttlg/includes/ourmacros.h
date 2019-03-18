#define EXPANDDIMS(funcname, param3, param4, param5, callstring) \
case 0: \
 CON(funcname,0) <<<param3, param4, param5>>> callstring;\
break;\
case 1: \
 CON(funcname,1) <<<param3, param4, param5>>> callstring;\
break;\
case 2:\
 CON(funcname,2) <<<param3, param4, param5>>> callstring;\
break;\
case 3:\
 CON(funcname,3) <<<param3, param4, param5>>> callstring;\
break;\
case 4:\
 CON(funcname,4) <<<param3, param4, param5>>> callstring;\
break;\
case 5:\
 CON(funcname,5) <<<param3, param4, param5>>> callstring;\
break;\
case 6:\
 CON(funcname,6) <<<param3, param4, param5>>> callstring;\
break;\
case 7:\
 CON(funcname,7) <<<param3, param4, param5>>> callstring;\
break;\
case 8:\
 CON(funcname,8) <<<param3, param4, param5>>> callstring;\
break;\
case 9:\
 CON(funcname,9) <<<param3, param4, param5>>> callstring;\
break;\
case 10:\
 CON(funcname,10) <<<param3, param4, param5>>> callstring;\
break;\
case 11:\
 CON(funcname,11) <<<param3, param4, param5>>> callstring;\
break;\
case 12:\
 CON(funcname,12) <<<param3, param4, param5>>> callstring;\
break;\
case 13:\
 CON(funcname,13) <<<param3, param4, param5>>> callstring;\
break;\
case 14:\
 CON(funcname,14) <<<param3, param4, param5>>> callstring;\
break;\
case 15:\
 CON(funcname,15) <<<param3, param4, param5>>> callstring;\
break;

#ifdef TCE_CUDA
#define SAFECUDAMALLOC(pointer, size) cudaMalloc(pointer, size) ;\
        {cudaError_t err = cudaGetLastError();\
                if(err != cudaSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}

#define SAFECUDAMEMCPY(dest, src, size, type) cudaMemcpy(dest, src, size, type) ;\
        {cudaError_t err = cudaGetLastError();\
                if(err != cudaSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}
#endif

#ifdef TCE_HIP
#define SAFECUDAMALLOC(pointer, size) hipMalloc(pointer, size) ;\
        {hipError_t err = hipGetLastError();\
                if(err != hipSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", hipGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}

#define SAFECUDAMEMCPY(dest, src, size, type) hipMemcpy(dest, src, size, type) ;\
        {hipError_t err = hipGetLastError();\
                if(err != hipSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", hipGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}
#endif
