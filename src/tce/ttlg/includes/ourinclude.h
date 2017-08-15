#ifndef type
#define type double
#endif
extern "C" void  
fvimatchl32_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm, type alpha, type beta);
extern "C" void
fvimatchg32_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm, type alpha, type beta);
extern "C" void
fvimatchg32_blocking_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm, type alpha, type beta);
//extern "C" void
//fvinomatchg32_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm);
extern "C" void
fvinomatchgeneral_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm, const int* perm, type alpha, type beta);
extern "C" void
fvigeneralolap_transpose_kernel(int ndim, const type *A, type *B, const int *lda, const int *ldb, const int* params, const int * rperm, const int* perm, type alpha, type beta);


extern "C"
double getEfficiency_nooverlap(int ilimit, int olimit, int asize, int bsize, int blockA, int blockB);
extern "C"
double getEfficiency_overlap(int ilimit, int olimit, int asize, int bsize, int blockA, int blockB);
extern "C"
double getEfficiency_nomatchg32(int ilimit, int olimit);
extern "C"
double getEfficiency_matchl32(int size0, int inp1, int out1, int block);


extern "C"
double getBW_nooverlap(double eff);
extern "C"
double getBW_overlap(double eff);
extern "C"
double getBW_nomatchg32(double eff);
extern "C"
double getBW_matchg32();
extern "C"
double getBW_matchl32(double eff, int tbsize);


extern "C"
double getTime(double, unsigned long);
#define SAFECUDAMALLOC(callstring) cudamalloc(callstring) ;\
	{cudaError_t err = cudaGetLastError();\
                if(err != cudaSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}

#define SAFECUDAMEMCPY(callstring) cudamemcpy(callstring) ;\
        {cudaError_t err = cudaGetLastError();\
                if(err != cudaSuccess){\
                        printf("\nKernel ERROR in hostCall: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);\
                        exit(-1);\
                }}


void transpose_check(int n, const type *A, type *B, const type alpha, const type beta,  const int *lda,  int *perm);
//template <typename tensortype>
//void ourproj_tranpose(int ndim, int *dims, int *perm, tensortype *input, tensortype *output, TensorType dataType = 1, tensortype alpha = 0.0, tensortype beta = 1.0);
extern "C"
void ttlg_transpose(int ndim, int *dims, int *perm, double *input, double *output, double alpha, double beta);

void transpose_init(type *A, int size);

int transpose_equal(const type *A, const type*B, int total_size);
