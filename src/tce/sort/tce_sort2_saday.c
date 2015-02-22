/* Index permutation routines in C for IA-32 processors
 * (interface compatible with IP routines in nwchem/tce)
 * Author: Qingda Lu (luq@cse.ohio-state.edu)
 * Time: 12/2006
 */


#include <stdio.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#if defined(__GNUC__)
static double buf[256] __attribute__ ((aligned(128)));
static double buf1[256] __attribute__ ((aligned(128)));
#elif defined(__INTEL_COMPILER) || defined(__ICC) 
__declspec(align(128)) static double buf[256],buf1[256];
#endif

#define linesize 64
#define tilesize 8

void copy(double* a,double* b,long size,double factor) {
  int i;
  for (i = 0; i < size; i++)
    b[i] = a[i]*factor;
}

void transpose_aligned(double *a, double *b, int N1, int N2, double factor) {

  int i,j,k,k1,it,jt,itt,jtt,conflict,tmp,tmpN;
  double *pA, *pB;


  register __m128d x, y, z, w,fac_vector; 

  fac_vector = _mm_load_sd(&factor);
  fac_vector = _mm_unpacklo_pd(fac_vector,fac_vector);

  for (it = 0; it < N1; it=it+tilesize) {      
    for (jt = 0; jt < N2; jt=jt+tilesize) {     

      k = 0;      
      for (j = jt; j < jt+tilesize; j=j+2) {
	for (i = it; i < it+tilesize; i=i+2) {
	  pA = a+i*N2+j;
	  x = _mm_load_pd(pA);    
	  y = _mm_load_pd(pA + N2);
	  x = _mm_mul_pd(x,fac_vector);  
	  y = _mm_mul_pd(y,fac_vector);
	  z = _mm_shuffle_pd( x, y, 0);
	  w = _mm_shuffle_pd( x, y, 3);
	  k = (j-jt)*tilesize + (i-it);
	  _mm_store_pd(buf + k,z);
	  _mm_store_pd(buf + k + tilesize,w);
	}
      }
        
      k = 0;
      k1 = 0;
      for (j = jt; j < jt+tilesize; j++) {
	pB = b+j*N1+it;
	k = (j-jt)*tilesize;
	x = _mm_load_pd(&buf[k]);
	y = _mm_load_pd(&buf[k]+2);
	z = _mm_load_pd(&buf[k]+2*2);
	w = _mm_load_pd(&buf[k]+3*2);
	_mm_stream_pd(pB,x);
	_mm_stream_pd(pB+2,y);
	_mm_stream_pd(pB+2*2,z);
	_mm_stream_pd(pB+3*2,w);
	
      } 
    }
  }
}


void transpose_misaligned(double *a, double *b, int N1, int N2, double factor) {

  int i,j,k,k1,it,jt,itt,jtt,it_bound,jt_bound,itt_bound,jtt_bound;
  int conflict,tmp,tmpN,offset,line_offset,setnum,set[8192/(4*sizeof(double))];
  double *pA, *pB;


  register __m128d x, y, z, w, t, t1,fac_vector; 
   
  fac_vector = _mm_load_sd(&factor);
  fac_vector = _mm_unpacklo_pd(fac_vector,fac_vector);

  itt_bound = (N1/tilesize)*tilesize;
  for (itt = 0; itt < itt_bound; itt=itt+5*tilesize) {
    jtt_bound =(N2/tilesize)*tilesize;
    for (jtt = 0; jtt < jtt_bound; jtt=jtt+5*tilesize) {
      it_bound = (itt+5*tilesize > itt_bound)?itt_bound:itt+5*tilesize; 
      for (it = itt; it < it_bound; it = it+tilesize) {
	jt_bound = (jtt+5*tilesize>itt_bound)?jtt_bound:jtt+5*tilesize;
	for (jt = jtt; jt < jt_bound; jt = jt+tilesize) { 
	  k = 0;      
	  for (j = jt; j < jt+tilesize; j=j+2) {
	    for (i = it; i < it+tilesize; i=i+2) {
	      pA = a+i*N2+j;
	      pB = b+j*N1+i;
	      x = _mm_loadu_pd(pA);
	      x = _mm_mul_pd(x,fac_vector);
	      y = _mm_loadu_pd(pA + N2);
	      y = _mm_mul_pd(y,fac_vector);
	      z = _mm_shuffle_pd( x, y, 0);
	      w = _mm_shuffle_pd( x, y, 3);
	      _mm_storeu_pd(pB,z);
	      _mm_storeu_pd(pB + N1,w);    
	    }
	  }
	}
      }
    }
    for (i = itt; i < itt+5*tilesize && i < itt_bound; i++) {
      for (j = jtt_bound; j < N2; j++) {
	b[j*N1+i] = factor * a[i*N2+j];
      }
    }
  }
  for (i = itt_bound; i < N1; i++) {
    for (j = 0; j < N2; j++) {
      b[j*N1+i] = factor * a[i*N2+j];
    }
  }
}

void tce_sort_2_(double* unsorted,double* sorted,
		 int* a_in, int* b_in, int* i_in, int* j_in, 
		 double* factor_in) {
  int a = *a_in;
  int b = *b_in;
  int i = *i_in;
  int j = *j_in;
  double factor = *factor_in;

  //check if trivial
  if (i==1 && j==2) {
    copy(unsorted,sorted,a*b,factor);
    return;
  }

  // check if aligned

  if (a <= 0 || b <= 0) {
    printf ("Dimensions must be larger than zero! \n");
  } else if (a%(linesize/sizeof(double))== 0 && b%(linesize/sizeof(double)) == 0 
	     && (long)unsorted % 64 == 0  && (long)sorted % 64 == 0) {
    //transpose_aligned, SSE2 used
    transpose_aligned(unsorted,sorted,a,b,factor);
  } else {
    // transpose_misaligned
    transpose_misaligned(unsorted,sorted,a,b,factor);
  }    
}
/* $Id$ */
