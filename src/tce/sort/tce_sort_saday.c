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

void tce_sort_0_(double* unsorted,double* sorted,double* factor) {
  *sorted = (*unsorted) * (*factor);
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



void tce_sort_4_scalar(double* unsorted,double* sorted,
		       int a, int b, int c, int d, 
		       int i, int j, int k, int l,
		       double factor) {
  int id[4],jd[4],ia,ib,j1,j2,j3,j4;
  int l1,l2,l3,l4;
  int ia1,ia2,ia3,ia4;
  int ib1,ib2,ib3,ib4;
  int rangea1,rangea2,rangea3,rangea4;
  int rangeb1,rangeb2,rangeb3,rangeb4;
  int range[4],order[4],order_r[4];
  int jj1,jj2,jj3,jj4;
  int jj1_bound,jj2_bound,jj3_bound,jj4_bound;
  int count,ir,jr,kr,lr,N1,N2;
  
  double *pA, *pB;

  jd[0] = a;
  jd[1] = b;
  jd[2] = c;
  jd[3] = d;

  // prefer writes
  
  range[0] = b*c*d;
  range[1] = c*d;
  range[2] = d;
  range[3] = 1;

  l1 = jd[i];
  l2 = jd[j];
  l3 = jd[k];
  l4 = jd[l];
      
  rangea1 = range[i];
  rangea2 = range[j];
  rangea3 = range[k];
  rangea4 = range[l];

  rangeb1 = l2*l3*l4;
  rangeb2 = l3*l4;
  rangeb3 = l4;
  rangeb4 = 1;
  

  if (l == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia = ia3 + j4*rangea4;                           
	    ib = ib3 + j4*rangeb4;
	    sorted[ib] = unsorted[ia] * factor;                       
	  }
	}
      }
    }
  }

  if (k == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3 += tilesize) {   
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj3_bound = (j3 + tilesize > l3)? l3 :j3+tilesize;
	    for (jj3 = j3; jj3 < jj3_bound; jj3++) {
	      ia3 = ia2 + jj3*rangea3;
	      ib3 = ib2 + jj3*rangeb3;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4++) {  
		ia = ia3 + jj4*rangea4;                           
		ib = ib3 + jj4*rangeb4;
		sorted[ib] = unsorted[ia] * factor;
	      }
	    }
	  }
	}
      }
    }
  }

  if (j == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;        
      for (j2 = 0; j2 < l2; j2 += tilesize) { 
	for (j3 = 0; j3 < l3; j3++) {    
	  ia3 = ia1 + j3*rangea3;
	  ib3 = ib1 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj2_bound = (j2 + tilesize > l2)? l2 :j2+tilesize;       
	    for (jj2 = j2; jj2 < jj2_bound; jj2++) {
	      ia2 = ia3 + jj2*rangea2;
	      ib2 = ib3 + jj2*rangeb2;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4++) {  
		ia = ia2 + jj4*rangea4;                           
		ib = ib2 + jj4*rangeb4;
		sorted[ib] = unsorted[ia] * factor;
	      }
	    }
	  }
	}
      }
    }
  }

  if (i == 3) {
    for (j1 = 0; j1 < l1; j1 += tilesize) { 
      for (j2 = 0; j2 < l2; j2++) {
	ia2 = j2*rangea2;
	ib2 = j2*rangeb2;        
	for (j3 = 0; j3 < l3; j3++) {    
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj1_bound = (j1 + tilesize > l1)? l1 :j1+tilesize;       
	    for (jj1 = j1; jj1 < jj1_bound; jj1++) {
	      ia1 = ia3 + jj1*rangea1;
	      ib1 = ib3 + jj1*rangeb1;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4++) {  
		ia = ia1 + jj4*rangea4;                           
		ib = ib1 + jj4*rangeb4;
		sorted[ib] = unsorted[ia] * factor;                 
	      }
	    }
	  }
	}
      }
    }
  }

}


void tce_sort_4_simd(double* unsorted,double* sorted,
		     int a, int b, int c, int d, 
		     int i, int j, int k, int l,
		     double factor) {
  int id[4],jd[4],ia,ib,j1,j2,j3,j4;
  int l1,l2,l3,l4;
  int ia1,ia2,ia3,ia4;
  int ib1,ib2,ib3,ib4;
  int rangea1,rangea2,rangea3,rangea4;
  int rangeb1,rangeb2,rangeb3,rangeb4;
  int range[4],order[4],order_r[4];
  int jj1,jj2,jj3,jj4;
  int jj1_bound,jj2_bound,jj3_bound,jj4_bound;
  int count,ir,jr,kr,lr,N1,N2;
  
  double *pA, *pB;
  register __m128d x, y, z, w, t, t1,fac_vector; 
   
  fac_vector = _mm_load_sd(&factor);
  fac_vector = _mm_unpacklo_pd(fac_vector,fac_vector);

  jd[0] = a;
  jd[1] = b;
  jd[2] = c;
  jd[3] = d;
  // prefer writes
  
  range[0] = b*c*d;
  range[1] = c*d;
  range[2] = d;
  range[3] = 1;

  l1 = jd[i];
  l2 = jd[j];
  l3 = jd[k];
  l4 = jd[l];
      
  rangea1 = range[i];
  rangea2 = range[j];
  rangea3 = range[k];
  rangea4 = range[l];

  rangeb1 = l2*l3*l4;
  rangeb2 = l3*l4;
  rangeb3 = l4;
  rangeb4 = 1;
  

  // here vectorization can rely on the compiler
  if (l == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia = ia3 + j4*rangea4;                           
	    ib = ib3 + j4*rangeb4;
	    sorted[ib] = unsorted[ia] * factor;                       
	  }
	}
      }
    }
  }

  if (k == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3 += tilesize) {   
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj3_bound = (j3 + tilesize > l3)? l3 :j3+tilesize;
	    for (jj3 = j3; jj3 < jj3_bound; jj3 += 2) {
	      ia3 = ia2 + jj3*rangea3;
	      ib3 = ib2 + jj3*rangeb3;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4 += 2) {  
		ia = ia3 + jj4*rangea4;                           
		ib = ib3 + jj4*rangeb4;
		N1 = rangeb3;
		N2 = rangea4;
		       
		pA = unsorted+ia;
		pB = sorted+ib;
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
    }
  }

  if (j == 3) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;        
      for (j2 = 0; j2 < l2; j2 += tilesize) { 
	for (j3 = 0; j3 < l3; j3++) {    
	  ia3 = ia1 + j3*rangea3;
	  ib3 = ib1 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj2_bound = (j2 + tilesize > l2)? l2 :j2+tilesize;       
	    for (jj2 = j2; jj2 < jj2_bound; jj2 += 2) {
	      ia2 = ia3 + jj2*rangea2;
	      ib2 = ib3 + jj2*rangeb2;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4 += 2) {  
		ia = ia2 + jj4*rangea4;                           
		ib = ib2 + jj4*rangeb4;
		N1 = rangeb2;
		N2 = rangea4;
		       
		pA = unsorted+ia;
		pB = sorted+ib;
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
    }
  }

  if (i == 3) {
    for (j1 = 0; j1 < l1; j1 += tilesize) { 
      for (j2 = 0; j2 < l2; j2++) {
	ia2 = j2*rangea2;
	ib2 = j2*rangeb2;        
	for (j3 = 0; j3 < l3; j3++) {    
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) {
	    jj1_bound = (j1 + tilesize > l1)? l1 :j1+tilesize;       
	    for (jj1 = j1; jj1 < jj1_bound; jj1 += 2) {
	      ia1 = ia3 + jj1*rangea1;
	      ib1 = ib3 + jj1*rangeb1;
	      jj4_bound = (j4 + tilesize > l4)? l4:j4+tilesize;
	      for (jj4 = j4; jj4 < jj4_bound; jj4 += 2) {  
		ia = ia1 + jj4*rangea4;                           
		ib = ib1 + jj4*rangeb4;
		N1 = rangeb1;
		N2 = rangea4;
		       
		pA = unsorted+ia;
		pB = sorted+ib;
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
    }
  }

}

void tce_sort_4_(double* unsorted,double* sorted,
		 int* a_in, int* b_in, int* c_in, int* d_in, 
		 int* i_in, int* j_in, int* k_in, int* l_in,
		 double* factor_in) {
  int a = *a_in;
  int b = *b_in;
  int c = *c_in;
  int d = *d_in;
  int i = *i_in - 1;
  int j = *j_in - 1;
  int k = *k_in - 1;
  int l = *l_in - 1;
  double factor = *factor_in;

  int dim[4] ;

  dim[0] = a;
  dim[1] = b;
  dim[2] = c;
  dim[3] = d;
  
  if (dim[l]%2 == 0 && dim[3]%2 == 0)
    tce_sort_4_simd(unsorted,sorted,a,b,c,d,i,j,k,l,factor);
  else
    tce_sort_4_scalar(unsorted,sorted,a,b,c,d,i,j,k,l,factor);
}


void tce_sort_6_simd(double* unsorted,double* sorted,
		     int a, int b, int c, int d, int e, int f,  
		     int i, int j, int k, int l, int m, int n,
		     double factor) {
  int id[6],jd[6],ia,ib,j1,j2,j3,j4,j5,j6;
  int l1,l2,l3,l4,l5,l6;
  int ia1,ia2,ia3,ia4,ia5,ia6;
  int ib1,ib2,ib3,ib4,ib5,ib6;
  int rangea1,rangea2,rangea3,rangea4,rangea5,rangea6;
  int rangeb1,rangeb2,rangeb3,rangeb4,rangeb5,rangeb6;
  int range[6],order[6],order_r[6];
  int jj1,jj2,jj3,jj4,jj5,jj6;
  int jj1_bound,jj2_bound,jj3_bound,jj4_bound,jj5_bound,jj6_bound;
  int N1,N2;
  
  double *pA, *pB;
  register __m128d x, y, z, w, p, q,fac_vector; 
   
  fac_vector = _mm_load_sd(&factor);
  fac_vector = _mm_unpacklo_pd(fac_vector,fac_vector);

  jd[0] = a;
  jd[1] = b;
  jd[2] = c;
  jd[3] = d;
  jd[4] = e;
  jd[5] = f;

  // prefer writes
  range[0] = b*c*d*e*f;
  range[1] = c*d*e*f;
  range[2] = d*e*f;
  range[3] = e*f;
  range[4] = f;
  range[5] = 1;

  l1 = jd[i];
  l2 = jd[j];
  l3 = jd[k];
  l4 = jd[l];
  l5 = jd[m];
  l6 = jd[n]; 

      
  rangea1 = range[i];
  rangea2 = range[j];
  rangea3 = range[k];
  rangea4 = range[l];
  rangea5 = range[m];
  rangea6 = range[n];


  rangeb1 = l2*l3*l4*l5*l6;
  rangeb2 = l3*l4*l5*l6;
  rangeb3 = l4*l5*l6;
  rangeb4 = l5*l6;
  rangeb5 = l6;
  rangeb6 = 1;
   
  // here vectorization can rely on the compiler
  if (n == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia4 = ia3 + j4*rangea4;                           
	    ib4 = ib3 + j4*rangeb4;
	    for (j5 = 0; j5 < l5; j5++) { 
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;
	      for (j6 = 0; j6 < l6; j6++) {   
		ia = ia5 + j6*rangea6;                           
		ib = ib5 + j6*rangeb6;
		sorted[ib] = unsorted[ia] * factor;  
	      }
	    }                     
	  }
	}
      }
    }
  }

  if (m == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia4 = ia3 + j4*rangea4;                           
	    ib4 = ib3 + j4*rangeb4;
	    for (j5 = 0; j5 < l5; j5 += tilesize) {   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj5_bound = (j5 + tilesize > l5)? l5 :j5+tilesize;
		for (jj5 = j5; jj5 < jj5_bound; jj5 += 2) {
		  ia5 = ia4 + jj5*rangea5;
		  ib5 = ib4 + jj5*rangeb5;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6 += 2) {  
		    ia = ia5 + jj6*rangea6;                           
		    ib = ib5 + jj6*rangeb6;
		    N1 = rangeb5;
		    N2 = rangea6;
       
		    pA = unsorted+ia;
		    pB = sorted+ib;
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
	}
      }
    }
  }

  if (l == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) { 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia3 + j5*rangea5;
	      ib5 = ib3 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj4_bound = (j4 + tilesize > l4)? l4 :j4+tilesize;
		for (jj4 = j4; jj4 < jj4_bound; jj4 += 2) {
		  ia4 = ia5 + jj4*rangea4;
		  ib4 = ib5 + jj4*rangeb4;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6 += 2) {  
		    ia = ia4 + jj6*rangea6;                           
		    ib = ib4 + jj6*rangeb6;
		    N1 = rangeb4;
		    N2 = rangea6;
       
		    pA = unsorted+ia;
		    pB = sorted+ib;
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
	}
      }
    }
  }

  if (k == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3 += tilesize) { 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia2 + j4*rangea4;
	    ib4 = ib2 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj3_bound = (j3 + tilesize > l3)? l3 :j3+tilesize;
		for (jj3 = j3; jj3 < jj3_bound; jj3 += 2) {
		  ia3 = ia5 + jj3*rangea3;
		  ib3 = ib5 + jj3*rangeb3;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6 += 2) {  
		    ia = ia3 + jj6*rangea6;                           
		    ib = ib3 + jj6*rangeb6;
		    N1 = rangeb3;
		    N2 = rangea6;
       
		    pA = unsorted+ia;
		    pB = sorted+ib;
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
	}
      }
    }
  }


  if (j == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2 += tilesize) {    
	for (j3 = 0; j3 < l3; j3++) {
	  ia3 = ia1 + j3*rangea3;
	  ib3 = ib1 + j3*rangeb3; 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia3 + j4*rangea4;
	    ib4 = ib3 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj2_bound = (j2 + tilesize > l2)? l2 :j2+tilesize;
		for (jj2 = j2; jj2 < jj2_bound; jj2 += 2) {
		  ia2 = ia5 + jj2*rangea2;
		  ib2 = ib5 + jj2*rangeb2;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6 += 2) {  
		    ia = ia2 + jj6*rangea6;                           
		    ib = ib2 + jj6*rangeb6;
		    N1 = rangeb2;
		    N2 = rangea6;
       
		    pA = unsorted+ia;
		    pB = sorted+ib;
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
	}
      }
    }
  }

  if (i == 5) {
    for (j1 = 0; j1 < l1; j1 += tilesize) {         
      for (j2 = 0; j2 < l2; j2++) {
	ia2 = j2*rangea2;
	ib2 = j2*rangeb2;      
	for (j3 = 0; j3 < l3; j3++) {
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3; 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia3 + j4*rangea4;
	    ib4 = ib3 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj1_bound = (j1 + tilesize > l1)? l1 :j1+tilesize;
		for (jj1 = j1; jj1 < jj1_bound; jj1 += 2) {
		  ia1 = ia5 + jj1*rangea1;
		  ib1 = ib5 + jj1*rangeb1;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6 += 2) {  
		    ia = ia1 + jj6*rangea6;                           
		    ib = ib1 + jj6*rangeb6;
		    N1 = rangeb1;
		    N2 = rangea6;
       
		    pA = unsorted+ia;
		    pB = sorted+ib;
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
	}
      }
    }
  }

}

void tce_sort_6_scalar(double* unsorted,double* sorted,
		       int a, int b, int c, int d, int e, int f,  
		       int i, int j, int k, int l, int m, int n,
		       double factor) {
  int id[6],jd[6],ia,ib,j1,j2,j3,j4,j5,j6;
  int l1,l2,l3,l4,l5,l6;
  int ia1,ia2,ia3,ia4,ia5,ia6;
  int ib1,ib2,ib3,ib4,ib5,ib6;
  int rangea1,rangea2,rangea3,rangea4,rangea5,rangea6;
  int rangeb1,rangeb2,rangeb3,rangeb4,rangeb5,rangeb6;
  int range[6],order[6],order_r[6];
  int jj1,jj2,jj3,jj4,jj5,jj6;
  int jj1_bound,jj2_bound,jj3_bound,jj4_bound,jj5_bound,jj6_bound;
  int N1,N2;
  

  jd[0] = a;
  jd[1] = b;
  jd[2] = c;
  jd[3] = d;
  jd[4] = e;
  jd[5] = f;

  // prefer writes
  range[0] = b*c*d*e*f;
  range[1] = c*d*e*f;
  range[2] = d*e*f;
  range[3] = e*f;
  range[4] = f;
  range[5] = 1;

  l1 = jd[i];
  l2 = jd[j];
  l3 = jd[k];
  l4 = jd[l];
  l5 = jd[m];
  l6 = jd[n]; 

      
  rangea1 = range[i];
  rangea2 = range[j];
  rangea3 = range[k];
  rangea4 = range[l];
  rangea5 = range[m];
  rangea6 = range[n];


  rangeb1 = l2*l3*l4*l5*l6;
  rangeb2 = l3*l4*l5*l6;
  rangeb3 = l4*l5*l6;
  rangeb4 = l5*l6;
  rangeb5 = l6;
  rangeb6 = 1;
   
  // here vectorization can rely on the compiler
  if (n == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia4 = ia3 + j4*rangea4;                           
	    ib4 = ib3 + j4*rangeb4;
	    for (j5 = 0; j5 < l5; j5++) { 
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;
	      for (j6 = 0; j6 < l6; j6++) {   
		ia = ia5 + j6*rangea6;                           
		ib = ib5 + j6*rangeb6;
		sorted[ib] = unsorted[ia] * factor;  
	      }
	    }                     
	  }
	}
      }
    }
  }

  if (m == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4++) {   
	    ia4 = ia3 + j4*rangea4;                           
	    ib4 = ib3 + j4*rangeb4;
	    for (j5 = 0; j5 < l5; j5 += tilesize) {   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj5_bound = (j5 + tilesize > l5)? l5 :j5+tilesize;
		for (jj5 = j5; jj5 < jj5_bound; jj5++) {
		  ia5 = ia4 + jj5*rangea5;
		  ib5 = ib4 + jj5*rangeb5;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6++) {  
		    ia = ia5 + jj6*rangea6;                           
		    ib = ib5 + jj6*rangeb6;
		    sorted[ib] = unsorted[ia] * factor;    
		  }
		}                     
	      }
	    }
	  }
	}
      }
    }
  }

  if (l == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3++) { 
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3;
	  for (j4 = 0; j4 < l4; j4 += tilesize) { 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia3 + j5*rangea5;
	      ib5 = ib3 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj4_bound = (j4 + tilesize > l4)? l4 :j4+tilesize;
		for (jj4 = j4; jj4 < jj4_bound; jj4++) {
		  ia4 = ia5 + jj4*rangea4;
		  ib4 = ib5 + jj4*rangeb4;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6++) {  
		    ia = ia4 + jj6*rangea6;                           
		    ib = ib4 + jj6*rangeb6;
		    sorted[ib] = unsorted[ia] * factor;
		  }
		}                     
	      }
	    }
	  }
	}
      }
    }
  }

  if (k == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2++) {    
	ia2 = ia1 + j2*rangea2;
	ib2 = ib1 + j2*rangeb2;
	for (j3 = 0; j3 < l3; j3 += tilesize) { 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia2 + j4*rangea4;
	    ib4 = ib2 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj3_bound = (j3 + tilesize > l3)? l3 :j3+tilesize;
		for (jj3 = j3; jj3 < jj3_bound; jj3++) {
		  ia3 = ia5 + jj3*rangea3;
		  ib3 = ib5 + jj3*rangeb3;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6++) {  
		    ia = ia3 + jj6*rangea6;                           
		    ib = ib3 + jj6*rangeb6;
		    sorted[ib] = unsorted[ia] * factor;
		  }
		}                     
	      }
	    }
	  }
	}
      }
    }
  }


  if (j == 5) {
    for (j1 = 0; j1 < l1; j1++) {
      ia1 = j1*rangea1;
      ib1 = j1*rangeb1;           
      for (j2 = 0; j2 < l2; j2 += tilesize) {    
	for (j3 = 0; j3 < l3; j3++) {
	  ia3 = ia1 + j3*rangea3;
	  ib3 = ib1 + j3*rangeb3; 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia3 + j4*rangea4;
	    ib4 = ib3 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj2_bound = (j2 + tilesize > l2)? l2 :j2+tilesize;
		for (jj2 = j2; jj2 < jj2_bound; jj2++) {
		  ia2 = ia5 + jj2*rangea2;
		  ib2 = ib5 + jj2*rangeb2;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6++) {  
		    ia = ia2 + jj6*rangea6;                           
		    ib = ib2 + jj6*rangeb6;
		    sorted[ib] = unsorted[ia] * factor;
		  }
		}                     
	      }
	    }
	  }
	}
      }
    }
  }

  if (i == 5) {
    for (j1 = 0; j1 < l1; j1 += tilesize) {         
      for (j2 = 0; j2 < l2; j2++) {
	ia2 = j2*rangea2;
	ib2 = j2*rangeb2;      
	for (j3 = 0; j3 < l3; j3++) {
	  ia3 = ia2 + j3*rangea3;
	  ib3 = ib2 + j3*rangeb3; 
	  for (j4 = 0; j4 < l4; j4++) {
	    ia4 = ia3 + j4*rangea4;
	    ib4 = ib3 + j4*rangeb4; 
	    for (j5 = 0; j5 < l5; j5++) {
	      ia5 = ia4 + j5*rangea5;
	      ib5 = ib4 + j5*rangeb5;   
	      for (j6 = 0; j6 < l6; j6 += tilesize) {
		jj1_bound = (j1 + tilesize > l1)? l1 :j1+tilesize;
		for (jj1 = j1; jj1 < jj1_bound; jj1++) {
		  ia1 = ia5 + jj1*rangea1;
		  ib1 = ib5 + jj1*rangeb1;
		  jj6_bound = (j6 + tilesize > l6)? l6:j6+tilesize;
		  for (jj6 = j6; jj6 < jj6_bound; jj6++) {  
		    ia = ia1 + jj6*rangea6;                           
		    ib = ib1 + jj6*rangeb6;
		    sorted[ib] = unsorted[ia] * factor;
		  }
		}                     
	      }
	    }
	  }
	}
      }
    }
  }

}



void tce_sort_6_(double* unsorted,double* sorted,
		 int* a_in, int* b_in, int* c_in, int* d_in,int* e_in, int* f_in, 
		 int* i_in, int* j_in, int* k_in, int* l_in,int* m_in, int* n_in,
		 double* factor_in) {
  int a = *a_in;
  int b = *b_in;
  int c = *c_in;
  int d = *d_in;
  int e = *e_in;
  int f = *f_in;
  int i = *i_in - 1;
  int j = *j_in - 1;
  int k = *k_in - 1;
  int l = *l_in - 1;
  int m = *m_in - 1;
  int n = *n_in - 1;
  double factor = *factor_in;

  int dim[6] ;

  dim[0] = a;
  dim[1] = b;
  dim[2] = c;
  dim[3] = d;
  dim[4] = e;
  dim[5] = f;
  
  if (dim[n]%2 == 0 && dim[5]%2 == 0)
    tce_sort_6_simd(unsorted,sorted,a,b,c,d,e,f,i,j,k,l,m,n,factor);
  else
    tce_sort_6_scalar(unsorted,sorted,a,b,c,d,e,f,i,j,k,l,m,n,factor);
}

/* $Id$ */
