/*
 $Id$
 */

#include <math.h>
#include <stdio.h>
#ifndef WIN32
#include <strings.h>
#endif
#include <stdlib.h>
#include "ga.h"

extern int fortchar_to_string();

/* Maximum number of open shells and corresponding maximum number 
   of spin functions for any possible multiplicity */
#define nsmax 14
#define nfmax 1001

/* #define nsmax 15 */
/* #define nfmax 2002 */

static void Error(char *string, int integer)
{
/*
  (void) fflush(stdout);
  (void) fprintf(stderr,"\n\nError was called.\n");
  (void) fprintf(stderr,string);
  (void) fprintf(stderr," %d (%#x).\n",integer,integer);
  exit(1); */
  GA_Error(string, (long) integer);
}


static void Dscal(int n, double s, double *u, int iu)
/*
  Scale double precision vector u by s
*/
{
  while (n--) {
    *u = *u * s;
    u += iu;
  }
}

static void Dfill(int n, double s, double *u, int iu)
/*
  Fill double precision vector u with s
*/
{
  while (n--) {
    *u = s;
    u += iu;
  }
}

static void Screen(int n, double s, double *u, int iu)
/*
  zero elements of double precision vector u if they are less
  in absolute magnitude than the threshold s
*/
{
  s = fabs(s);
  while (n--) {
    if (fabs(*u) < s)
      *u = 0.0;
    u += iu;
  }
}

static double Sparsity(int n, double *u, int iu)
/*
  Return (no. zero elements)/(total no. elements)
*/
{
  int zero=0, nn;

  nn = n;
  while (nn--) {
    if (*u == 0.0)
      ++zero;
    u += iu;
  }

  return (double) zero / (double) n;
}

static void Identity(double *u, int n)
/*
  Form the n*n double precision identity matrix
*/
{
  Dfill(n*n, 0.0, u, 1);

  Dfill(n, 1.0, u, n+1);
}

static void PrintMatrix(const double *u, int mcol, int mrow, 
			int ncol, int nrow)
/*
  Print a double precision matrix.
  mcol = skip between elements down a column
  mrow = skip between elements across a column
  ncol = no. of columns
  nrow = no. of rows
*/
{
  double *uu, *uuu;
  int todo, nnrow, nncol, ttodo, i;

  (void) printf("\n\n");
  nncol = ncol;
  while(nncol) {

    if (nncol < 6)
      todo = nncol;
    else
      todo = 6;
    
    (void) printf("  ");
    for (i=1; i<=todo; i++)
      (void) printf("%12d",ncol-nncol+i);
    (void) printf("\n\n");
    
    
    uu = (double *)u;
    nnrow = nrow;
    while(nnrow--) {
      
      (void) printf("%3d  ",nrow-nnrow);
      uuu = uu;
      ttodo = todo;
      while(ttodo--) {
	(void) printf("%12.6f",*uuu);
	uuu += mrow;
      }
      (void) printf("\n");
      uu += mcol;
    }
    
    (void) printf("\n\n");
    u += todo*mrow;
    nncol -= todo;
  }
}

static void MakeBranchingDiagram(int bd[nsmax+3][nsmax+3], 
				 int ns, int multi, int order, int info)
/*
  Form the branching diagram weights
  bd[m][n] = weight of node with multiplicty m and n electrons
  ns = maximum no. of open shells or unpaired electrons
  multi = state multiplicity
  order = 0 -> normal branching diagram
          1 -> reversed diagram (for use with dictionary order functions)
*/
{
  int i, j, wrote;
  
  for (i=0; i<ns+3; i++)
    for (j=0; j<ns+3; j++)
      bd[i][j] = 0;

  if (order == 0) {
    bd[1][0] = 1;                   
    
    for (i=1; i<=ns; i++) {
      for (j=i%2+1; j<=ns+1; j+=2) {
	if (abs(multi-j) > (ns-i))
	  continue;
	bd[j][i] = bd[j-1][i-1] + bd[j+1][i-1];
      }
    }
  }
  else {
    bd[multi][ns] = 1;
    
    for (i=ns-1; i>=0; i--) {
      for (j=i%2+1; j<=ns+1; j+=2) {
	if (j > i+1)
	  continue;
	bd[j][i] = bd[j-1][i+1] + bd[j+1][i+1];
      }
    }
  }
  
#ifndef DEBUG
  if (order == 1)
    return;
#endif

  if (info) {
      if (order == 0)
	  (void) printf("\nBranching diagram weights.\n\n");
      else
	  (void) printf("\nReversed branching diagram weights.\n\n");
      
      
      for (j=ns+1; j>=1; j--) {
	  wrote = 0;
	  for (i=0; i<=ns; i++)
	      if (bd[j][i])
		  wrote=1;
	  if (wrote) {
	      for (i=0; i<=ns; i++) {
		  if (bd[j][i])
		      (void) printf("%4d",bd[j][i]);
		  else
		      (void) printf("    ");
	      }
	      (void) printf("\n");
	  }
      }
      (void) printf("\n");
  }
}

static void MakeSpinFunctions(int bd[nsmax+3][nsmax+3], 
			      int ns, int multi, int f[nsmax+1][nfmax+1],
			      int nf, int order, int info)
/*
  Form the ordered spin functions as an array of 1s and 2s.
  bd = branching diagram (reversed if want dictionary ordered functions)
  ns = required no. of unpaired electrons
  multi = state multiplicity
  f[j][i] = occupation (1 or 2) of orbital j in function i
  nf = no. of functions
  order = 0 -> reversed dictionary, 1 -> dictionary

  if (reversed dictionary)
    weight of arc (m,k)->(m+1,k+1) = weight of node (m+2,k)
    and the graph has to be traversed right to left
  else
    weight of arc (m,k)->(m+1,k+1) = weight of node (m-1,k+1)
    and the graph has to be traversed left to right
*/
{
  int i, test, j, m;

  if (order ==0) {
    for (i=1; i<=nf; i++) {
      test = i-1;
      m = multi;
      for (j=ns; j>0; j--){
	
	if (bd[m-1][j-1] && (test-bd[m+1][j-1]) >= 0) {
	  test -= bd[m+1][j-1];
	  f[j][i] = 1;
	  m--;
	}
	else if (bd[m+1][j-1]){
	  f[j][i] = 2;
	  m++;
	}
	else
	  Error("mis-match walking spin-graph",j);
      }
      if (test != 0 || m != 1) {
	(void) printf("test=%d, m=%d\n",test,m);
	Error("invalid tail on spin graph",i);
      }
    }
  }
  else {
    for (i=1; i<=nf; i++) {
      test = i-1;
      m = 1;
      for (j=0; j<ns; j++){
	
	if (bd[m+1][j+1] && (test-bd[m-1][j+1]) >= 0) {
	  test -= bd[m-1][j+1];
	  f[j+1][i] = 1;
	  m++;
	}
	else if (bd[m-1][j+1]){
	  f[j+1][i] = 2;
	  m--;
	}
	else
	  Error("mis-match walking reversed spin-graph",j);
      }
      if (test != 0 || m != multi) {
	(void) printf("test=%d, m=%d\n",test,m);
	Error("invalid tail on reversed spin graph",i);
      }
    }
  }

  if (info) {
      (void) printf("\nSpin functions\n\n");
      
      test = 0;
      for (i=1; i<=nf; i++) {
	  test += printf("%4d [",i);
	  for (j=1; j<=ns; j++)
	      test += printf("%1d",f[j][i]);
	  if (test > (71-ns)) {
	      test = 0;
	      (void) printf("]\n");
	  }
	  else
	      test += printf("] ");
      }
      if (test)
	  (void) printf("\n");
  }
}


static void MakeAxialDistances(int ns, int f[nsmax+1][nfmax+1], 
			       int nf, int d[nsmax][nfmax+1], 
			       double p[nsmax][nfmax+1], double g[nsmax][nfmax+1])
/*
  Compute the axial distances for transpostions (j, j+1) for
  the Standard Young Tableaux represented by the branching diagram
  paths in f. (See Paunz, "Spin Eigenfunctions", p106.)
  ns = no. of electrons
  f[j][i] = 1 or 2, step j in the ith function
  d[j][i] = axial distance for (j,j+1) in tableau i
  p[j][i] = -1.0/d[j][i]
  g[j][i] = sqrt(1-p[j][i]^2)
*/
{
  int i, j, k, jb;

#ifdef DEBUG
  int test;
#endif

  for (i=1; i<=nf; i++) {
    for (j=1; j<ns; j++) {
      
      if (f[j][i] == f[j+1][i])
	d[j][i] = -1;
      else {
	jb=0;
	for (k=1; k<j; k++) {
	  jb += (f[k][i] == f[j][i]);
	  jb -= (f[k][i] == f[j+1][i]);
	}
	if (f[j][i] == 1)
	  d[j][i] = jb + 1;
	else
	  d[j][i] = jb - 1;
      }
    }
  }

  for (j=1; j<ns; j++) {
    for (i=1; i<=nf; i++) {
      p[j][i] = -1.0 / d[j][i];
      g[j][i] = sqrt(1.0 - p[j][i]*p[j][i]);
    }
  }
  
#ifdef DEBUG
  (void) printf("\nAxial distances for (j,j+1)\n\n");

  test = 0;
  for (i=1; i<=nf; i++) {
    test += printf("%4d [",i);
    for (j=1; j<ns; j++)
      test += printf("%4d",d[j][i]);
    if (test > (71-4*ns)) {
      test = 0;
      (void) printf("]\n");
    }
    else
      test += printf("] ");
  }
  if (test)
    (void) printf("\n");
#endif
}


static void MakePerm(int bd[nsmax+3][nsmax+3], int ns, int f[nsmax+1][nfmax+1], 
		     int nf, int perm[nsmax][nfmax+1], int order)
/*
  Compute the effect of the transpositions (j j+1) on the spin functions
  bd = branching diagram (reversed if dictionary order)
  ns = no. of electrons
  f[j][i] = 1 or 2, step j in ith function
  nf = no. of functions
  perm[j][i] = result of permutation (j j+1) on the ith tableau.
               is the no. of another standard tableau or zero
  order = 0 -> reversed dictionary order, 1 -> dictionary order
  See comments in MakeBranchingDiagram
*/
{
  int i, j, k, m, test, step;

  for (i=1; i<=nf; i++) {

    for (j=1; j<ns; j++) {

      test = 1;
      m = 1;
      for (k=1; k<=ns; k++) {

	if ( k == j )
	  step = f[j+1][i];
	else if ( k == j+1 )
	  step = f[j][i];
	else
	  step = f[k][i];

	if (step == 1) {
	  if (order == 0)
	    test = test + bd[m+2][k-1];
	  else
	    test = test + bd[m-1][k];
	  m++;
	}
	else
	  m--;

	if ( m<=0 || bd[m][k] == 0) {
	    test = 0;
	    break;
	  }
      }
      perm[j][i] = test;
    }
  }

#ifdef DEBUG
  (void) printf("\nResults of permutations for (j,j+1)\n\n");

  test = 0;
  for (i=1; i<=nf; i++) {
    test += printf("%4d [",i);
    for (j=1; j<ns; j++)
      test += printf("%4d",perm[j][i]);
    if (test > (71-4*ns)) {
      test = 0;
      (void) printf("]\n");
    }
    else
      test += printf("] ");
  }
  if (test)
    (void) printf("\n");
#endif
}


static void RightMultiply(int i, double *u, double *uu, int nf, 
			  double pp[nsmax][nfmax+1], double g[nsmax][nfmax+1], 
			  int perm[nsmax][nfmax+1])
/*
  uu = (-1) * u * U(i i+1), where u and uu are pointers to
  other represenation matrices.
  i = index of required permutation (i i+1)
  u = pointer to input matrix
  uu = pointer to result matrix
  nf = no. of functions = dimension of representation
  pp[i][j] = negative reciprocal axial distance of (i i+1) in tableau j
  gg[i][j] = sqrt(1-pp[i][j]^2)
  perm[i][j] = action of (i i+1) on tableau j
*/
{
  int p, q, m;

  for (p=0; p<nf; p++) {
    for (q=0; q<nf; q++) {
      
      *(uu+p+q*nf) = -(*(u+p+q*nf) * pp[i][q+1]);
      m = perm[i][q+1];
      if (m)
	*(uu+p+q*nf) -= *(u+p+(m-1)*nf)*g[i][q+1];
    }
  }
}
static void LeftMultiply(int i, double *u, double *uu, int nf, 
			 double pp[nsmax][nfmax+1], double g[nsmax][nfmax+1], 
			 int perm[nsmax][nfmax+1])
/*
  uu = (-1) * U(i i+1) * u, where u and uu are pointers to
  other represenation matrices.
  i = index of required permutation (i i+1)
  u = pointer to input matrix
  uu = pointer to result matrix
  nf = no. of functions = dimension of representation
  pp[i][j] = negative reciprocal axial distance of (i i+1) in tableau j
  gg[i][j] = sqrt(1-pp[i][j]^2)
  perm[i][j] = action of (i i+1) on tableau j
*/
{
  int p, q, m;

  for (p=0; p<nf; p++) {
    for (q=0; q<nf; q++) {
      
      *(uu+p+q*nf) = -(*(u+p+q*nf) * pp[i][p+1]);
      m = perm[i][p+1];
      if (m)
	*(uu+p+q*nf) -= *(u+m-1+q*nf)*g[i][p+1];
    }
  }
}


#ifdef OLDCODE
static void ReadArguments(argc, argv, ns, multi, print)
     int argc, *ns, *multi, *print;
     char **argv;
/*
  The syntax is

  command [-p] multiplicity [ nsmax ]

  If no arguments are specified the command syntax is returned.
  If only the multiplicity is specified then a sensible default
  for nsmax is used (this may be too small for large calculations).
  The -p option, which must precede the other arguments will cause
  more printout to be generated.
*/  
{
  *print = 0;

  if (argc==1 || (argc==2 && !strcmp(*(argv+1),"-p"))) {
    (void) fprintf(stderr,"Usage: %s [-p] multiplicity [ nsmax ]\n",*argv);
    exit(1);
  }

  argv++;
  argc--;

  if (!strcmp(*argv,"-p")) {
    *print = 1;
    argv++;
    argc--;
  }

  *multi = atoi(*(argv++));
  if (argc == 1)
    *ns = 8 + (*multi-1)%2;
  else
    *ns = atoi(*argv);
}

/*int main(argc, argv)
     int argc;
     char **argv;
     */
#endif

#if (defined(CRAY) ||defined(WIN32))&& !defined(__crayx1) && !defined(__MINGW32__)
void FATR SELCI_COUPLE(Integer *pmulti, Integer *pns, Integer *pprint, char *pfilename, int flen)
#else
void FATR selci_couple_(Integer *pmulti, Integer *pns, Integer *pprint, char *pfilename, int flen)
#endif
/*
  Generate the one particle coupling coefficients between two orbital
  occupancies including all spin couplings. The spin functions are
  the standard genealogical functions in the order specified by the
  variable order. If dictionary order is used (order=1) then the matrices
  for smaller no.s of open shells are just submatrices of the largest
  dimension matrix. Thus order is hardwired to 1 just below, but
  order = 0 (reverse dictionary) does work.

  n(i)=1
        W1(u,v,i) = <Iu|Eia|Jv> = (-1)^p U(P)uv , P=(i...ns)

  n(i)=2
        W2(u,v,i) = <Iu|Eia|Jv> = sqrt(2) (-1)^p U(P)uv , P=(1...ns)(1...i)

  i is the position of the orbital i in the open shell orbitals.
  ns is the no. of open shells, p is the parity of P, and a is an
  orbital higher in index (so coupled last) than any in |I>, which
  represents an orbital occupation of the form

  |I> = d1(1)d1(2)d2(3)d2(4) ... dn(2n-1)dn(2n)dn+1(2n+1)...dn+ns(2n+ns)

  That is the ordered doubly occupied orbitals are put first and the
  the singly occupied are ordered on the end.

  The matrices U(P) are the representation matrices for the permutation
  P in the genealogical basis.

  The input is just the two command arguments multiplicity and ns.

  The output file 'wmatrix' is a binary stream containing the following

  multi
  nsmax
  (nf(i),i=multi-1,nsmax,2)
  (((w1(u,v,i),u=1,nf(nsmax)),v=1,nf(nsmax)),i=nsmax,1)
  (((w2(u,v,i),u=1,nf(nsmax-2)),v=1,nf(nsmax)),i=1,nsmax-1)
  
  multi  = multiplicity
  nf(ns) = no. of spin functions for ns open shells
  w1,w2 are the desired matrices. note the odd i order for w1.

  Note that in dictionary order w*(*,*,i) for ns=nsmax-m is just
  w*(*,*,i+m).

*/
{
  int order = 1;
  static int bd[nsmax+3][nsmax+3], rbd[nsmax+3][nsmax+3];
  static int d[nsmax][nfmax+1], f[nsmax+1][nfmax+1];
  int ns, multi, nf, nf2, print, info;
  static int perm[nsmax][nfmax+1], i, j, k, kk;
  static double p[nsmax][nfmax+1], g[nsmax][nfmax+1], *u, *uu, *temp, *ttemp;
  FILE *wmatrix;
  char filename[255];

  /*ReadArguments(argc, argv, &ns, &multi, &print);*/

  ns = *pns;
  multi = *pmulti;
  if (!fortchar_to_string(pfilename, flen, filename, sizeof filename))
      Error("couple: failed to convert fortran filename", 0);
  info = *pprint;		/*  Information output */
  print = 0;			/*  Debug output */

  if (info) {
      (void) printf("\n");
      (void) printf("      Symmetric Group Coupling Coefficient Generation Program (8/9/89)\n");
      (void) printf("      ----------------------------------------------------------------\n\n");
      
      (void) printf("Maximum dimension of open shell  %2d\n",ns);
  }

  if (ns < 0 || ns > nsmax)
    Error("ns is invalid",ns);

  if (info) {
      (void) printf("Multiplicity of electronic state %2d\n",multi);
      
      if (order == 0)
	  (void) printf("Spin functions will be in reverse dictionary order\n");
      else
	  (void) printf("Spin functions will be in dictionary order\n");
  }
  
  if (print)
    (void) printf("High level print requested\n");

  if ( (multi<1) || (multi>ns+1) || (((multi+ns)%2) != 1) )
    Error("multiplicity invalid or not consistent with ns",multi);

  MakeBranchingDiagram(bd, ns, multi, 0, info);


  nf = bd[multi][ns];
  if (nf > nfmax)
    Error("number of spin functions exceeds maximum",nf);
  nf = nf*nf*sizeof(double);
  if ( (u = (double *) malloc((unsigned) nf)) == (double *) NULL )
    Error("failed to allocate memory for u",nf);
  if ( (uu = (double *) malloc((unsigned) nf)) == (double *) NULL )
    Error("failed to allocate memory for uu",nf);

/* open the file for coupling coefficients and write header info */

  if ( (wmatrix = fopen(filename, "w")) == (FILE *) NULL)
    Error("failed to open file wmatrix for write",0);
  
  (void) fprintf(wmatrix,"%d\n%d\n",multi,ns);
  for (i=ns%2; i<=ns; i+=2)
    (void) fprintf(wmatrix,"%d ",bd[multi][i]);
  (void) fprintf(wmatrix,"\n");

/* with dictionary order can get matrices for smaller values of ns
   from the largest value, so only make that. Using the
   commented for loop computes the matrices for all values of ns.
   for (i=ns%2; i<=ns; i+=2) {
*/
  for (i=ns; i<=ns; i+=2) {

    nf = bd[multi][i];
    if (nf == 0)
      continue;

    if (print) {
      (void) printf("\n\n---------------------\n");
      (void) printf("No. of open shells %d\n",i);
      (void) printf("No. of functions   %d\n",nf);
      (void) printf("---------------------\n\n");
    }

    if (order == 0)
      MakeSpinFunctions(bd, i, multi, f, nf, 0, info);
    else {
      MakeBranchingDiagram(rbd, i, multi, 1, info);
      MakeSpinFunctions(rbd, i, multi, f, nf, 1, info);
    }

    MakeAxialDistances(i, f, nf, d, p, g);

    if (order==0)
      MakePerm(bd, i, f, nf, perm,0);
    else
      MakePerm(rbd, i, f, nf, perm,1);

/* Make the matrices (-1)^p U(P), P=(j...m) */

    Identity(u, nf);
    if (print) {
      (void) printf("\n(-1)^p U(P), P=(%d...%d)\n",i,i);
      PrintMatrix(u, 1, nf, nf, nf);
    }

/*  Write matrix out to disk  */
    temp = u;
    for (k=0; k<nf*nf; k++)
      (void) fprintf(wmatrix,"%.14g\n",*temp++);

    for (j=i-1; j>=1; j--) {
      LeftMultiply(j, u, uu, nf, p, g, perm);
      temp = uu;
      uu = u;
      u = temp;
      
      Screen(nf*nf, 1.0e-8, u, 1);
      if (print) {
	(void) printf("\n(-1)^p U(P), P=(%d...%d)\n",j,i);
	(void) printf("Sparsity of the matrix is %4.2f\n",Sparsity(nf*nf,u,1));
	PrintMatrix(u,1,nf,nf,nf);
      }

/*  Write matrix out to disk  */
      temp = u;
      for (k=0; k<nf*nf; k++)
	(void) fprintf(wmatrix,"%.14g\n",*temp++);
      
    }
    
    if (i >= 2 && (nf2 = bd[multi][i-2]) ) {
      
/* Make the matrices (-1)^p sqrt(2) U(P), P=(1...m)(1...j)
   Note that in u we already have (-1)^m U(1...m) */

      Dscal(nf*nf, sqrt(2.0), u, 1);
      
      if (print) {
	(void) printf("\n(-1)^p sqrt(2) U(P), P=(1...%d)(1...%d)\n",i,1);
	(void) printf("Sparsity of the matrix is %4.2f\n",Sparsity(nf*nf,u,1));
	PrintMatrix(u, 1, nf, nf, nf2);
      }
      
/*  Write matrix out to disk  */
      temp = u;
      for (k=0; k<nf; k++) {
	ttemp = temp;
	for (kk=0; kk<nf2; kk++)
	  (void) fprintf(wmatrix,"%.14g\n",*ttemp++);
	temp += nf;
      } 
      
      for (j=1; j<i-1; j++) {
	RightMultiply(j, u, uu, nf, p, g, perm);
	temp = uu;
	uu = u;
	u = temp;
	
	Screen(nf*nf, 1.0e-8, u, 1);
	if (print) {
	  (void) printf("\n(-1)^p sqrt(2) U(P), P=(1...%d)(1...%d)\n",i,j+1);
	  (void) printf("Sparsity of the matrix is %4.2f\n",
			Sparsity(nf*nf,u,1));
	  PrintMatrix(u,1,nf,nf,nf2);
	}
/*  Write matrix out to disk  */
	temp = u;
	for (k=0; k<nf; k++) {
	  ttemp = temp;
	  for (kk=0; kk<nf2; kk++)
	    (void) fprintf(wmatrix,"%.14g\n",*ttemp++);
	  temp += nf;
	} 
	
      } 
    }
    
  }

  (void) free((char *) u);
  (void) free((char *) uu);
  (void) fclose(wmatrix);
}

  

