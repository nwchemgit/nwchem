#define max(a,b) (a > b?a:b)
#define min(a,b) (a < b?a:b)
#include<stdio.h>

#ifdef __cplusplus
extern "C" 
#endif
double getTime(double bandwidth, unsigned long int volume)//return time in s
{
	//printf("vol = %lu, bandwidth = %lf\t", volume, bandwidth);
	return (2*volume*8)/(bandwidth*1000000000);
}

#ifdef __cplusplus
extern "C" 
#endif
double getBW_nooverlap(double eff)
{
	if(eff == 0) return 0;
	
double intercept = 137.113;
double W = 65.348;
return intercept + W * eff;
}

#ifdef __cplusplus
extern "C" 
#endif
double getBW_overlap(double eff)
{
	if(eff == 0) return 0;
double intercept = 139.053;
double W = 55.536;
return intercept + W * eff;

}

#ifdef __cplusplus
extern "C" 
#endif
double getBW_nomatchg32(double eff)
{
	if(eff == 0) return 0;
double intercept = 66.560;
double W = 164.83;
return intercept + W * eff;
}

#ifdef __cplusplus
extern "C" 
#endif
double getBW_matchl32(double eff, unsigned tbsize)
{
	if(eff == 0) return 0;
double intercept = 119.2547;
double W = 74.8634;
double T = -1.4558;
return intercept + W * eff + T * tbsize;
}

#ifdef __cplusplus
extern "C" 
#endif
double getBW_matchg32()
{
double intercept = 196.82;
return intercept;
}

#ifdef __cplusplus
extern "C" 
#endif
double getEfficiency_nooverlap(int ilimit, int olimit, int asize, int bsize, int blockA, int blockB)
{
//	return 0;
	const int remainder1 = asize % blockA;
        const int remainder2 = bsize % blockB;
        const int ilimitr = ilimit * remainder1 / blockA;
        const int olimitr = olimit * remainder2 / blockB;
//printf("\tilimit=%d\tolimit=%d\t", ilimit, olimit);
//printf("\t%d\t%d\t%d\t%d\t", ilimit/32, ilimit%32, olimit/32,olimit%32 );
double f1, f2, f3, f4, f;
f1 =  ((ilimit/32) * (olimit/32) + (double)(ilimit/32) * (olimit%32) /32+ (double)(ilimit%32) * (olimit/32) /32 + (double)(ilimit%32) * (olimit%32) /(32*32) )/ (int)(((ilimit+31)/32) * ((olimit+31)/32));
f2 =  ((ilimitr/32) * (olimit/32) + (double)(ilimitr/32) * (olimit%32) /32+ (double)(ilimitr%32) * (olimit/32) /32 + (double)(ilimitr%32) * (olimit%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimit+31)/32)));
f3 =  ((ilimit/32) * (olimitr/32) + (double)(ilimit/32) * (olimitr%32) /32+ (double)(ilimit%32) * (olimitr/32) /32 + (double)(ilimit%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimit+31)/32) * ((olimitr+31)/32)));
f4 =  ((ilimitr/32) * (olimitr/32) + (double)(ilimitr/32) * (olimitr%32) /32+ (double)(ilimitr%32) * (olimitr/32) /32 + (double)(ilimitr%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimitr+31)/32)));
f = ((asize/blockA) * (bsize/blockB) *f1 + (double)((asize/blockA) * (bsize%blockB > 0) *f3)+ (double)((asize%blockA > 0) * (bsize/blockB)*f2)  + (double)((asize%blockA>0) * (bsize%blockB > 0) *f4) )/ (int)(((asize+blockA-1)/blockA) * ((bsize+blockB-1)/blockB));
//printf("\t%lf\t", f);
return f;
}

#ifdef __cplusplus
extern "C" 
#endif
double getEfficiency_overlap(int ilimit, int olimit, int asize, int bsize, int blockA, int blockB)
{
//	return 0;
	const int remainder1 = asize % blockA;
        const int remainder2 = bsize % blockB;
        const int ilimitr = ilimit * remainder1 / blockA;
        const int olimitr = olimit * remainder2 / blockB;
double f1, f2, f3, f4, f;
f1 =  ((ilimit/32) * (olimit/32) + (double)(ilimit/32) * (olimit%32) /32+ (double)(ilimit%32) * (olimit/32) /32 + (double)(ilimit%32) * (olimit%32) /(32*32) )/ (int)(((ilimit+31)/32) * ((olimit+31)/32));
f2 =  ((ilimitr/32) * (olimit/32) + (double)(ilimitr/32) * (olimit%32) /32+ (double)(ilimitr%32) * (olimit/32) /32 + (double)(ilimitr%32) * (olimit%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimit+31)/32)));
f3 =  ((ilimit/32) * (olimitr/32) + (double)(ilimit/32) * (olimitr%32) /32+ (double)(ilimit%32) * (olimitr/32) /32 + (double)(ilimit%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimit+31)/32) * ((olimitr+31)/32)));
f4 =  ((ilimitr/32) * (olimitr/32) + (double)(ilimitr/32) * (olimitr%32) /32+ (double)(ilimitr%32) * (olimitr/32) /32 + (double)(ilimitr%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((olimitr+31)/32) * ((olimitr+31)/32)));
int amax = blockA;
int bmax = blockB;
f = ((asize/amax) * (bsize/bmax) *f1 + (double)((asize/amax) * (bsize%bmax > 0) *f3)+ (double)((asize%amax > 0) * (bsize/bmax)*f2)  + (double)((asize%amax>0) * (bsize%bmax > 0) *f4) )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax));
return f;


}

#ifdef __cplusplus
extern "C" 
#endif
double getEfficiency_nomatchg32(int ilimit, int olimit)
{
	return -1;
	double f = ((ilimit/32) * (olimit/32) + (double)(ilimit/32) * (olimit%32) /32+ (double)(ilimit%32) * (olimit/32) /32 + (double)(ilimit%32) * (olimit%32) /(32*32) )/ (int)(((ilimit+31)/32) * ((olimit+31)/32));
	return f;
}

#ifdef __cplusplus
extern "C" 
#endif
double getEfficiency_matchg32(int ilimit, int olimit)
{
	return 1;
}

#ifdef __cplusplus
extern "C" 
#endif
double getEfficiency_matchl32(int size0, int asize, int bsize, int blockA)
{
double f1, f2, f3, f4, f;
	const int remainder1 = asize % blockA;
        const int remainder2 = bsize % blockA;
const int ilimit = remainder1 * size0;
        const int olimit = remainder2 * size0;
	const int plain = blockA * size0;
int minlimit = min(ilimit, olimit);

 //f1 =  ((plain/32)  + (double)(plain%32) /32)/ (int)((plain+31)/32);
 //f2 =  ((ilimit/32)  + (double)(ilimit%32) /32)/ (int)(max(1,(ilimit+31)/32));
 //f3 =  ((olimit/32)  + (double)(olimit%32) /32)/ (int)(max(1,(olimit+31)/32));
 //f4 =  ((minlimit/32)  + (double)(minlimit%32) /32)/ (int)(max(1,(minlimit+31)/32));
 f1 =  ((plain/32)  + (double)(plain%32) /32)/ (int)((plain+31)/32);
 f2 =  ((ilimit/32)  + (double)(ilimit%32) /32)/ (int)(max(1,(plain+31)/32));
 f3 =  ((olimit/32)  + (double)(olimit%32) /32)/ (int)(max(1,(plain+31)/32));
 f4 =  ((minlimit/32)  + (double)(minlimit%32) /32)/ (int)(max(1,(plain+31)/32));
int amax = blockA;
int bmax = blockA;
//printf("\t%lf %lf %lf %lf\t", f1, f2, f3, f4);
int n1, n2, n3, n4;
n1 = (asize/amax) * (bsize/bmax);
n2 = (asize%amax > 0 ) * (bsize/bmax);
n3 = (asize/amax) * (bsize%bmax>0);
n4 = (asize%amax > 0) * (bsize%bmax > 0);
//printf("\t%d %d %d %d\t", n1, n2, n3, n4);
//f = ((asize/amax) * (bsize/bmax) *f1 + (double)(asize/amax) * (bsize%bmax > 0) *f3+ (double)(asize%amax>0) * (bsize/bmax)*f2 + (double)(asize%amax > 0) * (bsize%bmax > 0) *f4 )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax));
f = (n1*f1 + n2*f2 + n3*f3 + n4*f4)/(n1+n2+n3+n4);
#ifdef MODEL
printf("\t%d\t%lf\t",blockA, f);
#endif
return f;
}
