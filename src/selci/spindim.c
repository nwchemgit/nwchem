/*
 $Id$
 */

#include <stdio.h>
#include <stdlib.h>

/*
  spindim multiplicity N

  Compute the dimension of the spin degeneracy for N open-shells
  coupled to the specified multiplicity
  */

double factorial(int a)
{
    double f=1.0;
    int n = a;


    while (a-- > 0)
	f *= (a+1);

    printf("%d! = %g\n",n, f);

    return f;
}

double C(int a, int b)
/*
  return aCb = a!/((a-b)! b!)
  */
{
    double c;

    if (a>=0 && b>=0)
	c = factorial(a)/(factorial(a-b)*factorial(b));
    else
	c = 0;
    printf("%dC%d = %g\n", a, b, c);
    return c;
}


int main(int argc, char **argv)
{
    int multi = atoi(argv[1]);
    int n = atoi(argv[2]);
    int dim = (int) (C(n,(n-multi+1)/2) - C(n,(n-multi-1)/2));

    printf("multi=%d nelec=%d dim=%d\n", multi, n, dim);

    return 0;
}
