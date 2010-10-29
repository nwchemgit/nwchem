/*
 $Id$
*/

/**********************************
 *				  *
 *            Norm		  *
 *				  *
 **********************************/

/* This routine calculates the Norm
   of a wavefunction assuming that
   the wavefunction decays like an
   exponential as r-> infinity

*/
double Norm(int M,double gamma,double *u)
{
    int	  i;
    double r0,sum,amesh,log_amesh,
    *r;

    amesh     = amesh_LogGrid();
    log_amesh = log_amesh_LogGrid();
    r         = r_LogGrid();

    /* Find Integral(u**2) */
    r0 = r[0]/sqrt(amesh);
    sum = pow(r0,(2.0*gamma+1.0))/(2.0*gamma+1.0);
    for (i=0; i<=(M-3); ++i)
        sum += log_amesh*r[i]*(u[i]*u[i]);

    sum += log_amesh*(  23.0*r[M-2]*(u[M-2]*u[M-2])
                        + 28.0*r[M-1]*(u[M-1]*u[M-1])
                        +  9.0*r[M]  *(u[M]  *u[M]))/24.0;

    return sum;

} /* Norm */

