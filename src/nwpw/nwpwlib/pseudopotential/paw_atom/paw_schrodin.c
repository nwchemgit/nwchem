/*
   $Id$
*/

/* paw_Schrodinger.c - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Schodinger equation.

*/
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "paw_my_constants.h"
#include "paw_loggrid.h"
#include "paw_pred_cor.h"
#include "paw_schrodin.h"

#define Max_Iterations		1000
#define tolerance 		1.0e-10
#define Corrector_Iterations	6

#define Min(x,y)	((x<y) ? x : y)
#define Max(x,y)	((x>y) ? x : y)


/**********************************
 *                                *
 *          paw_R_Schrodinger         *
 *                                *
 **********************************/

int paw_R_Schrodinger(int n,
                      int l,
                      double* v,
                      double *Eig,
                      double* u,
                      double* uprime)
{
    int     i,j,
    iteration,
    node,
    match,
    Ninf, Ngrid;

    double  E, de,
    Emax,
    Emin,
    log_amesh,
    log_amesh2,
    gamma,
    L2,
    sum, a, scale, m1scale,
    uout,upout,upin,
    *r,
    *f_upp,
    *r2,
    *upp;

    /* define eigenvalues */
    E     = *Eig;
    gamma = ((double) (l+1));
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = paw_N_LogGrid();
    log_amesh  = paw_log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = paw_r_LogGrid();
    r2    = paw_r2_LogGrid();
    f_upp = paw_alloc_LogGrid();
    upp   = paw_alloc_LogGrid();


    /* set up bounds for eigenvalue */
    Emin = 0.0;
    Emax = v[Ngrid-1]  + 0.5*L2/(r[Ngrid-1]*r[Ngrid-1]);

    for (i=0; i <= Ngrid-1; i++)
        Emin = Min(Emin, (v[i] + 0.5*L2/r2[i]));


    if (E > Emax) E = 1.25*Emax;
    if (E < Emin) E = 0.75*Emin;
    if (E > Emax) E = 0.5*(Emax+Emin);

    iteration = 0;
    while (iteration < Max_Iterations)
    {
        ++iteration;

        /* define f_upp */
        for (i=0; i<Ngrid; ++i)
            f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2[i]);

        /* classical turning point is used for matching  */
        match = Ngrid-1;
        while (f_upp[match-1]*f_upp[match] > 0.0)
        {
            match = match - 1;

            if (match < 2)
            {
                printf("Error in paw_R_Schrodinger: no turning point\n");
                return False;
            }
        }


        /* set the boundry condition near zero */
        m1scale = 1.0;
        for (i=0; i<(n-l-1); ++i)
            m1scale *= -1.0;
        for (i=0; i<=3; i++)
        {
            u[i]      = m1scale*pow(r[i],gamma);
            uprime[i] = log_amesh*gamma*u[i];
            upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
        }

        /* integrate from 0 to match */
        node = 0;
        for (i=3; i <= match-1; i++)
        {
            /* predictors */
            u[i+1]      = paw_Predictor_Out(i,u,uprime);
            uprime[i+1] = paw_Predictor_Out(i,uprime,upp);

            /* correctors */
            for (j=0; j<= Corrector_Iterations-1; j++)
            {
                upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
                uprime[i+1] =  paw_Corrector_Out(i,uprime,upp);
                u[i+1]      =  paw_Corrector_Out(i,u,uprime);
            }

            /* finding nodes */
            if (u[i+1]*u[i] <= 0) node = node + 1;
        }

        uout  = u[match];
        upout = uprime[match];

        /* not enough nodes in u */
        if ((node-n+l+1) < 0)
        {
            Emin = E;
            E    = 0.5*(Emin+Emax);
        }
        /* too many nodes in u */
        else if ((node-n+l+1) > 0)
        {
            Emax = E;
            E    = 0.5*(Emin+Emax);
        }
        /* number of nodes ok, start integration inward */
        else
        {

            /* find infinity

            */
            /* find infinity */
            Ninf = match + floor(2.3/log_amesh);
            if ((Ninf+5) > Ngrid) Ninf = Ngrid - 5;

            /* define boundry near infinity */
            a = sqrt( L2/(r[Ninf]*r[Ninf]) + 2.0*(v[Ninf]-E) );
            for (i=Ninf; i<=(Ninf+4); i++)
            {
                u[i]      = exp(-a*(r[i]-r[Ninf]));
                uprime[i] = -r[i]*log_amesh*a*u[i];
                upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
            }

            /* integrate from infinity to match */
            for (i=Ninf; i>=(match+1); --i)
            {
                /* predictors */
                u[i-1]      = paw_Predictor_In(i,u,uprime);
                uprime[i-1] = paw_Predictor_In(i,uprime,upp);

                /* Correctors */
                for (j=0; j<Corrector_Iterations; ++j)
                {
                    upp[i-1]    = log_amesh*uprime[i-1] + f_upp[i-1]*u[i-1];
                    uprime[i-1] =  paw_Corrector_In(i,uprime,upp);
                    u[i-1]      =  paw_Corrector_In(i,u,uprime);
                }
            }

            /* make the outside u, match the inside u */
            scale = uout/u[match];
            for (i=match; i<= Ninf; i++)
            {
                u[i]      = scale*u[i];
                uprime[i] = scale*uprime[i];
            }
            upin = uprime[match];

            /* Find Integral(u^2) */
            sum = paw_Norm_LogGrid(Ninf,gamma,u);


            sum = 1.0/sqrt(sum);
            uout  = sum*uout;
            upout = sum*upout;
            upin  = sum*upin;

            for (i=0; i<=Ninf; i++)
            {
                u[i]      = sum*u[i];
                uprime[i] = sum*uprime[i];
            }


            /* figure out new eigenvalue */
            de = 0.5*uout*(upout-upin)/(log_amesh*r[match]);

            /* eigenvalue is converged, exit */
            if (fabs(de) <  (Max(fabs(E),0.2)*tolerance))
            {
                *Eig  = E;

                for (i=Ninf+1; i<Ngrid; ++i)
                {

                    u[i]      = u[Ninf]*exp(-a*(r[i]-r[Ninf]));
                    uprime[i] = -r[i]*log_amesh*a*u[i];

                }

                /* dealloc memory */
                paw_dealloc_LogGrid(upp);
                paw_dealloc_LogGrid(f_upp);
                return True;
            }

            if (de > 0.0)
                Emin = E;
            else
                Emax = E;
            E = E + de;
            if (Emax<=Emin)
            {

                /* dealloc memory */
                paw_dealloc_LogGrid(upp);
                paw_dealloc_LogGrid(f_upp);
                *Eig = E;
                return False;

            }

            if ( (E > Emax) || (E < Emin))
                E = 0.5*(Emin+Emax);

        } /* nodes ok */
    } /* while */

    printf("Error paw_R_Schrodinger: More than %d iterations\n",Max_Iterations);
    *Eig = E;

    /* dealloc memory */
    paw_dealloc_LogGrid(upp);
    paw_dealloc_LogGrid(f_upp);

    return False;

} /* paw_R_Schrodinger */




/**********************************
 *                                *
 *     paw_R_Schrodinger_Fixed_E      *
 *                                *
 **********************************/

int paw_R_Schrodinger_Fixed_E(
    int l,
    double *v,
    int match,
    double E,
    double *u,
    double *uprime
)
{
    int     i,j,
    Ngrid;

    double  log_amesh,
    log_amesh2,
    gamma,
    L2,
    r2,
    *r,
    *f_upp,
    *upp;

    double max_u;


    /*set maximum allowed value for u*/
    max_u = 100000000;

    /* define power of the orbital near 0 */
    gamma = ((double) (l+1));

    /* define square of the angular momentum */
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = paw_N_LogGrid();
    log_amesh  = paw_log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    if (match > Ngrid-1)
    {
        printf("match point is outside the allowed region\n");
        printf("abborting the program\n");
        exit(1);
    }
    /* get pointer rgrid, and extra memory */
    r     = paw_r_LogGrid();
    f_upp = paw_alloc_LogGrid();
    upp   = paw_alloc_LogGrid();

    for (i=0; i<=match; i++)
    {
        u[i]      = 0.0;
        uprime[i] = 0.0;
    }

    /* define f_upp */
    for (i=0; i<match+1; ++i)
    {
        r2 = r[i]*r[i];
        f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
    }

    for (i=0; i<4; ++i)
    {
        u[i]      = pow(r[i],gamma);
        uprime[i] = log_amesh*gamma*u[i];
        upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
    }

    /* integrate from 0 to match */
    i=3;
    while (i<match && fabs(u[i])<max_u)
    {
        /* predictors */
        u[i+1]      = paw_Predictor_Out(i,u,uprime);
        uprime[i+1] = paw_Predictor_Out(i,uprime,upp);

        /* correctors */
        for (j=0; j<Corrector_Iterations; ++j)
        {
            upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
            uprime[i+1] =  paw_Corrector_Out(i,uprime,upp);
            u[i+1]      =  paw_Corrector_Out(i,u,uprime);
        }

        i = i+1;
    }


    /* dealloc memory */
    paw_dealloc_LogGrid(upp);
    paw_dealloc_LogGrid(f_upp);

    return True;

} /* paw_R_Schrodinger */


/****************************************
 Function name	  : paw_R_Schrodinger_Fixed_Logderiv
 Description	    :

 Return type		  : void

 Argument         : int n
 Argument         : int l
 Argument         : double *v
 Argument         : int match
 Argument         : double u_logderiv
 Argument         : double *Eig
 Argument         : double *u
 Argument         : double *uprime

 Author     		  : Eric Bylaska & Marat Valiev
 Date & Time		  : 1/8/99 2:04:28 PM
****************************************/
int paw_R_Schrodinger_Fixed_Logderiv(
    int n,
    int l,
    double *v,
    int match,
    double u_logderiv,
    double *Eig,
    double *u,
    double *uprime
)

{
    int     i;
    int     j;
    int     iteration;
    int     node;
    int     Ngrid;

    double  sum;
    double  E;
    double  de;
    double  Emax;
    double  Emin;
    double  log_amesh;
    double  log_amesh2;
    double  gamma;
    double  L2;
    double  r2;
    double  uout;
    double  upout;
    double  upin;
    double  *r;
    double  *f_upp;
    double  *upp;

    /* define eigenvalues */
    E     = *Eig;

    /* define power of the orbital near 0 */
    gamma = ((double) (l+1));

    /* define square of the angular momentum */
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = paw_N_LogGrid();
    log_amesh  = paw_log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = (double *) paw_r_LogGrid();
    f_upp = (double *) paw_alloc_LogGrid();
    upp   = (double *) paw_alloc_LogGrid();

    /* set up bounds for eigenvalue */
    Emax = E + 10.0;
    Emin = E - 10.0;


    iteration = 0;
    while (iteration < Max_Iterations)
    {
        ++iteration;
        /* define f_upp */
        for (i=0; i<Ngrid; ++i)
        {
            r2 = r[i]*r[i];
            f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
        }


        /* set the boundry condition near zero */
        for (i=0; i<4; ++i)
        {
            u[i]      = pow(r[i],gamma);
            uprime[i] = log_amesh*gamma*u[i];
            upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
        }


        /* integrate from 0 to match */
        node = 0;
        for (i=3; i<match; ++i)
        {
            /* predictors */
            u[i+1]      = paw_Predictor_Out(i,u,uprime);
            uprime[i+1] = paw_Predictor_Out(i,uprime,upp);

            /* correctors */
            for (j=0; j<Corrector_Iterations; ++j)
            {
                upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
                uprime[i+1] =  paw_Corrector_Out(i,uprime,upp);
                u[i+1]      =  paw_Corrector_Out(i,u,uprime);
            }

            /* finding nodes */
            if (u[i+1]*u[i] <= 0) node = node + 1;
        }


        uout  = u[match];
        upout = uprime[match];


        /* not enough nodes in u */
        if ((node-n+l+1) < 0)
        {
            Emin = E;
            E    = 0.5*(Emin+Emax);
        }

        /* too many nodes in u */
        else if ((node-n+l+1) > 0)
        {
            Emax = E;
            E    = 0.5*(Emin+Emax);
        }

        /* number of nodes ok, adjust e to log derivative */
        else
        {
            upin = u_logderiv*uout;

            /* Find Integral(u^2) */
            sum = paw_Norm_LogGrid(match,gamma,u);

            sum = 1.0/sqrt(sum);
            uout  = sum*uout;
            upout = sum*upout;
            upin  = sum*upin;

            for (i=0; i<=match; ++i)
            {
                u[i]      = sum*u[i];
                uprime[i] = sum*uprime[i];
            }

            for (i=match+1; i<Ngrid; ++i)
            {
                u[i]      = 0.0;
                uprime[i] = 0.0;
            }

            /* figure out new eigenvalue */
            de = 0.5*uout*(upout-upin)/(log_amesh*r[match]);

            /* eigenvalue is converged, exit */
            if (fabs(de) <  (Max(fabs(E),0.2)*tolerance))
            {
                *Eig  = E;

                /* dealloc memory */
                paw_dealloc_LogGrid(upp);
                paw_dealloc_LogGrid(f_upp);
                return True;
            }

            if (de > 0.0)
                Emin = E;
            else
                Emax = E;

            E = E + de;

            if ( (E > Emax) || (E < Emin))
                E = 0.5*(Emin+Emax);

        } /* nodes ok */
    } /* while */

    printf("paw_Error R_Schrodinger_Fixed_Logderiv: More than %d iterations\n",
           Max_Iterations);
    *Eig = E;


    /* dealloc memory */
    paw_dealloc_LogGrid(upp);
    paw_dealloc_LogGrid(f_upp);

    return False;

} /* paw_R_Schrodinger_Fixed_LogDeriv */


int paw_R_Schrodinger_Fixed_Logderiv1(
    int n,
    int l,
    double *v,
    int match,
    double u_logderiv,
    double *Eig,
    double *u,
    double *uprime
)

{
    int     i;
    int     j;
    int     iteration;
    int     Ngrid;

    double  sum;
    double  E;
    double  de;
    double  Emax;
    double  Emin;
    double  log_amesh;
    double  log_amesh2;
    double  gamma;
    double  L2;
    double  r2;
    double  uout;
    double  upout;
    double  upin;
    double  *r;
    double  *f_upp;
    double  *upp;

    /* define eigenvalues */
    E     = *Eig;

    /* define power of the orbital near 0 */
    gamma = ((double) (l+1));

    /* define square of the angular momentum */
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = paw_N_LogGrid();
    log_amesh  = paw_log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = (double *) paw_r_LogGrid();
    f_upp = (double *) paw_alloc_LogGrid();
    upp   = (double *) paw_alloc_LogGrid();

    /* set up bounds for eigenvalue */
    Emax = E + 10.0;
    Emin = E - 10.0;


    iteration = 0;
    while (iteration < Max_Iterations)
    {
        ++iteration;
        /* define f_upp */
        for (i=0; i<Ngrid; ++i)
        {
            r2 = r[i]*r[i];
            f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
        }


        /* set the boundry condition near zero */
        for (i=0; i<4; ++i)
        {
            u[i]      = pow(r[i],gamma);
            uprime[i] = log_amesh*gamma*u[i];
            upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
        }


        /* integrate from 0 to match */
        for (i=3; i<match; ++i)
        {
            /* predictors */
            u[i+1]      = paw_Predictor_Out(i,u,uprime);
            uprime[i+1] = paw_Predictor_Out(i,uprime,upp);

            /* correctors */
            for (j=0; j<Corrector_Iterations; ++j)
            {
                upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
                uprime[i+1] =  paw_Corrector_Out(i,uprime,upp);
                u[i+1]      =  paw_Corrector_Out(i,u,uprime);
            }

        }


        uout  = u[match];
        upout = uprime[match];


        upin = u_logderiv*uout;

        /* Find Integral(u^2) */
        sum = paw_Norm_LogGrid(match,gamma,u);

        sum = 1.0/sqrt(sum);
        uout  = sum*uout;
        upout = sum*upout;
        upin  = sum*upin;

        for (i=0; i<=match; ++i)
        {
            u[i]      = sum*u[i];
            uprime[i] = sum*uprime[i];
        }

        for (i=match+1; i<Ngrid; ++i)
        {
            u[i]      = 0.0;
            uprime[i] = 0.0;
        }

        /* figure out new eigenvalue */
        de = 0.5*uout*(upout-upin)/(log_amesh*r[match]);

        /* eigenvalue is converged, exit */
        if (fabs(de) <  (Max(fabs(E),0.2)*tolerance))
        {
            *Eig  = E;

            /* dealloc memory */
            paw_dealloc_LogGrid(upp);
            paw_dealloc_LogGrid(f_upp);
            return True;
        }

        if (de > 0.0)
            Emin = E;
        else
            Emax = E;

        E = E + de;

        if ( (E > Emax) || (E < Emin))
            E = 0.5*(Emin+Emax);

    } /* while */

    printf("paw_Error R_Schrodinger_Fixed_Logderiv: More than %d iterations\n",
           Max_Iterations);
    *Eig = E;


    /* dealloc memory */
    paw_dealloc_LogGrid(upp);
    paw_dealloc_LogGrid(f_upp);

    return False;

} /* paw_R_Schrodinger_Fixed_LogDeriv */


void paw_R_Schrodinger_Fixed_E1(
    int l,
    double *v,
    double *f,
    int match,
    double E,
    double *u,
    double *uprime
)
{
    int     i,j,
    Ngrid;

    double  log_amesh,
    log_amesh2,
    gamma,
    L2,
    r2,
    *r,
    *f_upp,
    *upp;

    double max_u;

    /*set maximum allowed value for u*/
    max_u = 10000000;

    /* define power of the orbital near 0 */
    gamma = ((double) (l+1));

    /* define square of the angular momentum */
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = paw_N_LogGrid();
    log_amesh  = paw_log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    if (match > Ngrid-1)
    {
        printf("match point is outside the allowed region\n");
        printf("abborting the program\n");
        exit(1);
    }
    /* get pointer rgrid, and extra memory */
    r     = (double *) paw_r_LogGrid();
    f_upp = (double *) paw_alloc_LogGrid();
    upp   = (double *) paw_alloc_LogGrid();

    for (i=0; i<=match; i++)
    {
        u[i]      = 0.0;
        uprime[i] = 0.0;
    }

    /* define f_upp */
    for (i=0; i<match+1; ++i)
    {
        r2 = r[i]*r[i];
        f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
    }

    for (i=0; i<4; ++i)
    {
        u[i]      = pow(r[i],gamma);
        uprime[i] = log_amesh*gamma*u[i];
        upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
    }

    /* integrate from 0 to match */
    i=3;
    while (i<match && fabs(u[i])<max_u)
    {
        /* predictors */
        u[i+1]      = paw_Predictor_Out(i,u,uprime);
        uprime[i+1] = paw_Predictor_Out(i,uprime,upp);

        /* correctors */
        for (j=0; j<Corrector_Iterations; ++j)
        {
            upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1]- 2.0*log_amesh2*r[i+1]*r[i+1]*f[i+1];
            uprime[i+1] =  paw_Corrector_Out(i,uprime,upp);
            u[i+1]      =  paw_Corrector_Out(i,u,uprime);
        }

        i = i+1;
    }


    /* dealloc memory */
    paw_dealloc_LogGrid(upp);
    paw_dealloc_LogGrid(f_upp);

    return;

} /* paw_R_Schrodinger */

