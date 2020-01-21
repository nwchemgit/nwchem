/*
 $Id$
   Schrodinger.c - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Schodinger equation.

*/

#include <stdio.h>

#include "loggrid.h"
#include "pred_cor.h"
#include "schrodin.h"

#define Max_Iterations		500
#define tolerance 		1.0e-10
#define Corrector_Iterations	6

#define Min(x,y)	((x<y) ? x : y)
#define Max(x,y)	((x>y) ? x : y)


/**********************************
 *                                *
 *          R_Schrodinger         *
 *                                *
 **********************************/

void R_Schrodinger(n,l,v,mch,Eig,u,uprime)
int    n,l;
double v[];
int    *mch;
double *Eig;
double u[],
uprime[];

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
    r2,
    sum, a, scale, m1scale,
    uout,upout,upin,
    *r,
    *f_upp,
    *upp;

    /* define eigenvalues */
    E     = *Eig;
    gamma = ((double) (l+1));
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = N_LogGrid();
    log_amesh  = log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = (double *) r_LogGrid();
    f_upp = (double *) alloc_LogGrid();
    upp   = (double *) alloc_LogGrid();

    /* set up bounds for eigenvalue */
    Emax = v[Ngrid-1]  + 0.5*L2/(r[Ngrid-1]*r[Ngrid-1]);
    Emin = 0.0;
    for (i=0; i<Ngrid; ++i)
    {
        r2 = r[i];
        r2 = r2*r2;
        Emin = Min(Emin, (v[i] + 0.5*L2/r2));
    }
    if (E > Emax) E = 1.25*Emax;
    if (E < Emin) E = 0.75*Emin;
    if (E > Emax) E = 0.5*(Emax+Emin);

    iteration = 0;
    while (iteration < Max_Iterations)
    {
        ++iteration;
        /* define f_upp */
        for (i=0; i<Ngrid; ++i)
        {
            r2 = r[i];
            r2 = r2*r2;
            f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
        }

        /* find the classical turning point, */
        /* which is used for matching        */
        match = Ngrid-1;
        while (f_upp[match-1]*f_upp[match] > 0.0)
        {
            match = match - 1;

            if (match < 2)
            {
                printf("Error in R_Schrodinger: no turning point\n");
                return;
            }
        }


        /* set the boundry condition near zero */
        m1scale = 1.0;
        for (i=0; i<(n-l-1); ++i)
            m1scale *= -1.0;
        for (i=0; i<4; ++i)
        {
            u[i]      = m1scale*pow(r[i],gamma);
            uprime[i] = log_amesh*gamma*u[i];
        }
        for (i=0; i<4; ++i)
            upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];

        /* integrate from 0 to match */
        node = 0;
        for (i=3; i<match; ++i)
        {
            /* predictors */
            u[i+1]      = Predictor_Out(i,u,uprime);
            uprime[i+1] = Predictor_Out(i,uprime,upp);

            /* correctors */
            for (j=0; j<Corrector_Iterations; ++j)
            {
                upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
                uprime[i+1] =  Corrector_Out(i,uprime,upp);
                u[i+1]      =  Corrector_Out(i,u,uprime);
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

            /* find infinity */
            Ninf = match + floor(2.3/log_amesh);
            if ((Ninf+5) > Ngrid) Ninf = Ngrid - 5;

            /* define boundry near infinity */
            a = sqrt( L2/(r[Ninf]*r[Ninf]) + 2.0*(v[Ninf]-E) );
            for (i=Ninf; i<=(Ninf+4); ++i)
            {
                u[i]      = exp(-a*(r[i]-r[Ninf]));
                uprime[i] = -r[i]*log_amesh*a*u[i];
                upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];
            }

            /* integrate from infinity to match */
            for (i=Ninf; i>=(match+1); --i)
            {
                /* predictors */
                u[i-1]      = Predictor_In(i,u,uprime);
                uprime[i-1] = Predictor_In(i,uprime,upp);

                /* Correctors */
                for (j=0; j<Corrector_Iterations; ++j)
                {
                    upp[i-1]    = log_amesh*uprime[i-1] + f_upp[i-1]*u[i-1];
                    uprime[i-1] =  Corrector_In(i,uprime,upp);
                    u[i-1]      =  Corrector_In(i,u,uprime);
                }
            }

            /* make the outside u, match the inside u */
            scale = uout/u[match];
            for (i=match; i<=Ninf; ++i)
            {
                u[i]      = scale*u[i];
                uprime[i] = scale*uprime[i];
            }
            upin = uprime[match];

            /* Find Integral(u^2) */
            sum = Norm_LogGrid(Ninf,gamma,u);



            sum = 1.0/sqrt(sum);
            uout  = sum*uout;
            upout = sum*upout;
            upin  = sum*upin;
            for (i=0; i<=Ninf; ++i)
            {
                u[i]      = sum*u[i];
                uprime[i] = sum*uprime[i];
            }
            for (i=Ninf+1; i<Ngrid; ++i)
            {
                u[i]      = 0.0;
                uprime[i] = 0.0;
            }

            /* figure out new eigenvalue */
            de = 0.5*uout*(upout-upin)/(log_amesh*r[match]);

            /* eigenvalue is converged, exit */
            if (fabs(de) <  (Max(fabs(E),0.2)*tolerance))
            {
                *mch  = match;
                *Eig  = E;

                /* dealloc memory */
                dealloc_LogGrid(upp);
                dealloc_LogGrid(f_upp);
                return;
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

    printf("R_Schrodinger: More than %d iterations, Eig=%le, de=%le\n",Max_Iterations,E,de);

    match = Ngrid-1;
    R_Schrodinger_Fixed_Logderiv(n,l,v,Ngrid-1,0.0,Eig,u,uprime);
    R_Schrodinger_Fixed_E(n,l,v,Ngrid-1,*Eig,u,uprime);
    printf("R_Schrodinger: Running R_Schrodginger_Fixed_E with eig=%le match=%d\n",*Eig,match);

    *mch = match;
    *Eig = E;


    /* dealloc memory */
    dealloc_LogGrid(upp);
    dealloc_LogGrid(f_upp);

    return;

} /* R_Schrodinger */




/**********************************
 *                                *
 *     R_Schrodinger_Fixed_E      *
 *                                *
 **********************************/

void R_Schrodinger_Fixed_E(n,l,v,match,E,u,uprime)
int    n,l;
double v[];
int    match;
double E;
double u[],
uprime[];

{
    int     i,j,
    Ngrid;

    double  log_amesh,
    log_amesh2,
    gamma,
    L2,
    r2,
    sum,
    *r,
    *f_upp,
    *upp;

    /* define eigenvalues */
    gamma = ((double) (l+1));
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = N_LogGrid();
    log_amesh  = log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = (double *) r_LogGrid();
    f_upp = (double *) alloc_LogGrid();
    upp   = (double *) alloc_LogGrid();

    /* define f_upp */
    for (i=0; i<Ngrid; ++i)
    {
        r2 = r[i];
        r2 = r2*r2;
        f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
    }

    /* set the boundry condition near zero */
    for (i=0; i<4; ++i)
    {
        u[i]      = pow(r[i],gamma);
        uprime[i] = log_amesh*gamma*u[i];
    }
    for (i=0; i<4; ++i)
        upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];

    /* integrate from 0 to match */
    for (i=3; i<match; ++i)
    {
        /* predictors */
        u[i+1]      = Predictor_Out(i,u,uprime);
        uprime[i+1] = Predictor_Out(i,uprime,upp);

        /* correctors */
        for (j=0; j<Corrector_Iterations; ++j)
        {
            upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
            uprime[i+1] =  Corrector_Out(i,uprime,upp);
            u[i+1]      =  Corrector_Out(i,u,uprime);
        }

    }

    /* Find Integral(u^2) */
    sum = Norm_LogGrid(match,gamma,u);
    sum = 1.0/sqrt(sum);

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

    /* dealloc memory */
    dealloc_LogGrid(upp);
    dealloc_LogGrid(f_upp);

    return;

} /* R_Schrodinger */


/**********************************
 *                                *
 *  R_Schrodinger_Fixed_Logderiv  *
 *                                *
 **********************************/

void R_Schrodinger_Fixed_Logderiv(n,l,v,match,u_logderiv,Eig,u,uprime)
int    n,l;
double v[];
int	match;
double u_logderiv;
double *Eig;
double u[],
uprime[];

{
    int     i,j,
    iteration,
    node,
    Ngrid;

    double  E, de,
    Emax,
    Emin,
    log_amesh,
    log_amesh2,
    gamma,
    L2,
    r2,
    sum,
    uout,upout,upin,
    *r,
    *f_upp,
    *upp;

    /* define eigenvalues */
    E     = *Eig;
    gamma = ((double) (l+1));
    L2    = ((double) (l*(l+1)));

    /* define log grid parameters */
    Ngrid      = N_LogGrid();
    log_amesh  = log_amesh_LogGrid();
    log_amesh2 = log_amesh*log_amesh;

    /* get pointer rgrid, and extra memory */
    r     = (double *) r_LogGrid();
    f_upp = (double *) alloc_LogGrid();
    upp   = (double *) alloc_LogGrid();

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
            r2 = r[i];
            r2 = r2*r2;
            f_upp[i] = log_amesh2*(L2 + 2.0*(v[i] - E)*r2);
        }


        /* set the boundry condition near zero */
        for (i=0; i<4; ++i)
        {
            u[i]      = pow(r[i],gamma);
            uprime[i] = log_amesh*gamma*u[i];
        }
        for (i=0; i<4; ++i)
            upp[i]    = log_amesh*uprime[i] + f_upp[i]*u[i];

        /* integrate from 0 to match */
        node = 0;
        for (i=3; i<match; ++i)
        {
            /* predictors */
            u[i+1]      = Predictor_Out(i,u,uprime);
            uprime[i+1] = Predictor_Out(i,uprime,upp);

            /* correctors */
            for (j=0; j<Corrector_Iterations; ++j)
            {
                upp[i+1]    = log_amesh*uprime[i+1] + f_upp[i+1]*u[i+1];
                uprime[i+1] =  Corrector_Out(i,uprime,upp);
                u[i+1]      =  Corrector_Out(i,u,uprime);
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
            sum = Norm_LogGrid(match,gamma,u);



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
                dealloc_LogGrid(upp);
                dealloc_LogGrid(f_upp);
                return;
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

    printf("Error R_Schrodinger_Fixed_Logderiv: More than %d iterations\n",
           Max_Iterations);
    *Eig = E;


    /* dealloc memory */
    dealloc_LogGrid(upp);
    dealloc_LogGrid(f_upp);

    return;

} /* R_Schrodinger_Fixed_LogDeriv */


