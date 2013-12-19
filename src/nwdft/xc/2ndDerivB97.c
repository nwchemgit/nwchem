//provided by 
//Tobias Risthaus
//tobias.risthaus@thch.uni-bonn.de
//Mulliken Center for Theoretical Chemistry
//Institut für Physikalische und Theoretische Chemie
//Universität Bonn, Beringstr. 4, D-53115 Bonn
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "2ndDerivB97.h"
void dft_xckernel_pwlda_(double* , double* , double* );

void dft_xckernel_pwlda(double  ra, double rb , double FCLDA[_FCLDA_ELEMENTS]){
  dft_xckernel_pwlda_(&ra, &rb, FCLDA);
};


static void gss3tar(double r, double g, double alpha, double *a, 
                    double *f, double *fr, double *fg, 
                    double *frr, double *frg, double *fgg)
{
  // B97 like 
  // S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
  // doi: 10.1002/jcc.20495
  // for now, this is only valid for 3 terms. 

  double r83,r113,r143,s2,s2p1,u;

  r83=pow(r,8.0/3.0);
  r113=r83*r;
  r143=r113*r;
  s2=g/r83;
  s2p1=1.0+alpha*s2;

  u=alpha*s2/s2p1;

  double D1uDs21= alpha/( s2p1 * s2p1 );                         // first  derivative of u wrt s2
  double D2uDs22= -2.0*alpha*alpha/( s2p1 * s2p1 * s2p1 );       // second derivative of u wrt s2

  double D1s2Dr1= ( -8.0/3.0)*g/r113;                            // first  derivative of s2 wrt r 
  double D2s2Dr2= (+88.0/9.0)*g/r143;

  double D1s2Dg1= 1.0/r83;                                       // first  derivative of s2 wrt g

  double D2s2Dr1g1= (-8.0/3.0)/r113;

  double D1uDr1=D1uDs21*D1s2Dr1;                                 // first  derivative of u wrt r
  double D1uDg1=D1uDs21*D1s2Dg1;

  double D2uDr2=D2uDs22*D1s2Dr1*D1s2Dr1 + D1uDs21*D2s2Dr2;
  double D2uDg2=D2uDs22*D1s2Dg1*D1s2Dg1; // + D1uDs21*D2s2Dg2; // second factor always zero, see above 

  double D2uDr1g1=D2uDs22*D1s2Dr1*D1s2Dg1 + D1uDs21*D2s2Dr1g1;   // mixed second derivative of u (wrt r,g)

  *f=a[0] + a[1]*u + a[2]*u*u;

  *fr=0.0 + a[1]*D1uDr1 + a[2]*2.0*u*D1uDr1;

  *fg=0.0 + a[1]*D1uDg1 + a[2]*2.0*u*D1uDg1;

  *frr=0.0 + a[1]*D2uDr2 + a[2]*(2.0*D1uDr1*D1uDr1 + 2.0*u*D2uDr2);

  *fgg=0.0 + a[1]*D2uDg2 + a[2]*(2.0*D1uDg1*D1uDg1 + 2.0*u*D2uDg2);

  *frg=0.0 + a[1]*D2uDr1g1 + a[2]*(2.0*D1uDr1*D1uDg1 + 2.0*u*D2uDr1g1);

  return;
};

static void gss5tar(double r, double g, double alpha, double *a, 
                    double *f, double *fr, double *fg, 
                    double *frr, double *frg, double *fgg)
{
  // B97 like 
  // S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
  // doi: 10.1002/jcc.20495
  // for now, this is only valid for 5 terms. 

  double r83,r113,r143,s2,s2p1,u;

  r83=pow(r,8.0/3.0);
  r113=r83*r;
  r143=r113*r;
  s2=g/r83;
  s2p1=1.0+alpha*s2;

  u=alpha*s2/s2p1;

  double D1uDs21= alpha/( s2p1 * s2p1 );
  double D2uDs22= -2.0*alpha*alpha/( s2p1 * s2p1 * s2p1 );

  double D1s2Dr1= ( -8.0/3.0)*g/r113;
  double D2s2Dr2= (+88.0/9.0)*g/r143;

  double D1s2Dg1= 1.0/r83;

  double D2s2Dr1g1= (-8.0/3.0)/r113;

  double D1uDr1=D1uDs21*D1s2Dr1;
  double D1uDg1=D1uDs21*D1s2Dg1;

  double D2uDr2=D2uDs22*D1s2Dr1*D1s2Dr1 + D1uDs21*D2s2Dr2;
  double D2uDg2=D2uDs22*D1s2Dg1*D1s2Dg1; // + D1uDs21*D2s2Dg2; // second factor always zero, see above 

  double D2uDr1g1=D2uDs22*D1s2Dr1*D1s2Dg1 + D1uDs21*D2s2Dr1g1; 

  *f=a[0] + a[1]*u + a[2]*u*u + a[3]*u*u*u + a[4]*u*u*u*u;

  *fr=0.0 + a[1]*D1uDr1 + a[2]*2.0*u*D1uDr1 + a[3]*3.0*u*u*D1uDr1 + a[4]*4.0*u*u*u*D1uDr1;

  *fg=0.0 + a[1]*D1uDg1 + a[2]*2.0*u*D1uDg1 + a[3]*3.0*u*u*D1uDg1 + a[4]*4.0*u*u*u*D1uDg1;

  *frr=0.0 + a[1]*D2uDr2 + a[2]*2.0*(D1uDr1*D1uDr1 + u*D2uDr2)
    + a[3]*3.0*(2.0*u*D1uDr1*D1uDr1 + u*u*D2uDr2)
    + a[4]*4.0*(3.0*u*u*D1uDr1*D1uDr1 + u*u*u*D2uDr2);

  *fgg=0.0 + a[1]*D2uDg2 + a[2]*2.0*(D1uDg1*D1uDg1 + u*D2uDg2)
    + a[3]*3.0*(2.0*u*D1uDg1*D1uDg1 + u*u*D2uDg2)
    + a[4]*4.0*(3.0*u*u*D1uDg1*D1uDg1 + u*u*u*D2uDg2);

  *frg=0.0 + a[1]*D2uDr1g1 + a[2]*2.0*(D1uDr1*D1uDg1 + u*D2uDr1g1)
    + a[3]*3.0*(2.0*u*D1uDr1*D1uDg1 + u*u*D2uDr1g1)
    + a[4]*4.0*(3.0*u*u*D1uDr1*D1uDg1 + u*u*u*D2uDr1g1);

  return;
};

static void gab3tar(double ra, double rb, double ga, double gb, double alpha, double *a, 
                    double  *f, double *fra, double *frb, double *fga, double *fgb,
                    double *frara, double *frarb, double *frbrb,
                    double *fraga, double *fragb, double *frbga, double *frbgb,
                    double *fgaga, double *fgagb, double *fgbgb)
{
  // B97 like 
  // S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
  // doi: 10.1002/jcc.20495
  // for now, this is only valid for 3 terms. 

  double ra83,rb83,ra113,rb113,ra143,rb143,sa2,sb2,sv2,sv2p1,u;

  ra83=pow(ra,8.0/3.0);
  ra113=ra83*ra;
  ra143=ra113*ra;
  rb83=pow(rb,8.0/3.0);
  rb113=rb83*rb;
  rb143=rb113*rb;

  sa2=ga/ra83;
  sb2=gb/rb83;
  sv2=0.5*(sa2 + sb2);                  // called s_av in some papers 
  sv2p1=1.0+alpha*sv2;

  u=alpha*sv2/sv2p1;
  *f=a[0] + a[1]*u + a[2]*u*u;

  double D1uDsv21= alpha/( sv2p1 * sv2p1 );
  double D2uDsv22= -2.0*alpha*alpha/( sv2p1 * sv2p1 * sv2p1 );


  double D1sv2Dra1= (-4.0/3.0)*ga/ra113;
  double D2sv2Dra2= (44.0/9.0)*ga/ra143;

  double D1sv2Drb1= (-4.0/3.0)*gb/rb113;
  double D2sv2Drb2= (44.0/9.0)*gb/rb143;

  double D2sv2Dra1ga1 = (-4.0/3.0)/ra113;

  double D1sv2Dga1= 0.5/ra83;
  double D2sv2Dga2= 0.0;

  double D1sv2Dgb1= 0.5/rb83;
  double D2sv2Dgb2= 0.0;

  double D2sv2Drb1gb1 = (-4.0/3.0)/rb113;

  double D2sv2Dra1rb1 = 0.0; 
  double D2sv2Dga1gb1 = 0.0; 
  double D2sv2Drb1ga1 = 0.0;
  double D2sv2Dra1gb1 = 0.0;

  double D1uDra1=D1uDsv21*D1sv2Dra1;

  double D2uDra2=D2uDsv22*D1sv2Dra1*D1sv2Dra1 + D1uDsv21*D2sv2Dra2;

  double D2uDra1ga1=D2uDsv22*D1sv2Dra1*D1sv2Dga1 + D1uDsv21*D2sv2Dra1ga1;
  double D2uDra1rb1=D2uDsv22*D1sv2Dra1*D1sv2Drb1 + D1uDsv21*D2sv2Dra1rb1;
  double D2uDra1gb1=D2uDsv22*D1sv2Dra1*D1sv2Dgb1 + D1uDsv21*D2sv2Dra1gb1;


  double D1uDrb1=D1uDsv21*D1sv2Drb1;

  double D2uDrb2=D2uDsv22*D1sv2Drb1*D1sv2Drb1 + D1uDsv21*D2sv2Drb2;
  double D2uDrb1gb1=D2uDsv22*D1sv2Drb1*D1sv2Dgb1 + D1uDsv21*D2sv2Drb1gb1;
  // D2uDrb1ra1 exists already as D2uDra1rb1
  double D2uDrb1ga1=D2uDsv22*D1sv2Drb1*D1sv2Dga1 + D1uDsv21*D2sv2Drb1ga1;     // mixed second derivative of u wrt to rho_beta and gaa  


  double D1uDga1=D1uDsv21*D1sv2Dga1;
  double D1uDgb1=D1uDsv21*D1sv2Dgb1;

  double D2uDga2=D2uDsv22*D1sv2Dga1*D1sv2Dga1 + D1uDsv21*D2sv2Dga2;
  double D2uDgb2=D2uDsv22*D1sv2Dgb1*D1sv2Dgb1 + D1uDsv21*D2sv2Dgb2;

  double D2uDga1gb1=D2uDsv22*D1sv2Dga1*D1sv2Dgb1 + D1uDsv21*D2sv2Dga1gb1;

  // D2uDga1gb1 exists

  *fra=0.0 + a[1]*D1uDra1 + a[2]*2.0*u*D1uDra1;
  *frb=0.0 + a[1]*D1uDrb1 + a[2]*2.0*u*D1uDrb1;

  *fga=0.0 + a[1]*D1uDga1 + a[2]*2.0*u*D1uDga1;

  *fgb=0.0 + a[1]*D1uDgb1 + a[2]*2.0*u*D1uDgb1;

  *frara=a[1]*D2uDra2 + 2.0*a[2]*(D1uDra1*D1uDra1 + u*D2uDra2);

  *fraga=a[1]*D2uDra1ga1 + 2.0*a[2]*(D1uDga1*D1uDra1 + u*D2uDra1ga1); 
  *frarb=a[1]*D2uDra1rb1 + 2.0*a[2]*(D1uDrb1*D1uDra1 + u*D2uDra1rb1);
  *fragb=a[1]*D2uDra1gb1 + 2.0*a[2]*(D1uDgb1*D1uDra1 + u*D2uDra1gb1);

  *frbrb=a[1]*D2uDrb2 + 2.0*a[2]*(D1uDrb1*D1uDrb1 + u*D2uDrb2);
  // frbra exists
  *frbga=a[1]*D2uDrb1ga1 + 2.0*a[2]*(D1uDga1*D1uDrb1 + u*D2uDrb1ga1);
  *frbgb=a[1]*D2uDrb1gb1 + 2.0*a[2]*(D1uDgb1*D1uDrb1 + u*D2uDrb1gb1);

  *fgaga=a[1]*D2uDga2 + 2.0*a[2]*(D1uDga1*D1uDga1 + u*D2uDga2);
  *fgbgb=a[1]*D2uDgb2 + 2.0*a[2]*(D1uDgb1*D1uDgb1 + u*D2uDgb2);

  *fgagb=a[1]*D2uDga1gb1 + 2.0*a[2]*(D1uDga1*D1uDgb1 + u*D2uDga1gb1);

  return;
};


static void gab5tar(double ra, double rb, double ga, double gb, double alpha, double *a, 
                    double *f, double *fra, double *frb, double *fga, double *fgb,
                    double *frara, double *frarb, double *frbrb,
                    double *fraga, double *fragb, double *frbga, double *frbgb,
                    double *fgaga, double *fgagb, double *fgbgb)
{
  // B97 like 
  // S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
  // doi: 10.1002/jcc.20495
  // for now, this is only valid for 5 terms. 

  double ra83,rb83,ra113,rb113,ra143,rb143,sa2,sb2,sv2,sv2p1,u;

  ra83=pow(ra,8.0/3.0);
  ra113=ra83*ra;
  ra143=ra113*ra;
  rb83=pow(rb,8.0/3.0);
  rb113=rb83*rb;
  rb143=rb113*rb;

  sa2=ga/ra83;
  sb2=gb/rb83;
  sv2=0.5*(sa2 + sb2);
  sv2p1=1.0+alpha*sv2;

  u=alpha*sv2/sv2p1;

  *f=a[0] + a[1]*u + a[2]*u*u + a[3]*u*u*u + a[4]*u*u*u*u;

  double D1uDsv21= alpha/( sv2p1 * sv2p1 );
  double D2uDsv22= -2.0*alpha*alpha/( sv2p1 * sv2p1 * sv2p1 );

  double D1sv2Dra1= (-4.0/3.0)*ga/ra113;
  double D2sv2Dra2= (44.0/9.0)*ga/ra143;

  double D1sv2Drb1= (-4.0/3.0)*gb/rb113;
  double D2sv2Drb2= (44.0/9.0)*gb/rb143;

  double D2sv2Dra1ga1 = (-4.0/3.0)/ra113;

  double D1sv2Dga1= 0.5/ra83;
  double D2sv2Dga2= 0.0;

  double D1sv2Dgb1= 0.5/rb83;
  double D2sv2Dgb2= 0.0;

  double D2sv2Drb1gb1 = (-4.0/3.0)/rb113;

  double D2sv2Dra1rb1 = 0.0; 
  double D2sv2Dga1gb1 = 0.0; 
  double D2sv2Drb1ga1 = 0.0;
  double D2sv2Dra1gb1 = 0.0;

  double D1uDra1=D1uDsv21*D1sv2Dra1;

  double D2uDra2=D2uDsv22*D1sv2Dra1*D1sv2Dra1 + D1uDsv21*D2sv2Dra2;

  double D2uDra1ga1=D2uDsv22*D1sv2Dra1*D1sv2Dga1 + D1uDsv21*D2sv2Dra1ga1;
  double D2uDra1rb1=D2uDsv22*D1sv2Dra1*D1sv2Drb1 + D1uDsv21*D2sv2Dra1rb1;
  double D2uDra1gb1=D2uDsv22*D1sv2Dra1*D1sv2Dgb1 + D1uDsv21*D2sv2Dra1gb1;


  double D1uDrb1=D1uDsv21*D1sv2Drb1;

  double D2uDrb2=D2uDsv22*D1sv2Drb1*D1sv2Drb1 + D1uDsv21*D2sv2Drb2;
  double D2uDrb1gb1=D2uDsv22*D1sv2Drb1*D1sv2Dgb1 + D1uDsv21*D2sv2Drb1gb1;
  // D2uDrb1ra1 exists already as D2uDra1rb1
  double D2uDrb1ga1=D2uDsv22*D1sv2Drb1*D1sv2Dga1 + D1uDsv21*D2sv2Drb1ga1;


  double D1uDga1=D1uDsv21*D1sv2Dga1;
  double D1uDgb1=D1uDsv21*D1sv2Dgb1;

  double D2uDga2=D2uDsv22*D1sv2Dga1*D1sv2Dga1 + D1uDsv21*D2sv2Dga2;
  double D2uDgb2=D2uDsv22*D1sv2Dgb1*D1sv2Dgb1 + D1uDsv21*D2sv2Dgb2;

  double D2uDga1gb1=D2uDsv22*D1sv2Dga1*D1sv2Dgb1 + D1uDsv21*D2sv2Dga1gb1;

  // D2uDga1gb1 exists

  *fra=0.0 + a[1]*D1uDra1 + a[2]*2.0*u*D1uDra1 + a[3]*3.0*u*u*D1uDra1 + a[4]*4.0*u*u*u*D1uDra1;
  *frb=0.0 + a[1]*D1uDrb1 + a[2]*2.0*u*D1uDrb1 + a[3]*3.0*u*u*D1uDrb1 + a[4]*4.0*u*u*u*D1uDrb1;

  *fga=0.0 + a[1]*D1uDga1 + a[2]*2.0*u*D1uDga1 + a[3]*3.0*u*u*D1uDga1 + a[4]*4.0*u*u*u*D1uDga1;
  *fgb=0.0 + a[1]*D1uDgb1 + a[2]*2.0*u*D1uDgb1 + a[3]*3.0*u*u*D1uDgb1 + a[4]*4.0*u*u*u*D1uDgb1;

  *frara=a[1]*D2uDra2 + 2.0*a[2]*(D1uDra1*D1uDra1 + u*D2uDra2) 
    + 3.0*a[3]*(2.0*u*D1uDra1*D1uDra1 + u*u*D2uDra2)
    + 4.0*a[4]*(3.0*u*u*D1uDra1*D1uDra1 + u*u*u*D2uDra2);

  *fraga=a[1]*D2uDra1ga1 + 2.0*a[2]*(D1uDga1*D1uDra1 + u*D2uDra1ga1) 
    + 3.0*a[3]*(2.0*u*D1uDra1*D1uDga1 + u*u*D2uDra1ga1)
    + 4.0*a[4]*(3.0*u*u*D1uDra1*D1uDga1 + u*u*u*D2uDra1ga1);

  *frarb=a[1]*D2uDra1rb1 + 2.0*a[2]*(D1uDrb1*D1uDra1 + u*D2uDra1rb1)
    + 3.0*a[3]*(2.0*u*D1uDra1*D1uDrb1 + u*u*D2uDra1rb1)
    + 4.0*a[4]*(3.0*u*u*D1uDra1*D1uDrb1 + u*u*u*D2uDra1rb1);

  *fragb=a[1]*D2uDra1gb1 + 2.0*a[2]*(D1uDgb1*D1uDra1 + u*D2uDra1gb1)
    + 3.0*a[3]*(2.0*u*D1uDra1*D1uDgb1 + u*u*D2uDra1gb1)
    + 4.0*a[4]*(3.0*u*u*D1uDra1*D1uDgb1 + u*u*u*D2uDra1gb1);

  *frbrb=a[1]*D2uDrb2 + 2.0*a[2]*(D1uDrb1*D1uDrb1 + u*D2uDrb2)
    + 3.0*a[3]*(2.0*u*D1uDrb1*D1uDrb1 + u*u*D2uDrb2)
    + 4.0*a[4]*(3.0*u*u*D1uDrb1*D1uDrb1 + u*u*u*D2uDrb2);

  // frbra exists
  *frbga=a[1]*D2uDrb1ga1 + 2.0*a[2]*(D1uDga1*D1uDrb1 + u*D2uDrb1ga1)
    + 3.0*a[3]*(2.0*u*D1uDrb1*D1uDga1 + u*u*D2uDrb1ga1)
    + 4.0*a[4]*(3.0*u*u*D1uDrb1*D1uDga1 + u*u*u*D2uDrb1ga1);

  *frbgb=a[1]*D2uDrb1gb1 + 2.0*a[2]*(D1uDgb1*D1uDrb1 + u*D2uDrb1gb1)
    + 3.0*a[3]*(2.0*u*D1uDrb1*D1uDgb1 + u*u*D2uDrb1gb1)
    + 4.0*a[4]*(3.0*u*u*D1uDrb1*D1uDgb1 + u*u*u*D2uDrb1gb1);

  *fgaga=a[1]*D2uDga2 + 2.0*a[2]*(D1uDga1*D1uDga1 + u*D2uDga2)
    + 3.0*a[3]*(2.0*u*D1uDga1*D1uDga1 + u*u*D2uDga2)
    + 4.0*a[4]*(3.0*u*u*D1uDga1*D1uDga1 + u*u*u*D2uDga2);

  *fgbgb=a[1]*D2uDgb2 + 2.0*a[2]*(D1uDgb1*D1uDgb1 + u*D2uDgb2)
    + 3.0*a[3]*(2.0*u*D1uDgb1*D1uDgb1 + u*u*D2uDgb2)
    + 4.0*a[4]*(3.0*u*u*D1uDgb1*D1uDgb1 + u*u*u*D2uDgb2);

  *fgagb=a[1]*D2uDga1gb1 + 2.0*a[2]*(D1uDga1*D1uDgb1 + u*D2uDga1gb1)
    + 3.0*a[3]*(2.0*u*D1uDga1*D1uDgb1 + u*u*D2uDga1gb1)
    + 4.0*a[4]*(3.0*u*u*D1uDga1*D1uDgb1 + u*u*u*D2uDga1gb1);

  return;
};

// ======================================================================================
//
//                      B97 exchange
//
// ======================================================================================



void dft_xckernel_xb97_(double *rho_a,double *rho_b, double *ScalGGAXin,
                        double *tol_rho, double *FX, long *max_pow_u, double *sol)
{
  double ra   = rho_a[0];
  double rb   = rho_b[0];
  double ScalGGAX = *ScalGGAXin;
  double gaa  =  rho_a[1] + rho_a[2] + rho_a[3]; // 1,2,3 are the x,y,z components of the gradient of the density, sqr is the square function
  double gbb  =  rho_b[1] + rho_b[2] + rho_b[3];

  double CX,CX43,CX49,ra43,ra13,rb43,rb13,raM23,rbM23;

  // enhancement factors and derivatives thereof 
  double exaa,exbb;
  double Gex_aa_ra,Gex_bb_rb,Gex_aa_gaa,Gex_bb_gbb;
  double Gex_aa_rara,Gex_bb_rbrb;
  double Gex_aa_gaagaa,Gex_bb_gbbgbb;
  double Gex_aa_ragaa,Gex_bb_rbgbb;   

  //  double xss[3] = {0.8094, 0.5073, 0.7481};
  double xss[5];
  double gamma_xss = 0.004 ;
  int i;

  for (i = 0; i <= *max_pow_u; i++) {
    xss[i]=sol[i*3];
  }

  //  CX    = 0.930525736349;
  CX = 0.9305257363490993;
  CX43  = (4.0/3.0)*CX;
  CX49  = (4.0/9.0)*CX;

  if(ra > *tol_rho ) {
    if(*max_pow_u == 2 ){
    gss3tar(ra,gaa,gamma_xss,xss,
            &exaa,&Gex_aa_ra,&Gex_aa_gaa,
            &Gex_aa_rara,&Gex_aa_ragaa,&Gex_aa_gaagaa);
    }else{
    gss5tar(ra,gaa,gamma_xss,xss,
            &exaa,&Gex_aa_ra,&Gex_aa_gaa,
            &Gex_aa_rara,&Gex_aa_ragaa,&Gex_aa_gaagaa);
    }
    //Slater-X and derivatives 
    ra43 = -CX  *pow(ra,(4.0/3.0));
    ra13 = -CX43*pow(ra,(1.0/3.0));
    raM23= -CX49*pow(ra,(-2.0/3.0));
  } else {
    ra43 = 0.0;
    ra13 = 0.0;
    raM23= 0.0;
    Gex_aa_ra     = 0.0; 
    Gex_aa_gaa    = 0.0;
    Gex_aa_rara   = 0.0;
    Gex_aa_ragaa  = 0.0;
    Gex_aa_gaagaa = 0.0;
  }

  if(rb > *tol_rho ) {
    if(*max_pow_u == 2 ){
    gss3tar(rb,gbb,gamma_xss,xss,
            &exbb,&Gex_bb_rb,&Gex_bb_gbb,
            &Gex_bb_rbrb,&Gex_bb_rbgbb,&Gex_bb_gbbgbb);
    }else{
    gss5tar(rb,gbb,gamma_xss,xss,
            &exbb,&Gex_bb_rb,&Gex_bb_gbb,
            &Gex_bb_rbrb,&Gex_bb_rbgbb,&Gex_bb_gbbgbb);
    }
    //Slater-X, more or less
    rb43 = -CX  *pow(rb,(4.0/3.0));
    rb13 = -CX43*pow(rb,(1.0/3.0));
    rbM23= -CX49*pow(rb,(-2.0/3.0));
  } else {
    rb43 = 0.0;
    rb13 = 0.0;
    rbM23 = 0.0;
    Gex_bb_rb     = 0.0; 
    Gex_bb_gbb    = 0.0;
    Gex_bb_rbrb   = 0.0;
    Gex_bb_rbgbb  = 0.0;
    Gex_bb_gbbgbb = 0.0;
  }

  // ScalGGAX is the ACM parameter for the hybrid
  // recall that B97 does not obey the UEG limit 

  FX[ _FXC_E      ] = ScalGGAX*(ra43*exaa + rb43*exbb);
  FX[ _FXC_RA     ] = ScalGGAX*(ra13*exaa + ra43*Gex_aa_ra);
  FX[ _FXC_RB     ] = ScalGGAX*(rb13*exbb + rb43*Gex_bb_rb);
  FX[ _FXC_GAA    ] = ScalGGAX*(ra43*Gex_aa_gaa);
  FX[ _FXC_GBB    ] = ScalGGAX*(rb43*Gex_bb_gbb);
  FX[ _FXC_RARA   ] = ScalGGAX*(raM23*exaa + 2.0*ra13*Gex_aa_ra + ra43*Gex_aa_rara);
  FX[ _FXC_RBRB   ] = ScalGGAX*(rbM23*exbb + 2.0*rb13*Gex_bb_rb + rb43*Gex_bb_rbrb);
  FX[ _FXC_GAAGAA ] = ScalGGAX*(ra43*Gex_aa_gaagaa);
  FX[ _FXC_GBBGBB ] = ScalGGAX*(rb43*Gex_bb_gbbgbb);
  FX[ _FXC_RAGAA  ] = ScalGGAX*(ra13*Gex_aa_gaa + ra43*Gex_aa_ragaa);
  FX[ _FXC_RBGBB  ] = ScalGGAX*(rb13*Gex_bb_gbb + rb43*Gex_bb_rbgbb);

  FX[ _FXC_GABGAB ] = 0.0;
  FX[ _FXC_GAB    ] = 0.0;  
  FX[ _FXC_RARB   ] = 0.0;
  FX[ _FXC_RAGAB  ] = 0.0;
  FX[ _FXC_RAGBB  ] = 0.0;
  FX[ _FXC_RBGAA  ] = 0.0;
  FX[ _FXC_RBGAB  ] = 0.0;
  FX[ _FXC_GAAGAB ] = 0.0;
  FX[ _FXC_GAAGBB ] = 0.0;
  FX[ _FXC_GBBGAB ] = 0.0;

};

// ======================================================================================
//
//                      B97 correlation 
//
// ======================================================================================
 void dft_xckernel_cb97_(double *rho_a,double *rho_b, double *ScalGGACin, 
			 double *tol_rho, double *FC, long *max_pow_u, double *sol)
{ 
  double ra   = rho_a[0];
  double rb   = rho_b[0];
  double ScalGGAC = *ScalGGACin;
  double gaa  =  rho_a[1] + rho_a[2] + rho_a[3];
  double gbb  =  rho_b[1] + rho_b[2] + rho_b[3];

  //LDA and derivatives
  double ec_a0,ec_b0;
  double Gec_a0=0.0,Gec_b0=0.0;
  double Gec_a0_ra,Gec_a0_rara;
  double Gec_b0_rb,Gec_b0_rbrb;

  double ec_ab;
  double Gec_ab_ra=0.0,Gec_ab_rb=0.0;
  double Gec_ab_rara,Gec_ab_rarb,Gec_ab_rbrb;

  // enhancement factors and derivatives
  // same spin
  double gcaa,gcbb;
  double Ggc_aa_ra,Ggc_bb_rb,Ggc_ab_ra,Ggc_ab_rb;
  double Ggc_aa_gaa,Ggc_bb_gbb;
  double Ggc_aa_rara,Ggc_aa_ragaa,Ggc_aa_gaagaa;
  double Ggc_bb_rbrb,Ggc_bb_rbgbb,Ggc_bb_gbbgbb;

  // opposite spin 
  double gcab;
  //  double *gcab=123.4567890123;
  double Ggc_ab_gaa,Ggc_ab_gbb;
  double Ggc_ab_rara,Ggc_ab_rarb,Ggc_ab_rbrb;
  double Ggc_ab_ragaa,Ggc_ab_ragbb,Ggc_ab_rbgaa,Ggc_ab_rbgbb;
  double Ggc_ab_gaagaa,Ggc_ab_gaagbb,Ggc_ab_gbbgbb;

  //old  double css[3] = {0.1737, 2.3487, -2.4868};
  //old  double cab[3] = {0.9454, 0.7471, -4.5961};
  double css[5];
  double cab[5];
  double gamma_css=0.2;
  double gamma_cab=0.006;
  int i;
  double FCLDA[_FCLDA_ELEMENTS];

  for (i = 0; i <= *max_pow_u; i++) {
    css[i]=sol[i*3+1];
    cab[i]=sol[i*3+2];
  }

  //normal
  if(ra > *tol_rho || rb > *tol_rho) {
    dft_xckernel_pwlda(ra, rb,
                       FCLDA);
    //LDA derivatives
    ec_ab      =FCLDA[_FXC_E]   ;
    Gec_ab_ra  =FCLDA[_FXC_RA]  ;
    Gec_ab_rb  =FCLDA[_FXC_RB]  ;
    Gec_ab_rara=FCLDA[_FXC_RARA];
    Gec_ab_rbrb=FCLDA[_FXC_RBRB];
    Gec_ab_rarb=FCLDA[_FXC_RARB];
  } else {
    ec_ab = 0.0;
    Gec_ab_ra = 0.0;
    Gec_ab_rb = 0.0;
    Gec_ab_rara=0.0;
    Gec_ab_rbrb=0.0;
    Gec_ab_rarb=0.0;
  }

  //alpha density only
  if(ra > *tol_rho) {
    dft_xckernel_pwlda(ra, 0.0,
                       FCLDA);
    //ec_ra0 *= ra;
    ec_a0      =FCLDA[_FXC_E]   ;
    Gec_a0_ra  =FCLDA[_FXC_RA]  ;
    //Gec_a0_rb  =0.0;//FCLDA[_FXC_RB]  ;
    Gec_a0_rara=FCLDA[_FXC_RARA];
    //Gec_a0_rbrb=0.0;//FCLDA[_FXC_RBRB];
    //Gec_a0_rarb=0.0;//FCLDA[_FXC_RARB];
  } else {
    ec_a0 = 0.0;
    Gec_a0 = 0.0;
    Gec_a0_rara = 0.0;
  }

  //beta density only
  if(rb > *tol_rho) {
    dft_xckernel_pwlda(0.0, rb,
                       FCLDA);
    //ec_rb0 *= rb;
    ec_b0      =FCLDA[_FXC_E]   ;
    //Gec_b0_ra  =0.0;//FCLDA[_FXC_RA]  ;
    Gec_b0_rb  =FCLDA[_FXC_RB]  ;
    //Gec_b0_rara=0.0;//FCLDA[_FXC_RARA];
    Gec_b0_rbrb=FCLDA[_FXC_RBRB];
    //Gec_b0_rarb=0.0;//FCLDA[_FXC_RARB];
  } else {
    ec_b0 = 0.0;
    Gec_b0 = 0.0;
    Gec_b0_rbrb = 0.0;
  }

  ec_ab = ec_ab - ec_a0 - ec_b0;
  Gec_ab_ra = Gec_ab_ra - Gec_a0_ra;
  Gec_ab_rb = Gec_ab_rb - Gec_b0_rb;
  Gec_ab_rara = Gec_ab_rara - Gec_a0_rara;
  Gec_ab_rbrb = Gec_ab_rbrb - Gec_b0_rbrb;

  if(ra > *tol_rho) {
    if(*max_pow_u == 2 ){
    gss3tar(ra,gaa,gamma_css,css,&gcaa,
            &Ggc_aa_ra,&Ggc_aa_gaa,
            &Ggc_aa_rara,&Ggc_aa_ragaa,&Ggc_aa_gaagaa);
    }else{
    gss5tar(ra,gaa,gamma_css,css,&gcaa,
            &Ggc_aa_ra,&Ggc_aa_gaa,
            &Ggc_aa_rara,&Ggc_aa_ragaa,&Ggc_aa_gaagaa);
    }
  } else {
    gcaa = css[0] + css[1] + css[2];
    Ggc_aa_ra = 0.0; 
    Ggc_aa_gaa = 0.0; 
    Ggc_aa_rara = 0.0;
    Ggc_aa_ragaa = 0.0;
    Ggc_aa_gaagaa = 0.0;
  }

  if(rb > *tol_rho) {
    if(*max_pow_u == 2 ){
    gss3tar(rb,gbb,gamma_css,css,&gcbb,
            &Ggc_bb_rb,&Ggc_bb_gbb,
            &Ggc_bb_rbrb,&Ggc_bb_rbgbb,&Ggc_bb_gbbgbb);
    }else{
    gss5tar(rb,gbb,gamma_css,css,&gcbb,
            &Ggc_bb_rb,&Ggc_bb_gbb,
            &Ggc_bb_rbrb,&Ggc_bb_rbgbb,&Ggc_bb_gbbgbb);
    }
  } else {
    gcbb = css[0] + css[1] + css[2];
    Ggc_bb_rb = 0.0; 
    Ggc_bb_gbb = 0.0; 
    Ggc_bb_rbrb = 0.0;
    Ggc_bb_rbgbb = 0.0; 
    Ggc_bb_gbbgbb = 0.0;
  }

  if(ra > *tol_rho && rb > *tol_rho) {
    if(*max_pow_u == 2 ){
    gab3tar(ra,rb,gaa,gbb,gamma_cab, cab,
            &gcab,&Ggc_ab_ra,&Ggc_ab_rb,&Ggc_ab_gaa,&Ggc_ab_gbb,
            &Ggc_ab_rara,&Ggc_ab_rarb,&Ggc_ab_rbrb,
            &Ggc_ab_ragaa,&Ggc_ab_ragbb,&Ggc_ab_rbgaa,&Ggc_ab_rbgbb,
            &Ggc_ab_gaagaa,&Ggc_ab_gaagbb,&Ggc_ab_gbbgbb);
    }else{
    gab5tar(ra,rb,gaa,gbb,gamma_cab, cab,
            &gcab,&Ggc_ab_ra,&Ggc_ab_rb,&Ggc_ab_gaa,&Ggc_ab_gbb,
            &Ggc_ab_rara,&Ggc_ab_rarb,&Ggc_ab_rbrb,
            &Ggc_ab_ragaa,&Ggc_ab_ragbb,&Ggc_ab_rbgaa,&Ggc_ab_rbgbb,
            &Ggc_ab_gaagaa,&Ggc_ab_gaagbb,&Ggc_ab_gbbgbb);
    }
    
  } else { 
    gcab = cab[0] + cab[1] + cab[2]; // u == 1 if s->inf. 
    // Not necessary for exchange, since there is no remapping going on.
    Ggc_ab_ra = 0.0;
    Ggc_ab_rb = 0.0;
    Ggc_ab_gaa = 0.0;
    Ggc_ab_gbb = 0.0;
    Ggc_ab_rara = 0.0;
    Ggc_ab_rarb = 0.0;
    Ggc_ab_rbrb = 0.0;
    Ggc_ab_ragaa = 0.0;
    Ggc_ab_ragbb = 0.0; 
    Ggc_ab_rbgaa = 0.0;
    Ggc_ab_rbgbb = 0.0;
    Ggc_ab_gaagaa = 0.0;
    Ggc_ab_gaagbb = 0.0; 
    Ggc_ab_gbbgbb = 0.0;
  }

  // ScalGGAC is the ACM D parameter as used in double hybrids 

   FC[ _FXC_E      ] = ScalGGAC*(ec_a0*gcaa + ec_b0*gcbb + ec_ab*gcab);

   FC[ _FXC_RA     ] = ScalGGAC*(Gec_a0_ra*gcaa + ec_a0*Ggc_aa_ra + 0.0 + Gec_ab_ra*gcab + ec_ab*Ggc_ab_ra);

   FC[ _FXC_RB     ] = ScalGGAC*(0.0 + Gec_b0_rb*gcbb + ec_b0*Ggc_bb_rb + Gec_ab_rb*gcab + ec_ab*Ggc_ab_rb);
   FC[ _FXC_GAA    ] = ScalGGAC*(0.0 + ec_a0*Ggc_aa_gaa + 0.0 + ec_ab*Ggc_ab_gaa);
   FC[ _FXC_GBB    ] = ScalGGAC*(0.0 + ec_b0*Ggc_bb_gbb + 0.0 + ec_ab*Ggc_ab_gbb);

   FC[ _FXC_RARA   ] = ScalGGAC*(Gec_a0_rara*gcaa + 2.0*Gec_a0_ra*Ggc_aa_ra + ec_a0*Ggc_aa_rara + Gec_ab_rara*gcab + 2.0*Gec_ab_ra*Ggc_ab_ra + ec_ab*Ggc_ab_rara);
   FC[ _FXC_RBRB   ] = ScalGGAC*(Gec_b0_rbrb*gcbb + 2.0*Gec_b0_rb*Ggc_bb_rb + ec_b0*Ggc_bb_rbrb + Gec_ab_rbrb*gcab + 2.0*Gec_ab_rb*Ggc_ab_rb + ec_ab*Ggc_ab_rbrb);
   FC[ _FXC_RARB   ] = ScalGGAC*(0.0 + 0.0 + Gec_ab_rarb*gcab + Gec_ab_ra*Ggc_ab_rb + Gec_ab_rb*Ggc_ab_ra + ec_ab*Ggc_ab_rarb);

   FC[ _FXC_GAAGAA ] = ScalGGAC*(0.0 + ec_a0*Ggc_aa_gaagaa + 0.0 + ec_ab*Ggc_ab_gaagaa);
   FC[ _FXC_GBBGBB ] = ScalGGAC*(0.0 + ec_b0*Ggc_bb_gbbgbb + 0.0 + ec_ab*Ggc_ab_gbbgbb);
   FC[ _FXC_GAAGBB ] = ScalGGAC*(ec_ab*Ggc_ab_gaagbb); 


   FC[ _FXC_RAGAA  ] = ScalGGAC*(0.0 + Gec_a0_ra*Ggc_aa_gaa + 0.0 + ec_a0*Ggc_aa_ragaa + 0.0 + Gec_ab_ra*Ggc_ab_gaa + 0.0 + ec_ab*Ggc_ab_ragaa);
     
   FC[ _FXC_RBGBB  ] = ScalGGAC*(0.0 + Gec_b0_rb*Ggc_bb_gbb + 0.0 + ec_b0*Ggc_bb_rbgbb + 0.0 + Gec_ab_rb*Ggc_ab_gbb + 0.0 + ec_ab*Ggc_ab_rbgbb);

   FC[ _FXC_RAGBB  ] = ScalGGAC*(0.0 + 0.0 + 0.0 + 0.0 + 0.0 + Gec_ab_ra*Ggc_ab_gbb + 0.0 + ec_ab*Ggc_ab_ragbb);

   FC[ _FXC_RBGAA  ] = ScalGGAC*(0.0 + 0.0 + 0.0 + 0.0 + 0.0 + Gec_ab_rb*Ggc_ab_gaa + 0.0 + ec_ab*Ggc_ab_rbgaa);

   // rest is zero, there is no GAB dependence
   FC[ _FXC_GABGAB ] = 0.0;
   FC[ _FXC_GAAGAB ] = 0.0;
   FC[ _FXC_GBBGAB ] = 0.0;
   FC[ _FXC_RAGAB  ] = 0.0;
   FC[ _FXC_RBGAB  ] = 0.0;
   FC[ _FXC_GAB    ] = 0.0;

}

/* $Id$ */
