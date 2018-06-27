/* dft.c -
   Author - Eric Bylaska
 $Id$

*/

#include	<stdio.h>
#include 	<string.h>
#include	"loggrid.h"
#include	"get_word.h"
#include	"dft.h"

/* Kawai-Weare default solver */
static int	Hartree_Type     = Hartree_On;
static int	Exchange_Type    = Exchange_Dirac;
static int	Correlation_Type = Correlation_Vosko;
static double   screening_cut = 0.0;
static double   blyp_screening_cut = 0.0;


void	init_DFT(char	*filename)
{
    FILE	*fp;
    char	*w;
    double alpha;

    /* set hartree type */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<hartree>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        w = get_word(fp);
        if (strcmp("on",w)==0)    Hartree_Type=Hartree_On;
        if (strcmp("off",w)  ==0) Hartree_Type=Hartree_Off;
    }
    fclose(fp);


    /* set exchange type */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<exchange>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        w = get_word(fp);
        if (strcmp("dirac",w)==0) Exchange_Type=Exchange_Dirac;
        if (strcmp("pbe96",w)==0) Exchange_Type=Exchange_PBE96;
        if (strcmp("becke",w)==0) Exchange_Type=Exchange_Becke;
        if (strcmp("revpbe",w)==0) Exchange_Type=Exchange_revPBE;
        if (strcmp("off",w)  ==0) Exchange_Type=Exchange_Off;
    }
    fclose(fp);

    /* set exchange alpha */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<alpha>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        fscanf(fp,"%lf",&alpha);
        set_Dirac_alpha(alpha);
    }
    fclose(fp);


    /* set correlation type */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<correlation>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        w = get_word(fp);
        if (strcmp("perdew-zunger",w)==0)
            Correlation_Type=Correlation_Perdew_Zunger;
        if (strcmp("vosko",w)  ==0)
            Correlation_Type=Correlation_Vosko;
        if (strcmp("pbe96",w)  ==0)
            Correlation_Type=Correlation_PBE96;
        if (strcmp("perdew-wang",w)  ==0)
            Correlation_Type=Correlation_Perdew_Wang;
        if (strcmp("lyp",w)  ==0)
            Correlation_Type=Correlation_LYP;
        if (strcmp("revpbe",w)  ==0)
            Correlation_Type=Correlation_revPBE;
        if (strcmp("off",w)  ==0)
            Correlation_Type=Correlation_Off;
    }
    fclose(fp);

   /* set screening_cut */
   screening_cut = 0.0;
   fp = fopen(filename,"r+");
   w = get_word(fp);
   while ((w!=NIL) && (strcmp("<screening_cut>",w)!=0))
      w = get_word(fp);
   if (w!=NIL)
   {
      fscanf(fp,"%lf",&screening_cut);
   }
   fclose(fp);

   /* set blyp_screening_cut */
   blyp_screening_cut = 0.0;
   fp = fopen(filename,"r+");
   w = get_word(fp);
   while ((w!=NIL) && (strcmp("<blyp_screening_cut>",w)!=0))
      w = get_word(fp);
   if (w!=NIL)
   {
      fscanf(fp,"%lf",&blyp_screening_cut);
   }
   fclose(fp);


}

void set_Exchange_DFT(int exchange)
{
    Exchange_Type = exchange;
}

void set_Correlation_DFT(int correlation)
{
    Correlation_Type = correlation;
}

void set_Hartree_DFT(int hartree)
{
    Hartree_Type = hartree;
}


void R_Screening_Cut(double * Vx)
{
   int k,NN,n0,n1=0;
   double r0,r1=0.0,v0,v1=0.0,m,b;
   double *r;
   if (screening_cut>0.0)
   {
      r = r_LogGrid();
      NN = index_r_LogGrid(screening_cut) + 5;
      for (k=0; k<NN; ++k)
        if (r[k] < screening_cut)
           { n0=n1; r0=r1; v0=v1; n1=k; r1=r[k]; v1=Vx[k]; }
      m = (v1-v0)/(r1-r0);
      b =  v1 - m*r1;
      for (k=0; k<n1; ++k)
         Vx[k] = m*r[k] + b;
   } 
   else if ((blyp_screening_cut>0.0) && (Exchange_Type==Exchange_Becke) && (Correlation_Type==Correlation_LYP)) 
   {
      r = r_LogGrid();
      NN = index_r_LogGrid(blyp_screening_cut) + 5;
      for (k=0; k<NN; ++k)
        if (r[k] < blyp_screening_cut)
           { n0=n1; r0=r1; v0=v1; n1=k; r1=r[k]; v1=Vx[k]; }
      m = (v1-v0)/(r1-r0);
      b =  v1 - m*r1;
      for (k=0; k<n1; ++k)
         Vx[k] = m*r[k] + b;
   } 

}




void R_Exchange_DFT(double * rho, double * Vx, double * Ex, double * Px)
{
    int k,Ngrid;

    if (Exchange_Type==Exchange_Dirac)
        R_Dirac_Exchange(rho,Vx,Ex,Px);
    else if (Exchange_Type==Exchange_PBE96)
        R_PBE96_Exchange(rho,Vx,Ex,Px);
    else if (Exchange_Type==Exchange_Becke)
        R_Becke_Exchange(rho,Vx,Ex,Px);
    else if (Exchange_Type==Exchange_revPBE)
        R_revPBE_Exchange(rho,Vx,Ex,Px);
    else
    {
        *Ex   = 0.0;
        *Px   = 0.0;
        Ngrid = N_LogGrid();
        for (k=0; k<Ngrid; ++k)
            Vx[k] = 0.0;
    }
}

void R_Correlation_DFT(double * rho, double * Vc, double * Ec, double * Pc)
{
    int k,Ngrid;

    if (Correlation_Type==Correlation_Vosko)
        R_Vosko(rho,Vc,Ec,Pc);
    else if (Correlation_Type==Correlation_Perdew_Zunger)
        R_Perdew_Zunger(rho,Vc,Ec,Pc);
    else if (Correlation_Type==Correlation_PBE96)
        R_PBE96_Correlation(rho,Vc,Ec,Pc);
    else if (Correlation_Type==Correlation_Perdew_Wang)
        R_Perdew_Wang(rho,Vc,Ec,Pc);
    else if (Correlation_Type==Correlation_LYP)
        R_LYP_Correlation(rho,Vc,Ec,Pc);
    else if (Correlation_Type==Correlation_revPBE)
        R_revPBE_Correlation(rho,Vc,Ec,Pc);
    else
    {
        *Ec   = 0.0;
        *Pc   = 0.0;
        Ngrid = N_LogGrid();
        for (k=0; k<Ngrid; ++k)
            Vc[k] = 0.0;
    }
}

double	R_Hartree_DFT(double * rho, double charge, double * Vh)
{
    int k,Ngrid;
    double ph;

    if (Hartree_Type==Hartree_On)
        ph = R_Hartree(rho,charge,Vh);
    else
    {
        ph   = 0.0;
        Ngrid = N_LogGrid();
        for (k=0; k<Ngrid; ++k)
            Vh[k] = 0.0;
    }

    return ph;
}




char	*hartree_Name_DFT()
{
    char *s;
    if (Hartree_Type==Hartree_On)
        s = "On";
    else
        s = "Off";

    return s;
}

char	*exchange_Name_DFT()
{
    char *s;
    if (Exchange_Type==Exchange_Dirac)
        s = "Dirac";
    else if (Exchange_Type==Exchange_PBE96)
        s = "PBE96 (Perdew, Burke, and Ernzerhof) parameterization";
    else if (Exchange_Type==Exchange_Becke)
        s = "Becke88 parameterization";
    else if (Exchange_Type==Exchange_revPBE)
        s = "revPBE (Norskov) parameterization";
    else
        s = "No Exchange";

    return s;
}

char	*correlation_Name_DFT()
{
    char *s;
    if (Correlation_Type==Correlation_Vosko)
        s = "Vosko parameterization";
    else if (Correlation_Type==Correlation_Perdew_Zunger)
        s = "Perdew and Zunger parameterization";
    else if (Correlation_Type==Correlation_PBE96)
        s = "PBE96 (Perdew, Burke, and Ernzerhof) parameterization";
    else if (Correlation_Type==Correlation_Perdew_Wang)
        s = "Perdew and Wang parameterization";
    else if (Correlation_Type==Correlation_LYP)
        s = "LYP (Lee, Yang, Parr) parameterization";
    else if (Correlation_Type==Correlation_revPBE)
        s = "revPBE (Norskov) parameterization";
    else
        s = "No Correlation";

    return s;
}

