/*
 $Id: dft.c,v 1.1 2001-08-30 16:58:35 bylaska Exp $
   dft.c - 
   Author - Eric Bylaska

*/

#include	<stdio.h>
#include	"loggrid.h"
#include	"get_word.h"
#include	"dft.h"

/* Kawai-Weare default solver */
static int	Hartree_Type     = Hartree_On;
static int	Exchange_Type    = Exchange_Dirac;
static int	Correlation_Type = Correlation_Vosko;

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
      
      if (strcmp("off",w)  ==0)       
	Correlation_Type=Correlation_Off;
   }
   fclose(fp);

}

void set_Exchange_DFT(exchange)
int exchange;
{
    Exchange_Type = exchange;
}
void set_Correlation_DFT(correlation)
int correlation;
{
    Correlation_Type = correlation;
}
void set_Hartree_DFT(hartree)
int hartree;
{
    Hartree_Type = hartree;
}

void	R_Exchange_DFT(rho,Vx,Ex,Px)

double	*rho;
double	*Vx;
double	*Ex;
double	*Px;
{
   int k,Ngrid;

   if (Exchange_Type==Exchange_Dirac)
      R_Dirac_Exchange(rho,Vx,Ex,Px);
   else if (Exchange_Type==Exchange_PBE96)
      R_PBE96_Exchange(rho,Vx,Ex,Px);
   else
   {
      *Ex   = 0.0;
      *Px   = 0.0;
      Ngrid = N_LogGrid();
      for (k=0; k<Ngrid; ++k)
         Vx[k] = 0.0;
   }
}
      
void	R_Correlation_DFT(rho,Vc,Ec,Pc)

double	*rho;
double	*Vc;
double	*Ec;
double	*Pc;
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
   else
   {
      *Ec   = 0.0;
      *Pc   = 0.0;
      Ngrid = N_LogGrid();
      for (k=0; k<Ngrid; ++k)
         Vc[k] = 0.0;
   }
}

double	R_Hartree_DFT(rho,charge,Vh)

double	*rho;
double	charge;
double	*Vh;
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
     s = "PBE96 (Perdew, Kurke, and Ernzerhof) parameterization";
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
       s = "PBE96 (Perdew, Kurke, and Ernzerhof) parameterization";
    else if (Correlation_Type==Correlation_Perdew_Wang)
       s = "Perdew and Wang parameterization";
    else
       s = "No Correlation";

   return s;
}

