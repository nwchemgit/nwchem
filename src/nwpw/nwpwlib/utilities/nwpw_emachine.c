/*
 $Id: nwpw_emachine.c 25745 2014-06-08 07:46:01Z d3y133 $
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include "typesf2c.h"

#include        <ctype.h>
#include        <stdio.h>



#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif


/*************************  parser routines *******************************/

static void nwpw_removespaces(int ns,char str[], char nstr[])
{
   int i,j;
   j = 0;
   for (i=0; i<ns; ++i)
   {
      if ((str[i]!=' ') && (str[i]!='\t'))
      {
        nstr[j] = tolower(str[i]);
        ++j;
      }
   }
   nstr[j] = '\0';
}


static void nwpw_findsestride(char str[], int *s, int *e, int *stride)
{
    int done;
    int i = 0;

    if (sscanf(str,"%d:%d:%d",s,e,stride)==3) return;

    if (isdigit(str[i])) sscanf(&str[i],"%d",s);
    if (str[i]=='n') *s = -9999;
    if (str[i]=='i') *s = -8888;
    if (str[i]=='j') *s = -7777;
    if (str[i]=='k') *s = -6666;
    done = 0; while (!done) { ++i; if (str[i]==':') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",e);
    if (str[i]=='n') *e = -9999;
    if (str[i]=='i') *e = -8888;
    if (str[i]=='j') *e = -7777;
    if (str[i]=='k') *e = -6666;
    done = 0; while (!done) { ++i; if (str[i]==':') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",stride);
    if (str[i]=='n') *stride = -9999;
    if (str[i]=='i') *stride = -8888;
    if (str[i]=='j') *stride = -7777;
    if (str[i]=='k') *stride = -6666;
}

static void nwpw_findrion(char str[], int *s, int *e) 
{
    int done;
    int i = 0;

    if (sscanf(str,"%d,%d",s,e)==2) return;

    if (isdigit(str[i])) sscanf(&str[i],"%d",s);
    if (str[i]=='n') *s = -9999;
    if (str[i]=='i') *s = -8888;
    if (str[i]=='j') *s = -7777;
    if (str[i]=='k') *s = -6666;
    done = 0; while (!done) { ++i; if (str[i]==',') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",e);
    if (str[i]=='n') *e = -9999;
    if (str[i]=='i') *e = -8888;
    if (str[i]=='j') *e = -7777;
    if (str[i]=='k') *e = -6666;
}



static int nwpw_findnum(char str[])
{
   int pstack=0;
   int ns=strlen(str);
   int i=-1;
   int done=0;

   while (!done)
   {
      ++i;

      if (isdigit(str[i]) && (pstack<1))
         done = 1;
      else if (str[i]=='(')
         ++pstack;
      else if (str[i]==')')
         --pstack;
     if (i>=ns)
     {
         i = -1;
         done = 1;
     }
   }
   return i;
}

static int nwpw_findtoken(char str[], char token[])
{
   int t;
   int pstack=0;
   int ns=strlen(str);
   int nt=strlen(token);
   int done=0;
   int i=-1;
   char oldch = 'a';

   while (!done)
   {
      ++i;
      if ((str[i]==token[0]) && (pstack<1) && (oldch!='e'))
      {
         done = 1;
         for (t=1; t<nt; ++t) if (str[i+t]!=token[t]) done = 0;
      }
      else if (str[i]=='(')
         ++pstack;
      else if (str[i]==')')
         --pstack;

     if (i>=ns)
     {
        i    = -1;
        done = 1;
     }
     oldch = str[i];
   }

   return i;
}

static int cln = 0;
static int nf = 0;

static void nwpw_gencode(Integer code[], double fconst[], int ln, char expr[])
{
   int  lnl,lnr,t,t2,nt,s,e,stride,i,ii;
   char lexpr[500],rexpr[500],tmps[10],tmpss[10];
   double f;

   //printf("EXPR: %s xxx\n",expr);

   if ((t=nwpw_findtoken(expr,"+")) >= 0)
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+1]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   +  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -1; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"-")) >= 0)
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+1]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   -  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -2; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if (((t=nwpw_findtoken(expr,"*")) >= 0) && (nwpw_findtoken(expr,"**") < 0))
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+1]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   *  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -3; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"/")) >= 0)
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+1]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   /  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -4; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"**")) >= 0)
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+2]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   **  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -5; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"^")) >= 0)
   {
      strncpy(lexpr,expr,t); lexpr[t] = '\0'; strcpy(rexpr,&expr[t+1]); ++cln; lnl = cln; ++cln; lnr = cln;
      //printf("code: ln=%d   ^  %d  %d\n",ln,lnl,lnr);
      code[5*ln] = -5; code[5*ln+1] = lnl; code[5*ln+2] = lnr; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnl,lexpr);
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"abs")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   abs  %d  \n",ln,lnr);
      code[5*ln] = -6; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"exp")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   exp  %d  \n",ln,lnr);
      code[5*ln] = -7; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"log")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   log  %d  \n",ln,lnr);
      code[5*ln] = -8; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sqrt")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   sqrt  %d  \n",ln,lnr);
      code[5*ln] = -9; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sinh")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   sinh  %d  \n",ln,lnr);
      code[5*ln] = -10; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"cosh")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   cosh  %d  \n",ln,lnr);
      code[5*ln] = -11; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"tanh")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   tanh  %d  \n",ln,lnr);
      code[5*ln] = -12; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"asin")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   asin  %d  \n",ln,lnr);
      code[5*ln] = -16; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"acos")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   acos  %d  \n",ln,lnr);
      code[5*ln] = -17; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"atan")) >= 0)
   {
      strcpy(rexpr,&expr[t+4]); ++cln; lnr = cln;
      //printf("code: ln=%d   atan  %d  \n",ln,lnr);
      code[5*ln] = -18; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sin")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   sin  %d  \n",ln,lnr);
      code[5*ln] = -13; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"cos")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   cos  %d  \n",ln,lnr);
      code[5*ln] = -14; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"tan")) >= 0)
   {
      strcpy(rexpr,&expr[t+3]); ++cln; lnr = cln;
      //printf("code: ln=%d   tan  %d  \n",ln,lnr);
      code[5*ln] = -15; code[5*ln+1] = lnr; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sumii")) >= 0)
   {
      nwpw_findsestride(&expr[t+6],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+6],"]");
      strcpy(rexpr,&expr[t+6+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumii  %d   ii=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -22; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
      
   }
   else if ((t=nwpw_findtoken(expr,"sumjj")) >= 0)
   {
      nwpw_findsestride(&expr[t+6],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+6],"]");
      strcpy(rexpr,&expr[t+7+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumjj  %d   jj=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -23; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sumkk")) >= 0)
   {
      nwpw_findsestride(&expr[t+6],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+6],"]");
      strcpy(rexpr,&expr[t+6+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumkk  %d   kk=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -24; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sumi")) >= 0)
   {
      nwpw_findsestride(&expr[t+5],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+5],"]");
      strcpy(rexpr,&expr[t+6+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumi  %d   i=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -19; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sumj")) >= 0)
   {
      nwpw_findsestride(&expr[t+5],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+5],"]");
      strcpy(rexpr,&expr[t+5+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumj  %d   j=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -20; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"sumk")) >= 0)
   {
      nwpw_findsestride(&expr[t+5],&s,&e,&stride);
      t2 = nwpw_findtoken(&expr[t+5],"]");
      strcpy(rexpr,&expr[t+5+t2]); ++cln; lnr = cln;
      //printf("code: ln=%d   sumk  %d   k=%d:%d:%d \n",ln,lnr,s,e,stride);
      code[5*ln] = -21; code[5*ln+1] = lnr; code[5*ln+2] = s; code[5*ln+3] = e; code[5*ln+4] = stride;
      nwpw_gencode(code,fconst,lnr,rexpr);
   }
   else if ((t=nwpw_findtoken(expr,"rion")) >= 0)
   {
      nwpw_findrion(&expr[t+5],&i,&ii);
      //printf("code: ln=%d   rion[%d,%d]\n",ln,i,ii);
      code[5*ln] = -30; code[5*ln+1] = 0; code[5*ln+2] = i; code[5*ln+3] = ii; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"i")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -40; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"j")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -41; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"k")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -42; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"ii")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -43; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"jj")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -44; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findtoken(expr,"kk")) >= 0)
   {
      //printf("code: ln=%d  i\n",ln);
      code[5*ln] = -45; code[5*ln+1] = 0; code[5*ln+2] = 0; code[5*ln+3] = 0; code[5*ln+4] = 0;
   }
   else if ((t=nwpw_findnum(expr)) >= 0)
   {
      sscanf(&expr[t],"%le",&f);
      fconst[nf] = f;
      //printf("code: ln=%d  pointer=%d  fconst  %le  \n",ln,nf,f);
      code[5*ln] = -25; code[5*ln+1] = 0; code[5*ln+2] = nf; code[5*ln+3] = 0; code[5*ln+4] = 0;
      ++nf;
   }
   else if ((t=nwpw_findtoken(expr,"(")) >= 0)
   {
      t2 = nwpw_findtoken(&expr[t+1],")");
      nt=t2-t;
      strncpy(lexpr,&expr[t+1],nt); lexpr[nt] = '\0'; lnl = ln;
      nwpw_gencode(code,fconst,lnl,lexpr);
   }
   else
   {
      f = 0.0;
      fconst[nf] = f;
      //printf("code: ln=%d pointer=%d  fconst  %le  \n",ln,nf,f);
      code[5*ln] = -25; code[5*ln+1] = 0; code[5*ln+2] = nf; code[5*ln+3] = 0; code[5*ln+4] = 0;
      ++nf;
   }
}
 





/*************************  parser routines *******************************/



/********************************************************
 *                                                      *
 *               nwpw_emachine_parse                    *
 *                                                      *
 ********************************************************/

#if (defined(CRAY) &&!defined(__crayx1)) || defined(CRAY_T3D) || defined(WIN32)
#define nwpw_emachine_parse_ nwpw_emachine_parse
#endif

void FATR nwpw_emachine_parse_
#if defined(USE_FCD)
( const _fcd fcd_eqn_string,
 Integer *n1,
 Integer *nc0,
 Integer code0[],
 Integer *nf0,
 double fconst0[])
{
    char *eqn_string = _fcdtocp(fcd_eqn_string);

#else
(eqn_string,n1,nc0,code0,nf0,fconst0)
char	eqn_string[];
Integer	*n1;
Integer	*nc0;
Integer code0[];
Integer	*nf0;
double fconst0[];
{

#endif

   int i;
   char nstr[500];
   int ns = *n1;

   cln = 0;
   nf  = 0;
   nwpw_removespaces(ns,eqn_string,nstr);

   //printf("eqnstring = %s xxx\n",eqn_string);
   //printf("nstr = %s xxx\n",nstr);
   nwpw_gencode(code0, fconst0,0,nstr); ++cln;


   *nc0 = cln;
   *nf0 = nf;

}

