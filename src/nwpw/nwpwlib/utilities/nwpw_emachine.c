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
#if !defined(__MINGW32__)
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif
#include "typesf2c.h"

#include        <ctype.h>
#include        <stdio.h>



#if defined(CRAY)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

extern void lattice_min_difference_();



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
    if (str[i+1]=='i') *s = -5555;
    if (str[i+1]=='j') *s = -4444;
    if (str[i+1]=='k') *s = -3333;
    done = 0; while (!done) { ++i; if (str[i]==':') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",e);
    if (str[i]=='n') *e = -9999;
    if (str[i]=='i') *e = -8888;
    if (str[i]=='j') *e = -7777;
    if (str[i]=='k') *e = -6666;
    if (str[i+1]=='i') *e = -5555;
    if (str[i+1]=='j') *e = -4444;
    if (str[i+1]=='k') *e = -3333;
    done = 0; while (!done) { ++i; if (str[i]==':') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",stride);
    if (str[i]=='n') *stride = -9999;
    if (str[i]=='i') *stride = -8888;
    if (str[i]=='j') *stride = -7777;
    if (str[i]=='k') *stride = -6666;
    if (str[i+1]=='i') *stride = -5555;
    if (str[i+1]=='j') *stride = -4444;
    if (str[i+1]=='k') *stride = -3333;
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
    if (str[i+1]=='i') *s = -5555;
    if (str[i+1]=='j') *s = -4444;
    if (str[i+1]=='k') *s = -3333;
    done = 0; while (!done) { ++i; if (str[i]==',') done = 1; }
    ++i;
    if (isdigit(str[i])) sscanf(&str[i],"%d",e);
    if (str[i]=='n') *e = -9999;
    if (str[i]=='i') *e = -8888;
    if (str[i]=='j') *e = -7777;
    if (str[i]=='k') *e = -6666;
    if (str[i+1]=='i') *e = -5555;
    if (str[i+1]=='j') *e = -4444;
    if (str[i+1]=='k') *e = -3333;
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
   int  lnl,lnr,t,t2,nt,s,e,stride,i,ii,jj;
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

   else if ((t=nwpw_findtoken(expr,"rion")) >= 0)
   {
      nwpw_findrion(&expr[t+5],&i,&ii);
      //printf("code: ln=%d   rion[%d,%d]\n",ln,i,ii);
      code[5*ln] = -30; code[5*ln+1] = 0; code[5*ln+2] = i; code[5*ln+3] = ii; code[5*ln+4] = 0;
   }

   else if ((t=nwpw_findtoken(expr,"d")) >= 0)
   {
      nwpw_findrion(&expr[t+2],&ii,&jj);
      //printf("code: ln=%d   d[%d,%d]\n",ln,ii,jj);
      code[5*ln] = -31; code[5*ln+1] = 0; code[5*ln+2] = ii; code[5*ln+3] = jj; code[5*ln+4] = 0;
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

/*************************  emachine routines *******************************/


/* registers */
static Integer i,j,k,ii,jj,kk;


/********************************************
 *                                          *
 *                nwpw_emachine             *
 *                                          *
 ********************************************/

double nwpw_emachine(int ln, Integer code[], double fconst[], Integer nion, double rion[])
{
   Integer op,arg1,arg2,arg3,arg4,ln1,ln2;
   double f,a,b,x,y,z;

   f = 0.0;
   op  = code[5*ln];
   ln1 = code[5*ln+1];
   ln2 = arg1 = code[5*ln+2];
   arg2       = code[5*ln+3];
   arg3       = code[5*ln+4];

   if (arg1==(-9999)) arg1 = nion;
   if (arg1==(-8888)) arg1 = i;
   if (arg1==(-7777)) arg1 = j;
   if (arg1==(-6666)) arg1 = k;
   if (arg1==(-5555)) arg1 = ii;
   if (arg1==(-4444)) arg1 = jj;
   if (arg1==(-3333)) arg1 = kk;
   if (arg2==(-9999)) arg2 = nion;
   if (arg2==(-8888)) arg2 = i;
   if (arg2==(-7777)) arg2 = j;
   if (arg2==(-6666)) arg2 = k;
   if (arg2==(-5555)) arg2 = ii;
   if (arg2==(-4444)) arg2 = jj;
   if (arg2==(-3333)) arg2 = kk;
   if (arg3==(-9999)) arg3 = nion;
   if (arg3==(-8888)) arg3 = i;
   if (arg3==(-7777)) arg3 = j;
   if (arg3==(-6666)) arg3 = k;
   if (arg3==(-5555)) arg3 = ii;
   if (arg3==(-4444)) arg3 = jj;
   if (arg3==(-3333)) arg3 = kk;
   
   f = 0.0;
   switch (op) {
    
      case -1: /* + */
         f = nwpw_emachine(ln1,code,fconst,nion,rion) + nwpw_emachine(ln2,code,fconst,nion,rion);
         break;

      case -2: /* - */
         f = nwpw_emachine(ln1,code,fconst,nion,rion) - nwpw_emachine(ln2,code,fconst,nion,rion);
         break;

      case -3: /* * */
         f = nwpw_emachine(ln1,code,fconst,nion,rion) * nwpw_emachine(ln2,code,fconst,nion,rion);
         break;

      case -4: /* / */
         f = nwpw_emachine(ln1,code,fconst,nion,rion) / nwpw_emachine(ln2,code,fconst,nion,rion);
         break;

      case -5: /* ** or ^ */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         b = nwpw_emachine(ln2,code,fconst,nion,rion);
         f = pow(a,b);
         break;

      case -6: /* abs */
         f = fabs(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -7: /* exp */
         f = exp(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -8: /* log */
         f = log(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -9: /* sqrt */
         f = sqrt(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -10: /* sinh */
         f = sinh(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -11: /* cosh */
         f = cosh(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -12: /* tanh */
         f = tanh(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -13: /* sin */
         f = sin(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -14: /* cos */
         f = cos(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -15: /* tan */
         f = tan(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -16: /* asin */
         f = asin(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -17: /* acos */
         f = acos(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -18: /* atan */
         f = atan(nwpw_emachine(ln1,code,fconst,nion,rion));
         break;

      case -19: /* sum(i=arg1,arg2:arg3) */
         for (i=arg1; i<=arg2; i+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -20: /* sum(j=arg1,arg2:arg3) */
         for (j=arg1; j<=arg2; j+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -21: /* sum(k=arg1,arg2:arg3) */
         for (k=arg1; k<=arg2; k+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -22: /* sum(ii=arg1,arg2:arg3) */
         for (ii=arg1; ii<=arg2; ii+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -23: /* sum(jj=arg1,arg2:arg3) */
         for (jj=arg1; jj<=arg2; jj+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -24: /* sum(kk=arg1,arg2:arg3) */
         for (kk=arg1; kk<=arg2; kk+=arg3) f += nwpw_emachine(ln1,code,fconst,nion,rion);
         break;

      case -25: /* fconst(arg1) */
         f = fconst[arg1];
         break;

      case -30: /* rion(arg1,arg2) */
         f = rion[(arg1-1)+3*(arg2-1)];

      case -31: /* d(arg1,arg2) */
         x = rion[0+3*(arg1-1)] - rion[0+3*(arg2-1)];
         y = rion[1+3*(arg1-1)] - rion[1+3*(arg2-1)];
         z = rion[2+3*(arg1-1)] - rion[2+3*(arg2-1)];
         lattice_min_difference_(&x,&y,&z);
         f = sqrt(x*x + y*y + z*z);
         break;


      case -40: /* i */
         f = ((double) i);
         break;

      case -41: /* j */
         f = ((double) j);
         break;

      case -42: /* k */
         f = ((double) k);
         break;

      case -43: /* ii */
         f = ((double) ii);
         break;

      case -44: /* jj */
         f = ((double) jj);
         break;

      case -45: /* kk */
         f = ((double) kk);
         break;

       default:
          f = 0.0;
          
   }

   return f;
}



/********************************************
 *                                          *
 *                nwpw_fmachine             *
 *                                          *
 ********************************************/

double nwpw_fmachine(Integer i0, Integer ii0, Integer ln, Integer code[], double fconst[], Integer nion, double rion[])
{
   Integer op,arg1,arg2,arg3,arg4,ln1,ln2;
   double df,a,b,c,d,bb,aa,fac1,fac2,x,y,z,xx[3];

   df = 0.0;
   op  = code[5*ln];
   ln1 = code[5*ln+1];
   ln2 = arg1 = code[5*ln+2];
   arg2       = code[5*ln+3];
   arg3       = code[5*ln+4];

   if (arg1==(-9999)) arg1 = nion;
   if (arg1==(-8888)) arg1 = i;
   if (arg1==(-7777)) arg1 = j;
   if (arg1==(-6666)) arg1 = k;
   if (arg1==(-5555)) arg1 = ii;
   if (arg1==(-4444)) arg1 = jj;
   if (arg1==(-3333)) arg1 = kk;
   if (arg2==(-9999)) arg2 = nion;
   if (arg2==(-8888)) arg2 = i;
   if (arg2==(-7777)) arg2 = j;
   if (arg2==(-6666)) arg2 = k;
   if (arg2==(-5555)) arg2 = ii;
   if (arg2==(-4444)) arg2 = jj;
   if (arg2==(-3333)) arg2 = kk;
   if (arg3==(-9999)) arg3 = nion;
   if (arg3==(-8888)) arg3 = i;
   if (arg3==(-7777)) arg3 = j;
   if (arg3==(-6666)) arg3 = k;
   if (arg3==(-5555)) arg3 = ii;
   if (arg3==(-4444)) arg3 = jj;
   if (arg3==(-3333)) arg3 = kk;

   df = 0.0;
   switch (op) {
    
      case -1: /* + */
         df = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion) + nwpw_fmachine(i0,ii0,ln2,code,fconst,nion,rion);
         break;

      case -2: /* - */
         df = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion) - nwpw_fmachine(i0,ii0,ln2,code,fconst,nion,rion);
         break;

      case -3: /* * */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         b = nwpw_emachine(ln2,code,fconst,nion,rion);
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         d = nwpw_fmachine(i0,ii0,ln2,code,fconst,nion,rion);
         df = c*b + a*d;
         break;

      case -4: /* / */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         b = nwpw_emachine(ln2,code,fconst,nion,rion);
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         d = nwpw_fmachine(i0,ii0,ln2,code,fconst,nion,rion);
         bb = b*b; if (bb<1.0e-8) bb = 1.0e-8;
         df = (b*c - a*d)/(bb);
         break;

      case -5: /* ** or ^ */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         b = nwpw_emachine(ln2,code,fconst,nion,rion);
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         d = nwpw_fmachine(i0,ii0,ln2,code,fconst,nion,rion);
         fac1 = pow(a,b-1.0);
         if (fabs(d)>1.0e-10)
            fac2 = (b*c + a*log(a)*d);
         else
            fac2 = b*c;
         df = fac1*fac2;
         break;

      case -6: /* abs */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         aa = fabs(a); if (aa<1.0e-8) aa = 1.0e-8;
         df = (a*c)/aa;
         break;

      case -7: /* exp */
         a = nwpw_emachine(ln1,code,fconst,nion,rion);
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = exp(a)*c;
         break;

      case -8: /* log */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         if (fabs(a)<1.0e-8)
         {
            if (a<0.0) 
               a = -1.0e-8;
            else
               a = 1.0e-8;
         }
         df = c/a;
         break;

      case -9: /* sqrt */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = 0.5*c/sqrt(a);
         break;

      case -10: /* sinh */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = c*cosh(a);
         break;

      case -11: /* cosh */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = c*sinh(a);
         break;

      case -12: /* tanh */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         b = cosh(a);
         if (fabs(b)<1.0e-8) b = 1.0e-8;
         df = c/(b*b);
         break;

      case -13: /* sin */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = c*cos(a);
         break;

      case -14: /* cos */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = -c*sin(a);
         break;

      case -15: /* tan */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         b = cos(a);
         if (fabs(b)<1.0e-8) b = 1.0e-8;
         df = c/(b*b);
         break;

      case -16: /* asin */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         b = 1.0-a*a;
         if (b<1.0e-8) b = 1.0e-8;
         b = sqrt(b);
         df = c/b;
         break;

      case -17: /* acos */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         b = 1.0-a*a;
         if (b<1.0e-8) b = 1.0e-8;
         b = sqrt(b);
         df = -c/b;
         break;

      case -18: /* atan */
         a = nwpw_emachine(ln1,code,fconst,nion,rion); 
         c = nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         df = c/(1+a*a);
         break;

      case -19: /* sum(i=arg1,arg2:arg3) */
         for (i=arg1; i<=arg2; i+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -20: /* sum(j=arg1,arg2:arg3) */
         for (j=arg1; j<=arg2; j+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -21: /* sum(k=arg1,arg2:arg3) */
         for (k=arg1; k<=arg2; k+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -22: /* sum(ii=arg1,arg2:arg3) */
         for (ii=arg1; ii<=arg2; ii+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -23: /* sum(jj=arg1,arg2:arg3) */
         for (jj=arg1; jj<=arg2; jj+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -24: /* sum(kk=arg1,arg2:arg3) */
         for (kk=arg1; kk<=arg2; kk+=arg3) df += nwpw_fmachine(i0,ii0,ln1,code,fconst,nion,rion);
         break;

      case -25: /* fconst(arg1) */
         df = 0.0;
         break;

      case -30: /* rion(arg1,arg2) */
         if ((arg1==i0) && (arg2==ii0))
            df = 1.0;
         else 
            df = 0.0;
         break;


      case -31: /* d(arg1,arg2) */
         if (arg1==ii0)
         {
            x = rion[0+3*(arg1-1)] - rion[0+3*(arg2-1)];
            y = rion[1+3*(arg1-1)] - rion[1+3*(arg2-1)];
            z = rion[2+3*(arg1-1)] - rion[2+3*(arg2-1)];
            lattice_min_difference_(&x,&y,&z);
            xx[0] = x; xx[1] = y; xx[2] = z;
            df = xx[i0]/sqrt(x*x + y*y + z*z);
         }
         else if (arg2==ii0)
         {
            x = rion[0+3*(arg2-1)] - rion[0+3*(arg1-1)];
            y = rion[1+3*(arg2-1)] - rion[1+3*(arg1-1)];
            z = rion[2+3*(arg2-1)] - rion[2+3*(arg1-1)];
            lattice_min_difference_(&x,&y,&z);
            xx[0] = x; xx[1] = y; xx[2] = z;
            df = xx[i0]/sqrt(x*x + y*y + z*z);
         }
         break;


      case -40: /* i */
         df = 0.0;
         break;

      case -41: /* j */
         df = 0.0;
         break;

      case -42: /* k */
         df = 0.0;
         break;

      case -43: /* ii */
         df = 0.0;
         break;
      
      case -44: /* jj */
         df = 0.0;
         break;

      case -45: /* kk */
         df = 0.0;
         break;


       default:
          df = 0.0;
          
   }

   return df;
}


/*************************  emachine routines *******************************/



/********************************************************
 *                                                      *
 *               nwpw_emachine_parse                    *
 *                                                      *
 ********************************************************/

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
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



/********************************************************
 *                                                      *
 *               nwpw_emachine_f                        *
 *                                                      *
 ********************************************************/

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define nwpw_emachine_f_ nwpw_emachine_f
#endif

double FATR nwpw_emachine_f_
#if defined(USE_FCD)
(
 Integer *nc0,
 Integer code0[],
 Integer *nf0,
 double fconst0[],
 Integer *nion0,
 double  *rion0[])
{

#else
(nc0,code0,nf0,fconst0,nion0,rion0)
Integer *nc0;
Integer code0[];
Integer *nf0;
double fconst0[];
Integer *nion0;
double rion0[];
{

#endif
   double f;

   f = nwpw_emachine(0,code0,fconst0,*nion0,rion0);

   return f;

}



/********************************************************
 *                                                      *
 *               nwpw_emachine_df                       *
 *                                                      *
 ********************************************************/

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define nwpw_emachine_df_ nwpw_emachine_df
#endif

double FATR nwpw_emachine_df_
#if defined(USE_FCD)
(
 Integer *i0,
 Integer *ii0,
 Integer *nc0,
 Integer code0[],
 Integer *nf0,
 double fconst0[],
 Integer *nion0,
 double  *rion0[])
{

#else
(i0,ii0,nc0,code0,nf0,fconst0,nion0,rion0)
Integer *i0;
Integer *ii0;
Integer *nc0;
Integer code0[];
Integer *nf0;
double fconst0[];
Integer *nion0;
double rion0[];
{

#endif
   double df;

   df = nwpw_fmachine(*i0,*ii0,0,code0,fconst0,*nion0,rion0);

   return df;

}


