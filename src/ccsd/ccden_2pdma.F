C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_2PDMoooo -- Form coupled cluster 2-particle density matrix
C
C REVISION
C     $Id: ccden_2pdma.F,v 2.1 1997-04-14 06:34:50 gg502 Exp $
C
C SYNOPSIS
      Subroutine CCDen_2PDMoooo(NOcc, NVir, g_Q)
      Implicit NONE
      Integer NOcc, NVir
      Integer g_Q
C
C ARGUMENTS
C DESCRIPTION
C     Compute the oooo part of the coupled cluster 2-PDM
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 46
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES

C
C     *********************************************************
C     * G(i,j,k,l) = 1/8 P(i,j,k,l) {1/2 [TauZ](i,k,j,l)      *
C     *              + 2 delta(k,l) ( Q(i,j) - delta(i,j) )   *
C     *              -   delta(j,l) ( Q(i,k) - delta(i,k) ) } *
C     * APR & TJL Eq. 46                                      *
C     *********************************************************
C
C     ****************************
C     * Form Q(p,q) - delta(p,q) *
C     ****************************
C
      Call GA_Duplicate(g_Q, g_Qm1, 'CCDen Modified 1PDM')
      Call GA_Copy(g_Q, g_Qm1)

C
      Return
      End
