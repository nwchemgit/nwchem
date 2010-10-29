C $Id$
************************************************************************
c:tex-\subsection{rel\_2e4cd\_sf}
c:tex-Wrapper routine for evaluation of the two-electron four-center
c:tex-spin-free relativistically-corrected derivative integrals. This
c:tex-only includes contributions from the Coulomb operator. Contributions 
c:tex-from the Gaunt operator could be added at a later stage.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_2e4cd_sf (
     &    Axyz,zeta_A,coefL_A,coefS_A,n_prim_A,n_cont_A,l_A,i_cent_A,
     &    Bxyz,zeta_B,coefL_B,coefS_B,n_prim_B,n_cont_B,l_B,i_cent_B,
     &    Cxyz,zeta_C,coefL_C,coefS_C,n_prim_C,n_cont_C,l_C,i_cent_C,
     &    Dxyz,zeta_D,coefL_D,coefS_D,n_prim_D,n_cont_D,l_D,i_cent_D,
     &    eri,leri,canAB,canCD,canPQ,DryRun,scr,lscr,
     &    AB_rel,CD_rel,SS_one_cent,do_SSSS,ibug)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "rel_consts.fh"
#include "errquit.fh"
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer n_prim_A ! [input] num. prims in shell A
      integer n_cont_A ! [input] num general conts in shell A
      integer l_A      ! [input] angular momentum of shell A
      integer i_cent_A ! [input] lexical atom index for shell A
      integer n_prim_B ! [input] num. prims in shell B
      integer n_cont_B ! [input] num general conts in shell B
      integer l_B      ! [input] angular momentum of shell B 
      integer i_cent_B ! [input] lexical atom index for shell B
      integer n_prim_C ! [input] num. prims in shell C
      integer n_cont_C ! [input] num general conts in shell C
      integer l_C      ! [input] angular momentum of shell C 
      integer i_cent_C ! [input] lexical atom index for shell C
      integer n_prim_D ! [input] num. prims in shell D
      integer n_cont_D ! [input] num general conts in shell D
      integer l_D      ! [input] angular momentum of shell D 
      integer i_cent_D ! [input] lexical atom index for shell D
      integer leri     ! length of integral buffer
      integer lscr     ! length of scratch array
      integer ibug     ! [input] debug variable
      double precision Axyz(3)          ! [input] position of center A
      double precision zeta_A(n_prim_A) ! [input] exponents of shell A
      double precision coefL_A(n_prim_A,n_cont_A) ! [input] A large coeffs
      double precision coefS_A(n_prim_A,n_cont_A) ! [input] A small coeffs
      double precision Bxyz(3)          ! [input] position of center B
      double precision zeta_B(n_prim_B) ! [input] exponents of shell B 
      double precision coefL_B(n_prim_B,n_cont_B)  ! [input] B large coeffs
      double precision coefS_B(n_prim_B,n_cont_B)  ! [input] B small coeffs
      double precision Cxyz(3)          ! [input] position of center C
      double precision zeta_C(n_prim_C) ! [input] exponents of shell C 
      double precision coefL_C(n_prim_C,n_cont_C)  ! [input] C large coeffs
      double precision coefS_C(n_prim_C,n_cont_C)  ! [input] C small coeffs
      double precision Dxyz(3)          ! [input] position of center D
      double precision zeta_D(n_prim_D) ! [input] exponents of shell D 
      double precision coefL_D(n_prim_D,n_cont_D)  ! [input] D large coeffs
      double precision coefS_D(n_prim_D,n_cont_D)  ! [input] D small coeffs
      double precision scr(lscr)    ! [scratch] scratch buffer
      double precision eri(leri,12)  ! [output] relativistic integrals
      logical canAB   ! [input] compute only canonical ints (false only)
      logical canCD   ! [input] compute only canonical ints (false only)
      logical canPQ   ! [input] compute only canonical ints (false only)
      logical DryRun  ! [input] true means only compute required memory
      logical AB_rel  ! [input] relativistic correction for A and B shells
      logical CD_rel  ! [input] relativistic correction for C and D shells
      logical SS_one_cent  ! [input] do relativistic correction only if
                           ! SS density is one-centre.
      logical do_SSSS ! [input] logical to do SSSS contributions
c:tex-\end{verbatim}
c:tex-
c:tex-{\it Subroutines called:} hf2, rel_LLSS, rel_SSLL, rel_SSSS
*
      integer iscr      ! address of free space in scr
      integer nscr      ! length of free space in scr
      integer neri      ! number of integrals per derivative
      integer neri12    ! number of derivative integrals
      integer max_mem   ! maximum memory used
      integer i         ! do loop index
      integer ibuf      ! integral buffer pointer
      logical do_CD     ! do relativistic CD corrections
      logical do_AB     ! do relativistic AB corrections
*
      max_mem = 0
      neri = n_cont_A*(l_A+1)*(l_A+2)/2
      neri = neri*n_cont_B*(l_B+1)*(l_B+2)/2
      neri = neri*n_cont_C*(l_C+1)*(l_C+2)/2
      neri = neri*n_cont_D*(l_D+1)*(l_D+2)/2
      if ((neri .gt. leri) .and. (.not.DryRun)) call errquit 
     &    ('rel_2e4cd_sf: eri buffer too small ',99, MEM_ERR)
      neri12 = neri*12
*
*   Do the LLLL integrals
*
      nscr = lscr
      call hf2d(
     &    Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,i_cent_A,
     &    Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,i_cent_B,
     &    Cxyz,zeta_C,coefL_C,n_prim_C,n_cont_C,l_C,i_cent_C,
     &    Dxyz,zeta_D,coefL_D,n_prim_D,n_cont_D,l_D,i_cent_D,
     &    eri,leri,canAB,canCD,canPQ,DryRun,scr,nscr)
      if (DryRun) max_mem = max(max_mem, nscr)
*
*  .. then the LLSS integrals
*
      iscr = neri12+1
      do_CD = CD_rel
      if (SS_one_cent) do_CD = do_CD .and. (i_cent_C .eq. i_cent_D)
      if (do_CD) then
        nscr = lscr-neri12
        call rel_LLSSd (
     &      Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,i_cent_A,
     &      Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,i_cent_B,
     &      Cxyz,zeta_C,coefS_C,n_prim_C,n_cont_C,l_C,i_cent_C,
     &      Dxyz,zeta_D,coefS_D,n_prim_D,n_cont_D,l_D,i_cent_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr)
        else
          ibuf = 1
          do i = 1,12
            call daxpy (neri,qalsq,scr(ibuf),1,eri(1,i),1)
            ibuf = ibuf+neri
          end do
        end if
      end if
*
*  .. and the SSLL integrals
*
      do_AB = AB_rel
      if (SS_one_cent) do_AB = do_AB .and. (i_cent_A .eq. i_cent_B)
      if (do_AB) then
        nscr = lscr-neri12
        call rel_SSLLd (
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,i_cent_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,i_cent_B,
     &      Cxyz,zeta_C,coefL_C,n_prim_C,n_cont_C,l_C,i_cent_C,
     &      Dxyz,zeta_D,coefL_D,n_prim_D,n_cont_D,l_D,i_cent_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr)
        else
          ibuf = 1
          do i = 1,12
            call daxpy (neri,qalsq,scr(ibuf),1,eri(1,i),1)
            ibuf = ibuf+neri
          end do
        end if
      end if
*
*  .. and finally the SSSS integrals
*
      if (do_AB .and. do_CD .and. do_SSSS) then
        nscr = lscr-neri12
        call rel_SSSSd (
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,i_cent_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,i_cent_B,
     &      Cxyz,zeta_C,coefS_C,n_prim_C,n_cont_C,l_C,i_cent_C,
     &      Dxyz,zeta_D,coefS_D,n_prim_D,n_cont_D,l_D,i_cent_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug,1,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr)
        else
          ibuf = 1
          do i = 1,12
            call daxpy (neri,salf,scr(ibuf),1,eri(1,i),1)
            ibuf = ibuf+neri
          end do
        end if
      end if
*
      if (DryRun) lscr = max_mem+neri
*
      return
      end
