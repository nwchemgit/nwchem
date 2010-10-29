C $Id$
************************************************************************
c:tex-\subsection{rel\_2e4c\_sf}
c:tex-Wrapper routine for evaluation of the two-electron four-center
c:tex-spin-free relativistically-corrected integrals. This only includes
c:tex-contributions from the Coulomb operator. Contributions from the 
c:tex-Gaunt operator could be added at a later stage.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_2e4c_sf (
     &    Axyz,zeta_A,coefL_A,coefS_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &    Bxyz,zeta_B,coefL_B,coefS_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &    Cxyz,zeta_C,coefL_C,coefS_C,n_prim_C,n_cont_C,l_C,ictr_C,
     &    Dxyz,zeta_D,coefL_D,coefS_D,n_prim_D,n_cont_D,l_D,ictr_D,
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
      integer ictr_A   ! [input] lexical atom index for shell A
      integer n_prim_B ! [input] num. prims in shell B
      integer n_cont_B ! [input] num general conts in shell B
      integer l_B      ! [input] angular momentum of shell B 
      integer ictr_B   ! [input] lexical atom index for shell B
      integer n_prim_C ! [input] num. prims in shell C
      integer n_cont_C ! [input] num general conts in shell C
      integer l_C      ! [input] angular momentum of shell C 
      integer ictr_C   ! [input] lexical atom index for shell C
      integer n_prim_D ! [input] num. prims in shell D
      integer n_cont_D ! [input] num general conts in shell D
      integer l_D      ! [input] angular momentum of shell D 
      integer ictr_D   ! [input] lexical atom index for shell D
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
      double precision eri(leri)    ! [output] relativistic integrals
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
      integer ibug_local ! local debug parameter
      integer iscr      ! address of free space in scr
      integer nscr      ! length of free space in scr
      integer neri      ! number of integrals
      integer nallAB    ! number of densities in bra
      integer nallCD    ! number of densities in ket
      integer max_mem   ! maximum memory used
      logical do_CD     ! do relativistic CD corrections
      logical do_AB     ! do relativistic AB corrections
*
      logical debug_gen       ! do general debug printing
      logical debug_addresses ! do address debug printing
      logical debug_arrays    ! do array debug printing
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      ibug_local = ibug/10
      if (debug_gen) write (LuOut,*) 'Entering rel_2e4c_sf'
*
      max_mem = 0
      nallAB = n_cont_A*(l_A+1)*(l_A+2)/2
      nallAB = nallAB*n_cont_B*(l_B+1)*(l_B+2)/2
      nallCD = n_cont_C*(l_C+1)*(l_C+2)/2
      nallCD = nallCD*n_cont_D*(l_D+1)*(l_D+2)/2
      neri = nallAB*nallCD
      if (debug_addresses) then
        write (LuOut,*) l_A,ictr_A,n_prim_A,n_cont_A
        write (LuOut,*) l_B,ictr_B,n_prim_B,n_cont_B
        write (LuOut,*) l_C,ictr_C,n_prim_C,n_cont_C
        write (LuOut,*) l_D,ictr_D,n_prim_D,n_cont_D
      end if
      if ((neri .gt. leri) .and. (.not.DryRun)) call errquit 
     &    ('rel_2e4c_sf: eri buffer too small ',99, MEM_ERR)
*
*   Do the LLLL integrals
*
      nscr = lscr
      if (debug_gen) then
        write (LuOut,*) 'Calling hf2'
        call util_flush(LuOut)
      end if
      call hf2(
     &    Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &    Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &    Cxyz,zeta_C,coefL_C,n_prim_C,n_cont_C,l_C,
     &    Dxyz,zeta_D,coefL_D,n_prim_D,n_cont_D,l_D,
     &    eri,leri,canAB,canCD,canPQ,DryRun,scr,nscr)
      if (DryRun) max_mem = max(max_mem, nscr)
      if (debug_arrays) call ecp_matpr(eri,1,nallAB,1,nallCD,
     &    1,nallAB,1,nallCD,'LLLL ints','E',120,6)
*
*  .. then the LLSS integrals
*
      iscr = neri+1
      do_CD = CD_rel
      if (SS_one_cent) do_CD = do_CD .and. (ictr_C .eq. ictr_D)
      if (do_CD) then
        nscr = lscr-neri
        if (debug_gen) then
          write (LuOut,*) 'Calling rel_LLSS'
          call util_flush(LuOut)
        end if
        call rel_LLSS (
     &      Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &      Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &      Cxyz,zeta_C,coefS_C,n_prim_C,n_cont_C,l_C,
     &      Dxyz,zeta_D,coefS_D,n_prim_D,n_cont_D,l_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug_local,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr+neri)
        else
          if (debug_arrays) call ecp_matpr(scr,1,nallAB,1,nallCD,
     &        1,nallAB,1,nallCD,'LLSS ints','E',120,6)
          call daxpy (neri,qalsq,scr,1,eri,1)
        end if
      end if
*
*  .. and the SSLL integrals
*
      do_AB = AB_rel
      if (SS_one_cent) do_AB = do_AB .and. (ictr_A .eq. ictr_B)
      if (do_AB) then
        nscr = lscr-neri
        if (debug_gen) then
          write (LuOut,*) 'Calling rel_SSLL'
          call util_flush(LuOut)
        end if
        call rel_SSLL (
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &      Cxyz,zeta_C,coefL_C,n_prim_C,n_cont_C,l_C,
     &      Dxyz,zeta_D,coefL_D,n_prim_D,n_cont_D,l_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug_local,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr+neri)
        else
          if (debug_arrays) call ecp_matpr(scr,1,nallAB,1,nallCD,
     &        1,nallAB,1,nallCD,'SSLL ints','E',120,6)
          call daxpy (neri,qalsq,scr,1,eri,1)
        end if
      end if
*
*  .. and finally the SSSS integrals
*
      if (do_AB .and. do_CD .and. do_SSSS) then
        nscr = lscr-neri
        if (debug_gen) then
          write (LuOut,*) 'Calling rel_SSSS'
          call util_flush(LuOut)
        end if
        call rel_SSSS (
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &      Cxyz,zeta_C,coefS_C,n_prim_C,n_cont_C,l_C,
     &      Dxyz,zeta_D,coefS_D,n_prim_D,n_cont_D,l_D,
     &      scr,neri,canAB,canCD,canPQ,
     &      DryRun,scr(iscr),nscr,ibug_local,1,1)
        if (DryRun) then
          max_mem = max(max_mem, nscr+neri)
        else
          if (debug_arrays) call ecp_matpr(scr,1,nallAB,1,nallCD,
     &        1,nallAB,1,nallCD,'SSSS ints','E',120,6)
          call daxpy (neri,salf,scr,1,eri,1)
        end if
      end if
*
      if (DryRun) lscr = max_mem+neri
      if (debug_gen) write (LuOut,*) 'Exiting rel_2e4c_sf'
*
      return
      end
