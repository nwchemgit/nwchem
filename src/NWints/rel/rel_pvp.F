C $Id$
************************************************************************
C> \ingroup nwint
C> @{
C>
C> \brief Generate the integral for various relativistic Hamiltonians
C>
C> This routine genrates the \f$pV\cdot p\f$ integrals needed for
C> various relativistic Hamiltonians.
C> \f{eqnarray*}{
C>   \langle a | pV\cdot p| b \rangle_{ab} &=& \nabla_A\cdot\nabla_B V_{ab}^{SS}
C> \f}
C>
C> Author: K. G. Dyall
C>
c:tex-\subsection{rel\_pvp}
c:tex-This routine generates the pV.p integrals needed for various
c:tex-relativistic Hamiltonians,
c:tex-\begin{equation}
c:tex- {\langle a | pV.p | b \rangle}_{ab} = 
c:tex-   \nabla_A\cdot\nabla_B V_{ab}^{SS}
c:tex- \nonumber 
c:tex-\end{equation}
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_pvp (
     &    Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &    Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &    Cxyz,zan,exinv,nat,pVp,lpvp,canAB,DryRun,scr,lscr,
     &    msg,ibug,ntyp)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "rel_consts.fh"
#include "errquit.fh"
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer n_prim_A !< [Input] num. prims in shell A
      integer n_cont_A !< [Input] num general conts in shell A
      integer l_A      !< [Input] angular momentum of shell A
      integer ictr_A   !< [Input] lexical atom index for shell A
      integer n_prim_B !< [Input] num. prims in shell B
      integer n_cont_B !< [Input] num general conts in shell B
      integer l_B      !< [Input] angular momentum of shell B 
      integer ictr_B   !< [Input] lexical atom index for shell B
      integer nat      !< [Input] number of atoms
      integer lscr     !< [Input] size of scratch array
      integer lpvp     !< [Input] size of integral buffer
      integer ibug     !< [Input] debug variable
      integer ntyp     !< [Input] potential energy integral type
      double precision Axyz(3)          !< [Input] position of center A
      double precision zeta_A(n_prim_A) !< [Input] exponents of shell A
      double precision coefS_A(n_prim_A,n_cont_A) !< [Input] A small coeffs
      double precision Bxyz(3)          !< [Input] position of center B
      double precision zeta_B(n_prim_B) !< [Input] exponents of shell B 
      double precision coefS_B(n_prim_B,n_cont_B) !< [Input] B small coeffs
      double precision Cxyz(3,nat)  !< [Input] all atom positions
      double precision zan(nat)     !< [Input] charges on all atoms
      double precision exinv(nat)   !< [Input] inverse nuclear exponents
      double precision scr(lscr)    !< [Scratch] scratch buffers
      double precision pVp(lpvp,ntyp) !< [Output] pV.p integrals
      logical canAB   !< [Input] compute only canonical ints (false only)
      logical DryRun  !< [Input] true means only compute required memory
      character*(*) msg !< [Input] calling func. identification message
c:tex-\end{verbatim}
c:tex-See rel_pot for a description of the allowed values of ibug and ntyp
c:tex-Note that in the current version of this routine, the call to rel_pot 
c:tex-uses a dummy ntyp=1. It is kept in the input so that in future, the 
c:tex-spin-orbit integrals can also be obtained with a call to this routine.
c:tex-
c:tex-{\it Subroutines called:} int\_hf1sp, rel\_pot, daxpy 
*
      integer n_cart_a  ! cartesian components of shell A
      integer n_cart_b  ! cartesian components of shell B
      integer n_cart_ab ! n_cart_a*n_cart_b
      integer n_cont_ab ! n_cont_a*n_cont_b
      integer n_all_b   ! n_cart_b*n_cont_b
      integer n_all_a   ! n_cart_a*n_cont_a
      integer n_ab      ! number of integrals
      integer n_cartp_a ! cartesian components for l_A+1
      integer n_cartp_b ! cartesian components for l_B+1
      integer n_cartm_a ! cartesian components for l_A-1
      integer n_cartm_b ! cartesian components for l_B-1
      integer n_intpp   ! number of integrals for l_A+1,l_B+1
      integer n_intpm   ! number of integrals for l_A-1,l_B+1
      integer n_intmp   ! number of integrals for l_A+1,l_B-1
      integer n_intmm   ! number of integrals for l_A-1,l_B-1
      integer i_xca     ! address in scr of exp*coef for shell A
      integer i_xcb     ! address in scr of exp*coef for shell B
      integer i_pp      ! address in scr of integrals for l_A+1,l_B+1
      integer i_pm      ! address in scr of integrals for l_A-1,l_B+1
      integer i_mp      ! address in scr of integrals for l_A+1,l_B-1
      integer i_mm      ! address in scr of integrals for l_A-1,l_B-1
      integer i_scr     ! address of free space in scr
      integer memscr    ! free space in scr
      integer max_mem   ! maximum memory used
      integer i,j,k     ! loop indices etc.
      double precision one ! Obvious!
      parameter (one = 1.0D0)
*
      logical debug_gen       ! do general debug printing
      logical debug_addresses ! do address debug printing
      logical debug_arrays    ! do array debug printing
      logical doS     ! compute overlap (True/False)
      logical doT     ! compute kinetic (True/False)
      logical doV     ! compute potential (True/False)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
      max_mem = 0
*
      if (debug_gen) then
        write (LuOut,*) 'l_A',l_A
        write (LuOut,*) 'n_prim_A',n_prim_A
        write (LuOut,*) 'n_cont_A',n_cont_A
        write (LuOut,*) 'ictr_A',ictr_A
        write (LuOut,*) 'l_B',l_B
        write (LuOut,*) 'n_prim_B',n_prim_B
        write (LuOut,*) 'n_cont_B',n_cont_B
        write (LuOut,*) 'ictr_B',ictr_B
        write (LuOut,*) 'msg ',msg
      end if
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = n_cont_a*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_all_b = n_cart_b*n_cont_b
      n_ab = n_cart_ab*n_cont_ab
      if (lpvp .lt. n_ab .and. .not.DryRun) call errquit (
     &      'Integral buffer length too small in rel_pvp',99, INT_ERR)
      if (debug_addresses) then
        write (LuOut,*) 'n_cart_a',n_cart_a
        write (LuOut,*) 'n_cart_b',n_cart_b
        write (LuOut,*) 'n_cart_ab',n_cart_ab
        write (LuOut,*) 'n_cont_ab',n_cont_ab
        write (LuOut,*) 'n_all_a',n_all_a
        write (LuOut,*) 'n_all_b',n_all_b
        write (LuOut,*) 'n_ab',n_ab
      end if
      if (debug_arrays) then
        call ecp_matpr (coefS_A,1,n_prim_a,1,n_cont_a,
     &      1,n_prim_a,1,n_cont_a,'S coef A','E',120,6)
        call ecp_matpr (coefS_B,1,n_prim_b,1,n_cont_b,
     &      1,n_prim_b,1,n_cont_b,'S coef B','E',120,6)
      end if
*
*   Generate small component potential arrays
*
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      n_cartp_a = n_cart_a+l_A+2
      n_cartp_b = n_cart_b+l_B+2
      n_cartm_a = n_cart_a-l_A-1
      n_cartm_b = n_cart_b-l_B-1
      n_intpp = n_cartp_a*n_cartp_b*n_cont_ab
      n_intpm = n_cartm_a*n_cartp_b*n_cont_ab
      n_intmp = n_cartp_a*n_cartm_b*n_cont_ab
      n_intmm = n_cartm_a*n_cartm_b*n_cont_ab
      i_xca = 1
      i_xcb = i_xca+n_prim_A*n_cont_A
      i_pp = i_xcb+n_prim_B*n_cont_B
      i_pm = i_pp+n_intpp
      i_mp = i_pm+n_intpm
      i_mm = i_mp+n_intmp
      i_scr = i_mm+n_intmm
*
      if (debug_addresses) then
        write (LuOut,*) 'n_cartp_a',n_cartp_a
        write (LuOut,*) 'n_cartp_b',n_cartp_b
        write (LuOut,*) 'n_cartm_a',n_cartm_a
        write (LuOut,*) 'n_cartm_b',n_cartm_b
        write (LuOut,*) 'n_intpp',n_intpp
        write (LuOut,*) 'n_intpm',n_intpm
        write (LuOut,*) 'n_intmp',n_intmp
        write (LuOut,*) 'n_intmm',n_intmm
        write (LuOut,*) 'i_xca,i_xcb',i_xca,i_xcb
        write (LuOut,*) 'i_pp,i_pm,i_mp,i_mm',i_pp,i_pm,i_mp,i_mm
        write (LuOut,*) 'i_scr',i_scr
      end if
*
*   Set up coefficients multiplied by exponents
*
      memscr = lscr-i_scr+1
      if (.not.DryRun) then
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in rel_pvp',99, MEM_ERR)
        k = i_xca-1
        do j = 1,n_cont_a
          do i = 1,n_prim_A
            scr(k+i) = zeta_A(i)*coefS_A(i,j)
          end do
          k = k+n_prim_A
        end do
        k = i_xcb-1
        do j = 1,n_cont_B
          do i = 1,n_prim_B
            scr(k+i) = zeta_B(i)*coefS_B(i,j)
          end do
          k = k+n_prim_A
        end do
      end if
      doS = .false.
      doT = .false.
      doV = .true.
*
*         Calculate integrals for l_A+1, l_B+1
*
      call int_hf1sp(
     &    Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,ictr_A,
     &    Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,ictr_B,
     &    Cxyz,zan,exinv,nat,scr,scr,scr(i_pp),n_intpp,doS,doT,doV,
     &    canAB,DryRun,scr(i_scr),memscr,msg)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        memscr = lscr-i_scr+1
      end if
*
*         Calculate integrals for l_A-1, l_B+1
*
      if (l_A .gt. 0) then
        call int_hf1sp(
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,ictr_A,
     &      Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,ictr_B,
     &      Cxyz,zan,exinv,nat,scr,scr,scr(i_pm),n_intpm,doS,doT,doV,
     &      canAB,DryRun,scr(i_scr),memscr,msg)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          memscr = lscr-i_scr+1
        end if
      end if
*
*         Calculate integrals for l_A+1, l_B-1
*
      if (l_B .gt. 0) then
        call int_hf1sp(
     &      Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,ictr_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,ictr_B,
     &      Cxyz,zan,exinv,nat,scr,scr,scr(i_mp),n_intmp,doS,doT,doV,
     &      canAB,DryRun,scr(i_scr),memscr,msg)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          memscr = lscr-i_scr+1
        end if
*     
*         Calculate integrals for l_A-1, l_B-1
*
        if (l_A .gt. 0) then
          call int_hf1sp(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,ictr_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,ictr_B,
     &        Cxyz,zan,exinv,nat,scr,scr,scr(i_mm),n_intmm,doS,doT,doV,
     &        canAB,DryRun,scr(i_scr),memscr,msg)
          if (DryRun) then
            max_mem = max(max_mem,i_scr+memscr-1)
            memscr = lscr-i_scr+1
          end if
        end if
      end if
*
*     Compute the relativistic potential energy integrals
*
      call rel_pot (scr(i_pp),scr(i_pm),scr(i_mp),scr(i_mm),
     &    pVp,lpvp,ntyp,
     &    l_A,n_cartp_a,n_cart_a,n_cartm_a,n_cont_A,
     &    l_B,n_cartp_b,n_cart_b,n_cartm_b,n_cont_B,
     &    DryRun,scr(i_scr),memscr,ibug/10)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
      else if (debug_arrays) then
        call ecp_matpr (pVp,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'pV.p integrals','E',120,6)
      end if
*
      return
      end
C>
C> @}
