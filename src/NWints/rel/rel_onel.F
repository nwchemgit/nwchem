C $Id: rel_onel.F,v 1.4 1999-04-20 17:39:19 mg201 Exp $
************************************************************************
c:tex-\subsection{rel\_onel}
c:tex-This routine generates the modified one-electron integrals for
c:tex-a relativistic basis set. These are the modified kinetic energy,
c:tex-the modified (scalar) potential energy and the modified overlap,
c:tex-\begin{eqnarray}
c:tex-&& \tilde{T}_{ab} = T_{ab}^{LS} + T_{ab}^{SL} - T_{ab}^{SS} 
c:tex- \nonumber \\ 
c:tex-&& \tilde{V}_{ab} = V_{ab}^{LL} + {{\alpha^2}\over{4}}
c:tex-   \nabla_A\cdot\nabla_B V_{ab}^{SS}
c:tex- \nonumber \\ 
c:tex-&& \tilde{S}_{ab} = S_{ab}^{LL} + {{\alpha^2}\over{2}} T_{ab}^{SS} 
c:tex- \nonumber 
c:tex-\end{eqnarray}
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_onel (
     &       Axyz,zeta_A,coefL_A,coefS_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &       Bxyz,zeta_B,coefL_B,coefS_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &       Cxyz,zan,nat,S,T,V,lstv,doStil,doTtil,doVtil,canAB,
     &       DryRun,scr,lscr,msg,ibug,ntyp)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "rel_consts.fh"
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
      integer nat      ! [input] number of atoms
      integer lscr     ! [input] size of scratch array
      integer lstv     ! [input] size of any integral buffer
      integer ibug    ! [input] debug variable
      integer ntyp    ! [input] potential energy integral type
      double precision Axyz(3)          ! [input] position of center A
      double precision zeta_A(n_prim_A) ! [input] exponents of shell A
      double precision coefL_A(n_prim_A,n_cont_A) ! [input] A large coeffs
      double precision coefS_A(n_prim_A,n_cont_A) ! [input] A small coeffs
      double precision Bxyz(3)          ! [input] position of center B
      double precision zeta_B(n_prim_B) ! [input] exponents of shell B 
      double precision coefL_B(n_prim_B,n_cont_B)  ! [input] B large coeffs
      double precision coefS_B(n_prim_B,n_cont_B)  ! [input] B small coeffs
      double precision Cxyz(3,nat)  ! [input] all atom positions
      double precision zan(nat)     ! [input] charges on all atoms
      double precision scr(lscr)    ! [scratch] scratch buffers
      double precision S(lstv)      ! [output] overlap integrals
      double precision T(lstv)      ! [output] kinetic energy integrals
      double precision V(lstv,ntyp) ! [output] potential integrals
      logical doStil     ! [input] compute modified overlap (True/False)
      logical doTtil     ! [input] compute modified kinetic (True/False)
      logical doVtil     ! [input] compute modified potential (True/False)
      logical canAB   ! [input] compute only canonical ints (false only)
      logical DryRun  ! [input] true means only compute required memory
      character*(*) msg ![input] calling func. identification message
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
      debug_arrays = mod(ibug,10)/2 .eq. 1
      max_mem = 0
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = n_cont_a*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_all_b = n_cart_b*n_cont_b
      n_ab = n_cart_ab*n_cont_ab
      if (debug_addresses) then
        write (LuOut,*) 'n_cart_a',n_cart_a
        write (LuOut,*) 'n_cart_b',n_cart_b
        write (LuOut,*) 'n_cart_ab',n_cart_ab
        write (LuOut,*) 'n_cont_ab',n_cont_ab
        write (LuOut,*) 'n_all_a',n_all_a
        write (LuOut,*) 'n_all_b',n_all_b
        write (LuOut,*) 'n_ab',n_ab
      end if
*
*   Calculate large component overlap and nuclear attraction integrals
*
      doS = doStil
      doT = .false.
      doV = doVtil
      memscr = lscr
      call int_hf1sp(
     &       Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &       Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &       Cxyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &       DryRun,scr,memscr,msg)
      if (DryRun) max_mem = max(max_mem,memscr)
*
*   Calculate kinetic energy integrals, correction to overlaps
*
      if (doTtil) then
        doS = .false.
        doT = .true.
        doV = .false.
        call int_hf1sp(
     &      Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &      Cxyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &      DryRun,scr,memscr,msg)
        if (DryRun) max_mem = max(max_mem,memscr)
        call int_hf1sp(
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &      Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &      Cxyz,zan,nat,S,scr,V,lstv,doS,doT,doV,canAB,
     &      DryRun,scr(n_ab+1),memscr,msg)
        if (DryRun) then
          max_mem = max(max_mem,memscr+n_ab)
        else
          call daxpy (n_ab,one,scr,1,T,1)
        end if
      end if
      if (doStil .or. doTtil) then
        doS = .false.
        doT = .true.
        doV = .false.
        call int_hf1sp(
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,ictr_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,ictr_B,
     &      Cxyz,zan,nat,S,scr,V,lstv,doS,doT,doV,canAB,
     &      DryRun,scr(n_ab+1),memscr,msg)
        if (DryRun) then
          max_mem = max(max_mem,memscr+n_ab)
        else
          if (doTtil) call daxpy (n_ab,-one,scr,1,T,1)
          if (doStil) call daxpy (n_ab,halsq,scr,1,S,1)
        end if
      end if
*
*   Generate small component potential arrays
*
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      if (.not.doVtil) return
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
      i_pp = max(i_xca+n_ab,i_xcb+n_prim_B*n_cont_B)
      i_pm = i_pp+n_intpp
      i_mp = i_pm+n_intpm
      i_mm = i_mp+n_intmp
      i_scr = i_mm+n_intmm
      memscr = lscr-i_scr+1

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
      if (.not.DryRun) then
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in rel_onel',99)
        k = i_xca-1
        do j = 1,n_cont_a
          do i = 1,n_prim_A
            scr(k+i) = zeta_A(i)*coefS_A(i,j)
          end do
          k = k+n_cont_A
        end do
        k = i_xcb-1
        do j = 1,n_cont_B
          do i = 1,n_prim_B
            scr(k+i) = zeta_B(i)*coefS_B(i,j)
          end do
          k = k+n_cont_A
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
     &    Cxyz,zan,nat,S,T,scr(i_pp),lstv,doS,doT,doV,canAB,
     &    DryRun,scr(i_scr),memscr,msg)
      if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
*
*         Calculate integrals for l_A-1, l_B+1
*
      if (l_A .gt. 0) then
        call int_hf1sp(
     &      Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,ictr_A,
     &      Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,ictr_B,
     &      Cxyz,zan,nat,S,T,scr(i_pm),lstv,doS,doT,doV,canAB,
     &      DryRun,scr(i_scr),memscr,msg)
        if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
      end if
*
*         Calculate integrals for l_A+1, l_B-1
*
      if (l_B .gt. 0) then
        call int_hf1sp(
     &      Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,ictr_A,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,ictr_B,
     &      Cxyz,zan,nat,S,T,scr(i_mp),lstv,doS,doT,doV,canAB,
     &      DryRun,scr(i_scr),memscr,msg)
        if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
*     
*         Calculate integrals for l_A-1, l_B-1
*
        if (l_A .gt. 0) then
          call int_hf1sp(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,ictr_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,ictr_B,
     &        Cxyz,zan,nat,S,T,scr(i_mm),lstv,doS,doT,doV,canAB,
     &        DryRun,scr(i_scr),memscr,msg)
          if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
        end if
      end if
*
*     Compute the relativistic potential energy integrals
*
      call rel_pot (scr(i_pp),scr(i_pm),scr(i_mp),scr(i_mm),
     &    scr,n_ab,ntyp,
     &    l_A,n_cartp_a,n_cart_a,n_cartm_a,n_cont_A,
     &    l_B,n_cartp_b,n_cart_b,n_cartm_b,n_cont_B,
     &    DryRun,scr(i_scr),memscr,ibug)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
      else
        i = 1
        do j = 1,ntyp
          call daxpy (n_ab,one,scr(i),1,V(1,j),1)
          i = i+n_ab
        end do
      end if
*
      return
      end

