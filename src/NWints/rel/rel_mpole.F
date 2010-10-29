C $Id$
************************************************************************
c:tex-\subsection{rel\_mpole}
c:tex-This routine generates the modified multipole integrals for
c:tex-a relativistic basis set. 
c:tex-\begin{equation}
c:tex- \tilde{L}^{sf}_{ab} = L_{ab}^{LL} + {{\alpha^2}\over{4}}
c:tex-   \nabla_A\cdot\nabla_B L_{ab}^{SS}
c:tex-\end{equation}
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_mpole (
     &    Axyz,zeta_A,coefL_A,coefS_A,n_prim_A,l_A,
     &    Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &    Bxyz,zeta_B,coefL_B,coefS_B,n_prim_B,l_B,
     &    mpints,lmpints,TriDiag,nonrel,DryRun,scr,lscr,ibug,ntyp)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "rel_consts.fh"
#include "errquit.fh"
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer n_prim_A ! [input] num. prims in shell A
      integer l_A      ! [input] angular momentum of shell A
      integer n_prim_C ! [input] num. prims in dummy shell C
      integer l_C      ! [input] angular momentum of shell C 
      integer n_prim_B ! [input] num. prims in shell B
      integer l_B      ! [input] angular momentum of shell B 
      integer lscr     ! [input] size of scratch array
      integer lmpints  ! [input] size of any integral buffer
      integer ibug    ! [input] debug variable
      integer ntyp    ! [input] integral type
      double precision Axyz(3)           ! [input] position of center A
      double precision zeta_A(n_prim_A)  ! [input] exponents of shell A
      double precision coefL_A(n_prim_A) ! [input] A large coeffs
      double precision coefS_A(n_prim_A) ! [input] A small coeffs
      double precision Bxyz(3)           ! [input] position of center B
      double precision zeta_B(n_prim_B)  ! [input] exponents of shell B 
      double precision coefL_B(n_prim_B) ! [input] B large coeffs
      double precision coefS_B(n_prim_B) ! [input] B small coeffs
      double precision Cxyz(3)           ! [input] position of center C
      double precision zeta_C(n_prim_C)  ! [input] exponents of dummy shell C
      double precision coef_C(n_prim_C)  ! [input] C dummy coeff
      double precision scr(lscr)         ! [scratch] scratch buffers
      double precision mpints(lmpints,ntyp)  ! [output] multipole integrals
      logical TriDiag ! [input] passed to hf3ois
      logical nonrel  ! [input] true if either centre is nonrelativistic
      logical DryRun  ! [input] true means only compute required memory
c:tex-\end{verbatim}
c:tex-See rel_pot for a description of the allowed values of ibug and ntyp
c:tex-
c:tex-{\it Subroutines called:} int\_hf1sp, rel\_pot, daxpy 
*
      integer n_cart_a  ! cartesian components of shell A
      integer n_cart_b  ! cartesian components of shell B
      integer n_cart_c  ! cartesian components of shell B
      integer n_cart_ab ! n_cart_a*n_cart_b
      integer n_abc     ! number of integrals
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
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
      max_mem = 0
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_c = (l_c+1)*(l_c+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_abc = n_cart_a*n_cart_b*n_cart_c
      if (lmpints .lt. n_abc .and. .not.DryRun) call errquit (
     &      'Integral buffer length too small in rel_onel',99, MEM_ERR)
      if (debug_addresses) then
        write (LuOut,*) 'n_cart_a',n_cart_a
        write (LuOut,*) 'n_cart_b',n_cart_b
        write (LuOut,*) 'n_abc',n_abc
      end if
*
*   Calculate large component integrals,
*
      memscr = lscr
      call hf3OIs(
     &    Axyz,zeta_A,coefL_A,n_prim_A,l_A,
     &    Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &    Bxyz,zeta_B,coefL_B,n_prim_B,l_B,
     &    mpints,n_abc,TriDiag,
     &    DryRun,scr,memscr)
      if (DryRun) max_mem = max(max_mem,memscr)
      if (nonrel) return
*
*   Generate small component integral arrays
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      n_cartp_a = n_cart_a+l_A+2
      n_cartp_b = n_cart_b+l_B+2
      n_cartm_a = n_cart_a-l_A-1
      n_cartm_b = n_cart_b-l_B-1
      n_intpp = n_cartp_a*n_cartp_b*n_cart_c
      n_intpm = n_cartm_a*n_cartp_b*n_cart_c
      n_intmp = n_cartp_a*n_cartm_b*n_cart_c
      n_intmm = n_cartm_a*n_cartm_b*n_cart_c
      i_xca = 1
      i_xcb = i_xca+n_prim_A
      i_pp = max(i_xca+n_abc*ntyp,i_xcb+n_prim_B)
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
     &      'Insufficient scratch memory in rel_onel',99, MEM_ERR)
        k = i_xca-1
        do i = 1,n_prim_A
          scr(k+i) = zeta_A(i)*coefS_A(i)
        end do
        k = i_xcb-1
        do i = 1,n_prim_B
          scr(k+i) = zeta_B(i)*coefS_B(i)
        end do
      end if
*
*         Calculate integrals for l_A+1, l_B+1
*
      call hf3OIs(
     &    Axyz,zeta_A,scr(i_xca),n_prim_A,l_A+1,
     &    Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &    Bxyz,zeta_B,scr(i_xcb),n_prim_B,l_B+1,
     &    scr(i_pp),n_intpp,TriDiag,
     &    DryRun,scr(i_scr),memscr)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        memscr = lscr-i_scr+1
      end if
*
*         Calculate integrals for l_A-1, l_B+1
*
      if (l_A .gt. 0) then
        call hf3OIs(
     &      Axyz,zeta_A,coefS_A,n_prim_A,l_A-1,
     &      Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &      Bxyz,zeta_B,scr(i_xcb),n_prim_B,l_B+1,
     &      scr(i_pm),n_intpm,TriDiag,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          memscr = lscr-i_scr+1
        end if
      end if
*
*         Calculate integrals for l_A+1, l_B-1
*
      if (l_B .gt. 0) then
        call hf3OIs(
     &      Axyz,zeta_A,scr(i_xca),n_prim_A,l_A+1,
     &      Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &      Bxyz,zeta_B,coefS_B,n_prim_B,l_B-1,
     &      scr(i_mp),n_intmp,TriDiag,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          memscr = lscr-i_scr+1
        end if
*     
*         Calculate integrals for l_A-1, l_B-1
*
        if (l_A .gt. 0) then
          call hf3OIs(
     &        Axyz,zeta_A,coefS_A,n_prim_A,l_A-1,
     &        Cxyz,zeta_C,coef_C,n_prim_C,l_C,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,l_B-1,
     &        scr(i_mm),n_intmm,TriDiag,
     &        DryRun,scr(i_scr),memscr)
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
     &    scr,n_abc,ntyp,
     &    l_A,n_cartp_a,n_cart_a,n_cartm_a,1,
     &    l_B,n_cartp_b,n_cart_b,n_cartm_b,n_cart_c,
     &    DryRun,scr(i_scr),memscr,ibug/10)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
      else
        i = 1
        do j = 1,ntyp
          call daxpy (n_abc,qalsq,scr(i),1,mpints(1,j),1)
          i = i+n_abc
        end do
      end if
*
      return
      end
