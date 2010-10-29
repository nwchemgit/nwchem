C $Id$
************************************************************************
c:tex-\subsection{rel\_LLSSd}
c:tex-This routine generates the LLSS-type two-electron derivative integrals
c:tex-for a relativistic basis set, within the modified Dirac formalism.
c:tex-The integrals returned are the derivatives of 
c:tex-$\nabla_C\cdot\nabla_D (ab|c^Sd^S)$,
c:tex-$\nabla_C\times\nabla_D (ab|c^Sd^S)$, or the 9 derivatives.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_LLSSd (
     &    Axyz,zeta_A,coef_A,NPA,NCA,l_A,ictra,
     &    Bxyz,zeta_B,coef_B,NPB,NCB,l_B,ictrb,
     &    Cxyz,zeta_C,coef_C,NPC,NCC,l_C,ictrc,
     &    Dxyz,zeta_D,coef_D,NPD,NCD,l_D,ictrd,
     &    rel_ints,n_ints,canAB,canCD,canPQ,
     &    DryRun,scr,lscr,ibug,ntyp)
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "rel_consts.fh"
#include "errquit.fh"
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer NPA     ! [input] num. prims in shell A
      integer NCA     ! [input] num general conts in shell A
      integer l_A     ! [input] angular momentum of shell A
      integer ictra   ! [input] lexical atom index of shell A
      integer NPB     ! [input] num. prims in shell B
      integer NCB     ! [input] num general conts in shell B
      integer l_B     ! [input] angular momentum of shell B 
      integer ictrb   ! [input] lexical atom index of shell B
      integer NPC     ! [input] num. prims in shell C
      integer NCC     ! [input] num general conts in shell C
      integer l_C     ! [input] angular momentum of shell C
      integer ictrc   ! [input] lexical atom index of shell C
      integer NPD     ! [input] num. prims in shell D
      integer NCD     ! [input] num general conts in shell D
      integer l_D     ! [input] angular momentum of shell D 
      integer ictrd   ! [input] lexical atom index of shell D
      integer lscr    ! [input] size of scratch array
      integer n_ints  ! [input] size of any integral buffer
      integer ibug    ! [input] debug variable
      integer ntyp    ! [input] potential energy integral type
      double precision Axyz(3)         ! [input] position of center A
      double precision zeta_A(NPA)     ! [input] exponents of shell A
      double precision coef_A(NPA,NCA) ! [input] A large coeffs
      double precision Bxyz(3)         ! [input] position of center B
      double precision zeta_B(NPB)     ! [input] exponents of shell B 
      double precision coef_B(NPB,NCB) ! [input] B large coeffs
      double precision Cxyz(3)         ! [input] position of center C
      double precision zeta_C(NPC)     ! [input] exponents of shell C
      double precision coef_C(NPC,NCC) ! [input] C small coeffs
      double precision Dxyz(3)         ! [input] position of center D
      double precision zeta_D(NPD)     ! [input] exponents of shell D 
      double precision coef_D(NPD,NCD) ! [input] D small coeffs
      double precision scr(lscr)       ! [scratch] scratch buffers
      double precision rel_ints(n_ints,12,ntyp)  ! [output] LLSS integrals
      logical canAB   ! [input] compute only canonical ints (false only)
      logical canCD   ! [input] compute only canonical ints (false only)
      logical canPQ   ! [input] compute only canonical ints (false only)
      logical DryRun  ! [input] true means only compute required memory
c:tex-\end{verbatim}
c:tex-See rel_pot for a description of the allowed values of ibug and ntyp
c:tex-Note that in the current version of this routine, the call to rel_pot 
c:tex-uses a dummy ntyp=1. It is kept in the input so that in future, the 
c:tex-spin-orbit integrals can also be obtained with a call to this routine.
c:tex-
c:tex-{\it Subroutines called:} hf2d, rel\_pot, daxpy 
*
      integer n_cart_a   ! cartesian components of shell A
      integer n_cart_b   ! cartesian components of shell B
      integer n_cart_c   ! cartesian components of shell C
      integer n_cart_d   ! cartesian components of shell D
      integer n_cart_cd  ! n_cart_c*n_cart_d
      integer n_cont_cd  ! NCC*NCD
      integer n_all_a    ! n_cart_a*NCA
      integer n_all_b    ! n_cart_b*NCB
      integer n_all_c    ! n_cart_c*NCC
      integer n_all_d    ! n_cart_d*NCD
      integer n_ab       ! number of ab densities
      integer n_cd       ! number of cd densities
      integer n_abcd     ! number of integrals
      integer n_cartp_c  ! cartesian components for l_C+1
      integer n_cartp_d  ! cartesian components for l_D+1
      integer n_cartm_c  ! cartesian components for l_C-1
      integer n_cartm_d  ! cartesian components for l_D-1
      integer n_intpp    ! number of integrals for l_C+1,l_D+1
      integer n_intpm    ! number of integrals for l_C-1,l_D+1
      integer n_intmp    ! number of integrals for l_C+1,l_D-1
      integer n_intmm    ! number of integrals for l_C-1,l_D-1
      integer i_xcc      ! address in scr of exp*coef for shell C
      integer i_xcd      ! address in scr of exp*coef for shell D
      integer i_pp       ! address in scr of integrals for l_A+1,l_B+1
      integer i_pm       ! address in scr of integrals for l_A-1,l_B+1
      integer i_mp       ! address in scr of integrals for l_A+1,l_B-1
      integer i_mm       ! address in scr of integrals for l_A-1,l_B-1
      integer i_scr      ! address of free space in scr
      integer memscr     ! free space in scr
      integer max_mem    ! maximum memory used
      integer i,j,k      ! loop indices etc.
      double precision one ! Obvious!
      parameter (one = 1.0D0)
*
      logical debug_gen       ! do general debug printing
      logical debug_addresses ! do address debug printing
      logical debug_arrays    ! do array debug printing
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      max_mem = 0
      if (debug_gen) write (LuOut,*) 'Entering rel_LLSSd'
*
      if ((ntyp .ne. 1) .and. (ntyp .ne. 3) .and. (ntyp .ne. 4) .and.
     &    (ntyp .ne. 9))
     &    call errquit('Invalid value of ntyp in rel_LLSSd',99,
     &       UNKNOWN_ERR)
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_c = (l_c+1)*(l_c+2)/2
      n_cart_d = (l_d+1)*(l_d+2)/2
      n_cart_cd = n_cart_c*n_cart_d
      n_cont_cd = NCC*NCD
      n_all_a = n_cart_a*NCA
      n_all_b = n_cart_b*NCB
      n_all_c = n_cart_c*NCC
      n_all_d = n_cart_d*NCD
      n_ab = n_all_a*n_all_b
      n_cd = n_all_c*n_all_d
      n_abcd = n_ab*n_cd
      if ((n_ints .lt. n_abcd) .and. (.not.DryRun)) call errquit (
     &    'Integral buffer n_ints too small in rel_LLSSd',99,
     &       MEM_ERR)
      if (debug_addresses) then
        write (LuOut,*) 'n_cart_a',n_cart_a
        write (LuOut,*) 'n_cart_b',n_cart_b
        write (LuOut,*) 'n_cart_c',n_cart_c
        write (LuOut,*) 'n_cart_d',n_cart_d
        write (LuOut,*) 'n_cart_cd',n_cart_cd
        write (LuOut,*) 'n_cont_cd',n_cont_cd
        write (LuOut,*) 'n_all_a',n_all_a
        write (LuOut,*) 'n_all_b',n_all_b
        write (LuOut,*) 'n_all_c',n_all_c
        write (LuOut,*) 'n_all_d',n_all_d
        write (LuOut,*) 'n_ab',n_ab
        write (LuOut,*) 'n_cd',n_cd
      end if
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      n_cartp_c = n_cart_c+l_C+2
      n_cartp_d = n_cart_d+l_D+2
      n_cartm_c = n_cart_c-l_C-1
      n_cartm_d = n_cart_d-l_D-1
      n_intpp = n_cartp_c*n_cartp_d*n_cont_cd*n_ab
      n_intpm = n_cartm_c*n_cartp_d*n_cont_cd*n_ab
      n_intmp = n_cartp_c*n_cartm_d*n_cont_cd*n_ab
      n_intmm = n_cartm_c*n_cartm_d*n_cont_cd*n_ab
      i_xcc = 1
      i_xcd = i_xcc+NPC*NCC
      i_pp = i_xcd+NPD*NCD
      i_pm = i_pp+n_intpp*12
      i_mp = i_pm+n_intpm*12
      i_mm = i_mp+n_intmp*12
      i_scr = i_mm+n_intmm*12
      memscr = lscr-i_scr+1

      if (debug_addresses) then
        write (LuOut,*) 'n_cartp_c',n_cartp_c
        write (LuOut,*) 'n_cartp_d',n_cartp_d
        write (LuOut,*) 'n_cartm_c',n_cartm_c
        write (LuOut,*) 'n_cartm_d',n_cartm_d
        write (LuOut,*) 'n_intpp',n_intpp
        write (LuOut,*) 'n_intpm',n_intpm
        write (LuOut,*) 'n_intmp',n_intmp
        write (LuOut,*) 'n_intmm',n_intmm
        write (LuOut,*) 'i_xcc,i_xcd',i_xcc,i_xcd
        write (LuOut,*) 'i_pp,i_pm,i_mp,i_mm',i_pp,i_pm,i_mp,i_mm
        write (LuOut,*) 'i_scr',i_scr
        write (LuOut,*) 'memscr,lscr',memscr,lscr
      end if
*
*   Set up coefficients multiplied by exponents
*
      if (.not.DryRun) then
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in rel_LLSSd',99, MEM_ERR)
        k = i_xcc-1
        do j = 1,NCC
          do i = 1,NPC
            scr(k+i) = zeta_C(i)*coef_C(i,j)
          end do
          k = k+NPC
        end do
        k = i_xcd-1
        do j = 1,NCD
          do i = 1,NPD
            scr(k+i) = zeta_D(i)*coef_D(i,j)
          end do
          k = k+NPD
        end do
      end if
*
*         Calculate integrals for l_C+1, l_D+1
*
      if (debug_gen) write (LuOut,*) 'calling hf2d ++'
      call hf2d(
     &    Axyz,zeta_A,coef_A,NPA,NCA,l_A,ictra,
     &    Bxyz,zeta_B,coef_B,NPB,NCB,l_B,ictrb,
     &    Cxyz,zeta_C,scr(i_xcc),NPC,NCC,l_C+1,ictrc,
     &    Dxyz,zeta_D,scr(i_xcd),NPD,NCD,l_D+1,ictrd,
     &    scr(i_pp),n_intpp,canAB,canCD,canPQ,
     &    DryRun,scr(i_scr),memscr)
      if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
*
*         Calculate integrals for l_C-1, l_D+1
*
      if (l_C .gt. 0) then
        memscr = lscr-i_scr+1
        if (debug_gen) write (LuOut,*) 'calling hf2d -+'
        call hf2d(
     &      Axyz,zeta_A,coef_A,NPA,NCA,l_A,ictra,
     &      Bxyz,zeta_B,coef_B,NPB,NCB,l_B,ictrb,
     &      Cxyz,zeta_C,coef_C,NPC,NCC,l_C-1,ictrc,
     &      Dxyz,zeta_D,scr(i_xcd),NPD,NCD,l_D+1,ictrd,
     &      scr(i_pm),n_intpm,canAB,canCD,canPQ,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
      end if
*
*         Calculate integrals for l_C+1, l_D-1
*
      if (l_D .gt. 0) then
        memscr = lscr-i_scr+1
        if (debug_gen) write (LuOut,*) 'calling hf2d +-'
        call hf2d(
     &      Axyz,zeta_A,coef_A,NPA,NCA,l_A,ictra,
     &      Bxyz,zeta_B,coef_B,NPB,NCB,l_B,ictrb,
     &      Cxyz,zeta_C,scr(i_xcc),NPC,NCC,l_C+1,ictrc,
     &      Dxyz,zeta_D,coef_D,NPD,NCD,l_D-1,ictrd,
     &      scr(i_mp),n_intmp,canAB,canCD,canPQ,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
*
*         Calculate integrals for l_C-1, l_D-1
*
        if (l_C .gt. 0) then
          memscr = lscr-i_scr+1
          if (debug_gen) write (LuOut,*) 'calling hf2d --'
          call hf2d(
     &        Axyz,zeta_A,coef_A,NPA,NCA,l_A,ictra,
     &        Bxyz,zeta_B,coef_B,NPB,NCB,l_B,ictrb,
     &        Cxyz,zeta_C,coef_C,NPC,NCC,l_C-1,ictrc,
     &        Dxyz,zeta_D,coef_D,NPD,NCD,l_D-1,ictrd,
     &        scr(i_mm),n_intmm,canAB,canCD,canPQ,
     &        DryRun,scr(i_scr),memscr)
          if (DryRun) max_mem = max(max_mem,i_scr+memscr-1)
        end if
      end if
*
*     Compute the relativistic integrals
*
      memscr = lscr-i_scr+1
      if (debug_gen) write (LuOut,*) 'calling rel_pot'
      call rel_pot (
     &    scr(i_pp),scr(i_pm),scr(i_mp),scr(i_mm),
     &    rel_ints,n_ints*12,ntyp,
     &    l_C,n_cartp_c,n_cart_c,n_cartm_c,NCC*n_ab*12,
     &    l_D,n_cartp_d,n_cart_d,n_cartm_d,NCD,
     &    DryRun,scr(i_scr),memscr,ibug)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
      end if
      if (debug_gen) write (LuOut,*) 'Exiting rel_LLSSd'
*
      return
      end
