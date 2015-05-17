C $Id$
************************************************************************
C> \ingroup nwint
C> @{
C>
C> \brief Generate the SSLL-type two-electron integrals
C>
C> This routine generates the SSLL-type two-electron integrals for
C> a relativistic basis set, within the modified Dirac formalism.
C> The integrals returned are \f$\nabla_A\cdot\nabla_B (a^Sb^S|cd)\f$,
C> \f$\nabla_A\times\nabla_B (a^Sb^S|cd)\f$, or the 9 derivatives.
C> 
C> Author: K. G. Dyall
C>
c:tex-\subsection{rel\_SSLL}
c:tex-This routine generates the SSLL-type two-electron integrals for
c:tex-a relativistic basis set, within the modified Dirac formalism.
c:tex-The integrals returned are $\nabla_A\cdot\nabla_B (a^Sb^S|cd)$,
c:tex-$\nabla_A\times\nabla_B (a^Sb^S|cd)$, or the 9 derivatives.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_SSLL (
     &    Axyz,zeta_A,coef_A,NPA,NCA,l_A,
     &    Bxyz,zeta_B,coef_B,NPB,NCB,l_B,
     &    Cxyz,zeta_C,coef_C,NPC,NCC,l_C,
     &    Dxyz,zeta_D,coef_D,NPD,NCD,l_D,
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
      integer NPA     !< [Input] num. prims in shell A
      integer NCA     !< [Input] num general conts in shell A
      integer l_A     !< [Input] angular momentum of shell A
      integer NPB     !< [Input] num. prims in shell B
      integer NCB     !< [Input] num general conts in shell B
      integer l_B     !< [Input] angular momentum of shell B 
      integer NPC     !< [Input] num. prims in shell C
      integer NCC     !< [Input] num general conts in shell C
      integer l_C     !< [Input] angular momentum of shell C
      integer NPD     !< [Input] num. prims in shell D
      integer NCD     !< [Input] num general conts in shell D
      integer l_D     !< [Input] angular momentum of shell D 
      integer lscr    !< [Input] size of scratch array
      integer n_ints  !< [Input] size of any integral buffer
      integer ibug    !< [Input] debug variable
      integer ntyp    !< [Input] potential energy integral type
      double precision Axyz(3)         !< [Input] position of center A
      double precision zeta_A(NPA)     !< [Input] exponents of shell A
      double precision coef_A(NPA,NCA) !< [Input] A small coeffs
      double precision Bxyz(3)         !< [Input] position of center B
      double precision zeta_B(NPB)     !< [Input] exponents of shell B 
      double precision coef_B(NPB,NCB) !< [Input] B small coeffs
      double precision Cxyz(3)         !< [Input] position of center C
      double precision zeta_C(NPC)     !< [Input] exponents of shell C
      double precision coef_C(NPC,NCC) !< [Input] C small coeffs
      double precision Dxyz(3)         !< [Input] position of center D
      double precision zeta_D(NPD)     !< [Input] exponents of shell D 
      double precision coef_D(NPD,NCD) !< [Input] D small coeffs
      double precision scr(lscr)       !< [Scratch] scratch buffers
      double precision rel_ints(n_ints,ntyp)  !< [Output] LLSS integrals
      logical canAB   !< [Input] compute only canonical ints (false only)
      logical canCD   !< [Input] compute only canonical ints (false only)
      logical canPQ   !< [Input] compute only canonical ints (false only)
      logical DryRun  !< [Input] true means only compute required memory
c:tex-\end{verbatim}
c:tex-See rel_pot for a description of the allowed values of ibug and ntyp
c:tex-Note that in the current version of this routine, the call to rel_pot 
c:tex-uses a dummy ntyp=1. It is kept in the input so that in future, the 
c:tex-spin-orbit integrals can also be obtained with a call to this routine.
c:tex-
c:tex-{\it Subroutines called:} hf2, rel\_pot, daxpy 
*
      integer n_cart_a   ! cartesian components of shell A
      integer n_cart_b   ! cartesian components of shell B
      integer n_cart_c   ! cartesian components of shell C
      integer n_cart_d   ! cartesian components of shell D
      integer n_cart_ab  ! n_cart_a*n_cart_b
      integer n_cont_ab  ! NCA*NCB
      integer n_all_a    ! n_cart_a*NCA
      integer n_all_b    ! n_cart_b*NCB
      integer n_all_c    ! n_cart_c*NCC
      integer n_all_d    ! n_cart_d*NCD
      integer n_ab       ! number of ab densities
      integer n_cd       ! number of cd densities
      integer n_abcd     ! number of integrals
      integer n_cartp_a  ! cartesian components for l_A+1
      integer n_cartp_b  ! cartesian components for l_B+1
      integer n_cartm_a  ! cartesian components for l_A-1
      integer n_cartm_b  ! cartesian components for l_B-1
      integer n_intpp    ! number of integrals for l_A+1,l_B+1
      integer n_intpm    ! number of integrals for l_A-1,l_B+1
      integer n_intmp    ! number of integrals for l_A+1,l_B-1
      integer n_intmm    ! number of integrals for l_A-1,l_B-1
      integer i_xca      ! address in scr of exp*coef for shell A
      integer i_xcb      ! address in scr of exp*coef for shell B
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
      if (debug_gen) write (LuOut,*) 'Entering rel_SSLL'
*
      if ((ntyp .ne. 1) .and. (ntyp .ne. 3) .and. (ntyp .ne. 4) .and.
     &    (ntyp .ne. 9))
     &    call errquit('Invalid value of ntyp in rel_SSLL',99,
     &       UNKNOWN_ERR)
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_c = (l_c+1)*(l_c+2)/2
      n_cart_d = (l_d+1)*(l_d+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = NCA*NCB
      n_all_a = n_cart_a*NCA
      n_all_b = n_cart_b*NCB
      n_all_c = n_cart_c*NCC
      n_all_d = n_cart_d*NCD
      n_ab = n_all_a*n_all_b
      n_cd = n_all_c*n_all_d
      n_abcd = n_ab*n_cd
      if ((n_ints .lt. n_abcd) .and. (.not.DryRun)) call errquit (
     &    'Integral buffer n_ints too small in rel_SSLL',99, MEM_ERR)
      if (debug_addresses) then
        write (LuOut,*) 'n_cart_a',n_cart_a
        write (LuOut,*) 'n_cart_b',n_cart_b
        write (LuOut,*) 'n_cart_c',n_cart_c
        write (LuOut,*) 'n_cart_d',n_cart_d
        write (LuOut,*) 'n_cart_ab',n_cart_ab
        write (LuOut,*) 'n_cont_ab',n_cont_ab
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
      n_cartp_a = n_cart_a+l_A+2
      n_cartp_b = n_cart_b+l_B+2
      n_cartm_a = n_cart_a-l_A-1
      n_cartm_b = n_cart_b-l_B-1
      n_intpp = n_cartp_a*n_cartp_b*n_cont_ab*n_cd
      n_intpm = n_cartm_a*n_cartp_b*n_cont_ab*n_cd
      n_intmp = n_cartp_a*n_cartm_b*n_cont_ab*n_cd
      n_intmm = n_cartm_a*n_cartm_b*n_cont_ab*n_cd
      i_xca = 1
      i_xcb = i_xca+NPA*NCA
      i_pp = i_xcb+NPB*NCB
      i_pm = i_pp+n_intpp
      i_mp = i_pm+n_intpm
      i_mm = i_mp+n_intmp
      i_scr = i_mm+n_intmm
      memscr = lscr-i_scr+1
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
      if (.not.DryRun) then
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in rel_SSLL',99, MEM_ERR)
        k = i_xca-1
        do j = 1,NCA
          do i = 1,NPA
            scr(k+i) = zeta_A(i)*coef_A(i,j)
          end do
          k = k+NPA
        end do
        k = i_xcb-1
        do j = 1,NCB
          do i = 1,NPB
            scr(k+i) = zeta_B(i)*coef_B(i,j)
          end do
          k = k+NPB
        end do
      end if
*
*         Calculate integrals for l_A+1, l_B+1
*
      call hf2(
     &    Axyz,zeta_A,scr(i_xca),NPA,NCA,l_A+1,
     &    Bxyz,zeta_B,scr(i_xcb),NPB,NCB,l_B+1,
     &    Cxyz,zeta_C,coef_C,NPC,NCC,l_C,
     &    Dxyz,zeta_D,coef_D,NPD,NCD,l_D,
     &    scr(i_pp),n_intpp,canAB,canCD,canPQ,
     &    DryRun,scr(i_scr),memscr)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
      else if (debug_arrays) then
        call ecp_matpr(scr(i_pp),1,n_intpp,1,1,
     &      1,n_intpp,1,1,'++ ints','E',120,6)
      end if
*
*         Calculate integrals for l_A-1, l_B+1
*
      if (l_A .gt. 0) then
        memscr = lscr-i_scr+1
        call hf2(
     &      Axyz,zeta_A,coef_A,NPA,NCA,l_A-1,
     &      Bxyz,zeta_B,scr(i_xcb),NPB,NCB,l_B+1,
     &      Cxyz,zeta_C,coef_C,NPC,NCC,l_C,
     &      Dxyz,zeta_D,coef_D,NPD,NCD,l_D,
     &      scr(i_pm),n_intpm,canAB,canCD,canPQ,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
        else if (debug_arrays) then
          call ecp_matpr(scr(i_pm),1,n_intpm,1,1,
     &        1,n_intpm,1,1,'+- ints','E',120,6)
        end if
      end if
*
*         Calculate integrals for l_A+1, l_B-1
*
      if (l_B .gt. 0) then
        memscr = lscr-i_scr+1
        call hf2(
     &      Axyz,zeta_A,scr(i_xca),NPA,NCA,l_A+1,
     &      Bxyz,zeta_B,coef_B,NPB,NCB,l_B-1,
     &      Cxyz,zeta_C,coef_C,NPC,NCC,l_C,
     &      Dxyz,zeta_D,coef_D,NPD,NCD,l_D,
     &      scr(i_mp),n_intmp,canAB,canCD,canPQ,
     &      DryRun,scr(i_scr),memscr)
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
        else if (debug_arrays) then
          call ecp_matpr(scr(i_mp),1,n_intmp,1,1,
     &        1,n_intmp,1,1,'-+ ints','E',120,6)
        end if
*
*         Calculate integrals for l_A-1, l_B-1
*
        if (l_A .gt. 0) then
          memscr = lscr-i_scr+1
          call hf2(
     &        Axyz,zeta_A,coef_A,NPA,NCA,l_A-1,
     &        Bxyz,zeta_B,coef_B,NPB,NCB,l_B-1,
     &        Cxyz,zeta_C,coef_C,NPC,NCC,l_C,
     &        Dxyz,zeta_D,coef_D,NPD,NCD,l_D,
     &        scr(i_mm),n_intmm,canAB,canCD,canPQ,
     &        DryRun,scr(i_scr),memscr)
          if (DryRun) then
            max_mem = max(max_mem,i_scr+memscr-1)
          else if (debug_arrays) then
            call ecp_matpr(scr(i_mm),1,n_intmm,1,1,
     &          1,n_intmm,1,1,'-- ints','E',120,6)
          end if
        end if
      end if
*
*     Compute the relativistic integrals
*
      memscr = lscr-i_scr+1
      call rel_pot2 (scr(i_pp),scr(i_pm),scr(i_mp),scr(i_mm),
     &    rel_ints,n_ints,ntyp,n_cd,
     &    l_A,n_cartp_a,n_cart_a,n_cartm_a,NCA,
     &    l_B,n_cartp_b,n_cart_b,n_cartm_b,NCB,
     &    DryRun,scr(i_scr),memscr,ibug)
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
      end if
*
      return
      end
C>
C> @}
