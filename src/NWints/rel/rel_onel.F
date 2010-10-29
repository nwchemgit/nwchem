C $Id$
************************************************************************
c:tex-\subsection{rel\_onel}
c:tex-This routine generates the modified one-electron integrals for
c:tex-a relativistic basis set. These are the modified kinetic energy,
c:tex-the modified potential energy and the modified overlap,
c:tex-\begin{eqnarray}
c:tex-&& \tilde{T}_{ab} = T_{ab}^{LS} + T_{ab}^{SL} - T_{ab}^{SS} 
c:tex- \nonumber \\ 
c:tex-&& \tilde{V}^{sf}_{ab} = V_{ab}^{LL} + {{\alpha^2}\over{4}}
c:tex-   \nabla_A\cdot\nabla_B V_{ab}^{SS}
c:tex-&& \tilde{V}^{so}_{ab} = V_{ab}^{LL} + {{\alpha^2}\over{4}}
c:tex-   \nabla_A\times\nabla_B V_{ab}^{SS}
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
     &    Axyz,zeta_A,coefL_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &    Bxyz,zeta_B,coefL_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &    Cxyz,zan,exinv,nat,S,T,V,lstv,doStil,doTtil,doVtil,
     &    canAB,do_nw,do_hnd,nonrel,DryRun,scr,lscr,ibug,ntyp)
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
      integer n_prim_B ! [input] num. prims in shell B
      integer n_cont_B ! [input] num general conts in shell B
      integer l_B      ! [input] angular momentum of shell B 
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
      double precision exinv(nat)   ! [input] inverse nuclear exponents
      double precision scr(lscr)    ! [scratch] scratch buffers
      double precision S(lstv)      ! [output] overlap integrals
      double precision T(lstv)      ! [output] kinetic energy integrals
      double precision V(lstv,ntyp) ! [output] potential integrals
      logical doStil  ! [input] compute modified overlap (True/False)
      logical doTtil  ! [input] compute modified kinetic (True/False)
      logical doVtil  ! [input] compute modified potential (True/False)
      logical canAB   ! [input] compute only canonical ints (false only)
      logical do_nw   ! [input] can do NW integrals
      logical do_hnd  ! [input] can do HONDO integrals
      logical nonrel  ! [input] true if either centre is nonrelativistic
      logical DryRun  ! [input] true means only compute required memory
c:tex-\end{verbatim}
c:tex-See rel_pot for a description of the allowed values of ibug and ntyp
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
      integer n_allp_b  ! n_cartp_b*n_cont_b
      integer n_allp_a  ! n_cartp_a*n_cont_a
      integer n_cartm_a ! cartesian components for l_A-1
      integer n_cartm_b ! cartesian components for l_B-1
      integer n_allm_b  ! n_cartm_b*n_cont_b
      integer n_allm_a  ! n_cartm_a*n_cont_a
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
      character*12 pot_type(4)
      data pot_type
     &    /'      Scalar','z spin-orbit','y spin-orbit','x spin-orbit'/
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
      max_mem = 0
*
      if (.not.(do_nw .or. do_hnd)) call errquit
     &    ('rel_onel: can''t do NW or HONDO integrals',99, INT_ERR)
*
      if (debug_gen) then
        write (LuOut,*)
        write (LuOut,*) 'rel_onel: ibug =',ibug
        write (LuOut,*) 'l_A,n_prim_A,n_cont_A',l_A,n_prim_A,n_cont_A
        write (LuOut,*) 'l_B,n_prim_B,n_cont_B',l_B,n_prim_B,n_cont_B
        if (do_hnd) then
          write (LuOut,*) 'Using HONDO integrals'
        else if (do_nw) then
          write (LuOut,*) 'Using NWChem integrals'
        end if
      end if
*
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = n_cont_a*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_all_b = n_cart_b*n_cont_b
      n_ab = n_cart_ab*n_cont_ab
      if (lstv .lt. n_ab .and. .not.DryRun) call errquit (
     &      'Integral buffer length too small in rel_onel',99,
     &       MEM_ERR)
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
        call ecp_matpr (coefL_A,1,n_prim_a,1,n_cont_a,
     &      1,n_prim_a,1,n_cont_a,'L coef A','E',120,6)
        call ecp_matpr (coefS_A,1,n_prim_a,1,n_cont_a,
     &      1,n_prim_a,1,n_cont_a,'S coef A','E',120,6)
        call ecp_matpr (coefL_B,1,n_prim_b,1,n_cont_b,
     &      1,n_prim_b,1,n_cont_b,'L coef B','E',120,6)
        call ecp_matpr (coefS_B,1,n_prim_b,1,n_cont_b,
     &      1,n_prim_b,1,n_cont_b,'S coef B','E',120,6)
      end if
*
*   Calculate large component overlap and nuclear attraction integrals,
*   and kinetic energy integrals if center is nonrelativistic.
*
      doS = doStil
      doT = doTtil.and.nonrel
      doV = doVtil
      memscr = lscr
      if (do_hnd) then
        call hnd_stvint(
     &      Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &      Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &      Cxyz,zan,nat,S,T,V,n_ab,doS,doT,doV,scr,memscr)
      else if (do_nw) then
        call hf1(
     &      Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &      Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &      Cxyz,zan,exinv,nat,S,T,V,n_ab,doS,doT,doV,canAB,
     &      DryRun,scr,memscr)
      end if
      if (debug_arrays) then
        if (doS) call ecp_matpr(S,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'LL overlap','E',120,6)
        if (doT) call ecp_matpr(T,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'LL kinetic','E',120,6)
        if (doV) call ecp_matpr(V,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'LL potential','E',120,6)
      end if
      if (DryRun) max_mem = max(max_mem,memscr)
      if (nonrel) return
*
*   Calculate kinetic energy integrals, correction to overlaps
*
      if (doTtil) then
        doS = .false.
        doT = .true.
        doV = .false.
        memscr = lscr-n_ab
        if (do_hnd) then
          call hnd_stvint(
     &        Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,nat,S,T,V,n_ab,doS,doT,doV,scr(n_ab+1),memscr)
        else if (do_nw) then
          call hf1(
     &        Axyz,zeta_A,coefL_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,exinv,nat,S,T,V,n_ab,doS,doT,doV,canAB,
     &        DryRun,scr(n_ab+1),memscr)
        end if
        if (DryRun) then
          max_mem = max(max_mem,memscr)
        else if (debug_arrays) then
          call ecp_matpr (T,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'LS kinetic','E',120,6)
        end if
        memscr = lscr-n_ab
        if (do_hnd) then
          call hnd_stvint(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,nat,S,scr,V,n_ab,doS,doT,doV,scr(n_ab+1),memscr)
        else if (do_nw) then
          call hf1(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefL_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,exinv,nat,S,scr,V,n_ab,doS,doT,doV,canAB,
     &        DryRun,scr(n_ab+1),memscr)
        end if
        if (DryRun) then
          max_mem = max(max_mem,memscr+n_ab)
        else
          if (debug_arrays) call ecp_matpr (scr,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'SL kinetic','E',120,6)
          call daxpy (n_ab,one,scr,1,T,1)
          if (debug_arrays) call ecp_matpr (T,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'LS+SL kinetic','E',120,6)
        end if
      end if
      if (doStil .or. doTtil) then
        doS = .false.
        doT = .true.
        doV = .false.
        memscr = lscr-n_ab
        if (do_hnd) then
          call hnd_stvint(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,nat,S,scr,V,n_ab,doS,doT,doV,scr(n_ab+1),memscr)
        else if (do_nw) then
          call hf1(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B,
     &        Cxyz,zan,exinv,nat,S,scr,V,n_ab,doS,doT,doV,canAB,
     &        DryRun,scr(n_ab+1),memscr)
        end if
        if (DryRun) then
          max_mem = max(max_mem,memscr+n_ab)
        else
          if (debug_arrays) call ecp_matpr (scr,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'SS kinetic','E',120,6)
          if (doTtil) then
            call daxpy (n_ab,-one,scr,1,T,1)
            if (debug_arrays) call ecp_matpr (T,1,n_all_b,1,n_all_a,
     &          1,n_all_b,1,n_all_a,'LS+SL-SS kinetic','E',120,6)
          end if
          if (doStil) then
            call daxpy (n_ab,halsq,scr,1,S,1)
            if (debug_arrays) call ecp_matpr (S,1,n_all_b,1,n_all_a,
     &          1,n_all_b,1,n_all_a,'LL+SS overlap','E',120,6)
          end if
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
      n_allp_a = n_cartp_a*n_cont_A
      n_allp_b = n_cartp_b*n_cont_B
      n_allm_a = n_cartm_a*n_cont_A
      n_allm_b = n_cartm_b*n_cont_B
      n_intpp = n_allp_a*n_allp_b
      n_intpm = n_allm_a*n_allp_b
      n_intmp = n_allp_a*n_allm_b
      n_intmm = n_allm_a*n_allm_b
      i_xca = 1
      i_xcb = i_xca+n_prim_A*n_cont_A
      i_pp = i_xcb+n_prim_B*n_cont_B
      i_pm = i_pp+n_intpp
      i_mp = i_pm+n_intpm
      i_mm = i_mp+n_intmp
      i_scr = max(i_xca+n_ab*ntyp,i_mm+n_intmm)
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
     &      'Insufficient scratch memory in rel_onel',99,
     &       MEM_ERR)
        k = i_xca-1
        do j = 1,n_cont_A
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
          k = k+n_prim_B
        end do
      end if
      doS = .false.
      doT = .false.
      doV = .true.
*
*         Calculate integrals for l_A+1, l_B+1
*
      if (do_hnd) then
        call hnd_stvint(
     &      Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,
     &      Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,
     &      Cxyz,zan,nat,S,T,scr(i_pp),n_intpp,doS,doT,doV,
     &      scr(i_scr),memscr)
      else if (do_nw) then
        call hf1(
     &      Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,
     &      Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,
     &      Cxyz,zan,exinv,nat,S,T,scr(i_pp),n_intpp,doS,doT,doV,
     &      canAB,DryRun,scr(i_scr),memscr)
      end if
      if (DryRun) then
        max_mem = max(max_mem,i_scr+memscr-1)
        if (debug_addresses) 
     &      write (luout,*) '++:memscr,max_mem = ',memscr,max_mem
        memscr = lscr-i_scr+1
      else if (debug_arrays) then
        call ecp_matpr (scr(i_pp),1,n_allp_b,1,n_allp_a,
     &      1,n_allp_b,1,n_allp_a,'V(la+1,lb+1)','E',120,6)
      end if
*
*         Calculate integrals for l_A-1, l_B+1
*
      if (l_A .gt. 0) then
        if (do_hnd) then
          call hnd_stvint(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,
     &        Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,
     &        Cxyz,zan,nat,S,T,scr(i_pm),n_intpm,doS,doT,doV,
     &        scr(i_scr),memscr)
        else if (do_nw) then
          call hf1(
     &        Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,
     &        Bxyz,zeta_B,scr(i_xcb),n_prim_B,n_cont_B,l_B+1,
     &        Cxyz,zan,exinv,nat,S,T,scr(i_pm),n_intpm,doS,doT,doV,
     &        canAB,DryRun,scr(i_scr),memscr)
        end if
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          if (debug_addresses) 
     &        write (luout,*) '-+:memscr,max_mem = ',memscr,max_mem
          memscr = lscr-i_scr+1
        else if (debug_arrays) then
          call ecp_matpr (scr(i_pm),1,n_allp_b,1,n_allm_a,
     &        1,n_allp_b,1,n_allp_a,'V(la-1,lb+1)','E',120,6)
        end if
      end if
*
*         Calculate integrals for l_A+1, l_B-1
*
      if (l_B .gt. 0) then
        if (do_hnd) then
          call hnd_stvint(
     &        Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,
     &        Cxyz,zan,nat,S,T,scr(i_mp),n_intmp,doS,doT,doV,
     &        scr(i_scr),memscr)
        else if (do_nw) then
          call hf1(
     &        Axyz,zeta_A,scr(i_xca),n_prim_A,n_cont_A,l_A+1,
     &        Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,
     &        Cxyz,zan,exinv,nat,S,T,scr(i_mp),n_intmp,doS,doT,doV,
     &        canAB,DryRun,scr(i_scr),memscr)
        end if
        if (DryRun) then
          max_mem = max(max_mem,i_scr+memscr-1)
          if (debug_addresses) 
     &        write (luout,*) '+-:memscr,max_mem = ',memscr,max_mem
          memscr = lscr-i_scr+1
        else if (debug_arrays) then
          call ecp_matpr (scr(i_mp),1,n_allm_b,1,n_allp_a,
     &        1,n_allm_b,1,n_allp_a,'V(la+1,lb-1)','E',120,6)
        end if
*     
*         Calculate integrals for l_A-1, l_B-1
*
        if (l_A .gt. 0) then
          if (do_hnd) then
            call hnd_stvint(
     &          Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,
     &          Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,
     &          Cxyz,zan,nat,S,T,scr(i_mm),n_intmm,doS,doT,doV,
     &          scr(i_scr),memscr)
          else if (do_nw) then
            call hf1(
     &          Axyz,zeta_A,coefS_A,n_prim_A,n_cont_A,l_A-1,
     &          Bxyz,zeta_B,coefS_B,n_prim_B,n_cont_B,l_B-1,
     &          Cxyz,zan,exinv,nat,S,T,scr(i_mm),n_intmm,doS,doT,doV,
     &          canAB,DryRun,scr(i_scr),memscr)
          end if
          if (DryRun) then
            max_mem = max(max_mem,i_scr+memscr-1)
            if (debug_addresses) 
     &          write (luout,*) '--:memscr,max_mem = ',memscr,max_mem
            memscr = lscr-i_scr+1
          else if (debug_arrays) then
            call ecp_matpr (scr(i_mm),1,n_allm_b,1,n_allm_a,
     &          1,n_allm_b,1,n_allm_a,'V(la+1,lb-1)','E',120,6)
          end if
        end if
      end if
*
*     Compute the relativistic potential energy integrals
*
      call rel_pot (scr(i_pp),scr(i_pm),scr(i_mp),scr(i_mm),
     &    scr,n_ab,ntyp,
     &    l_A,n_cartp_a,n_cart_a,n_cartm_a,n_cont_A,
     &    l_B,n_cartp_b,n_cart_b,n_cartm_b,n_cont_B,
     &    DryRun,scr(i_scr),memscr,ibug/10)
      if (DryRun) then
        if (debug_addresses) write (luout,*) 'rel_pot:max_mem = ',
     &      max_mem
        max_mem = max(max_mem,i_scr+memscr-1)
        lscr = max_mem
        if (debug_addresses) write (luout,*) 'max_mem = ',max_mem
      else
        i = 1
        do j = 1,ntyp
          if (debug_arrays) then
            write (LuOut,'(//2A)') pot_type(j),' potential'
            call ecp_matpr (scr(i),1,n_all_b,1,n_all_a,1,n_all_b,
     &          1,n_all_a,'SS potential','E',120,6)
          end if
          call daxpy (n_ab,qalsq,scr(i),1,V(1,j),1)
          i = i+n_ab
          if (debug_arrays) call ecp_matpr (V(1,j),1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'Relativistic potential','E',120,6)
        end do
      end if
*
      return
      end
