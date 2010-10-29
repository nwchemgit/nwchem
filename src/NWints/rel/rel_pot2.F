C $Id$
************************************************************************
c:tex-\subsection{rel\_pot2}
c:tex-This routine assembles the relativistic (small-small) potential 
c:tex-energy integrals from a set of potential energy integrals for 
c:tex-shifted angular momenta. The routine uses the modified Dirac 
c:tex-formalism, and returns integrals in various formats, depending
c:tex-on the input parameter ntyp, which also determines how many 
c:tex-integral blocks are returned. This routine generalizes rel\_pot
c:tex-to allow for two-electron integrals of any class. 
c:tex-
c:tex-NOTE: the factor of \alpha^2/4 is NOT included here.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_pot2 (Vpp,Vpm,Vmp,Vmm,Vrel,nints,ntyp,n_x,
     &    l_a,np_a,n_a,nm_a,n_ca,
     &    l_b,np_b,n_b,nm_b,n_cb,
     &    DryRun,scr,lscr,ibug)
c:tex-\end{verbatim}
************************************************************************
      implicit none
#include "stdio.fh"
#include "errquit.fh"
      double precision zero, one
      parameter (zero = 0.0d00, one = 1.0d00)
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer nints    ! [input] number of integrals of each type
      integer ntyp     ! [input] integral type
      integer l_A      ! [input] - angular momentum of shell A
      integer np_A     ! [input] number of cartesian components for l_A+1
      integer n_A      ! [input] number of cartesian components for l_A
      integer nm_A     ! [input] number of cartesian components for l_A-1
      integer n_cA     ! [input] number of contracted functions on centre A
      integer l_B      ! [input] angular momentum of functions on centre B
      integer np_B     ! [input] number of cartesian components for l_B+1
      integer n_B      ! [input] number of cartesian components for l_B
      integer nm_B     ! [input] number of cartesian components for l_B-1
      integer n_cB     ! [input] number of contracted functions on centre B
      integer n_x      ! number of functions preceding B
      integer lscr     ! [i/o] - length of scratch array/memory needed
      integer ibug     ! debug flag. 
      double precision Vpp(n_x,np_b,n_cb,np_a,n_ca) ! [input] integrals for l_B+1, l_A+1
      double precision Vpm(n_x,np_b,n_cb,nm_a,n_ca) ! [input] integrals for l_B+1, l_A-1
      double precision Vmp(n_x,nm_b,n_cb,np_a,n_ca) ! [input] integrals for l_B-1, l_A+1
      double precision Vmm(n_x,nm_b,n_cb,nm_a,n_ca) ! [input] integrals for l_B-1, l_A-1
      double precision Vrel(nints,ntyp) ! [output] final relativistic integrals
      double precision scr(lscr) ! scratch array
      logical DryRun   ! [input] logical for dry run
c:tex-\end{verbatim}
c:tex-The possible values for ntyp are
c:tex-{\leftskip 4pc\obeylines
c:tex-ntyp = 1 : scalar integrals only
c:tex-ntyp = 3 : spin-orbit integrals only
c:tex-ntyp = 4 : scalar and spin-orbit integrals
c:tex-ntyp = 9 : derivative integrals}
c:tex-The possible values for ibug are
c:tex-{\leftskip 4pc\obeylines
c:tex-ibug = 0 for no debug, 
c:tex-ibug = 1 for address printing, 
c:tex-ibug = 2 for array printing, 
c:tex-ibug = 3 for both.}
c:tex-
c:tex-Note that n_ca and n_cb are essentially only the second and
c:tex-fourth dimensions of the integral array; the code does not make use
c:tex-of their significance as number of contracted functions, so there
c:tex-is no reason why these cannot be used for the two-electron integrals.
*
      integer n_all_A  ! total basis functions on A
      integer n_all_B  ! total basis functions on B
      integer n_allp_B ! total basis functions on B for l_B+1
      integer n_allm_B ! total basis functions on B for l_B-1
      integer np  ! no. gradient integrals for l_B+1
      integer nm  ! no. gradient integrals for l_B-1
      integer nab ! size of integral block
      integer nn  ! no. second derivative integrals
      integer i_scr  ! address of scratch space
      integer i_grdp ! address of gradient integrals for l_B+1
      integer i_grdm ! address of gradient integrals for l_B+1
      integer i_free ! address of free space in scr
      integer k ! integer variable
      logical debug_gen       ! do general debug printing
      logical debug_addresses ! do address debug printing
      logical debug_arrays    ! do array debug printing
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering rel_pot ...'
      n_all_A = n_cA*n_A
      n_all_B = n_cB*n_B
      n_allp_B = n_cB*np_B*n_x
      n_allm_B = n_cB*nm_B*n_x
      np = 3*n_all_A*n_allp_B
      nm = 3*n_all_A*n_allm_B
      nab = n_all_A*n_all_B*n_x
      if ((nab .ne. nints) .and. (.not.DryRun)) call errquit
     &    ('Mismatch between nints and nab in rel_pot',99, INT_ERR)
      nn = 9*nab
*
      i_scr = 1
      i_grdp = i_scr+nn
      i_grdm = i_grdp+np
      i_free = i_grdm+nm
*
      if (DryRun) then
        lscr = i_free-1
      else
        if (lscr .lt. i_free-1) call errquit 
     &      ('Insufficient memory in rel_pot',99, MEM_ERR)
        call dcopy (np,zero,0,scr(i_grdp),1)
        call ecp_grad_a (l_A,np_A,n_A,nm_A,n_cA,Vpp,Vpm,
     &      scr(i_grdp),n_allp_B,ibug/10)
        call dcopy (nm,zero,0,scr(i_grdm),1)
        if (l_B .gt. 0) call ecp_grad_a (l_A,np_A,n_A,nm_A,n_ca,
     &      Vmp,Vmm,scr(i_grdm),n_allm_B,ibug/10)
        call dcopy (nn,zero,0,scr,1)
        call ecp_grad_a (l_b,np_b,n_b,nm_b,n_cb*n_all_A*3,
     &      scr(i_grdp),scr(i_grdm),scr,n_x,ibug/10)
        call dcopy (nints*ntyp,zero,0,Vrel,1)
*
*   Return with derivative integrals only
*
        if (ntyp .eq. 9) then
          call dcopy (nab*ntyp,scr(i_scr),1,Vrel,1)
          return
        end if 
*
*   Assemble scalar integrals
*
        if ((ntyp .eq. 1) .or. (ntyp .eq. 4)) then
          call daxpy (nab,one,scr(i_scr),1,Vrel,1)
          call daxpy (nab,one,scr(i_scr+4*nab),1,Vrel,1)
          call daxpy (nab,one,scr(i_scr+8*nab),1,Vrel,1)
          if (debug_arrays) call ecp_matpr (Vrel,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'Scalar integrals','E',78,4)
        end if
*
*   Assemble spin-orbit integrals
*
        k = ntyp-3
        if (k .ge. 0) then
          call daxpy (nab,one,scr(i_scr+5*nab),1,Vrel(1,k+3),1)
          call daxpy (nab,one,scr(i_scr+6*nab),1,Vrel(1,k+2),1)
          call daxpy (nab,one,scr(i_scr+nab),1,Vrel(1,k+1),1)
          call daxpy (nab,-one,scr(i_scr+7*nab),1,Vrel(1,k+3),1)
          call daxpy (nab,-one,scr(i_scr+2*nab),1,Vrel(1,k+2),1)
          call daxpy (nab,-one,scr(i_scr+3*nab),1,Vrel(1,k+1),1)
          if (debug_arrays) call ecp_matpr (Vrel(1,k+1),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'Z spin-orbit integrals','E',78,4)
          if (debug_arrays) call ecp_matpr (Vrel(1,k+2),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'Y spin-orbit integrals','E',78,4)
          if (debug_arrays) call ecp_matpr (Vrel(1,k+3),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'X spin-orbit integrals','E',78,4)
        end if
      end if
*
      if (debug_gen) write (LuOut,*) 'Exiting rel_pot'
*
      return
      end
