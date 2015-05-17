C $Id$
************************************************************************
C> \ingroup nwint
C> @{
C>
C> \brief Assemble the relativistic small-small potential energy integrals
C>
C> 
C> This routine assembles the relativistic (small-small) potential 
C> energy integrals from a set of potential energy integrals for 
C> shifted angular momenta. The routine uses the modified Dirac 
C> formalism, and returns integrals in various formats, depending
C> on the input parameter ntyp, which also determines how many 
C> integral blocks are returned. The factors of \f$\alpha^2/4\f$ are 
C> NOT included in the integrals, these must be added by the calling
C> routine.
C> 
C> Author: K. G. Dyall
C>
c:tex-\subsection{rel\_pot}
c:tex-This routine assembles the relativistic (small-small) potential 
c:tex-energy integrals from a set of potential energy integrals for 
c:tex-shifted angular momenta. The routine uses the modified Dirac 
c:tex-formalism, and returns integrals in various formats, depending
c:tex-on the input parameter ntyp, which also determines how many 
c:tex-integral blocks are returned. The factors of $\alpha^2/4$ are 
c:tex-NOT included in the integrals, these must be added by the calling
c:tex-routine.
c:tex-
c:tex-\noindent Author: K. G. Dyall
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine rel_pot (Vpp,Vpm,Vmp,Vmm,Vrel,nints,ntyp,
     &    l_a,np_a,n_a,nm_a,n_cont_A,
     &    l_b,np_b,n_b,nm_b,n_cont_B,
     &    DryRun,scr,lscr,ibug)
c:tex-\end{verbatim}
************************************************************************
      implicit none
#include "stdio.fh"
#include "errquit.fh"
      double precision zero,one
      parameter (zero = 0.0d00, one = 1.0d00)
*
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer nints    !< [Input] number of integrals of each type
      integer ntyp     !< [Input] integral type
      integer l_A      !< [Input] - angular momentum of shell A
      integer np_A     !< [Input] number of cartesian components for \f$l_A+1\f$
      integer n_A      !< [Input] number of cartesian components for \f$l_A\f$
      integer nm_A     !< [Input] number of cartesian components for \f$l_A-1\f$
      integer n_cont_A !< [Input] number of contracted functions on centre A
      integer l_B      !< [Input] angular momentum of functions on centre B
      integer np_B     !< [Input] number of cartesian components for \f$l_B+1\f$
      integer n_B      !< [Input] number of cartesian components for \f$l_B\f$
      integer nm_B     !< [Input] number of cartesian components for \f$l_B-1\f$
      integer n_cont_B !< [Input] number of contracted functions on centre B
      integer lscr     !< [In/Output] - length of scratch array/memory needed
      integer ibug     !< [Input] debug flag. 
      double precision Vpp(np_b,n_cont_b,np_a,n_cont_a) !< [Input] integrals for \f$l_B+1, l_A+1\f$
      double precision Vpm(np_b,n_cont_b,nm_a,n_cont_a) !< [Input] integrals for \f$l_B+1, l_A-1\f$
      double precision Vmp(nm_b,n_cont_b,np_a,n_cont_a) !< [Input] integrals for \f$l_B-1, l_A+1\f$
      double precision Vmm(nm_b,n_cont_b,nm_a,n_cont_a) !< [Input] integrals for \f$l_B-1, l_A-1\f$
      double precision Vrel(nints,ntyp) !< [Output] final relativistic integrals
      double precision scr(lscr) !< [Scratch] scratch array
      logical DryRun   !< [Input] logical for dry run
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
c:tex-Note that n_cont_A and n_cont_B are essentially only the second and
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
      integer j,k ! integer variables
      logical debug_gen       ! do general debug printing
      logical debug_addresses ! do address debug printing
      logical debug_arrays    ! do array debug printing
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) then
        write (LuOut,'(//A,/)') 'Entering rel_pot ...'
        write (LuOut,*) 'ibug = ',ibug
      end if
      n_all_A = n_cont_A*n_A
      n_all_B = n_cont_B*n_B
      n_allp_B = n_cont_B*np_B
      n_allm_B = n_cont_B*nm_B
      np = 3*n_all_A*n_allp_B
      nm = 3*n_all_A*n_allm_B
      nab = n_all_A*n_all_B
      if ((nab .gt. nints).and.(.not.Dryrun)) call errquit
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
        call ecp_grad_a (l_A,np_A,n_A,nm_A,n_cont_A,Vpp,Vpm,
     &      scr(i_grdp),n_allp_B,ibug/10)
        if (debug_arrays)  call ecp_matpr (scr(i_grdp),1,n_allp_b,
     &      1,n_all_a*3,1,n_allp_b,1,n_all_a*3,'Grad A +','E',120,6)
        call dcopy (nm,zero,0,scr(i_grdm),1)
        if (l_B .gt. 0) call ecp_grad_a (l_A,np_A,n_A,nm_A,n_cont_A,
     &      Vmp,Vmm,scr(i_grdm),n_allm_B,ibug/10)
        if (debug_arrays)  call ecp_matpr (scr(i_grdp),1,n_allp_b,
     &      1,n_all_a*3,1,n_allp_b,1,n_all_a*3,'Grad A -','E',120,6)
        call dcopy (nn,zero,0,scr,1)
        call ecp_grad_b (l_B,np_b,n_b,nm_b,n_cont_B,
     &      scr(i_grdp),scr(i_grdm),scr(i_scr),n_all_a*3,ibug/10)
        if (debug_arrays)  call ecp_matpr (scr(i_scr),1,n_all_b,
     &      1,n_all_a*9,1,n_all_b,1,n_all_a*9,'Gradients','E',120,6)
        call dcopy (nints*ntyp,zero,0,Vrel,1)
*
*   Return with derivative integrals only
*
        if (ntyp .eq. 9) then
          do j = 1,3
            do k = 0,6,3
              call dcopy (nab,scr(i_scr),1,Vrel(1,k+j),1)
              i_scr = i_scr+nab
            end do
          end do
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
C>
C> @}
