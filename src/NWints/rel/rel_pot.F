C $Id: rel_pot.F,v 1.1 1998-09-22 23:16:47 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine rel_pot (Vpp,Vpm,Vmp,Vmm,Vrel,ntyp,
     &    l_a,np_a,n_a,nm_a,n_cont_A,
     &    l_b,np_b,n_b,nm_b,n_cont_B,
     &    DryRun,scr,lscr,ibug)
*                                                                      *
*   Given the integrals for shifted angular momenta, assemble the      *
*   small component potential energy integrals in the modified Dirac   *
*   formalism. The output consists of the spin-free and spin-orbit     *
*   integrals or the 3x3 array of gradient integrals.                  *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   Vpp (inp) - integrals for l_B+1, l_A+1                             *
*   Vpm (inp) - integrals for l_B+1, l_A-1                             *
*   Vmp (inp) - integrals for l_B-1, l_A+1                             *
*   Vmm (inp) - integrals for l_B-1, l_A-1                             *
*   Vrel (out) - final integral array                                  *
*   ntyp (inp) - gives type of integral output, also acts as dimension *
*                ntyp = 1 : scalar integrals only                      *
*                ntyp = 3 : spin-orbit integrals only                  *
*                ntyp = 4 : scalar and spin-orbit integrals            *
*                ntyp = 9 : derivative integrals                       *
*   l_A (inp) - angular momentum of functions on centre A              *
*   np_A (inp) - number of cartesians for l_A+1                        *
*   n_A (inp) - number of cartesians for l_A                           *
*   nm_A (inp) - number of cartesians for l_A-1                        *
*   n_cont_A (inp) - number of contracted functions on centre A        *
*   l_B (inp) - angular momentum of functions on centre B              *
*   np_B (inp) - number of cartesians for l_B+1                        *
*   n_B (inp) - number of cartesians for l_B                           *
*   nm_B (inp) - number of cartesians for l_B-1                        *
*   n_cont_B (inp) - number of contracted functions on centre B        *
*   DryRun (inp) - logical for dry run. If true, routine only returns  *
*                  maximum scratch space needed, if false, integrals   *
*                  are returned.                                       *
*   scr (scr) - scratch array for work space                           *
*   lscr (i/o) - length of scratch array. Value returned if DryRun is  *
*                true, used as dimension if false.                     *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
      integer l_a,np_a,n_a,nm_a,l_b,np_b,n_b,nm_b,
     &    n_all_A,n_all_B,n_allp_B,n_allm_B,np,nm,nab,nn,
     &    i_scr,i_grdp,i_grdm,i_free
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision 
     &    Vpp(np_b,n_cont_b,np_a,n_cont_a),
     &    Vpm(np_b,n_cont_b,nm_a,n_cont_a),
     &    Vmp(nm_b,n_cont_b,np_a,n_cont_a),
     &    Vmm(nm_b,n_cont_b,nm_a,n_cont_a),
     &    Vrel(n_b,n_cont_b,n_a,n_cont_a,ntyp)
     &    scr(lscr)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering rel_pot ...'
      n_all_A = n_cont_A*n_A
      n_all_B = n_cont_B*n_B
      n_allp_B = n_cont_B*np_B
      n_allm_B = n_cont_B*nm_B
      np = 3*n_all_A*n_allp_B
      nab = n_all_A*n_all_B
      nn = 9*nab
      nm = 3*n_all_A*n_allm_B
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
     &      ('Insufficient memory in rel_pot',99)
        call dcopy (np,zero,0,scr(i_grdp),1)
        call ecp_grad_a (l_A,np_A,n_A,nm_A,n_cont_A,Vpp,Vpm,
     &      scr(i_grdp),n_allp_B)
        call dcopy (nm,zero,0,scr(i_grdm),1)
        if (l_B .gt. 0) call ecp_grad_a (l_A,np_A,n_A,nm_A,n_cont_A,
     &      Vmp,Vmm,scr,n_allm_B)
        call dcopy (nn,zero,0,scr(i_scr),1)
        call ecp_grad_b (l_B,np_b,n_b,nm_b,n_cont_B,
     &      scr(i_grdp),scr(i_grdm),scr(i_scr),n_all_a*3)
        call dcopy (nab*ntyp,zero,0,Vrel,1)
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
          call dcopy (nab,scr(i_scr),1,Vrel,1)
          call daxpy (nab,one,scr(i_scr+4*nab),1,Vrel,1)
          call daxpy (nab,one,scr(i_scr+8*nab),1,Vrel,1)
          if (debug_arrays) call ecp_matpr (Vrel,1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'Scalar integrals,'E',78,4)
        end if
*
*   Assemble spin-orbit integrals
*
        k = ntyp-3
        if (k .ge. 0) then
          call dcopy (nab,scr(i_scr+5*nab),1,Vrel(1,1,1,1,k+3),1)
          call dcopy (nab,scr(i_scr+6*nab),1,Vrel(1,1,1,1,k+2),1)
          call dcopy (nab,scr(i_scr+nab),1,Vrel(1,1,1,1,k+2),1)
          call daxpy (nab,-one,scr(i_scr+7*nab),1,Vrel(1,1,1,1,k+3),1)
          call daxpy (nab,-one,scr(i_scr+2*nab),1,Vrel(1,1,1,1,k+2),1)
          call daxpy (nab,-one,scr(i_scr+3*nab),1,Vrel(1,1,1,1,k+1),1)
          if (debug_arrays) call ecp_matpr (Vrel(1,1,1,1,k+1),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'Z spin-orbit integrals,'E',78,4)
          if (debug_arrays) call ecp_matpr (Vrel(1,1,1,1,k+2),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'Y spin-orbit integrals,'E',78,4)
          if (debug_arrays) call ecp_matpr (Vrel(1,1,1,1,k+3),
     &        1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &        'X spin-orbit integrals,'E',78,4)
        end if
      end if
*
      if (debug_gen) write (LuOut,*) 'Exiting rel_pot'
*
      return
      end
