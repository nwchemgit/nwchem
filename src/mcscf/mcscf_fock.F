c     
c     %A%
c     $Id$
c     
c     Compute the MCSCF MO-Fock matrices;
c     inactive, active and generalized
c     
c     
      subroutine mcscf_fock( geom, basis, nbf, nclosed, nact,
     $     oskel, tol2e, dm1, dm2, g_movecs,
     $     g_coul, eone, etwo, ecore,
     $     g_ifock, g_afock, g_gfock )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "mcscfprof.fh"
c     
c     
c     The definition for the *generalized* Fock
c     elements is defined as twice the normal Fock elements 
c     (see Siegbahn et al. [1]) and therefore the gradient 
c     is twice the generalized Fock element. This routine
c     follows this convention, i.e.,
c     
c     gen    gen
c     g   = 2(F    - F   )  = 4F
c     pq      pq     qp        pq
c     
c     
c     Note the asymmetry in the generalized Fock,
c     
c     
c     F   =  0      for p not closed-shell
c     pi
c     
c     F   =  0      for p not open-shell
c     pt
c     
c     [1] P.E.M. Siegbahn, J. Almlof, A. Heiberg, B.O. Roos,
c         "The complete active space SCF (CASSCF) method in a 
c          Newton-Raphson formulation with application to the HNO
c          molecule", J. Chem. Phys. 74 (1981),
c         pp 2384-2396, DOI:10.1063/1.441359. 
c         (Note this reference was not given originally in the comments
c          above. However, it seems to match the subject matter best of
c          whatever I could find. HvD)
c     
      integer geom, basis       ! [input] Geometry and basis handles
      integer nbf               ! [input] Number of basis functions
      integer nclosed           ! [input] Number of closed shells
      integer nact              ! [input] Number of open shells
      logical oskel             ! [input] Symmetry toggle
      double precision tol2e    ! [input] Integral tolerance
      double precision dm1(nact,nact) ! [input] Active space 1PDM
      double precision dm2(nact,nact,nact,nact) ! [input] Active space 2PDM
      integer g_movecs          ! [input] MO coefficients
      integer g_coul            ! [input] Coulomb integrals
      double precision eone, etwo ! [output] Energy components
      double precision ecore    ! [output] Inactive core energy
      integer g_ifock           ! [output] Inactive Fock
      integer g_afock           ! [output] Active Fock
      integer g_gfock           ! [output] Generalized Fock
c     
      integer g_cdens, g_adens, g_tmp, g_tmp2, g_tmp3
      double precision e2ii, e2ia, e2ai, e2aa, xx
      integer nset, aoff, aend, i
      integer iv_dens(2), iv_fock(2)
      double precision jfac(2), kfac(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision ga_trace_diag
      external ga_trace_diag
      data nset/2/
      data jfac/1.0d0, 1.d0/
      data kfac/-0.5d0, -0.5d0/
c     
c     
c     
      if (omcscfprof) call pstat_on(ps_gfock)
      aoff = nclosed + 1
      aend = nclosed + nact
      call ga_zero(g_ifock)
      call ga_zero(g_afock)
      g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'temp2')
      g_tmp3 = ga_create_atom_blocked(geom, basis, 'temp3')
      g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
      g_adens = ga_create_atom_blocked(geom, basis, 'active dens')
      call mcscf_occ2dens(nbf, nclosed, nact, dm1, g_movecs,
     $     g_cdens, g_adens )
c     
c     One-electron component
c     
      call ga_zero(g_tmp)
      call ga_zero(g_tmp2)
      call ga_zero(g_tmp3)
      call int_1e_ga( basis, basis, g_tmp, 'kinetic', oskel)
      call int_1e_ga( basis, basis, g_tmp, 'potential', oskel)
      if (oskel) call sym_symmetrize(geom, basis, .false., g_tmp)
      eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
c     
c     Two-electron component of the Fock matrices
c     
      iv_dens(1) = g_cdens
      iv_dens(2) = g_adens
      iv_fock(1) = g_tmp2
      iv_fock(2) = g_tmp3
      call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $     iv_dens, iv_fock, .false. )
c     
c     
c     Symmetrize Fock AO components
c     
      if (oskel) then
         call sym_symmetrize(geom, basis, .false., g_tmp2)
         call sym_symmetrize(geom, basis, .false., g_tmp3)
      endif
      e2ii = ga_ddot(g_tmp2,g_cdens)
      e2ia = ga_ddot(g_tmp2,g_adens)
      e2ai = ga_ddot(g_tmp3,g_cdens)
      call ga_dadd(1.d0, g_tmp, 1.d0, g_tmp2, g_tmp2 )
c     
c     Transform Fock to MO representation
c     
      call two_index_transf(g_tmp2, g_movecs, g_movecs,
     $     g_cdens, g_ifock )
      call two_index_transf(g_tmp3, g_movecs, g_movecs,
     $     g_cdens, g_afock )
c     
c     Generalized Fock matrix elements
c     I      A
c     F   =  2 (F   +  F  )
c     iq        iq     iq
c     
c     
      call ga_zero(g_gfock)
      call ga_dadd_patch( 2.d0, g_ifock, 1, nclosed, 1, nbf,
     $     2.d0, g_afock, 1, nclosed, 1, nbf,
     $     g_gfock, 1, nclosed, 1, nbf )
c     
c     I 
c     F   =  sum[ D  .F   ]  +  2.sum[ P    (qu|vx) ]
c     tq          tu  qu               tuvx
c     
c     
      call ga_zero(g_tmp)
      call ga_zero(g_cdens)
      call mcscf_denfock_trace(nclosed, nact, nbf, dm1,
     $     g_ifock, g_tmp)
      call mcscf_3index_trace( nbf, nclosed, nact, dm2,
     $     g_coul, g_cdens )
      e2aa = ga_trace_diag(g_cdens) ! <--- NB. Fii and Faa are zero in this
      call ga_dadd_patch( 1.d0, g_tmp, aoff, aend, 1, nbf, ! case, only Ftt is non-zero
     $     1.d0, g_cdens, aoff, aend, 1, nbf,
     $     g_gfock, aoff, aend, 1, nbf )
c     
c     Inactive core energy (when added to full CI energy = total energy)
c     
c     I       I
c     ecore = sum  F   -  e2    
c     i  ii          
c     
c     I
c     = sum  F   - sum  [ 2(ii|jj) - (ij|ij) ]
c     i  ii      ij
c     
      ecore = 0.d0
      do i=ga_nodeid()+1,nclosed,ga_nnodes()
         call ga_get(g_ifock,i,i,i,i,xx,1)
         ecore = ecore + xx*2.d0
      enddo
      call ga_sync()
      call ga_dgop(msg_mcscf_ifocktrace2,ecore,1,'+')
      ecore = ecore - e2ii*0.5d0
c     
c     Sum 2e-energy components
c     
      
      etwo = (e2ii + e2ia + e2ai + e2aa)*0.5d0

c$$$  WRITE(6,881) EONE,(E2II*0.5D0),E2IA,(E2AA*0.5d0)
c$$$  881   FORMAT(/,5X,'ECORE:',F20.14,
c$$$  $        /,5X,'EONE:',1X,F20.14,
c$$$  $        /,5X,'EII:',2X,F20.14,
c$$$  $        /,5X,'EIA:',2X,F20.14,
c$$$  $        /,5X,'EAA:',2X,F20.14)
c     
c     Clean up
c     
      if (.not.ga_destroy(g_tmp))
     $     call errquit('mcscf_fock: cannot destroy',0, GA_ERR)
      if (.not.ga_destroy(g_tmp2))
     $     call errquit('mcscf_fock: cannot destroy',0, GA_ERR)
      if (.not.ga_destroy(g_tmp3))
     $     call errquit('mcscf_fock: cannot destroy',0, GA_ERR)
      if (.not.ga_destroy(g_cdens))
     $     call errquit('mcscf_fock: cannot destroy',0, GA_ERR)
      if (.not.ga_destroy(g_adens))
     $     call errquit('mcscf_fock: cannot destroy',0, GA_ERR)
      if (omcscfprof) call pstat_off(ps_gfock)
      return
      end






       subroutine mcscf_denfock_trace(nclosed,nact,nbf,dm1,g_if,g_a)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       integer nclosed, nact, nbf
       double precision dm1(nact,nact)
       integer g_if, g_a
c
       integer l_ff, k_ff
       integer clo,chi,msize
       integer k_if,ld
       logical check_moblocked 
       external check_moblocked
c
       if (.not.check_moblocked(g_if,nbf,clo,chi))
     $      call errquit('mcscf_denfock_trace: wrong distrib.',0,
     &       INPUT_ERR)
       if ((clo.gt.0).and.(chi.gt.0)) then
         msize = (chi-clo+1)*nact
         if (.not.ma_push_get(MT_DBL,msize,'ff',l_ff,k_ff))
     $     call errquit('mcscf_denfock_trace: cannot allocate',0,
     &       MA_ERR)
         call dfill(msize,0.d0,dbl_mb(k_ff),1)
         call ga_access(g_if,1,nbf,clo,chi,k_if,ld)
         call mcscf_trace_xyz( nbf,nclosed,nact,clo,chi,dbl_mb(k_if),
     $                         dm1,dbl_mb(k_ff))
         call ga_release(g_if,1,nbf,clo,chi)
         call ga_acc(g_a,nclosed+1,nclosed+nact,clo,chi,dbl_mb(k_ff),
     $               nact,1.d0)
         if (.not.ma_pop_stack(l_ff))
     $     call errquit('mcscf_denfock_trace: cannot pop stack',0,
     &       MA_ERR)
       endif
       call ga_sync()
       return
       end









       subroutine mcscf_trace_xyz(nbf,nclosed,nact,clo,chi,inf,
     $                            dm1,ff)
       implicit none
       integer nbf,nclosed,nact
       integer clo,chi
       double precision inf(nbf,clo:chi)
       double precision dm1(nact,nact)
       double precision ff(nact,clo:chi)
       integer clen,offset
       
       clen = chi - clo + 1
       offset = nclosed + 1
       call dgemm('n','n',nact,clen,nact,1.d0,dm1,nact,inf(offset,clo),
     $            nbf,0.d0,ff,nact)
       return
       end







c
c  Take this trace
c
c      F    =  sum     Gamma    (qu|vx)
c       tq        uvx       tuvx
c
c                                vx
c           =  sum  sum  Gamma  J        (1st summation is local,
c                 vx   u      tu qu       2nd, 3rd via ga_acc)
c
c
       subroutine mcscf_3index_trace( nbf, nclosed, nact, dm2,
     $                                g_coul, g_a )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       integer nbf,nclosed,nact
       double precision dm2(nact,nact,nact,nact)
       integer g_coul
       integer g_a
c
       integer msize, nn, clo, chi, v, x, vx
       integer k_ff, l_ff, ld, k_c
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
c
       msize = nbf*nact
       nn = nbf*nbf
       if (.not.ma_push_get( MT_DBL, msize, 'ff', l_ff, k_ff ))
     $   call errquit('mcscf_3index_trace: cannot allocate',0, MA_ERR)
       call dfill(msize,0.d0,dbl_mb(k_ff),1)
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_3index_trace: wrong distrib operator',0,
     &       GA_ERR)
c
c   Loop over (vx), do summation if local
c
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access( g_coul, 1, nn, vx, vx, k_c, ld )
             call mcscf_trace_pqr( nbf, nclosed, nact, v, x,
     $                             dbl_mb(k_c), dm2, dbl_mb(k_ff) )
             call ga_release( g_coul, 1, nn, vx, vx )
           endif
         enddo
       enddo
       call ga_acc( g_a, nclosed+1, nclosed+nact, 1, nbf,
     $              dbl_mb(k_ff), nact, 1.d0 )
       if (.not.ma_pop_stack(l_ff))
     $   call errquit('mcscf_3index_trace: cannot pop stack',0, MA_ERR)
       call ga_sync()
       return
       end







       subroutine mcscf_trace_pqr( nbf, nclosed, nact, v, x, c,
     $                             dm2, ff )
       implicit none
       integer nbf, nclosed, nact, v, x
       double precision c(nbf,nbf)
       double precision dm2(nact,nact,nact,nact)
       double precision ff(nact,nbf)
       integer offset
       double precision scale

       offset = nclosed + 1
       scale = 2.d0
       if (v.eq.x) scale = 1.d0
       call dgemm( 'n', 'n', nact, nbf, nact, scale, dm2(1,1,x,v),
     $             nact, c(offset,1), nbf, 1.d0, ff, nact )
       call dgemm( 'n', 't', nact, nbf, nact, scale, dm2(1,1,v,x),
     $             nact, c(1,offset), nbf, 1.d0, ff, nact )
       return
       end






c
c
c   Return  (sum over active diagonals of above quantity)
c   This is the active-active energy component
c
c      sum    Gamma    (tu|vx)
c         tuvx     tuvx
c
c

       double precision function mcscf_act_ener( nbf, nclosed, nact,
     $                                           dm2, g_coul )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer nbf,nclosed,nact
       double precision dm2(nact,nact,nact,nact)
       integer g_coul
c
       integer nn, clo, chi, v, x, vx, ld, k_c
       double precision xx, s
       logical ga_check_JKblocked
       double precision mcscf_trace_pq
       external ga_check_JKblocked
       external mcscf_trace_pq
c
c
c
       nn = nbf*nbf
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_3index_trace: wrong distrib operator',0,
     &       GA_ERR)
c
c   Loop over (vx), do summation if local
c
       xx = 0.d0
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           s = 1.d0
           if (v.eq.x) s = 0.5d0
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access( g_coul, 1, nn, vx, vx, k_c, ld )
             xx = xx + s*mcscf_trace_pq( nbf, nclosed, nact,
     $                                   dbl_mb(k_c), dm2(1,1,v,x))
             call ga_release( g_coul, 1, nn, vx, vx )
           endif
         enddo
       enddo
       call ga_sync()
       call ga_dgop(msg_mcscf_actenergy,xx,1,'+')
       mcscf_act_ener = xx
       return
       end












      double precision function mcscf_trace_pq( nbf, nclo, nact, c, d )
      implicit none
      integer nbf, nclo, nact
      double precision c(nbf,nbf), d(nact,nact)
      double precision xx
      integer i,j

      xx = 0.d0
      do i=1,nact
        do j=1,nact
          xx = xx + c(nclo+j,nclo+i)*(d(j,i)+d(i,j))
        enddo
      enddo
      mcscf_trace_pq = xx
      return
      end


      subroutine mcscf_dump_for_grad(rtdb, geom, basis, nbf, nclosed,
     $     nact,oskel, tol2e, dm1, dm2, g_movecs,  g_coul)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "eaf.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
c
c     Dump out the information necessary for computing the
c     gradient w.r.t. nuclear coordinates.
c
c     (1) Total 1-pdm in AO basis
c     (2) Effective active space 2-pdm symmetrized
c     (3) Negative energy weighted density in the AO basis
c
c     NOTE: the dm2 is destroyed in this process
c
c     Also while we have the AO 1-pdm compute the multipoles
c
      integer rtdb
      integer geom, basis       ! [input] Geometry and basis handles
      integer nbf               ! [input] Number of basis functions
      integer nclosed           ! [input] Number of closed shells
      integer nact              ! [input] Number of open shells
      logical oskel             ! [input] Symmetry toggle
      double precision tol2e    ! [input] Integral tolerance
      double precision dm1(nact,nact) ! [input] Active space 1PDM
      double precision dm2(nact,nact,nact,nact) ! [input] Active space 2PDM
      integer g_movecs          ! [input] MO coefficients
      integer g_coul            ! [input] Coulomb integrals
c
      double precision offset   ! EAF file offset
      integer pdmfile           ! EAF handle for 2pdm file
      character*255 pdmfilename ! EAF 2pdm filename
c
      double precision eone, etwo ! Energy components
      double precision ecore    ! Inactive core energy
      integer g_ifock           ! Inactive Fock
      integer g_afock           ! Active Fock
      integer g_gfock           ! Generalized Fock
      integer g_over            ! Overlap matrix
c
      integer l_buf, k_buf, p, q, r, s, nocc
      logical odebug
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      nocc = nclosed+nact
      odebug = util_print('graddebug',print_debug)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'Act Fock',nbf,0,g_afock))
     $     call errquit('mcscf_orb: cannot allocate active Fock',0,
     &       GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'In Fock',nbf,0,g_ifock))
     $     call errquit('mcscf_orb: cannot allocate inactive Fock',0,
     &       GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'Gen Fock',nbf,0,g_gfock))
     $     call errquit('mcscf_orb: cannot allocate general Fock',0,
     &       GA_ERR)
c      
      call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $     oskel, tol2e, dm1, dm2, g_movecs,
     $     g_coul, eone, etwo, ecore,
     $     g_ifock, g_afock, g_gfock )
c
      if (.not. ma_push_get(mt_dbl,nbf*nbf,'buffer',l_buf,k_buf))
     $     call errquit('mcscf: grad dump ... memory problem',nbf**2,
     &       MA_ERR)
c
      if (odebug) then
         if (ga_nodeid().eq.0) then
            write(6,*) ' THE ACTIVE SPACE 1-PDM AT CONVERGENCE '
            call output(dm1, 1, nact, 1, nact, nact, nact, 1)
            write(6,*) ' THE LAGRANGIAN AT CONVERGENCE '
            call util_flush(6)
         endif
         call ga_print(g_gfock)
      endif
c
c     Open the file to hold the info
c
      if (ga_nodeid().eq.0) then
         call util_file_name('mcpdm',.true.,.false.,pdmfilename)
         if (eaf_open(pdmfilename, eaf_rw, pdmfile) .ne. 0)
     $        call errquit('mcscf: failed opening 2pdm file',0,
     &       DISK_ERR)
         offset = 0.0d0
      endif
c
c     Transform core+active Lagrangian to AO basis using afock
c     as workspace, putting AO result back into gfock.
c     Then write it out.
c
      call ga_dgemm('n','t',nocc,nbf,nocc,1.0d0,g_gfock,g_movecs,
     $     0.0d0,g_afock)
      call ga_dgemm('n','n',nbf,nbf,nocc,1.0d0,g_movecs,g_afock,
     $     0.0d0,g_gfock)
      if (ga_nodeid() .eq. 0) then
         call ga_get(g_gfock,1,nbf,1,nbf,dbl_mb(k_buf),nbf)
         if (eaf_write(pdmfile, offset, dbl_mb(k_buf), 8*nbf**2).ne.0)
     $        call errquit('mcscf: failed writing lagrangian',0, GA_ERR)
         offset = offset + 8*nbf**2
      endif
c
c     Form the full 1-particle density in the AO basis and output.
c     Do this by making full MO density in gfock and then transforming
c
      call ga_zero(g_gfock)
      if (ga_nodeid() .eq. 0) then
         call ga_put(g_gfock,
     $        nclosed+1,nclosed+nact,nclosed+1,nclosed+nact, dm1, nact)
         do p = 1, nclosed
            call ga_put(g_gfock,p,p,p,p,2.0d0,1)
         enddo
      endif
      if (odebug) then
         if (ga_nodeid().eq.0) then
            write(6,*) ' Full MO density '
            call util_flush(6)
         endif
         call ga_print(g_gfock)
      endif
      call ga_dgemm('n','t',nocc,nbf,nocc,1.0d0,g_gfock,g_movecs,
     $     0.0d0,g_afock)
      call ga_dgemm('n','n',nbf,nbf,nocc,1.0d0,g_movecs,g_afock,
     $     0.0d0,g_gfock)
c
      if (ga_nodeid() .eq. 0) then
         call ga_get(g_gfock,1,nbf,1,nbf,dbl_mb(k_buf),nbf)
         if (eaf_write(pdmfile, offset, dbl_mb(k_buf), 8*nbf**2).ne.0)
     $        call errquit('mcscf: failed writing lagrangian',0,
     &       DISK_ERR)
         offset = offset + 8*nbf**2
      endif
c
c     Peffective(p,q,r,s) = P(pq,rs) - (D(p,q)D(r,s) - 1/2 D(p,r)D(q,s))
c     symmetrized over r and s (or p and q)
c
c     Note that Adrian's 2-pdm includes the factor of 1/2 not
c     included in RJH's gradient derivation ... hence the factor
c     of 2 below
c
      do s = 1, nact
         do r = 1, nact
            do q = 1, nact
               do p = 1, nact
                  dm2(p,q,r,s) = 2.0d0*dm2(p,q,r,s) - 
     $                 (dm1(p,q)*dm1(r,s) - 0.5d0*dm1(p,r)*dm1(q,s))
               enddo
            enddo
            do p = 1, nact
               do q = 1, p-1
                  dm2(p,q,r,s) = 0.5d0*(dm2(p,q,r,s)+dm2(q,p,r,s))
                  dm2(q,p,r,s) = dm2(p,q,r,s)
               enddo
            enddo
         enddo
      enddo
c
      if (odebug .and. ga_nodeid().eq.0) then
         do s = 1, nact
            do r = 1, nact
               write(6,*) ' GAMMA *,*,',r,s
               call output(dm2(1,1,r,s),1,nact,1,nact,nact,nact,1)
            enddo
         enddo
         call util_flush(6)
      endif
c
      if (ga_nodeid() .eq. 0) then
         if (eaf_write(pdmfile, offset, dm2, 8*nact**4).ne.0)
     $        call errquit('mcscf: failed writing 2pdm file',0,
     &       DISK_ERR)
         if (eaf_close(pdmfile).ne.0)
     $        call errquit('mcscf: failed closing 2pdm file',0,
     &       DISK_ERR)
      endif
c
      if (ga_nodeid().eq.0 .and. 
     $     util_print('information',print_default)) then
         write(6,777) 
     $        pdmfilename(1:inp_strlen(pdmfilename))
 777     format(/' Wrote Lagrangian and 1- and 2-particle',
     $        ' density matrices to ',a/)
         call util_flush(6)
      endif
c
      call ga_sync()
c
c     Analyze the total density
c
      call mcscf_mpole(rtdb, geom, basis, g_gfock)
c
      g_over  = ga_create_atom_blocked(geom, basis, 'cas_mull:over')
c
      if (util_print('mulliken', print_low)) then
         call ga_zero(g_over)
         call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,
     $           'Mulliken analysis of the total density', 20,.true.)
         endif
         call mull_pop(geom, basis, g_gfock, g_over, 'total')
c     
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,
     $           'Mulliken analysis of the active space density',
     $           25,.true.)
         endif
c     
         call ga_zero(g_gfock)
         if (ga_nodeid() .eq. 0) then
            call ga_put(g_gfock,
     $           nclosed+1,nclosed+nact,nclosed+1,nclosed+nact,dm1,nact)
         endif
         call ga_dgemm('n','t',nocc,nbf,nocc,1.0d0,g_gfock,g_movecs,
     $        0.0d0,g_afock)
         call ga_dgemm('n','n',nbf,nbf,nocc,1.0d0,g_movecs,g_afock,
     $        0.0d0,g_gfock)
c     
         call mull_pop(geom, basis, g_gfock, g_over, 'open shell')
      endif
c
      if (.not. ma_pop_stack(l_buf)) call errquit('mcscf: dump MA?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_gfock)) call errquit('mcscf: dump GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_ifock)) call errquit('mcscf: dump GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_afock)) call errquit('mcscf: dump GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_over))  call errquit('mcscf: dump GA?',0,
     &       GA_ERR)
c
      end
      subroutine mcscf_mpole(rtdb, geom, basis, g_dens)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "rtdb.fh"
#include "bas.fh"
      integer rtdb, geom, basis, g_dens
c
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments) 
      double precision nuclmoments (len_moments)
      double precision center(3)
      logical oprintmultipole
      integer k, l, m, ltotal, ind
      data center/3*0.0d0/
c      
      if (.not. bas_cando_mpoles(basis)) return
c
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(basis, center, lmax, g_dens, 
     $     totalmoments, len_moments)
      do k = 1, len_moments
         totalmoments(k) = totalmoments(k) + nuclmoments(k)
      enddo
c     
      oprintmultipole = util_print('multipole', print_low)
      if (oprintmultipole .and. ga_nodeid().eq.0) then
         write(6,*)
         call util_print_centered(6,
     $        'Multipole analysis of the density', 20, .true.)
         write(6,*)
         write(6,*) '    L   x y z        total',
     $                 '        nuclear'
         write(6,*) '    -   - - -        -----',
     $        '        -------'
         ind = 0
         do ltotal = 0, lmax
            do k = ltotal, 0, -1
               do l = ltotal-k, 0, -1
                  m = ltotal - k - l
                  ind = ind + 1
                  write(6,12) ltotal, k, l, m, totalmoments(ind),
     $                 nuclmoments(ind)
 12               format(4x,i2,2x,3i2,2f14.6)
               enddo
            enddo
            write(6,*)
         enddo
      endif
      call ecce_print1('total dipole', mt_dbl, totalmoments(2), 3)
      call ecce_print1('nuclear dipole', mt_dbl, nuclmoments(2), 3)
      call ecce_print1('total quadrupole', mt_dbl, totalmoments(5), 6)
      call ecce_print1('nuclear quadrupole', mt_dbl, nuclmoments(5), 6)
c
c     Store the dipole moment
c
      if (.not. rtdb_put(rtdb, 'mcscf:dipole', mt_dbl, 
     $     3, totalmoments(2))) call errquit('multipole: rtdb?',0,
     &       RTDB_ERR)
c
      end
