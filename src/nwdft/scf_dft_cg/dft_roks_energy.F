      subroutine dft_roks_energy( g_vecs, eone, etwo, exc, enrep,
     $                            energy, edisp, g_grad, rtdb, nexc )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cosmo.fh"
c     
c     $Id$
c
      integer rtdb

      integer g_vecs
      double precision energy
      double precision edisp
      integer g_grad
c     
      integer nexc
      double precision eone, etwo, enrep, exc(2)
      integer gtype, grow, gcol
c
      logical odebug
      logical osome
      logical xc_gotxc
      external xc_gotxc
      double precision ecosmo,gcds
      logical do_cosmo_smd
c
      odebug=.false.
      osome=.false.
      osome=osome.or.odebug
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $    call errquit('dft_roks_energy: ROKS internal block invalid',0,
     &                 UNKNOWN_ERR)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.crohf_vlen).or.(gcol.ne.1))
     $     call errquit('dft_roks_energy: invalid vector length',grow,
     &                  UNKNOWN_ERR)
c     
      if(ga_nodeid().eq.0.and.osome) then
         write(6,*)'in dft_roks_energy calling dft_roks_fock ...'
         call util_flush(6)
      endif

      call dft_roks_fock(rtdb, geom, basis,
     $     nclosed, nopen, nmo, tol2e, g_vecs,
     $     eone, etwo, edisp, Exc, 
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp, nExc, oskel, iter )

      if(ga_nodeid().eq.0.and.osome) then
         write(6,*)'in dft_roks_energy returning from dft_roks_fock ...'
         call util_flush(6)
      endif
c
      if(ga_nodeid().eq.0.and.osome) then
         write(6,*) 'in  dft_roks_energy calling rohf_get_grad ...'
         call util_flush(6)
      endif

      call rohf_get_grad(g_grad)

      if(ga_nodeid().eq.0.and.osome) then
         write(6,*)'in dft_roks_energy returning from rohf_get_grad ...'
         call util_flush(6)
      endif
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('dft_roks_energy: no repulsion energy?', 0,
     &       UNKNOWN_ERR)
      energy = eone + etwo + enrep
      if (xc_gotxc()) then
        energy = energy + exc(1) + exc(2)
      endif
c
      if (ga_nodeid().eq.0 .and.
     $    util_print('roks energy',print_debug)) then
       write(6,33) energy, eone, etwo, enrep
33     format(' total energy = ',f20.12/
     $        ' 1-e energy   = ',f20.12/
     $        ' 2-e energy   = ',f20.12/
     $        ' nucl. energy = ',f20.12/)
       call util_flush(6)
      endif 
c
c     ----- cosmo and smd energy term -----
c
      odebug=.false.
      odebug=odebug.and.ga_nodeid().eq.0
c
      do_cosmo_smd = .false.
      gcds = 0.d0
      if(cosmo_on.and.cosmo_phase.eq.2) then
         if (.not. rtdb_get(
     $     rtdb,'cosmo:energy',mt_dbl,1,ecosmo))
     $     call errquit(
     $     'dft_roks_energy: rtdb get failed for ecosmo',
     $                  911, RTDB_ERR)
c
         if (.not. rtdb_get(
     $    rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $    call errquit(
     $    'dft_roks_energy: cannot get do_cosmo_smd from rtdb',
     $    0,rtdb_err)
c
          energy=energy+ecosmo  ! add cosmo energy
c
          if (do_cosmo_smd) then
            if (.not. rtdb_get
     $        (rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $      call errquit(
     $      'dft_fockbld: cannot get gcds from rtdb',
     $       0,rtdb_err)
c            
             energy=energy+gcds  ! add smd energy
c
          endif   ! do_cosmo_smd
c
         osome=.false.
         osome=osome.and.ga_nodeid().eq.0
         if(osome) then
             write(6,9999) (energy-ecosmo-gcds),ecosmo,energy
 9999 format(' energ0 = ',f20.10,/,' ecosmo = ',f20.10,/,
     $       ' energy = ',f20.10)
             call util_flush(6)
         endif  ! osome
c
      endif   ! cosmo_on.and.cosmo_phase.eq.2
c
      if(ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c
      end
