      block data crohf_data
#include "crohf.fh"
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c     
c     Invalidate data
c     
      data crohf_init_flag/-1/,
     $     crohf_vlen/-1/
     $     crohf_g_fcv/-1/,
     $     crohf_g_fpv/-1/,
     $     crohf_g_fcp/-1/
      end
      subroutine rohf_init( rtdb)
      implicit none
#include "mafdecls.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "global.fh"
#include "bas.fh"
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c     $Synopsis "Initialize the ROHF internal common block"
c     
      integer rtdb
      external crohf_data ! For T3D linker
c     
      if (crohf_init_flag.gt.0)
     $     call errquit('rohf internals already initialised?',0)
c     
      crohf_init_flag = 1
      crohf_vlen = (nmo-nclosed-nopen)*(nclosed+nopen)+(nclosed*nopen)
c
c     Special case for no degrees of freedom ... avoid breaking GA
c
      if (crohf_vlen .eq. 0) crohf_vlen = 1
c
c     Allocate persistent local and global arrays
c
c     FCV = closed-virtual fock matrix
c     FCV = closed-partial fock matrix
c     FPV = partial-virtual fock matrix
c
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_init: FCV',
     $     32, 32, crohf_g_fcv)) call errquit('rohf_init: FCV', 0)
c     
      if (nopen .gt. 0) then
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_init: FPV',
     $        32, 32, crohf_g_fpv)) call errquit('rohf_init: FPV', 0)
c     
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_init: FCP',
     $        32, 32, crohf_g_fcp)) call errquit('rohf_init: FCP', 0)
      endif
c     
      end

      subroutine rohf_tidy()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
      logical status
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c     $Synopsis "Clean up ROHF internal parameters"
c     
      if (crohf_init_flag.eq.1) then
         status = ga_destroy(crohf_g_fcv)
         if (nopen .gt. 0) then
            status = ga_destroy(crohf_g_fpv) .and. status
            status = ga_destroy(crohf_g_fcp)  .and. status
         endif
         if (.not.(status))
     $        call errquit('rohf_tidy: cannot deallocate handles',0)
c     
c     Invalidate all data
c     
         crohf_g_fcv = -1
         crohf_g_fpv = -1
         crohf_g_fcp = -1
         crohf_init_flag = -1
      endif
      end
      subroutine rohf_energy( g_vecs, eone, etwo, enrep, energy,
     $     g_grad, rtdb )
      implicit none
#include "mafdecls.fh"
#include "geom.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c
      integer rtdb

      integer g_vecs
      double precision energy
      integer g_grad
c     
      double precision eone, etwo, enrep
      integer gtype, grow, gcol
c
      logical odbug
      logical osome
      logical ocosmo
      double precision ecosmo
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_energy: ROHF internal block invalid',0)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.crohf_vlen).or.(gcol.ne.1))
     $     call errquit('rohf_energy: invalid vector length',grow)
c     
      call rohf_fock(rtdb, geom, basis,
     $     nclosed, nopen, nmo, tol2e, g_vecs,
     $     eone, etwo,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp, oskel )
c
*      call rohf_fock2grad( nmo, nclosed, nopen,
*     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
*     $     g_grad)
       call rohf_get_grad(g_grad)
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('rohf_energy: no repulsion energy?', 0)
      energy = eone + etwo + enrep
c
      if (ga_nodeid().eq.0 .and.
     $    util_print('rohf energy',print_debug)) then
       write(6,33) energy, eone, etwo, enrep
33     format(' total energy = ',f20.12/
     $        ' 1-e energy   = ',f20.12/
     $        ' 2-e energy   = ',f20.12/
     $        ' nucl. energy = ',f20.12/)
       call util_flush(6)
      endif 
c
c     ----- cosmo energy term ... -----
c
      odbug=.false.
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
         if(odbug) then
            write(6,*) '-cosmo- ... found in rohf_diis_solve',
     $                 ocosmo,ga_nodeid()
         endif
         if(ocosmo) then
            if (.not. rtdb_get(rtdb,'cosmo:energy',mt_dbl,1,ecosmo))
     $      call errquit('rohf_energy: rtdb get failed for ecosmo',911)
c
            energy=energy+ecosmo
c
            osome=.false.
            osome=osome.and.ga_nodeid().eq.0
            if(osome) then
               write(6,9999) (energy-ecosmo),ecosmo,energy
 9999 format(' energ0 = ',f20.10,/,' ecosmo = ',f20.10,/,
     $       ' energy = ',f20.10)
            endif
         else
            if(odbug) then
               write(6,*) '-cosmo- ... found but .false.',
     $                    ga_nodeid()
            endif
         endif
      else
         if(odbug) then
            write(6,*) '-cosmo- not found in rohf_diis_solve',
     $                 ga_nodeid()
         endif
      endif
      call ga_sync()
c
      end
      subroutine rohf_hessv( acc, g_x, g_ax )
      implicit none
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c
      integer g_x, g_ax
      double precision acc
c     
      integer gtype,grow,gcol,growp,gcolp
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_hessv: ROHF internal block invalid',0)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0)
      call ga_inquire(g_ax,gtype,growp,gcolp)
      if (growp.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0)
      if (gcol.ne.gcolp)
     $     call errquit('rohf_hessv: invalid no. of vectors',0)
c     
c     Call internal routine
c     
      call rohf_hessv_xx( basis, geom, nbf, nmo,
     $     nclosed, nopen,
     $     pflg, g_movecs, oskel, 
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     acc, lshift, g_x, g_ax )
c
c     Zap numbers much smaller than acc to ensure hard zeroes 
c     remain unpolluted ... cannot use a threshold larger than the
c     integral accuracy since can break symmetry in non-abelian groups
c     Also must ensure that the threshold tends to zero to permit
c     tight convergence.
c
*      call ga_screen(g_ax, acc*0.01d0)
      call ga_screen(g_ax, min(acc*acc,acc*0.01d0,1d-12))
c
      end

      subroutine rohf_hessv_precond( g_x, solveshift )
      implicit none
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id: rohf_wrap.F,v 1.16 2000-05-25 21:56:35 mdupuis Exp $
c
      integer g_x
      double precision solveshift
c     
      integer gtype, grow, gcol
      double precision crohf_lshift
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit(
     $     'rohf_hessv_precond: ROHF internal block invalid',0)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv_precond: invalid vector length',0)
c     
      crohf_lshift = lshift - solveshift
      call rohf_hdiag_scale(nmo, nclosed, nopen,
     $     lshift,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     g_x )
c
      end
