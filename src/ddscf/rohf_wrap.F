      block data crohf_data
#include "crohf.fh"
c     
c     $Id$
c     
c     Invalidate data
c     
      data crohf_init_flag/-1/,
     $     crohf_vlen/-1/
     $     crohf_g_fcv/-1/,
     $     crohf_g_fpv/-1/,
     $     crohf_g_fcp/-1/,
     $     noskew/.true./
      end
      subroutine rohf_init( rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "global.fh"
#include "bas.fh"
c     
c     $Id$
c     $Synopsis "Initialize the ROHF internal common block"
c     
      integer rtdb
      external crohf_data ! For T3D linker
c     
      if (crohf_init_flag.gt.0)
     $     call errquit('rohf internals already initialised?',0,
     &       UNKNOWN_ERR)
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
     $     32, 32, crohf_g_fcv)) call errquit('rohf_init: FCV', 0,
     &       GA_ERR)
c     
      if (nopen .gt. 0) then
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_init: FPV',
     $        32, 32, crohf_g_fpv)) call errquit('rohf_init: FPV', 0,
     &       GA_ERR)
c     
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_init: FCP',
     $        32, 32, crohf_g_fcp)) call errquit('rohf_init: FCP', 0,
     &       GA_ERR)
      endif
c     
      end

      subroutine rohf_tidy()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
      logical status
c     
c     $Id$
c     $Synopsis "Clean up ROHF internal parameters"
c     
      if (crohf_init_flag.eq.1) then
         status = ga_destroy(crohf_g_fcv)
         if (nopen .gt. 0) then
            status = ga_destroy(crohf_g_fpv) .and. status
            status = ga_destroy(crohf_g_fcp)  .and. status
         endif
         if (.not.(status))
     $        call errquit('rohf_tidy: cannot deallocate handles',0,
     &       GA_ERR)
c     
c     Invalidate all data
c     
         crohf_g_fcv = -1
         crohf_g_fpv = -1
         crohf_g_fcp = -1
         crohf_init_flag = -1
      endif
      end
      subroutine rohf_energy( g_vecs, eone, etwo, enrep, ecosmo, energy,
     $     g_grad, rtdb )
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
!
      integer g_vecs
      double precision energy
      integer g_grad
c     
      double precision eone, etwo, enrep
      integer gtype, grow, gcol
c
      logical odbug
      logical osome
      double precision ecosmo !< [Output] The COSMO solvation energy
c     DIM/QM JEM
      logical ldimqm
      double precision edimqm


      odbug=.false.
      osome=.false.
      osome=osome.or.odbug
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_energy: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.crohf_vlen).or.(gcol.ne.1))
     $     call errquit('rohf_energy: invalid vector length',grow,
     &       UNKNOWN_ERR)
c     
      if(ga_nodeid().eq.0.and.osome) then
         write(6,*) 'in  rohf_energy calling rohf_fock ...'
         call util_flush(6)
      endif

      call rohf_fock(rtdb, geom, basis,
     $     nclosed, nopen, nmo, tol2e, g_vecs,
     $     eone, etwo,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp, oskel, iter )

      if(ga_nodeid().eq.0.and.osome) then
         write(6,*) 'in  rohf_energy returning from rohf_fock ...'
         call util_flush(6)
      endif
c
*      call rohf_fock2grad( nmo, nclosed, nopen,
*     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
*     $     g_grad)

      if(ga_nodeid().eq.0.and.osome) then
         write(6,*) 'in  rohf_energy calling rohf_get_grad ...'
         call util_flush(6)
      endif

      call rohf_get_grad(g_grad)

      if(ga_nodeid().eq.0.and.osome) then
         write(6,*) 'in  rohf_energy returning from rohf_get_grad ...'
         call util_flush(6)
      endif
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('rohf_energy: no repulsion energy?', 0,
     &       UNKNOWN_ERR)
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
      odbug=odbug.and.ga_nodeid().eq.0
      if(cosmo_on.and.cosmo_phase.eq.2) then
         if (.not. rtdb_get(rtdb,'cosmo:energy',mt_dbl,1,ecosmo))
     $     call errquit('rohf_energy: rtdb get failed for ecosmo',911,
     &       RTDB_ERR)
c
         energy=energy+ecosmo
c
         osome=.false.
         osome=osome.and.ga_nodeid().eq.0
         if(osome) then
             write(6,9999) (energy-ecosmo),ecosmo,energy
 9999 format(' energ0 = ',f20.10,/,' ecosmo = ',f20.10,/,
     $       ' energy = ',f20.10)
             call util_flush(6)
         endif
      endif
c     
c     DIM/QM JEM
c     DIM/QM energy term
      if (.not. rtdb_get(rtdb, 'dimqm:ldimqm', MT_LOG, 1, ldimqm))
     $  ldimqm = .false. 
      if (ldimqm) then
        if (.not. rtdb_get(rtdb, 'dimqm:edimqm', MT_DBL, 1, edimqm))
     $    call errquit('rohf_energy: dimqm energy', 1, RTDB_ERR)
        energy = energy + edimqm
        if(ga_nodeid().eq.0.and.osome) then
          write(6,9998) (energy-edimqm), edimqm, energy
 9998 format(' energ0 = ',f20.10,/,' edimqm = ',f20.10,/,
     $       ' energy = ',f20.10)
          call util_flush(6)
        end if
      end if
      if(ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c
      end
      subroutine rohf_hessv( acc, g_x, g_ax )
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id$
c
      integer g_x, g_ax
      double precision acc
c     
      integer gtype,grow,gcol,growp,gcolp
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_hessv: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax,gtype,growp,gcolp)
      if (growp.ne.crohf_vlen)
     $     call errquit('rohf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (gcol.ne.gcolp)
     $     call errquit('rohf_hessv: invalid no. of vectors',0,
     &       UNKNOWN_ERR)
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
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id$
c
      integer g_x
      double precision solveshift
c     
      integer gtype, grow, gcol
      double precision crohf_lshift
c
      logical oprint, olprint
      integer ilo(2), ihi(2)
      integer icol
      double precision dnrm
c     
c     Check
c
      oprint = util_print('precond',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
      if (crohf_init_flag.ne.1)
     $     call errquit(
     $     'rohf_hessv_precond: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv_precond: invalid vector length',0,
     &       UNKNOWN_ERR)
c
      if (oprint) then
        do icol = 1, gcol
          ilo(1) = 1
          ilo(2) = icol
          ihi(1) = grow
          ihi(2) = icol
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"rohf_wrap: in g_x = ",i4,f24.8)')
     +      icol,dnrm
          endif
        enddo
      endif
c     
      crohf_lshift = lshift - solveshift
      call rohf_hdiag_scale(nmo, nclosed, nopen,
     $     lshift,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     g_x)
c
      if (oprint) then
        do icol = 1, gcol
          ilo(1) = 1
          ilo(2) = icol
          ihi(1) = grow
          ihi(2) = icol
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"rohf_wrap: out g_x = ",i4,f24.8)')
     +      icol,dnrm
          endif
        enddo
      endif
c
      end

      subroutine rohf_hessv_precond3( g_x, omega)
c
c ... jochen: modified version of rohf_hessv_precond to go along with
c     cphf_solve3. This preconditioner PRESUMABLY sets up the rhs of
c     the CPKS problem which is where we should add and subtract
c     the frequency from the energy denominator
c
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id$
c
      integer g_x               ! [in/ou] rhs vector of eq. system
      double precision omega    ! [input] frequency
c     
      integer gtype, grow, gcol
      double precision crohf_lshift
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit(
     $     'rohf_hessv_precond: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv_precond: invalid vector length',0,
     &       UNKNOWN_ERR)
c     
c ... jochen: in contrast to rohf_hessv_precond we add
c     omega here and use the result in the call to
c     rohf_hdiag_scale. Note: in rohf_hdiag_scale there is a division
c     by [ 4 (epsilon_a - epsilon_i) + omega ] so we should use
c     four times omega here to be consistent with the TDDFT formalism.
c
c ... jochen: Note:    
c     level shifting might have some strange effects here 
c     
      crohf_lshift = lshift + 4d0 * omega
c
      call rohf_hdiag_scale(nmo, nclosed, nopen,
     $   crohf_lshift,
     $   crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $   g_x )
c     
      end

c     =========================================================================
c     =========================================================================

      subroutine rohf_hessv_precond3_damp(g_x, g_x_im, omega, gamwidth)
c
c ... jochen: modified version of rohf_hessv_precond to go along with
c     cphf_solve3. This preconditioner PRESUMABLY sets up the rhs of
c     the CPKS problem which is where we should add and subtract
c     the frequency from the energy denominator.
c     NOTE: this is a modified version of rohf_hessv_precond3 for
c     calculations with finite-lifetime damping
c
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id$
c
      integer g_x               ! [in/ou] Re rhs vector of eq. system
      integer g_x_im            ! [in/ou] Im rhs vector of eq. system
      double precision omega    ! [input] frequency
      double precision gamwidth ! [input] damping constant
c     
      integer gtype, grow, gcol
      double precision crohf_lshift, crohf_gamwidth
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit(
     $     'rohf_hessv_precond: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv_precond: invalid vector length',0,
     &       UNKNOWN_ERR)
c     
c ... jochen: in contrast to rohf_hessv_precond we add
c     omega here and use the result in the call to
c     rohf_hdiag_scale. Note: in rohf_hdiag_scale there is a division
c     by [ 4 (epsilon_a - epsilon_i) + omega ] so we should use
c     four times omega here to be consistent with the TDDFT formalism.
c
c ... jochen: Note:    
c     level shifting might have some strange effects here

c     jochen 02/2016: un-did some of the changes introduced by Lasse's
c     student.  Fredy applied those during the spin-unrestricted
c     extension, the code lines are marked by FA-03-12-14.  There are
c     matching changes in rohf_hdiag_scale_damp
c     
      crohf_lshift = lshift + 4d0 * omega
c      crohf_lshift = lshift + omega ! FA-03-12-14
c
c     use 4 * lifetime parameter to be consistent with
c     using 4 * omega
      crohf_gamwidth = 4d0 * gamwidth
c      crohf_gamwidth = gamwidth ! FA-03-12-14
c
      call rohf_hdiag_scale_damp(nmo, nclosed, nopen,
     $   crohf_lshift, crohf_gamwidth, 
     $   crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $   g_x, g_x_im)
c
      end
