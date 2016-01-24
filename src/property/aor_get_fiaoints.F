      subroutine aor_get_fiaoints(
     &  g_dipel,                ! out: el dipole 
     &  g_dipmag,               ! out: mag dipole
     &  g_smat0,                ! out: basis overlap
     &  g_quadel,               ! out: el-quadrupole
     &  rtdb,                   ! in : rtdb  handle
     &  basis,                  ! in : basis handle
     &  geom,                   ! in : geom  handle
     &  lvelocity,              ! in : velocity flag
     &  lgiao,                  ! in : giao flag  
     &  lbeta,                  ! in : orbeta flag
     &  lquad,                  ! in : quad flag
     &  nbf)                    ! in : nr. basis functions
      
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c     updated by JA 01/16

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "stdio.fh"
#include "prop.fh"
#include "util.fh"
c     calculate required one-electron integrals in AO basis.
c     There is some redundancy the the AO matrix calculations
c     that we might get rid of later (e.g. dipole integrals 
c     are also used in fiao_e1_movecs). 
      integer rtdb,basis,geom
      integer g_dipel, ! out: electron dipole moment AO  integral
     &        g_dipmag,! out: magnetic dipole moment AO  integral
     &        g_smat0, ! out: zeroth-order basis overlap integral
     &        g_smat1, ! out: 1st-order B-field perturb  integral
     &        g_sket1, ! out: 1st-order S-ket matrix
     &        g_quadel ! out: e-quadrupole    moment AO  integral

      logical lvelocity,
     &        lbeta,
     &        lgiao,
     &        lquad      
      integer alo(3),ahi(3),
     &        blo(3),bhi(3)
      integer nbf,ndir
      double precision dummy(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      logical oskel,status,debug
      integer sh_atom
      parameter (sh_atom=1)
      integer nat
      parameter (nat=1)
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb (oskel = .false.)
      oskel=.false.
      ndir=3 ! nr directions (x,y,z)
c     Integral initialization 
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      status = rtdb_parallel(.true.)

      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
c     electron dipole moment integrals:
      if (.not.nga_create(MT_DBL,3,ahi,'e-dipole',
     &                    alo,g_dipel)) call 
     &    errquit('aor_driver: nga_create failed g_dipel',0,GA_ERR)
      call ga_zero(g_dipel)
c     note: origin has been set to (0,0,0) for multipole integs.
      if (lvelocity) then

        call int_giao_1ega(basis,basis,
     &                     g_dipel,'velocity', ! out
     &                     origin,1,oskel)

c       invert sign of velocity integrals for compatibility 
c       w/ dipole-length optical rotation
        call ga_scale (g_dipel, -1d0) 

      else

        call int_mpole_1ega(basis,basis,
     &                      g_dipel,'dipole',  ! out
     &                      origin,oskel)

      endif

c     zeroth-order basis-overlap integrals:
      if (.not.nga_create(MT_DBL,2,ahi,'smat0',
     &                    alo,g_smat0)) call 
     &    errquit('aor_driver: nga_create failed g_smat0',0,GA_ERR)
      call ga_zero(g_smat0)

      call int_1e_ga(basis,basis,
     &               g_smat0,'overlap', ! out
     &  oskel)

c     magnetic dipole moment integrals (r x nabla), times 0.5:
      
      if (.not.lbeta) then ! TEMP .and. .not.lgiao) then
 
        if (.not.nga_create(MT_DBL,3,ahi,'mag-dipole',
     &                      alo,g_dipmag)) call errquit
     &     ('aor_driver: nga_create failed g_dipmag',0,GA_ERR)
        call ga_zero(g_dipmag)   

          call int_giao_1ega(basis,basis,
     &                       g_dipmag,'angmom', ! out
     &                       origin,1,oskel)
          call ga_scale (g_dipmag, 0.5d0)

        end if                  ! lgiao



c     electron quadrupole moment integrals
      if (lquad) then
        alo(1) = nbf
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbf
        ahi(2) = nbf
        ahi(3) = 6 
        if (.not.nga_create(MT_DBL,3,ahi,'e-quad',
     &                      alo,g_quadel)) call 
     &     errquit('aor_driver: nga_create failed g_quadel',0,GA_ERR)
        call ga_zero(g_quadel)

        call int_mpole_1ega(basis,basis,
     &                      g_quadel,'quadrupole',
     &                      origin,oskel)

      end if
      call schwarz_tidy()
      call int_terminate()
      return
      end
c $Id$
