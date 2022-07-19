      Subroutine dftg_gridpart(d_qwght,qwght, qxyz, xyz, xyzm,rq,  
     .     ictr, natoms, iandex,nqpts, 
     ,      whichfunc,
     ,     rad,expo,basis)
c
C$Id$
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
      integer natoms  ![input/out] return atoms with significant contrib
      integer  nqpts  ![input]
      integer basis
c
      double precision qxyz(3,*)! grid points [input]
      double precision xyz(3,*)  ! atom coordinates [input] 
      double precision xyzm(3,*)  ! atom coordinates [out] 
      integer ictr,ictr_scr
      character*4 whichfunc ! [input]
      double precision qwght(*) ! weights [output]
      double precision d_qwght(3,nqpts,*) ! weight deriv [output]
      double precision rq(*)  ! sum of p(n) [ output]
      integer iandex(*)
      double precision expo(*) ! [local]
      double precision rad ! [in]
c
      integer nq_orig,natoms_l,ii
      integer k_rij,l_rij
      integer idp, ldp, 
     &        izeta, lzeta, idzeta, ldzeta,
     &        lbp,ibp,nq3
c
c     works for new grid
c     combine xc_signf
c     
      nq_orig=nqpts

      if(lssw) then
         nq3=3*natoms
      else
         nq3= nqpts*3*natoms
      endif
      if (.not.ma_push_get(mt_dbl, nq3, 'dp',
     &     ldp, idp))call errquit
     &     (' exiting in dftg_gridpart: insuff stack',133,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl, nq3, 
     &     'dzeta', ldzeta, idzeta))call errquit
     &     (' exiting in dftg_gridpart: insuff stack',133, MA_ERR)
      if (.not.ma_push_get(mt_dbl, nqpts, 'zeta',
     &     lzeta, izeta))call errquit
     &     (' exiting in dftg_gridpart: insuff stack',133, MA_ERR)
      if(lssw) then
         nq3=natoms
      else
         nq3= nqpts*natoms
      endif
      if (.not.ma_push_get(mt_dbl,nq3,'bp',lbp,ibp))
     &     call errquit
     &     (' exiting in dftg_quadv0: insuff stack',13, MA_ERR)
      natoms_l=0
c
c     no screening now
c     
      do ii = 1, natoms
            natoms_l=natoms_l+1
            xyzm(1,natoms_l) = xyz(1,ii)
            xyzm(2,natoms_l) = xyz(2,ii)
            xyzm(3,natoms_l) = xyz(3,ii)
            iandex(natoms_l) = ii
      enddo
      if (.not.ma_push_get(mt_dbl,natoms*natoms,'rij',
     ,     l_rij,k_rij)) call errquit
     &     (' exiting in dftg_grssw: insuff stack',13, MA_ERR)

      if(lssw) then
         ictr_scr=ictr
         if(trunc_neigh)
     .        call dftg_signf(basis,   natoms_l, 
     &        xyz,xyzm,ictr_scr,
     &        rad, expo, iandex,iatype_pt_chg)

       if(whichfunc.eq.'erf1') then
          call dftg_gridssw(.true.,d_qwght,qwght, qxyz, xyzm, 
     .        dbl_mb(k_Rij),  rq, dbl_mb(ibp),
     &        dbl_mb(idzeta),
     &        dbl_mb(idp),ictr_scr, natoms_l, 
     .        nq_orig, nqpts,
     ,        .true.)
       else
         call dftg_gridsswold(.true.,d_qwght,qwght, qxyz, xyzm, 
     .        dbl_mb(k_Rij),  rq, dbl_mb(ibp),
     &         dbl_mb(idzeta),
     &        dbl_mb(idp),ictr_scr, natoms_l, 
     .        nq_orig, nqpts,
     ,        whichfunc,.true.)
      endif
      else
         call a_dist(xyz, dbl_mb(k_Rij), natoms,.false.)
         call dftg_beckew(.true.,d_qwght, qwght, qxyz, xyzm,
     &        dbl_mb(k_Rij), rq, dbl_mb(ibp),
     &        dbl_mb(izeta), dbl_mb(idzeta),
     &        dbl_mb(idp), ictr, natoms_l,
     &        nq_orig,nqpts, iandex)
      endif

      if(nq_orig.ne.nqpts) then
         write(6,'(" dftg_part: nc ",4I5)') 
     .        natoms,natoms_l,nq_orig,nqpts
         call errquit ( ' missing ',11, UNKNOWN_ERR)
      endif
cedo      if(natoms.ne.natoms_l) then
cedo      do ii=1,natoms_l
cedo         write(6,'(2I5,3f20.12)')ii,iandex(ii),
cedo     ,        d_qwght(1,1,ii),d_qwght(2,1,ii),
cedo     ,        d_qwght(3,1,ii)
cedo      enddo
cedo      endif
      if (.not.ma_chop_stack(ldp))call errquit
     &     ('dftg_gridparta: cannot pop stack',0, MA_ERR)

      natoms=natoms_l

      return
      end
