      Subroutine dftg_gridssw(d_qwght,qwght, qxyz, xyz, xyzm,rq,  
     .     ictr, natoms, iandex,nqpts, 
     ,      whichfunc,
     ,     rad,expo)
c
C$Id: dftg_gridssw.F,v 1.2 2000-04-15 01:24:22 edo Exp $
c
      implicit none
#include "mafdecls.fh"
#include "cdft.fh"
      integer natoms  ![input/out] return atoms with significant contrib
      integer  nqpts  ![input]
c
      double precision qxyz(3,*)! grid points [input]
      double precision xyz(3,*)  ! atom coordinates [input] 
      double precision xyzm(3,*)  ! atom coordinates [out] 
      integer ictr
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
     &        lbp,ibp,lcmp,icmp
c
c     works for new grid
c     combine xc_signf
c     
      nq_orig=nqpts

      if (.not.ma_push_get(mt_dbl, nqpts*3*natoms, 'dp',
     &     ldp, idp))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl, nqpts, 'zeta',
     &     lzeta, izeta))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl, 3*nqpts*natoms, 
     &     'dzeta', ldzeta, idzeta))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl,nqpts*natoms,'bp',lbp,ibp))
     &     call errquit
     &     (' exiting in dftg_quadv0: insuff stack',13)
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
     &     (' exiting in dftg_grssw: insuff stack',13)
      call a_dist(xyz, dbl_mb(k_Rij), natoms)

      
      
      if (.not.MA_Push_Get(mt_int,nqpts,'compr scr',lcmp,icmp))
     &     call errquit('dftg_grssw: cannot allocate compr',0)
      if(lssw) then
         call xc_signf(AO_bas_han,  nbf_ao, iandex,natoms_l, 
     &        xyz,xyzm,ictr,
     &        rad,
     &        expo, -99d0,iatype_pt_chg)
         call dftg_ssw(.true.,d_qwght,qwght, qxyz, xyzm, 
     .        dbl_mb(k_Rij),  rq, dbl_mb(ibp),
     &                    dbl_mb(izeta), dbl_mb(idzeta),
     &        dbl_mb(idp),ictr, natoms_l, 
     .        nq_orig, nqpts,iandex,
     ,        int_mb(icmp),tol_rho,whichfunc)
      else
         call dftg_beckew(.true.,d_qwght, qwght, qxyz, xyzm,
     &        dbl_mb(k_Rij), rq, dbl_mb(ibp),
     &        dbl_mb(izeta), dbl_mb(idzeta),
     &        dbl_mb(idp), ictr, natoms_l,
     &        nq_orig,nqpts, iandex)
      endif

      if(nq_orig.ne.nqpts) then
         write(6,'(" nc ",4I5)') natoms,natoms_l,nq_orig,nqpts
         call errquit ( ' missing ',11)
      endif
cedo      if(natoms.ne.natoms_l) then
cedo      do ii=1,natoms_l
cedo         write(6,'(2I5,3f20.12)')ii,iandex(ii),
cedo     ,        d_qwght(1,1,ii),d_qwght(2,1,ii),
cedo     ,        d_qwght(3,1,ii)
cedo      enddo
cedo      endif
      if (.not.ma_chop_stack(ldp))call errquit
     &     ('dftg_gridsswa: cannot pop stack',0)

      natoms=natoms_l

      return
      end
