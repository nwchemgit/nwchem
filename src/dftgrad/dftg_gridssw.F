      Subroutine dftg_gridssw(d_qwght,qwght, qxyz, xyz, rq,  
     .     ictr, ncenters, iandex,nqpts, 
     ,     tol_rho, whichfunc)
c
C$Id: dftg_gridssw.F,v 1.1 2000-04-14 01:54:17 edo Exp $
c
      implicit none
#include "mafdecls.fh"
      integer ncenters, nqpts  ![input]
c
      double precision qxyz(3,*)! grid points [input]
      double precision xyz(3,*)  ! atom coordinates [input] 
      integer ictr
      character*4 whichfunc ! [input]
      double precision qwght(*) ! weights [output]
      double precision d_qwght(*) ! weight deriv [output]
      double precision rq(*)  ! sum of p(n) [ output]
      integer iandex(*)
      double precision tol_rho ! [in]
c
      integer nq_orig,ncenters_l,ii
      integer k_rij,l_rij
      integer idp, ldp, 
     &        izeta, lzeta, idzeta, ldzeta,
     &        lbp,ibp,lcmp,icmp
c
c     works for new grid
c     combine xc_signf
c     
      nq_orig=nqpts

      if (.not.ma_push_get(mt_dbl, nqpts*3*ncenters, 'dp',
     &     ldp, idp))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl, nqpts, 'zeta',
     &     lzeta, izeta))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl, 3*nqpts*ncenters, 
     &     'dzeta', ldzeta, idzeta))call errquit
     &     (' exiting in dftg_gridssw: insuff stack',133)
      if (.not.ma_push_get(mt_dbl,nqpts*ncenters,'bp',lbp,ibp))
     &     call errquit
     &     (' exiting in dftg_quadv0: insuff stack',13)
      ncenters_l=ncenters
c
c     no screening now
c     
      do ii = 1, ncenters
!         xyzm(1,ii) = xyz(1,ii)
!         xyzm(2,ii) = xyz(2,ii)
!         xyzm(3,ii) = xyz(3,ii)
         iandex(ii) = ii
      enddo
      if(ncenters.ne.ncenters_l.or.nq_orig.ne.nqpts) then
      write(6,'(" nc ",4I5)') ncenters,ncenters_l,nq_orig,nqpts
      call errquit ( ' missing ',11)
      endif
      if (.not.ma_push_get(mt_dbl,ncenters*ncenters,'rij',
     ,     l_rij,k_rij)) call errquit
     &     (' exiting in dftg_grssw: insuff stack',13)
      call a_dist(xyz, dbl_mb(k_Rij), ncenters)
      if (.not.MA_Push_Get(mt_int,nqpts,'compr scr',lcmp,icmp))
     &     call errquit('dftg_grssw: cannot allocate compr',0)
      call dftg_ssw(.true.,d_qwght,qwght, qxyz, xyz, 
     .                    dbl_mb(k_Rij),  rq, dbl_mb(ibp),
     &                    dbl_mb(izeta), dbl_mb(idzeta),
     &                    dbl_mb(idp),ictr, ncenters_l, 
     .                    nq_orig, nqpts,iandex,
     ,                    int_mb(icmp),tol_rho,whichfunc)
      if (.not.ma_chop_stack(ldp))call errquit
     &     ('dftg_gridsswa: cannot pop stack',0)
      return
      end
