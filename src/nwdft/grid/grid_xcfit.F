      subroutine grid_xcfit(xc_bas_han,nbf_xc,tmat,grad,
     ,     iniz,ifin,ibf,coef,expo,nscr,iscr,hess,
     ,     amat,cmat,rq,rq0,nbf_xc_mxprim, ncoef_max,
     ,     xyz,qxyz,nqpts,natoms,ipol,ictr,rad,
     ,     acc_xc_gauss)
      implicit none
c
c$Id: grid_xcfit.F,v 1.2 2000-04-21 19:16:47 edo Exp $
c
#include "mafdecls.fh"
      logical grad
      integer xc_bas_han ! [in]
      integer nbf_xc ! [in]
      integer ipol ! [in]
      integer nqpts,natoms
      integer ictr ! [in]
      double precision rad ! [in]
      double precision tmat(nbf_xc,*)
      integer iniz(natoms),ifin(natoms),ibf(nbf_xc)
      double precision coef(*),expo(*)
      integer nscr ! [in]
      double precision amat(nqpts,ipol)! [in]
      double precision Cmat(nqpts,3,*)
      double precision rq(*),rq0(*)
      double precision xyz(*),qxyz(*)
      double precision acc_xc_gauss ! [in]
      double precision iscr(nscr)
      double precision hess(*)
      integer nbf_xc_mxprim, ncoef_max
c  
      integer k,mbf_xc,itpmat,ltpmat
      integer lchi_xc, ichi_xc, 
     &        ldelchi_xc, idelchi_xc,l3nqscr, i3nqscr
      integer itpp
c
      mbf_xc = nbf_xc
c     
      call ifill(mbf_xc, 0, ibf, 1)
      call nbf_to_mbf( XC_bas_han,  mbf_xc,
     &     natoms, 
     &     ictr,rad,xyz,
     &     ibf,iniz, ifin,
     &     expo,acc_XC_gauss)
      if (mbf_xc.eq.0)return
c
      if (.not.MA_Push_Get(MT_Dbl, mbf_xc*ipol, 'Tpmat', lTpmat,
     &     iTpmat))
     &     call errquit('grid_xcfit: cannot allocate Tpmat',0)
      call dfill(mbf_xc*ipol,0.D0,dbl_mb(iTpmat),1)
      if(grad) then
         if (.not.MA_Push_Get(MT_Dbl,3*nqpts,
     &     'nqsmallx3 scr',l3nqscr,i3nqscr))call errquit
     &     ('grid_xcfit: cannot allocate nqsmallx3 scr',0)
      endif
c
      if (.not.MA_Push_Get(MT_Dbl, nqpts*mbf_xc, 'chi_xc', 
     &     lchi_xc, ichi_xc))
     &     call errquit('grid_xcfit: cannot allocate chi_xc',0)
c     
      if (.not.MA_Push_Get(MT_Dbl, 3*nqpts*mbf_xc, 'delchi_xc',
     &     ldelchi_xc, idelchi_xc))
     &     call errquit('grid_xcfit: cannot allocate delchi_xc',0)

      call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
      call new_eval_gbsets(xc_bas_han, grad, .false.,
     &     Dbl_MB(ichi_xc), Dbl_MB(idelchi_xc), 
     &     hess, ibf, 0, rq0, 
     .     rq, qxyz, xyz,nqpts, nqpts, 
     ,     mbf_xc, natoms,iscr,nscr,
     ,     expo,nbf_xc_mxprim, coef, ncoef_max,   
     &              iniz, ifin)
c     
c     
c     Gather T 
c     
      do  k = 1,ipol
         call dft_gather(mbf_xc, 
     &        Dbl_MB(iTpmat+(k-1)*mbf_xc),
     &        Tmat(1,k), ibf, 0)
      enddo
      call dgemv('T', nqpts, mbf_xc, 1.D0, 
     &     DBL_MB(ichi_xc), nqpts, Amat, 1, 
     &     1.D0, dbl_mb(iTpmat), 1)
      if (ipol.eq.2)then
         iTpp = iTpmat + mbf_xc
         call dgemv('T', nqpts, mbf_xc, 1.D0, 
     &        DBL_MB(ichi_xc), nqpts, Amat, 1, 
     &        1.D0, dbl_mb(iTpp), 1)
      endif
c     
c     C*D product ... see tabcd
c     
      if (GRAD)then
c     
c     get scratch vector to gather the 3 components for the 
c     dgemv
c     
         call dcopy(nqpts,cmat(1,1,1),1, dbl_mb(i3nqscr),1)
         call dcopy(nqpts,cmat(1,2,1),1,
     &        dbl_mb(i3nqscr+nqpts),1)
         call dcopy(nqpts,cmat(1,3,1),1,
     &        dbl_mb(i3nqscr+2*nqpts),1)
         call dgemv('T', 3*nqpts, mbf_xc, 1.D0, 
     &        DBL_MB(idelchi_xc), 3*nqpts, 
     &        dbl_mb(i3nqscr), 1, 1.D0, 
     &        dbl_mb(iTpmat), 1)
c     
         if (ipol.eq.2)then
            call dcopy(nqpts,cmat(1,1,2),1,dbl_mb(i3nqscr),1)
            call dcopy(nqpts,cmat(1,2,2),1,
     &           dbl_mb(i3nqscr+nqpts),1)
            call dcopy(nqpts,cmat(1,3,2),1,
     &           dbl_mb(i3nqscr+2*nqpts),1)
            iTpp = iTpmat + mbf_xc
            call dgemv('T', 3*nqpts, mbf_xc, 1.D0, 
     &           DBL_MB(idelchi_xc), 3*nqpts, 
     &           dbl_mb(i3nqscr), 1, 1.D0, 
     &           dbl_mb(iTpp), 1)
         endif
      endif
      do  k = 1,ipol
         call dft_scatter(mbf_xc, Tmat(1,k), ibf,
     &        Dbl_MB(iTpmat+(k-1)*mbf_xc))
      enddo
      
 2010 continue
      if (.not.ma_chop_stack(ltpmat))
     &     call errquit('grid_fitxc: cannot chop stack',19)
            
      return
      end
