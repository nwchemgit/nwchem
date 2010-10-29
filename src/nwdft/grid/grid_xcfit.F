      subroutine grid_xcfit(xc_bas_han,nbf_xc,tmat,grad,
     ,     iniz,ifin,ibf,coef,expo,nscr,iscr,hess,
     ,     amat,cmat,rq,rq0,nbf_xc_mxprim, ncoef_max,
     ,     xyz,qxyz,nqpts,natoms,ipol,ictr,rad,
     .     bas_cent_xinfo,bas_cset_xinfo,docsetx,
     ,     acc_xc_gauss)
      implicit none
#include "errquit.fh"
c
c$Id$
c
#include "mafdecls.fh"
#include "bas.fh"
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
      integer nbf_xc_mxprim, ncoef_max,maxdbas,i
      integer bas_cent_xinfo(3,natoms),bas_cset_xinfo(6,*)
      logical docsetx(*)
c  
      integer k,mbf_xc,itpmat,ltpmat
      integer lchi_xc, ichi_xc, 
     &        ldelchi_xc, idelchi_xc,l3nqscr, i3nqscr
      integer itpp,ncontrsetx
c
      mbf_xc = nbf_xc
      if (grad) then
         maxdbas = 1
      else
         maxdbas = 0
      endif
c     
      call ifill(mbf_xc, 0, ibf, 1)
      mbf_xc=nbf_xc
      if (.not.bas_numcont(xc_bas_han, ncontrsetx))
     &     call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
      if(natoms.gt.1) then
      call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
      Call oldnBF_to_mBF(xc_bas_han,rq0,
     &                mbf_xc,natoms,
     &                iniz,ifin,ibf,
     &                Expo,nbf_xc_mxprim,
     .        docsetx,acc_xc_gauss)
         if (mbf_xc.eq.0) return
      else
         iniz(1)=1
         ifin(1)=nbf_xc
         do i=1,ncontrsetx
            docsetx(i)=.true.
         enddo
         do i=1,nbf_xc
            ibf(i)=i
         enddo
      endif
      if (mbf_xc.eq.0)return
c
      if (.not.MA_Push_Get(MT_Dbl, mbf_xc*ipol, 'Tpmat', lTpmat,
     &     iTpmat))
     &     call errquit('grid_xcfit: cannot allocate Tpmat',0, MA_ERR)
      call dfill(mbf_xc*ipol,0.D0,dbl_mb(iTpmat),1)
      if(grad) then
         if (.not.MA_Push_Get(MT_Dbl,3*nqpts,
     &     'nqsmallx3 scr',l3nqscr,i3nqscr))call errquit
     &     ('grid_xcfit: cannot allocate nqsmallx3 scr',0, MA_ERR)
      endif
c
      if (.not.MA_Push_Get(MT_Dbl, nqpts*mbf_xc, 'chi_xc', 
     &     lchi_xc, ichi_xc))
     &     call errquit('grid_xcfit: cannot allocate chi_xc',0, MA_ERR)
c     
      if (.not.MA_Push_Get(MT_Dbl, 3*nqpts*mbf_xc, 'delchi_xc',
     &     ldelchi_xc, idelchi_xc))
     &     call errquit('grid_xcfit: cannot allocate delchi_xc',0,
     &       MA_ERR)

      call xc_eval_basis(xc_bas_han, maxdbas, dbl_mb(ichi_xc),
     &     dbl_mb(idelchi_xc), 0, 0,
     ,     rq, qxyz,xyz,
     ,     nqpts, natoms, iniz, docsetx, 
     .     bas_cent_xinfo, bas_cset_xinfo)
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
     &     call errquit('grid_fitxc: cannot chop stack',19, MA_ERR)
            
      return
      end
      Subroutine oldnbf_to_mbf(basis_hand, rq0, mbf, mcenters, 
     &                      iniz, ifin, ibf, zprim, nzprim, 
     .     docset, acc_gauss)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand
      integer mcenters ! [input]
      integer nzprim ! [input]
      double precision acc_gauss ! [input]
      logical docset(*)
      integer mbf  ! [input/output]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision rq0(mcenters)
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(nzprim)
      integer ncontrset, n1, icset, ictr, nprimo, ncontr, 
     &        isphere, nshbf, l, nang, iprimo
      integer n,npt,m
      integer ibf(*) ! [output]
      double precision zmin
      double precision alpha,logeps,acc_loc
      double precision gaussian_range,r_arg
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      acc_loc=-acc_gauss
      
      if (.not.bas_numcont(basis_hand, ncontrset))
     &   call errquit('Exiting in nbf_to_mbf.',1, BASIS_ERR)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
c
      n1 = 0
      npt=0
c
      do 400 icset = 1,ncontrset
         docset(icset) = .false.
         if (.not.bas_cn2ce(basis_hand, icset, ictr))
     &      call errquit('Exiting in nbf_to_mbf.',2, BASIS_ERR)
c         
c        get info about current contraction set
c         
         if (.not.bas_continfo(basis_hand, icset,
     &      l, nprimo, ncontr, isphere))
     &      call errquit('Exiting in nbf_to_mbf.',3, BASIS_ERR)
c         
         if (isphere .eq. 0)then !  cartesian set
            nang = (l+1)*(l+2)/2
         elseif (isphere .eq. 1)then !  spherical harmonic
            nang = 2*l+1
         endif 
c         
c        get exponents and contraction coefficients for this contraction set
c         
         if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &      call errquit('Exiting in nbf_to_mbf.',7, BASIS_ERR)
c
c        Determine the minimum Gaussian exponent.
c
         zmin = 1.D+06
         do 300 iprimo = 1,nprimo
            zmin = min(zprim(iprimo),zmin)
  300    continue
c
c        Only include those basis functions that are "non-zero" for at least one 
c        point in the sampling set.
c
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          else
            nshbf=ncontr*(((l+1)*(l+2))/2)
          endif
          r_arg=sqrt(rq0(ictr))
        if (r_arg.lt.gaussian_range(l,acc_loc,zmin)) then

          do m=1,nshbf
            ibf(n1+m ) = npt+m
          enddo
          docset(icset) = .true.
            if (iniz(ictr).eq.0) iniz(ictr)=n1+1
            n1 = n1 + nshbf
#ifdef DEBUG
         else
            write(0,*) ' got no '
#endif
         endif
c
         ifin(ictr)= n1
c
         npt=npt+nshbf
  400 continue
c
      mbf = n1
      return 
      end
c
