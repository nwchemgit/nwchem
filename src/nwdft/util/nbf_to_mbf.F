      Subroutine nbf_to_mbf(geom_hand, basis_hand, rq0, mbf, mcenters, 
     &                      iniz, ifin, zprim, nzprim, acc_gauss)
c
C$Id: nbf_to_mbf.F,v 1.11 1998-12-19 01:29:56 edo Exp $
c
      implicit none
c      
      integer basis_hand, geom_hand
      integer mcenters ! [input]
      integer nzprim ! [input]
      double precision acc_gauss ! [input]
      integer mbf  ! [input/output]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision rq0(mcenters)
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(nzprim)
      integer ncontrset, n1, icset, ictr, itype, nprimo, ncontr, 
     &        isphere, ifirst, ilast, nshbf, l, nang, iprimo
      double precision zmin
c
      if (.not.bas_numcont(basis_hand, ncontrset))
     &   call errquit('Exiting in nbf_to_mbf.',1)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
c
      n1 = 0
c
      do 400 icset = 1,ncontrset
         if (.not.bas_cn2ce(basis_hand, icset, ictr))
     &      call errquit('Exiting in nbf_to_mbf.',2)
c         
c        get info about current contraction set
c         
         if (.not.bas_continfo(basis_hand, icset,
     &      itype, nprimo, ncontr, isphere))
     &      call errquit('Exiting in nbf_to_mbf.',3)
         if (.not. bas_cn2bfr(basis_hand, icset, ifirst, ilast))
     &      call errquit('Exiting in nbf_to_mbf.',4)
c
         nshbf = ilast - ifirst + 1
c         
c        angular momentum
c         
         if (itype .lt. 0)then
            call errquit('nbf_to_mbf: sp-type orbital not coded', 5)
         else 
            l = itype
         endif 
c         
c        cartesian/spherical harmonic
c         
         if (isphere .eq. 0)then !  cartesian set
            nang = (l+1)*(l+2)/2
         elseif (isphere .eq. 1)then !  spherical harmonic
            nang = 2*l+1
         else 
            call errquit('nbf_to_mbf: illegal isphere value', 6)
         endif 
c         
c        get exponents and contraction coefficients for this contraction set
c         
         if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &      call errquit('Exiting in nbf_to_mbf.',7)
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
         if (zmin*rq0(ictr).gt.acc_gauss)then
            go to 400
         else
            if (iniz(ictr).eq.0) iniz(ictr)=n1+1
         endif
         if(isphere.eq.1) then
           n1=n1+ncontr*(2*l+1)
         else
           n1=n1 + ncontr*(((l+1)*(l+2))/2)
         endif
c
         ifin(ictr)= n1
c
  400 continue
c
      mbf = n1
c
      return
      end
