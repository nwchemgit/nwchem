      Subroutine nbf_to_mbf(basis_hand,  mbf, mcenters, 
     &     ctr_quad, rad_quad, xyz,
     &                  ibf,iniz, ifin, zprim,  acc_gauss)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand
      integer mcenters ! [input]
      double precision acc_gauss ! [input]
      double precision xyz(3,*) ! [input]
      integer ctr_quad ! grid center  [input]
      integer mbf  ! [input/output]
      integer ibf(*) ! [output]
      double precision rad_quad ! [input]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(*)
      integer ncontrset, n1, icset, ictr,  nprim, ncontr, 
     &        isphere,  l, iprimo,npt,nshbf
      double precision zmin,acc_loc,r2,x,y,z
      integer n,ic1,ic2,m
      double precision alpha,logeps
      double precision gaussian_range,r_arg
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      if (.not.bas_numcont(basis_hand, ncontrset))
     &   call errquit('Exiting in nbf_to_mbf.',1, BASIS_ERR)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
      acc_loc=-acc_gauss
c
      n1 = 0
      npt=0
      do 400 ictr=1,mcenters
        if (.not.bas_ce2cnr(basis_hand,ictr,ic1,ic2))
     &       call errquit('Exiting in xc_signf.',11, BASIS_ERR)
        r2=0d0
        if (ictr.ne.ctr_quad) then
          x = xyz(1,ctr_quad) - xyz(1,ictr)
          y = xyz(2,ctr_quad) - xyz(2,ictr)
          z = xyz(3,ctr_quad) - xyz(3,ictr)
          r2 = sqrt(x*x + y*y + z*z)
        endif
        r_arg=0d0
        if (rad_quad.lt.r2) r_arg = (r2-rad_quad)
        do icset = ic1,ic2
c         
c       get info about current contraction set
c       
          if (.not. bas_continfo(basis_hand, icset,  l ,nprim,
     &         ncontr, isphere))
     &         call errquit('Exiting in xc_signf.',4, BASIS_ERR)
c         
c       
c       get exponents and contraction coefficients for this contraction set
c       
        if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &       call errquit('Exiting in xc_signf.',7, BASIS_ERR)
c       
c       Determine the minimum Gaussian exponent.
c       

        zmin = 1.D+06
        do iprimo = 1,nprim
          zmin = min(zprim(iprimo),zmin)
        enddo
c       
c       Only include those basis functions that are "non-zero" for at least 
c       one  point in the sampling set.
c       
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          else
            nshbf=ncontr*(((l+1)*(l+2))/2)
          endif
        if (r_arg.lt.gaussian_range(l,acc_loc,zmin)) then
          do m=1,nshbf
            ibf(n1+m ) = npt+m
          enddo
          if (iniz(ictr).eq.0) iniz(ictr)=n1+1
          n1=n1+nshbf
        endif
          
        npt=npt+nshbf
        enddo
c           
        ifin(ictr)= n1
c
  400 continue
c
      mbf = n1
c
      return
      end
