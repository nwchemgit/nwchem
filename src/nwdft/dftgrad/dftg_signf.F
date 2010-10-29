      subroutine dftg_signf(basis,  mcenters, 
     &                      xyz,xyzm,ctr_quad, rad_quad,
     &                       zprim, iandex,iatype_pt_chg)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
      integer basis
      integer ctr_quad ! atom where the grid is centered [input]
      integer mcenters ! [input]
      double precision xyz(3,*) ! [input]
      double precision xyzm(3,*) ! [input/output]
      double precision rad_quad ! [input]
      logical iatype_pt_chg(*) ! [input]
      integer iandex(*) ! [scratch]
    
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision zprim(*),acc_sigf
      integer icset, ictr, 
     &        ifirst, ilast, nprim,  iprimo
      double precision zmin,r2,x,y,z,r_arg
      integer mcenters_scr,ctr_out
      integer ic1,ic2,dum
      integer n,l
      double precision alpha,logeps
      double precision gaussian_range
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)

c
c
c
      acc_sigf=log(1d-10)
      mcenters_scr = 0
c
      do ictr=1,mcenters
        if(iatype_pt_chg(ictr)) goto 2001
        if (ictr.eq.ctr_quad) then
          mcenters_scr=mcenters_scr+1
         iandex(mcenters_scr)=ictr
          goto 2001
        endif
        if (.not.bas_ce2cnr(basis,ictr,ic1,ic2))
     &       call errquit('Exiting in xc_signf.',11, BASIS_ERR)
c
c     no bf on center ie but we keep grid
c
        if(ic1.eq.0) then
           mcenters_scr=mcenters_scr+1
           iandex(mcenters_scr)=ictr
           goto 2001
        endif

        x = xyz(1,ctr_quad) - xyz(1,ictr)
        y = xyz(2,ctr_quad) - xyz(2,ictr)
        z = xyz(3,ctr_quad) - xyz(3,ictr)
        r2 = sqrt(x*x + y*y + z*z)
        r_arg=0d0
        if (rad_quad.lt.r2) r_arg = (r2-rad_quad)
        do icset = ic1,ic2
c         
c       get info about current contraction set
c       
        if (.not. bas_cn2bfr(basis, icset, ifirst, ilast))
     &       call errquit('Exiting in xc_signf.',4, BASIS_ERR)
            if (.not. bas_continfo(basis, icset,  l ,nprim,
     &                             dum, dum))
     &       call errquit('Exiting in xc_signf.',4, BASIS_ERR)
c       
c       
c       get exponents and contraction coefficients for this contraction set
c       
        if (.not.bas_get_exponent(basis, icset, zprim))
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
c        if (zmin*r_arg.lt.acc_sigf) then
        if (r_arg.lt.gaussian_range(l,acc_sigf,zmin)) then
          mcenters_scr=mcenters_scr+1
          iandex(mcenters_scr)=ictr
          goto 2001
        endif

      enddo
 2001 continue
      enddo
      ctr_out = 1                ! take care of compiler warnings
      do ictr = 1, mcenters_scr
        ic1=iandex(ictr)
        xyzm(1,ictr) = xyzm(1,ic1)
        xyzm(2,ictr) = xyzm(2,ic1)
        xyzm(3,ictr) = xyzm(3,ic1)
        if(ic1.eq.ctr_quad) ctr_out=ictr
      enddo
      mcenters=mcenters_scr
      ctr_quad=ctr_out
      return
      end
