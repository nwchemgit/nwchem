      subroutine hnd_jkbcdf(b00,b01,b10,c00,d00,f00,
     1                  abv,cv,rwv,numg)
c
c $Id$
c
      implicit double precision (a-h,o-z)
c
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_rys.fh"
#include "case.fh"
c
      dimension b00(numg,nroots),b01(numg,nroots),b10(numg,nroots)
      dimension c00(numg,nroots,3)
      dimension d00(numg,nroots,3)
      dimension f00(numg,nroots)
      dimension abv(6,numg),cv(12,numg)
      dimension rwv(2,numg,nroots)
      data pt5,one /0.5d+00,1.0d+00/
c
      do 40 nr=1,nroots
      do 30 ng=1,numg
      aa =abv(1,ng)
      bb =abv(2,ng)
      rho=abv(3,ng)
      qab=abv(4,ng)
      uu =rho*rwv(1,ng,nr)
      ww =    rwv(2,ng,nr)
      f00(ng,nr)=ww*qab
      if (doscreen) then
         co2 = cam_omega * cam_omega
         rho0 = abv(6,ng)
         f00(ng,nr) = f00(ng,nr)*sqrt(co2/(rho0+co2))
         t2=uu/(uu+rho)
         t2a=t2*rho/aa
         t2b=t2*rho/bb
         b00(ng,nr)=pt5/(aa+bb)*t2*rho/rho0
         b01(ng,nr)=pt5/bb*(one-t2b)
         b10(ng,nr)=pt5/aa*(one-t2a)
         ainv=one/aa
         binv=one/bb
         if(mmaxs) go to 10
         d00(ng,nr,1)=cv( 2,ng)*binv+t2b*(cv( 1,ng)*ainv-cv( 2,ng)*binv)
         d00(ng,nr,2)=cv( 4,ng)*binv+t2b*(cv( 3,ng)*ainv-cv( 4,ng)*binv)
         d00(ng,nr,3)=cv( 6,ng)*binv+t2b*(cv( 5,ng)*ainv-cv( 6,ng)*binv)
   10    if(nmaxs) go to 20
         c00(ng,nr,1)=cv( 8,ng)*binv-t2a*(cv( 7,ng)*ainv-cv( 8,ng)*binv)
         c00(ng,nr,2)=cv(10,ng)*binv-t2a*(cv( 9,ng)*ainv-cv(10,ng)*binv)
         c00(ng,nr,3)=cv(12,ng)*binv-t2a*(cv(11,ng)*ainv-cv(12,ng)*binv)
   20    continue
      else
         aauu=aa+uu
         bbuu=bb+uu
         dum2=one/(aa*bb+uu*(aa+bb))
         audum=aauu*dum2
         budum=bbuu*dum2
         b00(ng,nr)= udum
         b01(ng,nr)=audum
         b10(ng,nr)=budum
         if(mmaxs) go to 11
         d00(ng,nr,1)= udum*cv( 1,ng) + audum*cv( 2,ng)
         d00(ng,nr,2)= udum*cv( 3,ng) + audum*cv( 4,ng)
         d00(ng,nr,3)= udum*cv( 5,ng) + audum*cv( 6,ng)
   11    if(nmaxs) go to 21
         c00(ng,nr,1)= udum*cv( 8,ng) + budum*cv( 7,ng)
         c00(ng,nr,2)= udum*cv(10,ng) + budum*cv( 9,ng)
         c00(ng,nr,3)= udum*cv(12,ng) + budum*cv(11,ng)
   21    continue
      endif

   30 continue
   40 continue
c
      return
      end
