      subroutine hnd_jkbcdf(b00,b01,b10,c00,d00,f00,dij,dkl,
     1                  abv,cv,rwv,numg)
c
c $Id: hnd_jkbcdf.F,v 1.1 2004-07-21 16:24:35 bert Exp $
c
      implicit double precision (a-h,o-z)
c
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_rys.fh"
c
      dimension b00(numg,nroots,3),b01(numg,nroots,3),b10(numg,nroots,3)
      dimension c00(numg,nroots,3)
      dimension d00(numg,nroots,3)
      dimension f00(numg,nroots,3),dij(numg,nroots,3),dkl(numg,nroots,3)
      dimension abv(5,numg),cv(18,numg)
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
      aauu=aa+uu
      bbuu=bb+uu
      f00(ng,nr,1)=ww*qab
      f00(ng,nr,2)=one
      f00(ng,nr,3)=one
      dum2=pt5/(aa*bb+uu*(aa+bb))
      audum=aauu*dum2
      budum=bbuu*dum2
       udum=  uu*dum2
      b00(ng,nr,1)= udum
      b00(ng,nr,2)= udum
      b00(ng,nr,3)= udum
      b01(ng,nr,1)=audum
      b01(ng,nr,2)=audum
      b01(ng,nr,3)=audum
      b10(ng,nr,1)=budum
      b10(ng,nr,2)=budum
      b10(ng,nr,3)=budum
       udum= udum+ udum
      if(mmaxs) go to 10
      audum=audum+audum
      d00(ng,nr,1)= udum*cv( 1,ng) + audum*cv( 2,ng)
      d00(ng,nr,2)= udum*cv( 3,ng) + audum*cv( 4,ng)
      d00(ng,nr,3)= udum*cv( 5,ng) + audum*cv( 6,ng)
   10 if(nmaxs) go to 20
      budum=budum+budum
      c00(ng,nr,1)= udum*cv( 8,ng) + budum*cv( 7,ng)
      c00(ng,nr,2)= udum*cv(10,ng) + budum*cv( 9,ng)
      c00(ng,nr,3)= udum*cv(12,ng) + budum*cv(11,ng)
   20 continue
c
   30 continue
   40 continue
c
      do 60 nr=1,nroots
      do 50 ng=1,numg
      dij(ng,nr,1)=cv(13,ng)
      dij(ng,nr,2)=cv(14,ng)
      dij(ng,nr,3)=cv(15,ng)
      dkl(ng,nr,1)=cv(16,ng)
      dkl(ng,nr,2)=cv(17,ng)
      dkl(ng,nr,3)=cv(18,ng)
   50 continue
   60 continue
c
      return
      end
