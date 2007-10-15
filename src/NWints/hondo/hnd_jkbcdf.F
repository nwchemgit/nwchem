      subroutine hnd_jkbcdf(b00,b01,b10,c00,d00,f00,
     1                  abv,cv,rwv,numg)
c
c $Id: hnd_jkbcdf.F,v 1.2 2007-10-15 23:37:24 bert Exp $
c
      implicit double precision (a-h,o-z)
c
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_rys.fh"
c
      dimension b00(numg,nroots),b01(numg,nroots),b10(numg,nroots)
      dimension c00(numg,nroots,3)
      dimension d00(numg,nroots,3)
      dimension f00(numg,nroots)
      dimension abv(5,numg),cv(12,numg)
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
      f00(ng,nr)=ww*qab
      dum2=pt5/(aa*bb+uu*(aa+bb))
      audum=aauu*dum2
      budum=bbuu*dum2
       udum=  uu*dum2
      b00(ng,nr)= udum
      b01(ng,nr)=audum
      b10(ng,nr)=budum
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
      return
      end
