      subroutine esp_rtdbin
c
c $Id: esp_rtdbin.F,v 1.2 1997-10-07 18:09:56 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
      logical int_normalize
      external int_normalize
c
c     get input from rtdb
c
c
c     load geometry
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('esp: geom_create failed',me)
c
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('esp: geom_rtdb_load failed',me)
c
      if(.not.geom_ncent(igeom,natoms))
     + call errquit('esp: geom_ncent failed',me)
c
      if(.not.geom_nuc_charge(igeom,chgnuc))
     + call errquit('esp: geom_nuc_charge failed',me)
c
c     load the basis set
c
      if(.not.bas_create(ibasis,'ao basis'))
     + call errquit('esp: bas_create failed',me)
c
      if(.not.bas_rtdb_load(irtdb,igeom,ibasis,'ao basis'))
     + call errquit('esp: bas_rtdb_load failed',me)
c
      if(.not.int_normalize(ibasis))
     + call errquit('esp: int_normalize failed',me)
c
      if(.not.bas_numbf(ibasis,nbf))
     + call errquit('esp: bas_numbf failed',me)
c     
c     get cutoff radius and gridspacing
c
      if(.not.rtdb_get(irtdb,'esp:rcut',mt_dbl,1,rcut)) rcut=0.5d0
      if(.not.rtdb_get(irtdb,'esp:spac',mt_dbl,1,spac)) spac=0.1d0
      if(.not.rtdb_get(irtdb,'esp:radius',mt_dbl,107,radius)) then
      radius(  1)=0.145d0
      radius(  2)=0.200d0
      radius(  3)=0.200d0
      radius(  4)=0.200d0
      radius(  5)=0.200d0
      radius(  6)=0.170d0
      radius(  7)=0.180d0
      radius(  8)=0.190d0
      radius(  9)=0.200d0
      radius( 10)=0.200d0
      radius( 11)=0.200d0
      radius( 12)=0.200d0
      radius( 13)=0.200d0
      radius( 14)=0.200d0
      radius( 15)=0.200d0
      radius( 16)=0.200d0
      radius( 17)=0.200d0
      radius( 18)=0.200d0
      radius( 19)=0.200d0
      radius( 20)=0.200d0
      radius( 21)=0.200d0
      radius( 22)=0.200d0
      radius( 23)=0.200d0
      radius( 24)=0.200d0
      radius( 25)=0.200d0
      radius( 26)=0.200d0
      radius( 27)=0.200d0
      radius( 28)=0.200d0
      radius( 29)=0.200d0
      radius( 30)=0.200d0
      radius( 31)=0.200d0
      radius( 32)=0.200d0
      radius( 33)=0.200d0
      radius( 34)=0.200d0
      radius( 35)=0.200d0
      radius( 36)=0.200d0
      radius( 37)=0.200d0
      radius( 38)=0.200d0
      radius( 39)=0.200d0
      radius( 40)=0.200d0
      radius( 41)=0.200d0
      radius( 42)=0.200d0
      radius( 43)=0.200d0
      radius( 44)=0.200d0
      radius( 45)=0.200d0
      radius( 46)=0.200d0
      radius( 47)=0.200d0
      radius( 48)=0.200d0
      radius( 49)=0.200d0
      radius( 50)=0.200d0
      radius( 51)=0.200d0
      radius( 52)=0.200d0
      radius( 53)=0.200d0
      radius( 54)=0.200d0
      radius( 55)=0.200d0
      radius( 56)=0.200d0
      radius( 57)=0.200d0
      radius( 58)=0.200d0
      radius( 59)=0.200d0
      radius( 60)=0.200d0
      radius( 61)=0.200d0
      radius( 62)=0.200d0
      radius( 63)=0.200d0
      radius( 64)=0.200d0
      radius( 65)=0.200d0
      radius( 66)=0.200d0
      radius( 67)=0.200d0
      radius( 68)=0.200d0
      radius( 69)=0.200d0
      radius( 70)=0.200d0
      radius( 71)=0.200d0
      radius( 72)=0.200d0
      radius( 73)=0.200d0
      radius( 74)=0.200d0
      radius( 75)=0.200d0
      radius( 76)=0.200d0
      radius( 77)=0.200d0
      radius( 78)=0.200d0
      radius( 79)=0.200d0
      radius( 80)=0.200d0
      radius( 81)=0.200d0
      radius( 82)=0.200d0
      radius( 83)=0.200d0
      radius( 84)=0.200d0
      radius( 85)=0.200d0
      radius( 86)=0.200d0
      radius( 87)=0.200d0
      radius( 88)=0.200d0
      radius( 89)=0.200d0
      radius( 90)=0.200d0
      radius( 91)=0.200d0
      radius( 92)=0.200d0
      radius( 93)=0.200d0
      radius( 94)=0.200d0
      radius( 95)=0.200d0
      radius( 96)=0.200d0
      radius( 97)=0.200d0
      radius( 98)=0.200d0
      radius( 99)=0.200d0
      radius(100)=0.200d0
      radius(101)=0.200d0
      radius(102)=0.200d0
      radius(103)=0.200d0
      radius(104)=0.200d0
      radius(105)=0.200d0
      radius(106)=0.200d0
      radius(107)=0.200d0
      endif
      if(.not.rtdb_get(irtdb,'esp:nconst',mt_int,1,nconst)) nconst=0
      if(nconst.gt.0) then
      if(.not.rtdb_get(irtdb,'esp:ncon',mt_int,1,ncon))
     + call errquit('Error reading ncon from rtdb',ncon)
      if(.not.rtdb_get(irtdb,'esp:icon',mt_int,ncon,icon))
     + call errquit('Error reading icon from rtdb',ncon)
      if(.not.rtdb_get(irtdb,'esp:rcon',mt_dbl,nconst,rcon))
     + call errquit('Error reading rcon from rtdb',nconst)
      endif
      if(.not.rtdb_get(irtdb,'esp:irest',mt_int,1,irest)) irest=0
      if(irest.gt.0) then
      if(.not.rtdb_get(irtdb,'esp:resa',mt_dbl,1,resa))
     + call errquit('Error reading resa from rtdb',irest)
      if(.not.rtdb_get(irtdb,'esp:hfree',mt_int,1,ihfree))
     + call errquit('Error reading ihfree from rtdb',irest)
      if(irest.gt.1) then
      if(.not.rtdb_get(irtdb,'esp:resb',mt_dbl,1,resb))
     + call errquit('Error reading resb from rtdb',irest)
      if(.not.rtdb_get(irtdb,'esp:maxit',mt_int,1,maxit))
     + call errquit('Error reading maxit from rtdb',irest)
      if(.not.rtdb_get(irtdb,'esp:toler',mt_dbl,1,toler))
     + call errquit('Error reading toler from rtdb',irest)
      endif
      endif
      if(.not.rtdb_get(irtdb,'esp:dump',mt_int,1,idump)) idump=0
c
c     convert rcut and spac from nm to au
c
      rcut=cnm2au*rcut
      spac=cnm2au*spac
c
      nconst=nconst+1
c
      return
      end
