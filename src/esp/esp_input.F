      subroutine esp_input(irtdb)
c
c $Id: esp_input.F,v 1.1 1997-10-06 23:00:34 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item
      integer i,me,np,icon(2500),ncon,nconst,jcon,irest,maxit,ihfree
      real*8 rcut,spac,radius(107),rcon(500),resa,resb
      real*8 toler
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      if(me.eq.0) then
      write(6,*)
      write(6,*)
      call util_print_centered(6,'ESP Input Module',40,.true.)
      write(6,*)
      write(6,*)
      endif
c
c     set the defaults
c
      ncon=0
      nconst=0
      irest=0
      ihfree=0
      toler=1.0d-4
      maxit=25
      resa=0.0d0
      resb=0.0d0
      rcut=0.5d0
      spac=0.1d0
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
c
c     read the nwchem input file
c
    1 continue
c
      if(.not.inp_read()) call errquit('esp_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
      if(inp_compare(.false.,'range',item)) then
      if(.not.inp_f(rcut)) call errquit('esp_input: range error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'spacing',item)) then
      if(.not.inp_f(spac)) call errquit('esp_input: spacing error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'radius',item)) then
      if(.not.inp_i(i)) call errquit('esp_input: radius error',0)
      if(.not.inp_f(radius(i)))
     + call errquit('esp_input: radius error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'constrain',item)) then
      nconst=nconst+1
      if(.not.inp_f(rcon(nconst))) 
     + call errquit('esp_input: constrain error',0)
      ncon=ncon+1
      jcon=ncon
      icon(jcon)=0
    4 continue
      if(inp_i(i)) then
      icon(jcon)=icon(jcon)+1
      ncon=ncon+1
      icon(ncon)=i
      goto 4
      endif
      goto 2
      endif
c
      if(inp_compare(.false.,'restrain',item)) then
      irest=2
      resa=0.001d0
      resb=0.1d0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'hfree',item)) then
      ihfree=1
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'harmonic',item)) then
      irest=1
      resb=0.0d0
      if(.not.inp_f(resa)) resa=0.005d0
      goto 2
      endif
      if(inp_compare(.false.,'hyperbolic',item)) then
      irest=2
      if(.not.inp_f(resa)) resa=0.001d0
      if(.not.inp_f(resb)) resb=0.1d0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'maxiter',item)) then
      if(.not.inp_i(maxit)) call errquit('esp_input: maxit error',0)
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'tolerance',item)) then
      if(.not.inp_f(toler)) call errquit('esp_input: toler error',0)
      if(.not.inp_a(item)) goto 1
      endif
      endif
      goto 3
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      if(me.eq.0) write(6,9999) item
 9999 format('arg_input: unknown field ',a)
      if(np.gt.1) call ga_sync()
      call errquit('arg_input: error reading input file',0)
      endif
c
      if(.not.rtdb_put(irtdb,'esp:rcut',mt_dbl,1,rcut))
     + call errquit('esp: rtdb_put rcut failed',0)
      if(.not.rtdb_put(irtdb,'esp:spac',mt_dbl,1,spac))
     + call errquit('esp: rtdb_put spac failed',0)
      if(.not.rtdb_put(irtdb,'esp:radius',mt_dbl,107,radius))
     + call errquit('esp: rtdb_put radius failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:nconst',mt_int,1,nconst))
     + call errquit('esp: rtdb_put nconst failed',0)
      if(nconst.gt.0) then
      if(.not.rtdb_put(irtdb,'esp:ncon',mt_int,1,ncon))
     + call errquit('esp: rtdb_put ncon failed',0)
      if(.not.rtdb_put(irtdb,'esp:icon',mt_int,ncon,icon))
     + call errquit('esp: rtdb_put icon failed',0)
      if(.not.rtdb_put(irtdb,'esp:rcon',mt_dbl,nconst,rcon))
     + call errquit('esp: rtdb_put rcon failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'esp:irest',mt_int,1,irest))
     + call errquit('esp: rtdb_put irest failed',0)
      if(.not.rtdb_put(irtdb,'esp:resa',mt_dbl,1,resa))
     + call errquit('esp: rtdb_put resa failed',0)
      if(.not.rtdb_put(irtdb,'esp:resb',mt_dbl,1,resb))
     + call errquit('esp: rtdb_put resb failed',0)
      if(.not.rtdb_put(irtdb,'esp:hfree',mt_int,1,ihfree))
     + call errquit('esp: rtdb_put hfree failed',0)
      if(.not.rtdb_put(irtdb,'esp:maxit',mt_int,1,maxit))
     + call errquit('esp: rtdb_put maxit failed',0)
      if(.not.rtdb_put(irtdb,'esp:toler',mt_dbl,1,toler))
     + call errquit('esp: rtdb_put toler failed',0)
c
      return
      end
