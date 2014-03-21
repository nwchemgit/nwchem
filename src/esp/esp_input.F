      subroutine esp_input(irtdb)
c
c $Id$
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer maxcon
      parameter(maxcon=5000)
c
      integer irtdb
      character*255 item
      integer i,me,np,icon(maxcon),ncon,nconst,jcon,irest,maxit,ihfree
      integer idump,ifrom,j,lscrn,newgrd,k,l,m,denmat
      real*8 rcut,spac,radius(107),rcon(500),resa,resb,resc,resd,qh
      real*8 toler,factor,probe,scrtol
      logical lstage2
c
      me=ga_nodeid()
      np=ga_nnodes()
c
c      if(me.eq.0) then
c      write(6,*)
c      write(6,*)
c      call util_print_centered(6,'ESP Input Module',40,.true.)
c      write(6,*)
c      write(6,*)
c      endif
c
c     set the defaults
c
      ncon=0
      nconst=0
      irest=0
      ihfree=0
      idump=0
      lscrn=0
      newgrd=0
      denmat=0
      scrtol=1.0d-9
      toler=1.0d-4
      maxit=25
      resa=0.0d0
      resb=0.0d0
      resc=0.0d0
      resd=0.0d0
      rcut=0.3d0
      spac=0.05d0
      qh=0.0d0
      radius(  1)=0.030d0
      radius(  2)=0.122d0
      radius(  3)=0.123d0
      radius(  4)=0.089d0
      radius(  5)=0.088d0
      radius(  6)=0.077d0
      radius(  7)=0.070d0
      radius(  8)=0.066d0
      radius(  9)=0.058d0
      radius( 10)=0.160d0
      radius( 11)=0.140d0
      radius( 12)=0.136d0
      radius( 13)=0.125d0
      radius( 14)=0.117d0
      radius( 15)=0.110d0
      radius( 16)=0.104d0
      radius( 17)=0.099d0
      radius( 18)=0.191d0
      radius( 19)=0.203d0
      radius( 20)=0.174d0
      radius( 21)=0.144d0
      radius( 22)=0.132d0
      radius( 23)=0.122d0
      radius( 24)=0.119d0
      radius( 25)=0.117d0
      radius( 26)=0.117d0
      radius( 27)=0.116d0
      radius( 28)=0.115d0
      radius( 29)=0.117d0
      radius( 30)=0.125d0
      radius( 31)=0.125d0
      radius( 32)=0.122d0
      radius( 33)=0.121d0
      radius( 34)=0.117d0
      radius( 35)=0.114d0
      radius( 36)=0.198d0
      radius( 37)=0.222d0
      radius( 38)=0.192d0
      radius( 39)=0.162d0
      radius( 40)=0.145d0
      radius( 41)=0.134d0
      radius( 42)=0.129d0
      radius( 43)=0.127d0
      radius( 44)=0.124d0
      radius( 45)=0.125d0
      radius( 46)=0.128d0
      radius( 47)=0.134d0
      radius( 48)=0.141d0
      radius( 49)=0.150d0
      radius( 50)=0.140d0
      radius( 51)=0.141d0
      radius( 52)=0.137d0
      radius( 53)=0.133d0
      radius( 54)=0.209d0
      radius( 55)=0.235d0
      radius( 56)=0.198d0
      radius( 57)=0.169d0
      radius( 58)=0.165d0
      radius( 59)=0.165d0
      radius( 60)=0.164d0
      radius( 61)=0.165d0
      radius( 62)=0.166d0
      radius( 63)=0.165d0
      radius( 64)=0.161d0
      radius( 65)=0.159d0
      radius( 66)=0.159d0
      radius( 67)=0.158d0
      radius( 68)=0.157d0
      radius( 69)=0.156d0
      radius( 70)=0.156d0
      radius( 71)=0.156d0
      radius( 72)=0.144d0
      radius( 73)=0.134d0
      radius( 74)=0.130d0
      radius( 75)=0.128d0
      radius( 76)=0.126d0
      radius( 77)=0.126d0
      radius( 78)=0.129d0
      radius( 79)=0.134d0
      radius( 80)=0.144d0
      radius( 81)=0.155d0
      radius( 82)=0.154d0
      radius( 83)=0.152d0
      radius( 84)=0.153d0
      radius( 85)=0.150d0
      radius( 86)=0.220d0
      radius( 87)=0.324d0
      radius( 88)=0.268d0
      radius( 89)=0.225d0
      radius( 90)=0.216d0
      radius( 91)=0.193d0
      radius( 92)=0.166d0
      radius( 93)=0.157d0
      radius( 94)=0.181d0
      radius( 95)=0.221d0
      radius( 96)=0.143d0
      radius( 97)=0.142d0
      radius( 98)=0.140d0
      radius( 99)=0.139d0
      radius(100)=0.138d0
      radius(101)=0.137d0
      radius(102)=0.136d0
      radius(103)=0.134d0
      radius(104)=0.130d0
      radius(105)=0.130d0
      radius(106)=0.130d0
      radius(107)=0.130d0
c
      probe=0.07d0
      factor=1.00d0
c
c     read the nwchem input file
c
    1 continue
c
      if(.not.inp_read()) call md_abort('esp_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
c     print & noprint
c     ---------------
c
      if(inp_compare(.false.,'print',item).or.
     + inp_compare(.false.,'noprint',item)) then
      call util_print_input(irtdb,'esp')
      goto 1
      endif
c
      if(inp_compare(.false.,'range',item)) then
      if(.not.inp_f(rcut)) call md_abort('esp_input: range error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'spacing',item)) then
      if(.not.inp_f(spac)) call md_abort('esp_input: spacing error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'radius',item)) then
      if(.not.inp_i(i)) call md_abort('esp_input: radius error',0)
      if(.not.inp_f(radius(i)))
     + call md_abort('esp_input: radius error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'probe',item)) then
      if(.not.inp_f(probe)) call md_abort('esp_input: probe error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'factor',item)) then
      if(.not.inp_f(factor)) call md_abort('esp_input: factor error',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'screen',item)) then
      lscrn=1
      if(.not.inp_f(scrtol)) scrtol=1.0d-5
      goto 2
      endif
c
      if(inp_compare(.false.,'qh',item)) then
        if(.not.inp_f(qh)) qh=0.0d0
        goto 2
      endif
c
      if(inp_compare(.false.,'recalculate',item)) then
        newgrd=1
        goto 2
      endif
c
      if(inp_compare(.false.,'denmat',item)) then
       denmat=1
       goto 2
      endif
c
      if(inp_compare(.false.,'constrain',item)) then
      ifrom=0
      nconst=nconst+1
      if(.not.inp_f(rcon(nconst))) then
      if(.not.inp_a(item))
     + call md_abort('esp_input: constrain error',0)
      if(inp_compare(.false.,'stage2',item)) then
        if(.not.inp_a(item))
     +   call md_abort('esp_input: constrain error',0)
        lstage2=.true.
      else
        lstage2=.false.
      endif
      if(inp_compare(.false.,'equal',item)) then
        if(.not.inp_i(i))
     +   call md_abort('esp_input: constrain error',0)
        if(.not.inp_i(j))
     +   call md_abort('esp_input: constrain error',0)
   41   continue
        ncon=ncon+1
        if(ncon+2.gt.maxcon)
     +     call md_abort('Increase dimension maxcon',ncon)
        icon(ncon)=2
        if(lstage2) icon(ncon)=-2
        icon(ncon+1)=-i
        icon(ncon+2)=j
        rcon(nconst)=0.0d0
        ncon=ncon+2
        if(ncon+2.gt.maxcon)
     +     call md_abort('Increase dimension maxcon',ncon)
        if(inp_i(j)) then
          nconst=nconst+1
          goto 41
        endif
        goto 2
      elseif(inp_compare(.false.,'group',item)) then
        if(.not.inp_i(i))
     +     call md_abort('esp_input: constrain error',0)
        if(.not.inp_i(j))
     +     call md_abort('esp_input: constrain error',0)
        if(.not.inp_a(item))
     +     call md_abort('esp_input: constrain error',0)
        if(.not.inp_compare(.false.,'to',item))
     +     call md_abort('esp_input: constrain error',0)
        if(.not.inp_i(k))
     +     call md_abort('esp_input: constrain error',0)
        if(.not.inp_i(l))
     +     call md_abort('esp_input: constrain error',0)
        if(j.lt.i.or.l.lt.k.or.j-i.ne.l-k)
     +     call md_abort('esp_input: constrain error',0)
        nconst=nconst-1
        do 42 m=0,j-i
          nconst=nconst+1
          ncon=ncon+1
          if(ncon+2.gt.maxcon)
     +       call md_abort('Increase dimension maxcon',ncon)
          icon(ncon)=2
          icon(ncon+1)=-(i+m)
          icon(ncon+2)=k+m
          rcon(nconst)=0.0d0
          ncon=ncon+2
          if(ncon+2.gt.maxcon)
     +       call md_abort('Increase dimension maxcon',ncon)
   42   continue
        goto 2
      elseif(inp_compare(.false.,'xhn',item)) then
        ncon=ncon+1
        if(ncon+2.gt.maxcon)
     +     call md_abort('Increase dimension maxcon',ncon)
        jcon=ncon
        icon(jcon)=0
        rcon(nconst)=0.0d0
   43   continue
        if(inp_i(i)) then
          icon(jcon)=icon(jcon)+1
          ncon=ncon+1
          if(ncon+2.gt.maxcon)
     +       call md_abort('Increase dimension maxcon',ncon)
          icon(ncon)=i
          goto 43
        elseif(icon(jcon).gt.2) then
          do 44 j=1,icon(jcon)-2
            nconst=nconst+1
            rcon(nconst)=0.0d0
            ncon=ncon+1
            if(ncon+2.gt.maxcon)
     +         call md_abort('Increase dimension maxcon',ncon)
            icon(ncon)=2
            icon(ncon+1)=icon(jcon+j+1)
            icon(ncon+2)=-icon(jcon+j+2)
            ncon=ncon+2
            if(ncon+2.gt.maxcon)
     +         call md_abort('Increase dimension maxcon',ncon)
   44     continue
        endif
        goto 2
      else
        call md_abort('esp_input: constrain error',0)
      endif
      endif
      ncon=ncon+1
      if(ncon+2.gt.maxcon)
     + call md_abort('Increase dimension maxcon',ncon)
      jcon=ncon
      icon(jcon)=0
    4 continue
      if(inp_i(i)) then
      icon(jcon)=icon(jcon)+1
      ncon=ncon+1
      if(ncon+2.gt.maxcon)
     + call md_abort('Increase dimension maxcon',ncon)
      icon(ncon)=i
      ifrom=iabs(i)
      goto 4
      else
      if(.not.inp_a(item)) goto 2
      if(.not.inp_compare(.false.,'through',item))
     + call md_abort('esp_input: constrain error',0)
      if(.not.inp_i(i))
     + call md_abort('esp_input: constrain error',0)
      if(ifrom.eq.0) 
     + call md_abort('esp_input: constrain error',0)
      do 40 j=ifrom+1,iabs(i)
      icon(jcon)=icon(jcon)+1
      ncon=ncon+1
      if(ncon+2.gt.maxcon)
     + call md_abort('Increase dimension maxcon',ncon)
      icon(ncon)=j
   40 continue
      ifrom=0
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
      if(.not.inp_f(resc)) resc=0.010d0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'maxiter',item)) then
      if(.not.inp_i(maxit)) call md_abort('esp_input: maxit error',0)
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'tolerance',item)) then
      if(.not.inp_f(toler)) call md_abort('esp_input: toler error',0)
      if(.not.inp_a(item)) goto 1
      endif
      goto 2
      endif
      if(inp_compare(.false.,'hyperbolic',item)) then
      irest=2
      if(.not.inp_f(resa)) resa=0.001d0
      if(.not.inp_f(resb)) resb=0.1d0
      if(.not.inp_f(resc)) resc=0.002d0
      if(.not.inp_f(resd)) resd=0.1d0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'maxiter',item)) then
      if(.not.inp_i(maxit)) call md_abort('esp_input: maxit error',0)
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'tolerance',item)) then
      if(.not.inp_f(toler)) call md_abort('esp_input: toler error',0)
      if(.not.inp_a(item)) goto 1
      endif
      endif
      goto 3
      endif
c
      if(inp_compare(.false.,'dump',item)) then
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'grid',item)) then
      idump=1
      goto 2
      endif
      goto 3
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      if(me.eq.0) write(6,9999) item
 9999 format('esp_input: unknown field ',a)
      if(np.gt.1) call ga_sync()
      call md_abort('esp_input: error reading input file',0)
      endif
c
      if(factor.le.0.0d0) call md_abort('illegal value for factor',0)
      if(probe.le.0.0d0) call md_abort('illegal value for probe',0)
c
      if(.not.rtdb_put(irtdb,'esp:probe',mt_dbl,1,probe))
     + call md_abort('esp: rtdb_put rcut failed',0)
      if(.not.rtdb_put(irtdb,'esp:factor',mt_dbl,1,factor))
     + call md_abort('esp: rtdb_put rcut failed',0)
      if(.not.rtdb_put(irtdb,'esp:rcut',mt_dbl,1,rcut))
     + call md_abort('esp: rtdb_put rcut failed',0)
      if(.not.rtdb_put(irtdb,'esp:spac',mt_dbl,1,spac))
     + call md_abort('esp: rtdb_put spac failed',0)
      if(.not.rtdb_put(irtdb,'esp:radius',mt_dbl,107,radius))
     + call md_abort('esp: rtdb_put radius failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:nconst',mt_int,1,nconst))
     + call md_abort('esp: rtdb_put nconst failed',0)
      if(nconst.gt.0) then
      if(.not.rtdb_put(irtdb,'esp:ncon',mt_int,1,ncon))
     + call md_abort('esp: rtdb_put ncon failed',0)
      if(.not.rtdb_put(irtdb,'esp:icon',mt_int,ncon,icon))
     + call md_abort('esp: rtdb_put icon failed',0)
      if(.not.rtdb_put(irtdb,'esp:rcon',mt_dbl,nconst,rcon))
     + call md_abort('esp: rtdb_put rcon failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'esp:irest',mt_int,1,irest))
     + call md_abort('esp: rtdb_put irest failed',0)
      if(.not.rtdb_put(irtdb,'esp:resa',mt_dbl,1,resa))
     + call md_abort('esp: rtdb_put resa failed',0)
      if(.not.rtdb_put(irtdb,'esp:resb',mt_dbl,1,resb))
     + call md_abort('esp: rtdb_put resb failed',0)
      if(.not.rtdb_put(irtdb,'esp:resc',mt_dbl,1,resc))
     + call md_abort('esp: rtdb_put resc failed',0)
      if(.not.rtdb_put(irtdb,'esp:resd',mt_dbl,1,resd))
     + call md_abort('esp: rtdb_put resd failed',0)
      if(.not.rtdb_put(irtdb,'esp:hfree',mt_int,1,ihfree))
     + call md_abort('esp: rtdb_put hfree failed',0)
      if(.not.rtdb_put(irtdb,'esp:maxit',mt_int,1,maxit))
     + call md_abort('esp: rtdb_put maxit failed',0)
      if(.not.rtdb_put(irtdb,'esp:toler',mt_dbl,1,toler))
     + call md_abort('esp: rtdb_put toler failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:dump',mt_int,1,idump))
     + call md_abort('esp: rtdb_put dump failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:screen',mt_int,1,lscrn))
     + call md_abort('esp: rtdb_put screen failed',0)
      if(.not.rtdb_put(irtdb,'esp:scrtol',mt_dbl,1,scrtol))
     + call md_abort('esp: rtdb_put scrtol failed',0)
      if(.not.rtdb_put(irtdb,'esp:qh',mt_dbl,1,qh))
     + call md_abort('esp: rtdb_put qh failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:newgrd',mt_int,1,newgrd))
     + call md_abort('esp: rtdb_put newgrd failed',0)
c
      if(.not.rtdb_put(irtdb,'esp:denmat',mt_int,1,denmat))
     + call md_abort('esp: rtdb_put denmat failed',0)
c
      return
      end
