      subroutine argos_prop_proper(mdstep,stime,eww,esw,
     + ess,fss,esk,epme,
     + epot,epotw,epotsw,epots,vol,dwr,dsr,ekin,etot,
     + npolit,gsm,esa,box,xsm)
c
      implicit none
c
#include "argos_prop.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      external timer_wall,timer_wall_total
      real*8 timer_wall,timer_wall_total
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme(3)
      real*8 fss(msf,msf,3,2)
      real*8 epot,epotw,epotsw,epots,vol,ekin,etot,epmec,ubias
      real*8 eq,el,ep,ek,et,ewp,box(3)
      real*8 tempw,temps,stime,esk(msf),ewk,ep2,ep3,edrs,epmf
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits,npolit
      real*8 dwr,dsr(msm),gsm(msm,4,2),esa(nsa,2)
      real*8 prest(3,3),virt(3,3),ekct(3,3),ep2m,ep3m,epmec2,epmec3
      real*8 xsm(msm,3),dx(3),upmf(100)
c
      integer i,j,k,it,lenp
c
      call timer_start(54)
c
      if(ltwin) then
      do 1 i=1,mpe
      eww(i,1)=eww(i,1)+eww(i,2)
    1 continue
      do 2 i=1,msf
      do 3 j=1,8
      esw(i,j,1)=esw(i,j,1)+esw(i,j,2)
    3 continue
    2 continue
      do 4 i=1,msf
      do 5 j=1,msf
      do 6 k=1,8
      ess(i,j,k,1)=ess(i,j,k,1)+ess(i,j,k,2)
    6 continue
      fss(i,j,1,1)=fss(i,j,1,1)+fss(i,j,1,2)
      fss(i,j,2,1)=fss(i,j,2,1)+fss(i,j,2,2)
      fss(i,j,3,1)=fss(i,j,3,1)+fss(i,j,3,2)
    5 continue
    4 continue
      endif
c
      el=eww(5,1)+eww(7,1)
      eq=eww(6,1)+eww(8,1)
      ep=zero
      ewp=zero
      do 9 i=1,mpe
      ep=ep+eww(i,1)
      ewp=ewp+eww(i,1)
      do 10 j=1,msf
      ep=ep+esw(j,i,1)
      do 11 k=1,msf
      ep=ep+ess(k,j,i,1)
   11 continue
   10 continue
    9 continue
      epmec=eww(9,1)
      epmec2=eww(9,1)
      epmec3=eww(9,1)
      do 12 j=1,msf
      el=el+esw(j,5,1)
      eq=eq+esw(j,6,1)
      do 13 k=1,msf
      el=el+ess(k,j,5,1)+ess(k,j,7,1)
      eq=eq+ess(k,j,6,1)+ess(k,j,8,1)
      epmec=epmec+ess(k,j,9,1)
      epmec2=epmec2+ess(k,j,10,1)
      epmec3=epmec3+ess(k,j,11,1)
   13 continue
   12 continue
c
c     if using cafe get properties from it
c
      if(lcafe) then
      call argos_cafe_proper(volume,temp,tempw,temps,pres,tmpscl,tmsscl,
     + prsscl,
     + ewk,nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits,
     + ep2,ep3,ep2m,ep3m,edrs,epmf,virial,prest,virt,ekct,ubias,upmf)
      endif
c
c     if using space get properties from it
c
      if(lspac) then
      endif
c
      if(me.eq.0) ep=ep+edrs
      ep=ep+epmf
c
      if(lpme) then
      ep=ep+epme(iset)
      eq=eq+epme(iset)
      ep2=ep2+epme(2)-epme(1)+epmec2-epmec
      ep3=ep3+epme(3)-epme(1)+epmec3-epmec
      endif
c
      ek=ewk
      do 111 i=1,msf
      ek=ek+esk(i)
  111 continue
      et=ep+ek
c
c     fill the property vector
c
      do 7 i=1,maxpro
      p(i)=zero
    7 continue
c
      if(me.eq.0) then
c
      if(ntype.eq.0) call argos_space_gettp(temp,pres)
c
      p(1)=dble(mdstep)
      p(2)=stime
      endif
      p(3)=dble(nwws)
      p(4)=dble(nwwl)
      p(5)=dble(nsws)
      p(6)=dble(nswl)
      p(7)=dble(nsss)
      p(8)=dble(nssl)
      p(9)=dble(nshitw)
      p(10)=dble(nshits)
      p(11)=dble(npolit)
      p(33)=zero
      p(34)=zero
      if(me.eq.0) then
      p(12)=volume
      p(13)=1.6605655d0*wbox/volume
      p(14)=pres
      p(15)=prsscl
      p(16)=temp
      p(17)=tempw
      p(18)=temps
      p(21)=tmpscl
      p(22)=tmsscl
      p(27)=ek
      p(32)=ek
      p(33)=virial
      p(34)=pres*volume
      p(64)=ewk
      p(67)=ewk
      endif
      p(24)=eq
      p(25)=el
      p(26)=ep
      p(32)=p(32)+ep
      p(36)=ep2
      p(37)=ep3
      p(38)=ep2+ep2m
      p(39)=ep3+ep3m
      if(me.eq.0) p(40)=dfree
      p(50)=ubias
      p(52)=eww(6,1)
      p(53)=eww(5,1)
      p(54)=eww(8,1)
      p(55)=eww(7,1)
      p(56)=eww(1,1)
      p(57)=eww(2,1)
      p(58)=eww(13,1)
      p(59)=eww(3,1)
      p(60)=eww(4,1)
      p(67)=ewp
      p(68)=p(68)+ewp
      p(70)=dwr
      if(lnoe) p(75)=edrs
      if(lpmf) p(76)=epmf
      if(ntype.ne.3) then
      p(77)=epme(iset)
      p(80)=epmec
      else
      p(77)=epme(1)
      p(78)=epme(2)-p(77)
      p(79)=epme(3)-p(77)
      p(80)=epme(1)+epmec
      p(81)=epme(2)+epmec2-p(80)
      p(82)=epme(3)+epmec3-p(80)
      endif
      p(84)=timer_wall(202)
      p(85)=timer_wall(203)
      if(me.eq.0) p(86)=timer_wall(203)
      if(me.eq.0) then
      p(87)=virt(1,1)
      p(88)=virt(1,2)
      p(89)=virt(1,3)
      p(90)=virt(2,1)
      p(91)=virt(2,2)
      p(92)=virt(2,3)
      p(93)=virt(3,1)
      p(94)=virt(3,2)
      p(95)=virt(3,3)
      p(96)=prest(1,1)
      p(97)=prest(1,2)
      p(98)=prest(1,3)
      p(99)=prest(2,1)
      p(100)=prest(2,2)
      p(101)=prest(2,3)
      p(102)=prest(3,1)
      p(103)=prest(3,2)
      p(104)=prest(3,3)
      p(105)=ekct(1,1)
      p(106)=ekct(1,2)
      p(107)=ekct(1,3)
      p(108)=ekct(2,1)
      p(109)=ekct(2,2)
      p(110)=ekct(2,3)
      p(111)=ekct(3,1)
      p(112)=ekct(3,2)
      p(113)=ekct(3,3)
      p(114)=box(1)
      p(115)=box(2)
      p(116)=box(3)
      endif
c
      maxp=isprop
      it=isprop
c
      if(nsf.gt.0) then
      do 33 i=1,nsf
      it=isprop+(i-1)*30
      p(it+2)=ess(i,i,6,1)
      p(it+3)=ess(i,i,5,1)
      p(it+5)=esw(i,6,1)
      p(it+6)=esw(i,5,1)
      p(it+7)=ess(i,i,1,1)
      p(it+8)=ess(i,i,2,1)
      p(it+9)=ess(i,i,13,1)
      p(it+10)=ess(i,i,3,1)
      p(it+11)=ess(i,i,4,1)
      p(it+14)=ess(i,i,7,1)
      p(it+15)=ess(i,i,8,1)
      p(it+16)=ess(i,i,5,1)
      p(it+17)=ess(i,i,6,1)
      if(me.eq.0) then
      p(it+19)=esk(i)
      endif
      p(it+20)=ess(i,i,5,1)+ess(i,i,6,1)+ess(i,i,7,1)+ess(i,i,8,1)+
     + ess(i,i,1,1)+ess(i,i,2,1)+
     + ess(i,i,3,1)+ess(i,i,4,1)+half*(esw(i,6,1)+esw(i,5,1))
      do 34 j=i+1,nsf
      p(it+20)=p(it+18)+half*(ess(i,j,6,1)+ess(j,i,6,1)+
     + ess(i,j,5,1)+ess(j,i,5,1))
   34 continue
   33 continue
      it=it+30
      if(nsf.gt.1) then
      do 35 i=1,nsf-1
      do 36 j=i+1,nsf
      it=it+1
      it=it+1
      p(it)=ess(i,j,6,1)+ess(j,i,6,1)
      it=it+1
      p(it)=ess(i,j,5,1)+ess(j,i,5,1)
      it=it+1
      it=it+1
      p(it)=fss(i,j,1,1)-fss(j,i,1,1)
      it=it+1
      p(it)=fss(i,j,2,1)-fss(j,i,2,1)
      it=it+1
      p(it)=fss(i,j,3,1)-fss(j,i,3,1)
      it=it+1
      p(it)=sqrt((fss(i,j,1,1)-fss(j,i,1,1))**2+
     + (fss(i,j,2,1)-fss(j,i,2,1))**2+(fss(i,j,3,1)-fss(j,i,3,1))**2)
   36 continue
   35 continue
      endif
      maxp=it
c
      if(nsm.gt.0.and.npstat.gt.0) then
      it=maxp
      do 37 i=1,nsm
      it=it+1
      if(me.eq.0) p(it)=gsm(i,iset,1)
      if(p(it).lt.tiny) lp(it)=.false.
   37 continue
      do 55 i=1,nsm
      it=it+1
      if(me.eq.0) p(it)=gsm(i,4,1)
      if(p(it).lt.tiny) lp(it)=.false.
   55 continue
      do 155 i=1,nsm-1
      do 156 j=i+1,nsm
      dx(1)=xsm(i,1)-xsm(j,1)
      dx(2)=xsm(i,2)-xsm(j,2)
      dx(3)=xsm(i,3)-xsm(j,3)
      if(lpbc) call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      it=it+1
      if(me.eq.0) p(it)=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
  156 continue
  155 continue
      endif
c
      endif
c
      if(npmfi.gt.0) then
      do 356 i=1,npmfi
      it=it+1
      p(it)=upmf(i)
  356 continue
      endif
c
      call timer_stop(54)
c
      if(iprof.eq.1) then
      do 56 i=1,55
      it=it+1
      p(it)=timer_wall_total(i)
   56 continue
      endif
c
      maxp=it
c
      if(npener.gt.0) call argos_cafe_add_esa(esa)
      if(np.gt.0) then
      call ga_dgop(mrg_d44,p(9),3,'max')
      if(me.ne.0) then
      p(9)=zero
      p(10)=zero
      p(11)=zero
      endif
      if(maxp+mpe*(1+msf*(1+msf)).gt.maxpro) then
      call ga_dgop(mrg_d45,p,maxp,'+')
      call ga_dgop(mrg_d40,eww,mpe,'+')
      call ga_dgop(mrg_d41,esw,msf*mpe,'+')
      call ga_dgop(mrg_d42,ess,msf*msf*mpe,'+')
      else
      lenp=maxp
      do 501 i=1,mpe
      p(lenp+i)=eww(i,1)
  501 continue
      lenp=lenp+mpe
      do 502 i=1,mpe
      do 503 j=1,msf
      lenp=lenp+1
      p(lenp)=esw(j,i,1)
  503 continue
  502 continue
      do 504 i=1,mpe
      do 505 j=1,msf
      do 506 k=1,msf
      lenp=lenp+1
      p(lenp)=ess(k,j,i,1)
  506 continue
  505 continue
  504 continue
      call ga_dgop(mrg_d45,p,lenp,'+')
      lenp=maxp
      do 511 i=1,mpe
      eww(i,1)=p(lenp+i)
  511 continue
      lenp=lenp+mpe
      do 512 i=1,mpe
      do 513 j=1,msf
      lenp=lenp+1
      esw(j,i,1)=p(lenp)
  513 continue
  512 continue
      do 514 i=1,mpe
      do 515 j=1,msf
      do 516 k=1,msf
      lenp=lenp+1
      ess(k,j,i,1)=p(lenp)
  516 continue
  515 continue
  514 continue
      endif
cxxxxx      call ga_dgop(mrg_d43,epme,3,'+')
      if(npener.gt.0) call ga_dgop(mrg_d47,esa,2*nsa,'+')
      endif
c
      if(p(85).ne.0.0d0) then
      p(85)=(p(85)-p(84))/p(85)
      endif
c
      if(temp.eq.zero.or.me.ne.0) then
      p(117)=zero
      p(118)=zero
      p(119)=zero
      p(120)=zero
      p(51)=zero
      else
      p(117)=exp(-p(36)/(rgas*temp))
      p(118)=exp(-p(37)/(rgas*temp))
      p(119)=exp(-p(38)/(rgas*temp))
      p(120)=exp(-p(39)/(rgas*temp))
      p(51)=exp(-p(50)/(rgas*temp))
      if(nbias.gt.0) p(41)=p(40)*p(51)
      endif
c
      if(nwm.gt.0) then
      dwr=p(70)/dble(nwm)
      p(70)=dwr
      if(stime.ne.0) then
      p(71)=1.0d-6*dwr/(6.0d0*stime)
      else
      p(71)=zero
      endif
      endif
c
      epot=p(26)
      epots=zero
      epotsw=zero
      epotw=eww(5,1)+eww(6,1)+eww(8,1)
      do 40 i=1,msf
      epotsw=epotsw+esw(i,5,1)+esw(i,6,1)+esw(i,8,1)
      do 41 j=1,msf
      epots=epots+ess(i,j,1,1)+ess(i,j,2,1)+ess(i,j,3,1)+ess(i,j,4,1)
      epots=epots+ess(i,j,5,1)+ess(i,j,6,1)+ess(i,j,8,1)
   41 continue
   40 continue
c
      epots=zero
      epotsw=zero
      epotw=zero
      do 42 i=1,mpe
      epotw=epotw+eww(i,1)
      do 38 j=1,msf
      epotsw=epotsw+esw(j,i,1)
      do 39 k=1,msf
      epots=epots+ess(k,j,i,1)
   39 continue
   38 continue
   42 continue
c
      vol=volume
      ekin=p(27)
      etot=ekin+epot
c
      return
      end
c $Id$
