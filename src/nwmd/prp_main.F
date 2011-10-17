      subroutine prp_step(mdstep,stime,eww,esw,ess,fss,esk,
     + epme,uqmd,uqmmm)
c
c $Id$
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical frequency
      external frequency
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme(3)
      real*8 fss(msf,msf,3,2)
      real*8 stime,esk(msf),uqmd,uqmmm
c
      if(me.eq.0) call prp_stat(mdstep,stime,eww,esw,ess,esk,epme(iset))
c
      if(frequency(mdstep,npener)) then
      call cf_print_energy(lfnout)
      endif
c
      if(frequency(mdstep,nfprop)) call prp_record()
c
      return
      end
      subroutine prp_proper(mdstep,stime,eww,esw,ess,fss,esk,epme,
     + uqmd,uqmmm,epot,epotw,epotsw,epots,vol,dwr,dsr,ekin,etot,
     + npolit,gsm,esa,box,xsm)
c
      implicit none
c
#include "prp_common.fh"
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
      real*8 eq,el,ep,ek,et,ewp,uqmd,uqmmm,box(3)
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
      if(me.eq.0) ep=ep+uqmd+uqmmm
c
c     if using cafe get properties from it
c
      if(lcafe) then
      call cf_proper(volume,temp,tempw,temps,pres,tmpscl,tmsscl,prsscl,
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
      if(ntype.eq.0) call sp_gettp(temp,pres)
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
      if(lpbc) call cf_pbc(1,dx,1,dx,1,0,1,1)
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
      if(npener.gt.0) call cf_add_esa(esa)
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
      subroutine prp_stat(mdstep,stime,eww,esw,ess,esk,epme)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical frequency
      external frequency
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme
      real*8 stime,esk(msf)
c
      integer i,j
      character*10 pdate,ptime
      real*8 facs,fact,tfacs,tfact,rt
c
      nsum=nsum+1
      nsumt=nsumt+1
      nsump=nsump+1
      do 8 i=1,maxpro
      if(abs(p(i)).lt.tiny) p(i)=zero
      psum(i)=psum(i)+p(i)
      p2sum(i)=p2sum(i)+p(i)*p(i)
      pslop(i)=pslop(i)+stime*p(i)
      psumt(i)=psumt(i)+p(i)
      p2sumt(i)=p2sumt(i)+p(i)*p(i)
      pslopt(i)=pslopt(i)+stime*p(i)
      psump(i)=psump(i)+p(i)
    8 continue
      tsum=tsum+stime
      t2sum=t2sum+stime*stime
      tsumt=tsumt+stime
      t2sumt=t2sumt+stime*stime
c
      if(frequency(mdstep,nfoutp)) then
      if(.not.lhdr) then
      call swatch(pdate,ptime)
      write(lfnout,1000) pdate,ptime
 1000 format(/,' MOLECULAR DYNAMICS TIME STEP INFORMATION',T110,2A10,//,
     + '    Time      Temp    Pres      Volume Tscalw Tscals  Pscal ',
     + '   U(ele)    U(vdW)      U(pot)     U(kin)     U(tot) ',/,
     + '     ps         K      Pa        nm**3                      ',
     + '   kJ/mol    kJ/mol      kJ/mol     kJ/mol     kJ/mol ',/)
      if(lpstep) write(lfnout,1010)
 1010 format(14X,
     + '   U(bnd)     U(ang)     U(dih)     U(imp)  ',
     + '   Ui(3rd)    Ui(non) ',
     + '   Uw(ele)    Uw(vdW)    U(kin)     U(pot)     U(tot)',/,14X,
     + '   kJ/mol     kJ/mol     kJ/mol     kJ/mol  ',
     + '   kJ/mol     kJ/mol  ',
     + '   kJ/mol     kJ/mol     kJ/mol     kJ/mol     kJ/mol  ',/)
      lhdr=.true.
      endif
      write(lfnout,1001) stime,temp,pres,volume,p(21),p(22),p(15),
     + p(24),p(25),p(26),p(27),p(32)
 1001 format(1x,f10.5,0pf8.2,1pe9.2,0pf10.3,3f7.4,5(1pe11.4),i5,i7)
      if(lpstep) then
      if(nwm.gt.0) then
      rt=one/nwm
      write(lfnout,1011) rt*p(56),rt*p(57),rt*p(58),rt*p(59),rt*p(52),
     + rt*p(53),rt*p(54),rt*p(55),rt*p(64),rt*p(66),rt*p(67)
 1011 format('   solvent  ',11f11.2)
      endif
      do 1 i=1,nsf
      j=isprop+(i-1)*30
      rt=p(j+7)+p(j+8)+p(j+9)+p(j+10)+p(j+2)+p(j+3)+p(j+5)+p(j+6)
      write(lfnout,1012) i,p(j+7),p(j+8),p(j+9),p(j+10),
     + p(j+14)+p(j+15),p(j+2)+p(j+3),p(j+5),p(j+6),p(j+19),rt,rt+p(j+19)
 1012 format('   solute',i3,11f11.2) 
    1 continue
      endif
      endif
c
      if(frequency(mdstep,nfstat)) then
      call swatch(pdate,ptime)
      write(lfnout,2000) pdate,ptime,nsum,nsumt
 2000 format(/,' MOLECULAR DYNAMICS STATISTICAL INFORMATION',t110,2a10,
     + //,t41,2(3X,'Statistics over last ',I8,' steps',2X),/,
     + t41,2(3X,'Average',5X,'RMS fluct',5X,'Drift/ps',3X),/)
      facs=one/dble(nsum)
      fact=one/dble(nsumt)
      tfacs=one/(t2sum-facs*tsum*tsum)
      tfact=one/(t2sumt-fact*tsumt*tsumt)
      do 15 i=1,nprop
      j=ixp(i)
      if(lp(j)) write(lfnout,2001) pronam(j)(1:39),
     + psum(j)*facs,sqrt(abs((p2sum(j)-psum(j)*psum(j)*facs)*facs)),
     + (pslop(j)-facs*psum(j)*tsum)*tfacs,
     + psumt(j)*fact,sqrt(abs((p2sumt(j)-psumt(j)*psumt(j)*fact)*fact)),
     + (pslopt(j)-fact*psumt(j)*tsumt)*tfact,
     + pronam(j)(40:50)
 2001 format(1x,a39,t41,2(3(1pe12.5,1x),1x),a11)
      psum(j)=zero
      p2sum(j)=zero
      pslop(j)=zero
   15 continue
      tsum=zero
      t2sum=zero
      nsum=0
      lhdr=.false.
      endif
c
      return
      end
      subroutine prp_print()
c
      implicit none
c
#include "prp_common.fh"
c
      character*10 pdate,ptime
      integer i,j
c
      if(me.ne.0) return
c
      call swatch(pdate,ptime)
      write(lfnout,1000) pdate,ptime
 1000 format(/,' SINGLE POINT PROPERTIES',t110,2a10,/)
      do 1  i=1,nprop
      j=ixp(i)
      if(lp(j)) write(lfnout,1001) pronam(j)(1:39),p(j),pronam(j)(40:50)
 1001 format(1x,a39,t41,1pe18.9,1x,a11)
    1 continue
c
      return
      end
      subroutine prp_record()
c
      implicit none
c
#include "prp_common.fh"
c
      integer i
      character*10 pdate,ptime
c
      if(me.ne.0) return
c
      if(.not.lfhdr) then
      call swatch(pdate,ptime)
      write(lfnprp,3000) nprop,pdate,ptime,np,npfft
 3000 format(i7,1x,2a10,2i5)
      write(lfnprp,3001) (pronam(ixp(i)),i=1,nprop)
 3001 format(a50)
      lfhdr=.true.
      endif
      write(lfnprp,3002)
 3002 format('frame')
      if(iprop.eq.0) then
      write(lfnprp,3003) (p(ixp(i)),i=1,nprop)
      else
      write(lfnprp,3003) (psump(ixp(i))/dble(nsump),i=1,nprop)
 3003 format(4(1pe12.5))
      endif
      nsump=0
      do 1 i=1,nprop
      psump(ixp(i))=zero
    1 continue
c
      return
      end
      subroutine prp_header()
c
      implicit none
c
#include "prp_common.fh"
c
      lfhdr=.false.
c
      return
      end
      logical function prp_mcti_step(ida,lda)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical prp_mcti_acc
      external prp_mcti_acc
c
      integer ida,lda
c
      real*8 fdata(28)
      real*8 aver,drift,stderr,corerr,ratio
      logical done
c
      lerror=ida.gt.lda
c
      call cf_mcti(fdata)
c
      call ga_dgop(mrg_d44,fdata,28,'+')
c
      if(me.eq.0) then
      done=prp_mcti_acc(ida,dbl_mb(i_dfr),dbl_mb(i_dfrm),fdata,
     + aver,drift,stderr,corerr,ratio)
      endif
c
      if(np.gt.1) then
      call ga_brdcst(mrg_d46,done,ma_sizeof(mt_log,1,mt_byte),0)
      endif
c
      prp_mcti_step=done
c
      return
      end
      logical function prp_mcti_acc(ida,dfr,dfrm,fdata,
     + aver,drift,stderr,corerr,ratio)
c
      implicit none
c
#include "prp_common.fh"
c
      integer ida
      real*8 dfr(mda),dfrm(mda),fdata(28)
      real*8 aver,drift,stderr,corerr,ratio,cerror
c
      integer i
      logical done
      real*8 dfrnom
c
      dfree=zero
      dfrnom=zero
      do 1 i=1,24
      dfree=dfree+fdata(i)
      if(i.ne.1.and.i.ne.13) dfrnom=dfrnom+fdata(i)
      deriv(i)=deriv(i)+fdata(i)
    1 continue
      nderiv=nderiv+1
c
      dfr(ida)=dfree
      dfrm(ida)=dfrnom
      if(dfree.ne.dfrnom) lfreem=.true.
      nda=ida
c
    2 continue
c
      if(lerror) then
      call error(lauto,lappr,1000,dfr,ida,
     + aver,drift,stderr,corerr,ratio)
      cerror=corerr
      if(.not.lauto) cerror=samrat*corerr 
      if(.not.lauto.and.cerror.lt.edacq) then
      lauto=.true.
      lappr=.true.
      goto 2
      else
      if(lauto) samrat=ratio
      lauto=.false.
      lappr=.false.
      endif
      done=cerror.lt.edacq.and.drift.lt.ddacq
      else
      done=.false.
      endif
c
      prp_mcti_acc=done
      return
      end
      subroutine prp_mcti_run(rlambd,dlambd,ndec)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      real*8 rlambd,dlambd
      integer ndec
c
      call prp_mcti_r(rlambd,dlambd,ndec,dbl_mb(i_dfr),dbl_mb(i_dfrm),
     + psumt(16),psumt(117),psumt(118),psumt(119),psumt(120),
     + psumt(51),psumt(41))
c
      return
      end
      subroutine prp_mcti_r(rlambd,dlambd,ndec,dfr,dfrm,taver,
     + ep2ave,ep3ave,ep2avm,ep3avm,ebias,dfbias)
c
      implicit none
c
#include "prp_common.fh"
c
      real*8 rlambd,dlambd,dfr(mda),dfrm(mda)
      real*8 taver,ep2ave,ep3ave,ebias,dfbias,ep2avm,ep3avm
      integer i,ndec
c
      if(me.eq.0) then
      write(lfngib,1000) nderiv,nda,rlambd,dlambd,ndec,nsa,
     + nbias,ebias/dble(nsumt),lfreem
 1000 format(2i7,2f12.6,2i8,i4,e20.12,4x,l1)
      write(lfngib,1001) deriv
 1001 format(4e20.12)
      write(lfngib,1002) (dfr(i),i=1,nda)
      write(lfngib,1002) (dfrm(i),i=1,nda)
 1002 format(4e20.12)
      write(lfngib,1003) nsumt,taver/dble(nsumt),
     + ep2ave/dble(nsumt),ep3ave/dble(nsumt),dfbias/ebias,
     + ep2avm/dble(nsumt),ep3avm/dble(nsumt)
 1003 format(i10,/,4e20.12,/,2e20.12)
      endif
      if(ndec.gt.0) call cf_wrtgib(lfngib)
c
      return
      end
      subroutine prp_mcti(n,filnam)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer n,i_dec,l_dec,ibl,i,igp_handle
      character*80 filnam,string
      character*5 gid
c
      if(npg.gt.1) then
      if(me.eq.0) then
      close(unit=lfngib,status='keep')
      endif
c
c     now go global 
c
      igp_handle=ga_pgroup_get_default()
      call ga_pgroup_set_default(ga_pgroup_get_world())
      call ga_sync()
      if(ga_nodeid().eq.0) then
c
      ibl=index(filnam,' ')-1
      string=filnam(1:ibl)//'.gib'
      open(unit=lfngib,file=string,status='unknown',form='formatted')
      rewind(lfngib)
      do 1 i=1,npg
      write(gid,3300) i-1
 3300 format(i5.5)
      string=filnam(1:ibl)//gid//'.gib'
      open(unit=75,file=string,status='old',form='formatted')
      rewind(unit=75)
    2 continue
      read(75,3301,end=99,err=99) string
      write(lfngib,3301) string
 3301 format(a)
      goto 2
   99 continue
      close(unit=75)
    1 continue
      close(unit=lfngib,status='keep')
c
      endif
c
c     go back to groups
c
      call ga_sync()
      call ga_pgroup_set_default(igp_handle)
c
      else if(ipg.eq.npg-1.or.npg.eq.1) then
      if(me.eq.0) then
      close(unit=lfngib,status='keep')
      endif
      endif
c
      if(ipg.eq.npg-1.or.npg.eq.1) then
c
      if(me.eq.0) then
      ibl=index(filnam,' ')-1
      string=filnam(1:ibl)//'.gib'
      open(unit=lfngib,file=string,status='old',form='formatted')
      rewind(lfngib)
      endif
c
      npgdec=n
c
      if(npgdec.gt.0) then
      if(.not.ma_push_get(mt_dbl,6*nsa,'dec',l_dec,i_dec))
     + call md_abort('Failed to allocate dec',0)
      else
      if(.not.ma_push_get(mt_dbl,1,'dec',l_dec,i_dec))
     + call md_abort('Failed to allocate dec',0)
      endif
c
      call prp_mcti_s(dbl_mb(i_dfr),dbl_mb(i_dfrm),dbl_mb(i_dec))
c
      if(.not.ma_pop_stack(l_dec))
     + call md_abort('Failed to deallocate dec',0)
c
      endif
c      endif
c
      return
      end
      subroutine prp_mcti_s(dfr,dfrm,dec)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      real*8 dfr(mda),dfrm(mda),dec(6,nsa)
c
      character*10 pdate,ptime
      integer i,j,k,number,ndec
      real*8 rlambd,dlambd,ddrft,dsterr,dcerr,ratio
      real*8 freeti,errti,drftti,taver,ep2ave,ep3ave,ebias,dfbias
      real*8 ep2avm,ep3avm,dtmp
      real*8 epr,epf,etp,slambd,rnum,fterm(24),gbias,freeb,freem
      logical lfrm
c
      gbias=zero
c
      if(me.eq.0) then
c
      do 1 k=1,8
c
      if(k.le.4.and.npgdec.eq.0) goto 1
c
      call swatch(pdate,ptime)
      if(k.eq.1) then
      write(lfnout,1001)
 1001 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION ',
     + 'DECOMPOSITION',//,' Solvent derivatives in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      elseif(k.eq.2) then
      write(lfnout,1002)
 1002 format(//,' Solute derivatives in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      elseif(k.eq.3) then
      write(lfnout,1003)
 1003 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION ',
     + 'DECOMPOSITION',//,' Solvent contributions in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      elseif(k.eq.4) then
      write(lfnout,1004)
 1004 format(//,' Solute contributions in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      elseif(k.eq.5) then
      if(lfreem) then
      write(lfnout,1005) pdate,ptime
 1005 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION ',
     + 'EXCLUDING MASS CONTRIBUTIONS',t110,2a10,//,
     + ' Run Lambda   Size',
     + '  Derivative  Derivative  Derivative Lambda',
     + ' Free Energy Free Energy Free Energy Free Energy  Sampling',/,
     + 18x,
     + '     Average       Error       Drift       ',
     + ' Accumulated       Error       Drift   Corrected     Ratio',/,
     + 18x,
     + '      kJ/mol      kJ/mol   kJ/mol ps       ',
     + '      kJ/mol      kJ/mol   kJ/mol ps      kJ/mol',/)
      else
      write(lfnout,1006) pdate,ptime
 1006 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION',
     + t110,2a10,//,
     + ' Run Lambda   Size',
     + '  Derivative  Derivative  Derivative Lambda',
     + ' Free Energy Free Energy Free Energy Free Energy  Sampling',/,
     + 18x,
     + '     Average       Error       Drift       ',
     + ' Accumulated       Error       Drift   Corrected     Ratio',/,
     + 18x,
     + '      kJ/mol      kJ/mol   kJ/mol ps       ',
     + '      kJ/mol      kJ/mol   kJ/mol ps      kJ/mol',/)
      endif
      elseif(k.eq.6) then
      if(lfreem) then
      write(lfnout,1007) pdate,ptime
 1007 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION ',
     + 'INCLUDING MASS CONTRIBUTIONS',t110,2a10,//,
     + ' Run Lambda   Size',
     + '  Derivative  Derivative  Derivative Lambda',
     + ' Free Energy Free Energy Free Energy Free Energy  Sampling',/,
     + 18x,
     + '     Average       Error       Drift       ',
     + ' Accumulated       Error       Drift   Corrected     Ratio',/,
     + 18x,
     + '      kJ/mol      kJ/mol   kJ/mol ps       ',
     + '      kJ/mol      kJ/mol   kJ/mol ps      kJ/mol',/)
      endif
      elseif(k.eq.7) then
      if(lfreem) then
      write(lfnout,1008) pdate,ptime
 1008 format(//' MULTISTEP THERMODYNAMIC PERTURBATION ',
     + 'EXCLUDING MASS CONTRIBUTIONS',t110,2a10,//,
     + ' Run Lambda   Size',
     + ' Temperature     Reverse     Forward    Ensemble Lambda',
     + ' Accumulated        Bias   Corrected',/,18x,
     + '           K      kJ/mol      kJ/mol      kJ/mol       ',
     + '      kJ/mol      kJ/mol      kJ/mol',/)
      else
      write(lfnout,1009) pdate,ptime
 1009 format(//' MULTISTEP THERMODYNAMIC PERTURBATION',
     + t110,2a10,//,
     + ' Run Lambda   Size',
     + ' Temperature     Reverse     Forward    Ensemble Lambda',
     + ' Accumulated        Bias   Corrected',/,18x,
     + '           K      kJ/mol      kJ/mol      kJ/mol       ',
     + '      kJ/mol      kJ/mol      kJ/mol',/)
      endif
      else
      if(lfreem) then
      write(lfnout,1010) pdate,ptime
 1010 format(//' MULTISTEP THERMODYNAMIC PERTURBATION ',
     + 'INCLUDING MASS CONTRIBUTIONS',t110,2a10,//,
     + ' Run Lambda   Size',
     + ' Temperature     Reverse     Forward    Ensemble Lambda',
     + ' Accumulated        Bias   Corrected',/,18x,
     + '           K      kJ/mol      kJ/mol      kJ/mol       ',
     + '      kJ/mol      kJ/mol      kJ/mol',/)
      endif
      endif
c
      rewind(lfngib)
c
      freeti=zero
      freeb=zero
      freem=zero
      errti=zero
      drftti=zero
      etp=zero
      slambd=zero
      gbias=zero
      do 2 i=1,24
      fterm(i)=zero
    2 continue
      if(npgdec.gt.0) then
      do 22 i=1,nsa
      dec(1,i)=zero
      dec(2,i)=zero
      dec(3,i)=zero
      dec(4,i)=zero
      dec(5,i)=zero
      dec(6,i)=zero
   22 continue
      endif
      do 3 i=1,mrun
      read(lfngib,2000) nderiv,nda,rlambd,dlambd,ndec,nsa,nbias,ebias,
     + lfrm
      if(lfrm) lfreem=.true.
 2000 format(2i7,2f12.6,2i8,i4,e20.12,4x,l1)
      read(lfngib,2001) deriv
 2001 format(4e20.12)
      if(mda.lt.nda) then
      read(lfngib,2002) (dtmp,j=1,nda)
      read(lfngib,2002) (dtmp,j=1,nda)
      else
      read(lfngib,2002) (dfr(j),j=1,nda)
      read(lfngib,2002) (dfrm(j),j=1,nda)
      endif
 2002 format(4e20.12)
      read(lfngib,2003) number,taver,ep2ave,ep3ave,dfbias,
     + ep2avm,ep3avm
 2003 format(i10,/,4e20.12,/,2e20.12)
      if(nbias.gt.0) gbias=-rgas*taver*log(ebias)
      if(ndec.gt.0) call cf_rdgib(lfngib,dec,dlambd/dble(ndec))
c
      slambd=slambd+dlambd
      rnum=dlambd/dble(nderiv)
c
      do 4 j=1,24
      fterm(j)=fterm(j)+rnum*deriv(j)
    4 continue
c
      freem=freem+rnum*(deriv(1)+deriv(13))
c
      if(k.eq.1) then
      rnum=one/dble(nderiv)
      write(lfnout,1101) i,slambd,(rnum*deriv(j),j=1,12)
 1101 format(i4,f7.3,12f10.3)
      elseif(k.eq.2) then
      rnum=one/dble(nderiv)
      write(lfnout,1102) i,slambd,(rnum*deriv(j),j=13,24)
 1102 format(i4,f7.3,12f10.3)
      elseif(k.eq.3) then
      write(lfnout,1103) i,slambd,(fterm(j),j=1,12)
 1103 format(i4,f7.3,12f10.5)
      elseif(k.eq.4) then
      write(lfnout,1104) i,slambd,(fterm(j),j=13,24)
 1104 format(i4,f7.3,12f10.5)
      elseif(k.eq.5) then
      call error(.true.,.true.,1000,dfrm,nda,
     + dfree,ddrft,dsterr,dcerr,ratio)
c
      freeti=freeti+dfree*dlambd
      freeb=freeb+dfbias*dlambd
      errti=errti+(dcerr*dlambd)**2
      drftti=drftti+ddrft*dlambd
c
      write(lfnout,1105) i,rlambd,nderiv,dfree,dcerr,ddrft/tstep,
     + slambd,freeti,sqrt(errti),drftti/tstep,freeb,ratio
 1105 format(i4,f7.3,i7,f12.3,2f12.6,f7.3,f12.3,2f12.6,f12.3,f10.5)
c
      elseif(k.eq.6) then
      if(lfreem) then
      call error(.true.,.true.,1000,dfr,nda,
     + dfree,ddrft,dsterr,dcerr,ratio)
c
      freeti=freeti+dfree*dlambd
      freeb=freeb+dfbias*dlambd
      errti=errti+(dcerr*dlambd)**2
      drftti=drftti+ddrft*dlambd
c
      write(lfnout,1106) i,rlambd,nderiv,dfree,dcerr,ddrft/tstep,
     + slambd,freeti,sqrt(errti),drftti/tstep,freeb,ratio
 1106 format(i4,f7.3,i7,f12.3,2f12.6,f7.3,f12.3,2f12.6,f12.3,f10.5)
      endif
c
      elseif(k.eq.7) then
      epr=-rgas*taver*log(ep2ave)
      epf=-rgas*taver*log(ep3ave)
      etp=etp-epr+epf
      write(lfnout,1107) i,rlambd,number,taver,epr,epf,epf-epr,
     + slambd,etp,-gbias,etp-gbias
 1107 format(i4,f7.3,i7,4f12.3,f7.3,3f12.3)
      else
      if(lfreem) then
      epr=-rgas*taver*log(ep2avm)
      epf=-rgas*taver*log(ep3avm)
      etp=etp-epr+epf
      write(lfnout,1108) i,rlambd,number,taver,epr,epf,epf-epr,
     + slambd,etp,-gbias,etp-gbias
 1108 format(i4,f7.3,i7,4f12.3,f7.3,3f12.3)
      endif
      endif
c
    3 continue
c
      if(ndec.gt.0.and.k.eq.4) call cf_print_deco(lfnout,dec)
c
    1 continue
c
      close(unit=lfngib)
c
      endif
c
      return
      end
      subroutine prp_wrtmro(lfnmro,ndec)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      integer lfnmro,ndec
c
      call prp_wtmro(lfnmro,ndec,dbl_mb(i_dfr),dbl_mb(i_dfrm))
c
      return
      end
      logical function prp_rdmri(lfnmri,ndec,mropt)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      logical prp_rmri
      external prp_rmri
c
      integer lfnmri,ndec,mropt
c
      prp_rdmri=prp_rmri(lfnmri,ndec,mropt,dbl_mb(i_dfr),dbl_mb(i_dfrm))
c
      return
      end
      subroutine prp_wtmro(lfnmro,ndec,dfr,dfrm)
c
      implicit none
c
#include "prp_common.fh"
c
      integer lfnmro,ndec
      real*8 dfr(mda),dfrm(mda)
c
      integer i
c
      if(me.eq.0) then
      write(lfnmro) nderiv,nda,nprop,nsum,nsumt,ndec,maxp
      write(lfnmro) deriv
      write(lfnmro) (dfr(i),i=1,nda)
      write(lfnmro) (dfrm(i),i=1,nda)
      write(lfnmro) tsum,t2sum,tsumt,t2sumt
      write(lfnmro) (psum(i),i=1,maxp)
      write(lfnmro) (p2sum(i),i=1,maxp)
      write(lfnmro) (pslop(i),i=1,maxp)
      write(lfnmro) (psumt(i),i=1,maxp)
      write(lfnmro) (p2sumt(i),i=1,maxp)
      write(lfnmro) (pslopt(i),i=1,maxp)
      endif
      if(ndec.gt.0) call cf_wrtmro(lfnmro)
c
      return
      end
      logical function prp_rmri(lfnmri,ndec,mropt,dfr,dfrm)
c
      implicit none
c
#include "prp_common.fh"
c
      logical cf_rdmri
      external cf_rdmri
c
      integer lfnmri,ndec,mropt
      real*8 dfr(mda),dfrm(mda),dtmp
c
      integer i,nprp,nmxp
c
      read(lfnmri,err=9,end=9) nderiv,nda,nprp,nsum,nsumt,ndec,nmxp
      if(mropt.ne.2) then
      if(nprop.ne.nprp) then
      call md_abort('Number of properties nprop changed',0)
      endif
      if(maxp.ne.nmxp) then
      call md_abort('Number of properties maxp changed',0)
      endif
      endif
      read(lfnmri,err=9,end=9) deriv
      if(mda.lt.nda) then
      read(lfnmri,err=9,end=9) (dtmp,i=1,nda)
      read(lfnmri,err=9,end=9) (dtmp,i=1,nda)
      else
      read(lfnmri,err=9,end=9) (dfr(i),i=1,nda)
      read(lfnmri,err=9,end=9) (dfrm(i),i=1,nda)
      endif
      read(lfnmri,err=9,end=9) tsum,t2sum,tsumt,t2sumt
      if(mropt.ne.2) then
      read(lfnmri,err=9,end=9) (psum(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (p2sum(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (pslop(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (psumt(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (p2sumt(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (pslopt(i),i=1,maxp)
      else
      do 1 i=1,6
      read(lfnmri,err=9,end=9) dtmp
    1 continue
      endif
      if(ndec.gt.0) then
      prp_rmri=cf_rdmri(lfnmri)
      else
      prp_rmri=.true.
      endif
      return
c
    9 continue
      prp_rmri=.false.
      return
      end
      integer function prp_dfr(ncopy)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      integer prp_dfr_copy
      external prp_dfr_copy
c
      integer ncopy
c
      prp_dfr=prp_dfr_copy(ncopy,dbl_mb(i_dfr),dbl_mb(i_dfrm))
c
      return
      end
      integer function prp_dfr_copy(ncopy,dfr,dfrm)
c
      implicit none
c
#include "prp_common.fh"
c
      integer ncopy
      real*8 dfr(mda),dfrm(mda)
c
      integer i,ioff
c
      if(ncopy.lt.nda) then
      ioff=nda-ncopy
      do 1 i=1,ncopy
      dfr(i)=dfr(i+ioff)
      dfrm(i)=dfrm(i+ioff)
    1 continue
      nda=ncopy
      endif
c
      prp_dfr_copy=nda
      return
      end
