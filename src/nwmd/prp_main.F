      subroutine prp_step(mdstep,stime,eww,esw,ess,epme)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme
      real*8 epot,stime
c
      call prp_proper(mdstep,stime,eww,esw,ess,epme,epot)
      if(me.eq.0) call prp_stat(mdstep,stime,eww,esw,ess,epme)
c
      return
      end
      subroutine prp_proper(mdstep,stime,eww,esw,ess,epme,epot)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme,epot
      real*8 eq,el,ep,ek,et,ewp
      real*8 tempw,temps,stime,esk,ewk
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits
c
      integer i,j,k,it
      character*10 pdate,ptime
      real*8 facs,fact,tfacs,tfact
c
      if(ltwin) then
      do 1 i=1,mpe
      eww(i,1)=eww(i,1)+eww(i,2)
    1 continue
      do 2 i=1,msf
      do 3 j=1,mpe
      esw(i,j,1)=esw(i,j,1)+esw(i,j,2)
    3 continue
    2 continue
      do 4 i=1,msf
      do 5 j=1,msf
      do 6 k=1,mpe
      ess(i,j,k,1)=ess(i,j,k,1)+ess(i,j,k,2)
    6 continue
    5 continue
    4 continue
      endif
c
      el=eww(5,1)
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
      do 12 j=1,msf
      el=el+esw(j,5,1)
      eq=eq+esw(j,6,1)+esw(j,8,1)
      do 13 k=1,msf
      el=el+ess(k,j,5,1)
      eq=eq+ess(k,j,6,1)+ess(k,j,8,1)
   13 continue
   12 continue
c
c     if using cafe get properties from it
c
      if(lcafe) then
      call cf_proper(volume,temp,tempw,temps,pres,tmpscl,tmsscl,prsscl,
     + ewk,esk,nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits)
      endif
c
c     if using space get properties from it
c
      if(lspac) then
      endif
c
      ep=ep+epme
      eq=eq+epme
      ek=ewk+esk
      et=ep+ek
      epot=ep
c
c     fill the property vector
c
      do 7 i=1,250
      p(i)=zero
    7 continue
c
      if(me.eq.0) then
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
      if(me.eq.0) then
      p(12)=volume
      p(14)=pres
      p(15)=prsscl
      p(16)=temp
      p(17)=tempw
      p(18)=temps
      p(21)=tmpscl
      p(22)=tmsscl
      p(24)=eq
      p(25)=el
      p(26)=ep
      p(27)=ek
      p(32)=ep+ek
      p(54)=eww(6,1)+eww(8,1)
      p(55)=eww(5,1)
      p(56)=eww(1,1)
      p(57)=eww(2,1)
      p(58)=eww(3,1)
      p(59)=eww(4,1)
      p(64)=ewk
      p(66)=ewp
      p(67)=ewk+ewp
      endif
c
      isprop=96
      if(nsf.gt.0) then
      do 33 i=1,nsf
      it=isprop+(i-1)*27
      p(it+5)=esw(i,6,1)
      p(it+6)=esw(i,5,1)
      p(it+7)=ess(i,i,1,1)
      p(it+8)=ess(i,i,2,1)
      p(it+9)=ess(i,i,3,1)
      p(it+10)=ess(i,i,4,1)
   33 continue
      it=it+27
      if(nsf.gt.1) then
      do 35 i=1,nsf-1
      do 36 j=i+1,nsf
      it=it+1
      it=it+1
      it=it+1
      it=it+1
   36 continue
   35 continue
      endif
      endif
c
      if(np.gt.0) then
      call ga_dgop(mrg_d44,p(9),3,'max')
      if(me.ne.0) then
      p(9)=zero
      p(10)=zero
      p(11)=zero
      endif
      call ga_dgop(mrg_d43,p,22,'+')
      call ga_dgop(mrg_d40,eww,mpe,'+')
      call ga_dgop(mrg_d41,esw,msf*mpe,'+')
      call ga_dgop(mrg_d42,ess,msf*msf*mpe,'+')
      call ga_dgop(mrg_d43,epme,1,'+')
      endif
c
      return
      end
      subroutine prp_stat(mdstep,stime,eww,esw,ess,epme)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme
      real*8 eq,el,ep,ek,et,ewp
      real*8 tempw,temps,stime,esk,ewk
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits
c
      integer i,j,k,it
      character*10 pdate,ptime
      real*8 facs,fact,tfacs,tfact
c
      nsum=nsum+1
      nsumt=nsumt+1
      do 8 i=1,250
      if(abs(p(i)).lt.tiny) p(i)=zero
      psum(i)=psum(i)+p(i)
      p2sum(i)=p2sum(i)+p(i)*p(i)
      pslop(i)=pslop(i)+stime*p(i)
      psumt(i)=psumt(i)+p(i)
      p2sumt(i)=p2sumt(i)+p(i)*p(i)
      pslopt(i)=pslopt(i)+stime*p(i)
    8 continue
      tsum=tsum+stime
      t2sum=t2sum+stime*stime
      tsumt=tsumt+stime
      t2sumt=t2sumt+stime*stime
c
      if(mod(mdstep,nfoutp).eq.0) then
      if(.not.lhdr) then
      call swatch(pdate,ptime)
      write(lfnout,1000) pdate,ptime
 1000 format(/,' MOLECULAR DYNAMICS TIME STEP INFORMATION',T110,2A10,//,
     + '    Time      Temp    Pres      Volume Tscalw Tscals  Pscal ',
     + '   U(Cou)    U(L-J)      U(pot)     U(kin)     U(tot) ',/,
     + '    ps         K      Pa        nm**3              ',
     + '   kJ/mol    kJ/mol      kJ/mol     kJ/mol     kJ/mol ',/)
      lhdr=.true.
      endif
      write(lfnout,1001) stime,temp,pres,volume,p(21),p(22),p(15),
     + p(24),p(25),p(26),p(27),p(32)
 1001 format(1x,f10.5,0pf8.2,1pe9.2,0pf10.3,3f7.4,5(1pe11.4),i5,i7)
      endif
c
      if(mod(mdstep,nfstat).eq.0) then
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
      write(lfnout,2001) pronam(j)(1:39),
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
      if(mod(mdstep,nfprop).eq.0) then
      write(lfnprp,3000) (p(ixp(i)),i=1,nprop)
 3000 format(4e12.5)
      endif
c
      return
      end
