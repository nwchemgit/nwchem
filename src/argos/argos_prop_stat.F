      subroutine argos_prop_stat(mdstep,stime,eww,esw,ess,esk,epme)
c
      implicit none
c
#include "argos_prop.fh"
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
     + '   Ui(ele)    Ui(vdW) ',
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
      j=isprop+(i-1)*27
      rt=p(j+7)+p(j+8)+p(j+9)+p(j+10)+p(j+2)+p(j+3)+p(j+5)+p(j+6)
      write(lfnout,1012) i,p(j+7),p(j+8),p(j+9),p(j+10),
     + p(j+2),p(j+3),p(j+5),p(j+6),p(j+16),rt,rt+p(j+16)
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
c $Id$
