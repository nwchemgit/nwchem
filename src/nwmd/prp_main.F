      subroutine prp_step(mdstep,stime,eww,esw,ess,esk,epme,uqmd,uqmmm)
c
c $Id: prp_main.F,v 1.36 2000-05-08 19:29:25 d3j191 Exp $
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
      real*8 stime,esk(msf),uqmd,uqmmm
c
      if(me.eq.0) call prp_stat(mdstep,stime,eww,esw,ess,esk,epme(iset))
c
      if(frequency(mdstep,nfprop)) call prp_record()
c
      return
      end
      subroutine prp_proper(mdstep,stime,eww,esw,ess,esk,epme,
     + uqmd,uqmmm,epot,epotw,epotsw,epots,vol,dwr,dsr,ekin,etot,
     + npolit,gsm)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      external timer_wall
      real*8 timer_wall
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme(3)
      real*8 epot,epotw,epotsw,epots,vol,ekin,etot
      real*8 eq,el,ep,ek,et,ewp,uqmd,uqmmm
      real*8 tempw,temps,stime,esk(msf),ewk,ep2,ep3,edrs,epmf
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits,npolit
      real*8 dwr,dsr(msm),gsm(msm,4)
c
      integer i,j,k,it
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
      do 12 j=1,msf
      el=el+esw(j,5,1)
      eq=eq+esw(j,6,1)
      do 13 k=1,msf
      el=el+ess(k,j,5,1)
      eq=eq+ess(k,j,6,1)
   13 continue
   12 continue
c
      if(me.eq.0) ep=ep+uqmd+uqmmm
c
c     if using cafe get properties from it
c
      if(lcafe) then
      call cf_proper(volume,temp,tempw,temps,pres,tmpscl,tmsscl,prsscl,
     + ewk,nwwl,nwws,nswl,nsws,nssl,nsss,nshitw,nshits,ep2,ep3,
     + edrs,epmf,virial)
      endif
c
c     if using space get properties from it
c
      if(lspac) then
      endif
c
      if(me.eq.0) ep=ep+edrs+epmf
c
      ep=ep+epme(iset)
      eq=eq+epme(iset)
      ep2=ep2+epme(2)-epme(1)
      ep3=ep3+epme(3)-epme(1)
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
      if(me.eq.0) p(40)=dfree
      p(52)=eww(6,1)
      p(53)=eww(5,1)
      p(54)=eww(8,1)
      p(55)=eww(7,1)
      p(56)=eww(1,1)
      p(57)=eww(2,1)
      p(58)=eww(3,1)
      p(59)=eww(4,1)
      p(66)=ewp
      p(67)=p(67)+ewp
      p(69)=dwr
      p(76)=epme(iset)
      p(97)=timer_wall(2)
      p(98)=timer_wall(3)
      if(me.eq.0) p(99)=timer_wall(3)
c
      maxp=isprop
      if(nsf.gt.0) then
      it=0
      do 33 i=1,nsf
      it=isprop+(i-1)*27
      p(it+2)=ess(i,i,6,1)
      p(it+3)=ess(i,i,5,1)
      p(it+5)=esw(i,6,1)
      p(it+6)=esw(i,5,1)
      p(it+7)=ess(i,i,1,1)
      p(it+8)=ess(i,i,2,1)
      p(it+9)=ess(i,i,3,1)
      p(it+10)=ess(i,i,4,1)
      if(me.eq.0) then
      p(it+16)=esk(i)
      endif
   33 continue
      it=it+27
      if(nsf.gt.1) then
      do 35 i=1,nsf-1
      do 36 j=i+1,nsf
      it=it+1
      it=it+1
      p(it)=ess(i,j,6,1)+ess(j,i,6,1)
      it=it+1
      p(it)=ess(i,j,5,1)+ess(j,i,5,1)
      it=it+1
   36 continue
   35 continue
      endif
      maxp=it
c
      if(nsm.gt.0) then
      it=maxp
      do 37 i=1,nsm
      it=it+1
      p(it)=gsm(i,iset)
      if(p(it).lt.tiny) lp(it)=.false.
   37 continue
      nprop=0
      do 137 i=1,maxpro
      if(lp(i)) then
      nprop=nprop+1
      ixp(nprop)=i
      endif
  137 continue
      maxp=it
      endif
c
      endif
c
      if(np.gt.0) then
      call ga_dgop(mrg_d44,p(9),3,'max')
      if(me.ne.0) then
      p(9)=zero
      p(10)=zero
      p(11)=zero
      endif
      call ga_dgop(mrg_d45,p,maxp,'+')
      call ga_dgop(mrg_d40,eww,mpe,'+')
      call ga_dgop(mrg_d41,esw,msf*mpe,'+')
      call ga_dgop(mrg_d42,ess,msf*msf*mpe,'+')
      call ga_dgop(mrg_d43,epme,3,'+')
      endif
c
      p(98)=(p(98)-p(97))/p(98)
c
      p(38)=exp(-p(36)/(rgas*temp))
      p(39)=exp(-p(37)/(rgas*temp))
c
      dwr=p(69)/dble(6*nwm)
      if(nwm.gt.0) then
      p(69)=dwr
      p(70)=1.0d-6*dwr/stime
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
      do 8 i=1,maxpro
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
      rt=one/nwm
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
      write(lfnout,1001) pronam(j)(1:39),p(j),pronam(j)(40:50)
 1001 format(1x,a39,t41,1pe12.5,1x,a11)
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
      write(lfnprp,3000) nprop,pdate,ptime
 3000 format(i7,1x,2a10)
      write(lfnprp,3001) (pronam(ixp(i)),i=1,nprop)
 3001 format(a50)
      lfhdr=.true.
      endif
      write(lfnprp,3002) (p(ixp(i)),i=1,nprop)
 3002 format(4e12.5)
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
      real*8 fdata(26)
      real*8 aver,drift,stderr,corerr,ratio
      logical done
c
      lerror=ida.gt.lda
c
      call cf_mcti(fdata)
c
      call ga_dgop(mrg_d44,fdata,26,'+')
c
      if(me.eq.0) then
      done=prp_mcti_acc(ida,dbl_mb(i_dfr),fdata,
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
      logical function prp_mcti_acc(ida,dfr,fdata,
     + aver,drift,stderr,corerr,ratio)
c
      implicit none
c
#include "prp_common.fh"
c
      integer ida
      real*8 dfr(mda),fdata(26)
      real*8 aver,drift,stderr,corerr,ratio,cerror
c
      integer i
      logical done
c
      dfree=zero
      do 1 i=1,24
      dfree=dfree+fdata(i)
      deriv(i)=deriv(i)+fdata(i)
    1 continue
      nderiv=nderiv+1
c
      dfr(ida)=dfree
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
      subroutine prp_mcti_run(rlambd,dlambd)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      real*8 rlambd,dlambd
c
      call prp_mcti_r(rlambd,dlambd,dbl_mb(i_dfr),
     + psumt(16),psumt(38),psumt(39))
c
      return
      end
      subroutine prp_mcti_r(rlambd,dlambd,dfr,taver,ep2ave,ep3ave)
c
      implicit none
c
#include "prp_common.fh"
c
      real*8 rlambd,dlambd,dfr(mda),taver,ep2ave,ep3ave
      integer i
c
      if(me.eq.0) then
      write(lfngib,1000) nderiv,nda,rlambd,dlambd
 1000 format(2i7,2f12.6)
      write(lfngib,1001) deriv
 1001 format(4e20.12)
      write(lfngib,1002) (dfr(i),i=1,nda)
 1002 format(4e20.12)
      write(lfngib,1003) nsumt,taver/dble(nsumt),
     + ep2ave/dble(nsumt),ep3ave/dble(nsumt)
 1003 format(i10,3e20.12)
      endif
c
      return
      end
      subroutine prp_mcti()
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      call prp_mcti_s(dbl_mb(i_dfr))
c
      return
      end
      subroutine prp_mcti_s(dfr)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      real*8 dfr(mda)
c
      character*10 pdate,ptime
      integer i,j,k,number
      real*8 rlambd,dlambd,ddrft,dsterr,dcerr,ratio
      real*8 freeti,errti,drftti,taver,ep2ave,ep3ave
      real*8 epr,epf,etp,slambd,rnum,fterm(24)
c
      if(me.eq.0) then
c
      do 1 k=1,4
c
      call swatch(pdate,ptime)
      if(k.eq.1) then
      write(lfnout,1000) pdate,ptime
 1000 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION',
     + t110,2a10,//,
     + ' Run Lambda   Size',
     + '  Derivative  Derivative  Derivative Lambda',
     + ' Free Energy Free Energy Free Energy  Sampling',/,18x,
     + '     Average       Error       Drift       ',
     + ' Accumulated       Error       Drift     Ratio',/,18x,
     + '      kJ/mol      kJ/mol   kJ/mol ps       ',
     + '      kJ/mol      kJ/mol   kJ/mol ps',/)
      elseif(k.eq.2) then
      write(lfnout,1001)
 1001 format(//,' MULTICONFIGURATION THERMODYNAMIC INTEGRATION ',
     + 'DECOMPOSITION',//,' Solvent contributions in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      elseif(k.eq.3) then
      write(lfnout,1002)
 1002 format(//,' Solute contributions in kJ/mol',//,
     + ' Run Lambda',
     + '      Mass Solvnt LJ Solute LJ Solvnt el Solute el     Bonds',
     + '  Constrts    Angles Dihedrals Impropers Slvnt pol Solut pol',/)
      else
      write(lfnout,1003) pdate,ptime
 1003 format(//' MULTISTEP THERMODYNAMIC PERTURBATION',
     + t110,2a10,//,
     + ' Run Lambda   Size',
     + '     Reverse     Forward    Ensemble Lambda Accumulated ',/,18x,
     + '      kJ/mol      kJ/mol      kJ/mol             kJ/mol',/)
      endif
c
      rewind(lfngib)
c
      freeti=zero
      errti=zero
      drftti=zero
      etp=zero
      slambd=zero
      do 2 i=1,24
      fterm(i)=zero
    2 continue
      do 3 i=1,mrun
      read(lfngib,2000) nderiv,nda,rlambd,dlambd
 2000 format(2i7,2f12.6)
      read(lfngib,2001) deriv
 2001 format(4e20.12)
      read(lfngib,2002) (dfr(j),j=1,nda)
 2002 format(4e20.12)
      read(lfngib,2003) number,taver,ep2ave,ep3ave
 2003 format(i10,3e20.12)
c
      slambd=slambd+dlambd
      rnum=dlambd/dble(nderiv)
c
      do 4 j=1,24
      fterm(j)=fterm(j)+rnum*deriv(j)
    4 continue
c
      if(k.eq.1) then
      call error(.true.,.true.,1000,dfr,nda,
     + dfree,ddrft,dsterr,dcerr,ratio)
c
      freeti=freeti+dfree*dlambd
      errti=errti+(dcerr*dlambd)**2
      drftti=drftti+ddrft*dlambd
c
      write(lfnout,1004) i,rlambd,nderiv,dfree,dcerr,ddrft/tstep,
     + slambd,freeti,sqrt(errti),drftti/tstep,ratio
 1004 format(i4,f7.3,i7,f12.3,2f12.6,f7.3,f12.3,2f12.6,f10.5)
c
      elseif(k.eq.2) then
      write(lfnout,1005) i,slambd,(fterm(j),j=1,12)
 1005 format(i4,f7.3,12f10.5)
      elseif(k.eq.3) then
      write(lfnout,1006) i,slambd,(fterm(j),j=13,24)
 1006 format(i4,f7.3,12f10.5)
      else
      epr=-rgas*taver*log(ep2ave)
      epf=-rgas*taver*log(ep3ave)
      etp=etp-epr+epf
      write(lfnout,1007) i,rlambd,number,epr,epf,epf-epr,slambd,etp
 1007 format(i4,f7.3,i7,3f12.3,f7.3,f12.3)
      endif
c
    3 continue
    1 continue
c
      close(unit=lfngib)
c
      endif
c
      return
      end
      subroutine prp_wrtmro(lfnmro)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      integer lfnmro
c
      call prp_wtmro(lfnmro,dbl_mb(i_dfr))
c
      return
      end
      subroutine prp_rdmri(lfnmri)
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      integer lfnmri
c
      call prp_rmri(lfnmri,dbl_mb(i_dfr))
c
      return
      end
      subroutine prp_wtmro(lfnmro,dfr)
c
      implicit none
c
#include "prp_common.fh"
c
      integer lfnmro
      real*8 dfr(mda)
c
      integer i
c
      write(lfnmro) nderiv,nda,nprop,nsum,nsumt
      write(lfnmro) deriv
      write(lfnmro) (dfr(i),i=1,nda)
      write(lfnmro) tsum,t2sum,tsumt,t2sumt
      write(lfnmro) (psum(i),i=1,nprop)
      write(lfnmro) (p2sum(i),i=1,nprop)
      write(lfnmro) (pslop(i),i=1,nprop)
      write(lfnmro) (psumt(i),i=1,nprop)
      write(lfnmro) (p2sumt(i),i=1,nprop)
      write(lfnmro) (pslopt(i),i=1,nprop)
c
      return
      end
      subroutine prp_rmri(lfnmri,dfr)
c
      implicit none
c
#include "prp_common.fh"
c
      integer lfnmri
      real*8 dfr(mda)
c
      integer i,nprp
c
      read(lfnmri) nderiv,nda,nprp,nsum,nsumt
      if(nprop.ne.nprp) call errquit('Number of properties changed',0)
      read(lfnmri) deriv
      read(lfnmri) (dfr(i),i=1,nda)
      read(lfnmri) tsum,t2sum,tsumt,t2sumt
      read(lfnmri) (psum(i),i=1,nprop)
      read(lfnmri) (p2sum(i),i=1,nprop)
      read(lfnmri) (pslop(i),i=1,nprop)
      read(lfnmri) (psumt(i),i=1,nprop)
      read(lfnmri) (p2sumt(i),i=1,nprop)
      read(lfnmri) (pslopt(i),i=1,nprop)
c
      return
      end
