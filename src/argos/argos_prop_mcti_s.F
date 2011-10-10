      subroutine argos_prop_mcti_s(dfr,dfrm,dec)
c
      implicit none
c
#include "argos_prop.fh"
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
      if(ndec.gt.0) call argos_cafe_rdgib(lfngib,dec,dlambd/dble(ndec))
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
      if(ndec.gt.0.and.k.eq.4) call argos_cafe_print_deco(lfnout,dec)
c
    1 continue
c
      close(unit=lfngib)
c
      endif
c
      return
      end
c $Id$
