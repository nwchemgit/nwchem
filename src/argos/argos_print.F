      subroutine argos_print()
c
      implicit none
c
#include "argos.fh"
c
      character*22 ctype
c
      if(me.eq.0) then
      if(lfnout.ne.6)
     + open(unit=lfnout,file=filout(1:index(filout,' ')-1),
     + form='formatted',status='unknown')
c
      if(ntype.eq.0) then
      ctype='ENERGY EVALUATION     '
      elseif(ntype.eq.1) then
      ctype='ENERGY MINIMIZATION   '
      elseif(ntype.eq.2) then
      ctype='MOLECULAR DYNAMICS    '
      elseif(ntype.eq.3) then
      ctype='FREE ENERGY EVALUATION'
      else
      call md_abort('Unknown calculation type',ntype)
      endif
c
      call swatch(today,now)
      write(lfnout,1000) ctype,today,now,npw,npg,np,ncp
 1000 format(' ARGOS MOLECULAR DYNAMICS SIMULATION PACKAGE',//,
     + ' MODULE: ',a22,t110,2a10,//,
     + ' Number of processes',t30,i6,/,
     + ' Number of process groups',t30,i6,/,
     + ' Number of processes in group',t30,i6,/,
     + ' Number of cores per cell',t30,i6,/)
c
      write(lfnout,1001) titinp(1),datinp,timinp,titinp(2),titinp(3)
 1001 format(/,' Title ',t10,a,t110,2a10,/,t10,a,/,t10,a)
c
      write(lfnout,1002) filnam(1:index(filnam,' ')-1)
 1002 format(/,' System ',a)
c
      if(ntype.le.2) then
      write(lfnout,1003) iset
 1003 format(/,' Force field parameter set ',i4)
      endif
c
      if(nserie.eq.0) then
      write(lfnout,1026)
 1026 format(/,' Initial simulation',/)
      elseif(nserie.eq.1) then
      write(lfnout,1027)
 1027 format(/,' Restart simulation',/)
      else
      write(lfnout,1028)
 1028 format(/,' Continuation simulation',/)
      endif
c
      if(ntype.eq.3) then
      if(mropt.eq.0) then
      write(lfnout,1032)
 1032 format(' Multiple run initial simulation')
      elseif(mropt.eq.1) then
      write(lfnout,1033)
 1033 format(' Multiple run simulation with initial conditions')
      else
      write(lfnout,1034)
 1034 format(' Multiple run extension simulation')
      endif
      write(lfnout,1025) mrun
 1025 format(' Number of runs ',i7)
      endif
      if(ntype.ge.2) then
      write(lfnout,1024) mequi,mdacq
 1024 format(' Number of equilibration steps ',i7,/
     + ' Number of data gathering steps',i7)
      endif
c
      if(ntype.eq.3.and.issscl.gt.0) then
      write(lfnout,1029) ssshft
 1029 format(/,' Separation shifted scaling, delta ',f12.6,' nm**2')
      endif
c
      if(ipme.gt.0) then
      write(lfnout,1030) morder,ngx,ngy,ngz,nodpme,pmetol
 1030 format(/,' Particle-mesh Ewald, spline to order ',i5,/,
     + ' Grid size ',i5,'x',i5,'x',i5,/,
     + ' Number of processors for p-FFT ',i5,/,
     + ' Tolerance at cutoff ',1pe12.5)
      if(imfft.le.1) then
      write(lfnout,1036)
 1036 format(' Using NWChem 3D-pFFT')
      else
      write(lfnout,1037)
#if defined(ESSL)
 1037 format(' Using PESSL 3D-pFFT')
#else
 1037 format(' Using system specific 3D-pFFT')
#endif
      endif
      endif
c
      if(ipolt.eq.1) then
      write(lfnout,1031) mpolit,ptol
 1031 format(/,' Polarization model, maximum iterations ',i5,/,
     + ' Tolerance ',1pe12.5)
      endif
c
      if(.not.ltwin) then
      write(lfnout,1004) rshort
 1004 format(/,' Cutoff radius ',f12.6, ' nm')
      else
      write(lfnout,1005) rshort,rlong
 1005 format(/,' Cutoff radius short range ',f12.6, ' nm',
     +       /,' Cutoff radius long range  ',f12.6, ' nm')
      endif
c
      if(ifield.gt.0) then
      write(lfnout,1035) xfield,xffreq,xfvect
 1035 format(/,' External electrostatic field strength ',f12.6,/,
     + 25x,'frequency ',f12.6,/,
     + 25x,'field vector ',3f12.6)
      endif
c
      if(ntype.ge.2) then
      if(ivreas.ne.0) then
      if(ivreas.gt.0) then
      write(lfnout,1022) tgauss,ivreas
 1022 format(/,' Velocity reassignment temperature ',f12.6,' K',/,
     + ' Velocity reassignment frequency   ',i5,/)
      else
      write(lfnout,1122) tgauss
 1122 format(/,' Velocity reassignment temperature ',f12.6,' K',/,
     + ' Velocity reassignment in first step only',/)
      endif
      endif
      if(itscal.gt.0) then
      write(lfnout,1021) tmpext1,tmpext2,tmprlx,tmsrlx,tann1,tann2
 1021 format(/,' Isothermal ensemble external temperature ',
     + 2f12.6,' K',/,
     + ' Temperature relaxation time solvent      ',f12.6,' ps',/,
     + ' Temperature relaxation time solute       ',f12.6,' ps',/,
     + ' Temperature annealing between times      ',2f12.6,' ps')
      endif
      if(ipscal.ne.0) then
      write(lfnout,1023) prsext,prsrlx,compr
 1023 format(' Isobaric ensemble external pressure ',1pe12.5,' Pa',/,
     + ' Pressure relaxation time ',0pf12.6,' ps',/,
     + ' Compressebility ',1pe12.5,' m**2/N')
      if(ipopt.eq.12) write(lfnout,2023)
 2023 format(' Pressure scaling in x and y only')
      if(ipopt.eq.3) write(lfnout,3023)
 3023 format(' Pressure scaling in z only')
      if(ipopt.eq.123) write(lfnout,3024)
 3024 format(' Pressure scaling in x and y separate from z')
      endif
      if(iguide.gt.0) then
      write(lfnout,1038) fguide,tguide
 1038 format(/,' Self-guided molecular dynamics with force coupling of',
     + f8.5,' and relaxation time of',f8.5,' ps')
      endif
      if(mdalgo.eq.1) then
      write(lfnout,1039)
 1039 format(/,' Leap-frog integration')
      else
      write(lfnout,1040)
 1040 format(/,' Leap-frog integration (Brown-Clarke)')
      endif
      endif
c
      if(scaleq.ge.zero) then
      write(lfnout,1041) scaleq
 1041 format(/,' Solute charges scaled by factor ',f8.5)
      endif
c
      if(ntype.eq.1) then
      if(msdit.gt.0) then
      write(lfnout,1006) msdit,dx0sd,dxmsd,dxsdmx
 1006 format(/,' Maximum number of steepest descent steps ',i5,/,
     + ' Initial step size ',f12.6,' nm',/,
     + ' Maximum step size ',f12.6,' nm',/,
     + ' Minimum step size ',f12.6,' nm')
      endif
      if(mcgit.gt.0) then
      write(lfnout,1007) mcgit,ncgcy,dx0cg,dxsdmx
 1007 format(/,' Maximum number of conjugate gradient steps ',i5,/,
     + ' Number of conjugate gradient steps per cycle ',i5,/,
     + ' Initial interval size ',f12.6,' nm',/,
     + ' Minimum step size ',f12.6,' nm')
      endif
      endif
c
      write(lfnout,1008) mshitw,tlwsha,mshits,tlssha
 1008 format(/,' Maximum number of solvent SHAKE iterations ',i5,
     + ', solvent SHAKE tolerance ',f12.6,' nm',/,
     + ' Maximum number of solute SHAKE iterations  ',i5,
     + ', solute SHAKE tolerance  ',f12.6,' nm',/)
c
      if(ntype.ge.2) then
      write(lfnout,1019) nfpair
 1019 format(' Frequency update pair lists ',t55,i5)
      if(ltwin) write(lfnout,1020) nflong
 1020 format(' Frequency update long range forces ',t55,i5)
      write(lfnout,1017) nfslow
 1017 format(' Frequency removal overall motion ',t55,i5)
      write(lfnout,1018) nfcntr
 1018 format(' Frequency solute centering ',t55,i5)
      if(icentr.eq.1) write(lfnout,1818)
 1818 format(' Centering in z-direction only')
      if(icentr.eq.2) write(lfnout,1819)
 1819 format(' Centering in xyplane only')
      write(lfnout,1009) nfoutp
 1009 format(' Frequency printing step information ',t55,i5)
      write(lfnout,1010) nfstat
 1010 format(' Frequency printing statistical information ',t55,i5)
      endif
      if(ntype.eq.1) then
      write(lfnout,1011) nfqrs
 1011 format(' Frequency recording minimum energy structure ',t55,i5)
      endif
      if(ntype.ge.2) then
      write(lfnout,1012) nfrest
 1012 format(' Frequency recording restart file ',t55,i5)
      write(lfnout,1013) nfcoor
 1013 format(' Frequency recording system coordinates ',t55,i5)
      write(lfnout,1014) nfscoo
 1014 format(' Frequency recording solute coordinates ',t55,i5)
      write(lfnout,1015) nfvelo
 1015 format(' Frequency recording system velocities ',t55,i5)
      write(lfnout,1016) nfsvel
 1016 format(' Frequency recording solute velocities ',t55,i5)
      write(lfnout,1115) nfforc
 1115 format(' Frequency recording system forces ',t55,i5)
      write(lfnout,1116) nfsfor
 1116 format(' Frequency recording solute forces ',t55,i5)
      endif
c
      write(lfnout,2000)
 2000 format(/,' LOAD BALANCING',/)
c
 
      if(loadb.eq.0) write(lfnout,2001)
 2001 format(' None')
      if(loadb.eq.2) write(lfnout,2002)
 2002 format(' Redistribution of inter-processor box pairs')
      if(loadb.eq.1) write(lfnout,2003) factld
 2003 format(' Resizing of processor domains: smallest scaling ',
     + t55,f12.6)
      if(loadb.eq.3) write(lfnout,2004) lbpair,factld
 2004 format(' Redistribution of inter-processor box pairs: retry ',
     + t55,i5,/,' Resizing of processor domains: smallest scaling ',
     + t55,f12.6)
      if(isload.gt.0) write(lfnout,2104)
 2104 format(' Resizing in z-dimension only')
c
      if(ireset.gt.0) write(lfnout,2005)
 2005 format(/,' Load balancing reset')
c
      write(lfnout,2012) nfload
 2012 format(/,' Load balancing frequency ',i5)
c
      if(itload.eq.0) write(lfnout,2006)
 2006 format(/,' Load balancing based on last synchronization time')
      if(itload.eq.1) write(lfnout,2007)
 2007 format(/,' Load balancing based on minimum synchronization time')
      if(itload.eq.2) write(lfnout,2008)
 2008 format(/,' Load balancing based on average synchronization time')
      if(itload.eq.3) write(lfnout,2009)
 2009 format(/,' Load balancing based on average synchronization time',
     + ' with minimum synchronization time on processor 0')
c
      if(ioload.eq.1) write(lfnout,2010)
 2010 format(/,' Time for I/O counted at each step')
      if(ioload.eq.1) write(lfnout,2011)
 2011 format(/,' Experimental load balancing')
c
      if(nfhop.gt.0) then
c
      write(lfnout,3000)
 3000 format(/,' QHOP PROTON HOPPING',/)
c
      write(lfnout,3001) nfhop
 3001 format(' Frequency hopping attempts ',t55,i5)
      write(lfnout,3002) rhop
 3002 format(' Donor-acceptor cutoff distance',t55,f12.6,' nm')
      write(lfnout,3003) thop
 3003 format(' Minimum time before backhop',t55,f12.6,' ps')
c
      endif
c
      endif
c
      return
      end
c $Id$
