      subroutine prp_start(nseri,ntyp,nftri,lfn,fil,iout,iprp,igib,
     + nf1,nf2,nf3,ipropi,
     + lcf,lsp,ltwn,lpol,lpmei,llpstep,llpener,npstti,
     + nwmtot,msfi,nsfi,mpei,mdacq,mrn,ipset,isetp1,isetp2,tstp,
     + msmi,nsmi,nsai,ddi,edi,iprf,npmf,npenrg,ndistr,lpbci,nbiasi,
     + nodpme,nipmf,lguide,lqmmm,lqmd,iradgi,idifci,nbget,ipgi,npgi)
c
c $Id$
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nseri,ntyp,lfn,iout,iprp,igib,iprf,ndistr,nbiasi,nftri
      integer nwmtot,msfi,nsfi,nsai,mpei,nf1,nf2,nf3,npenrg,npstti
      integer mdacq,mrn,ipset,isetp1,isetp2,msmi,nsmi,npmf,nodpme
      integer ipropi,nipmf,iradgi,idifci,nbget,npgi,ipgi
      logical lcf,lsp,ltwn,lpol,lpmei,llpener,llpstep,lpbci
      real*8 tstp,ddi,edi
      character*255 fil
c
      integer i,j,it
      logical lpert2,lpert3,lguide,lqmmm,lqmd
c
      me=ga_nodeid()
      np=ga_nnodes()
      npfft=nodpme
      ipg=ipgi
      npg=npgi
c
      ntype=ntyp
      nserie=nseri
c
      limport=nserie.eq.0.and.nftri.gt.0
c
      lfnout=iout
      lfnprp=iprp
      lfngib=igib
      lcafe=lcf
      lspac=lsp
      ltwin=ltwn
      lpener=llpener
      lpstep=llpstep
      msf=msfi
      nsf=nsfi
      mpe=mpei
      nfoutp=nf1
      nfstat=nf2
      nfprop=nf3
      mrun=mrn
      nwm=nwmtot
      msm=msmi
      nsm=nsmi
      nsa=nsai
      nbias=nbiasi
      npstat=npstti
      npmfi=nipmf
      lpme=lpmei
c
      iset=ipset
      ipset1=isetp1
      ipset2=isetp2
c
      iprof=iprf
      iprop=ipropi
      iradgy=iradgi
      idifco=idifci
c
      npener=npenrg
      tstep=tstp
      ddacq=ddi
      edacq=edi
c
      lauto=.false.
      lappr=.false.
      lpmf=npmf.gt.0
      lnoe=ndistr.gt.0
      lpbc=lpbci
c
      lpert2=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.2.or.isetp2.eq.2))
      lpert3=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.3.or.isetp2.eq.3))
c
      do 1 i=1,maxpro
      lp(i)=.false.
    1 continue
c
      pronam(  1)='Step                                              '
      lp(1)=ntype.gt.1.or.limport
      pronam(  2)='Time                                            ps'
      lp(2)=ntype.gt.1.or.limport
      pronam(  3)='Number solvent-solvent interactions SR            '
      lp(3)=nwm.gt.0
      pronam(  4)='Number solvent-solvent interactions LR            '
      lp(4)=nwm.gt.0.and.ltwin
      pronam(  5)='Number solvent-solute interactions  SR            '
      lp(5)=nwm.gt.0.and.nsf.gt.0
      pronam(  6)='Number solvent-solute interactions  LR            '
      lp(6)=nwm.gt.0.and.nsf.gt.0.and.ltwin
      pronam(  7)='Number solute-solute interactions   SR            '
      lp(7)=nsf.gt.0
      pronam(  8)='Number solute-solute interactions   LR            '
      lp(8)=nsf.gt.0.and.ltwin
      pronam(  9)='Number of solvent SHAKE iterations                '
      lp(9)=nwm.gt.0
      pronam( 10)='Number of solute SHAKE iterations                 '
      lp(10)=nsf.gt.0
      pronam( 11)='Number of polarization SCF iterations             '
      lp(11)=lpol
      pronam( 12)='Volume                                       nm**3'
      lp(12)=.true.
      pronam( 13)='Mass density                               kg/m**3'
      lp(13)=.true.
      pronam( 14)='Pressure                                        Pa'
      lp(14)=.true.
      pronam( 15)='Pressure scaling                                  '
      lp(15)=.true.
      pronam( 16)='Temperature                                      K'
      lp(16)=.true.
      pronam( 17)='Temperature solvent                              K'
      lp(17)=nwm.gt.0
      pronam( 18)='Temperature solute                               K'
      lp(18)=nsf.gt.0
      pronam( 19)='Temperature translations                         K'
      pronam( 20)='Temperature rotations                            K'
      pronam( 21)='Temperature scaling solvent                       '
      lp(21)=nwm.gt.0
      pronam( 22)='Temperature scaling solute                        '
      lp(22)=nsf.gt.0
      pronam( 23)='Energy polarization                         kJ/mol'
      lp(23)=lpol
      pronam( 24)='Energy electrostatic                        kJ/mol'
      lp(24)=.true.
      pronam( 25)='Energy van der Waals                        kJ/mol'
      lp(25)=.true.
      pronam( 26)='Energy potential                            kJ/mol'
      lp(26)=.true.
      pronam( 27)='Energy kinetic                              kJ/mol'
      lp(27)=.true.
      pronam( 28)='Energy kinetic translations                 kJ/mol'
      pronam( 29)='Energy kinetic rotations                    kJ/mol'
      pronam( 30)='Energy kinetic translations * virial   (kJ/mol)**2'
      pronam( 31)='Energy kinetic rotations * virial      (kJ/mol)**2'
      pronam( 32)='Energy total                                kJ/mol'
      lp(32)=.true.
      pronam( 33)='Virial total                                kJ/mol'
      pronam( 34)='pV total                                    kJ/mol'
      pronam( 35)='Enthalpy                                    kJ/mol'
      pronam( 36)='Potential energy perturbation 1             kJ/mol'
      lp(36)=lpert2
      pronam( 37)='Potential energy perturbation 2             kJ/mol'
      lp(37)=lpert3
      pronam( 38)='Total energy perturbation 1                 kJ/mol'
      lp(38)=lpert2
      pronam( 39)='Total energy perturbation 2                 kJ/mol'
      lp(39)=lpert3
      pronam( 40)='Derivative free energy wrt lambda           kJ/mol'
      lp(40)=ntype.eq.3
      pronam( 41)='Bias corrected derivative                   kJ/mol'
      lp(41)=(ntype.eq.3).and.(nbias.gt.0)
      pronam( 42)='Boltzmann probability factor                      '
      pronam( 43)='System dipole moment                       e30 C m'
      pronam( 44)='Displacement by pressure scaling                nm'
      pronam( 45)='Finite system g-factor                            '
      pronam( 46)='Infinite system g-factor                          '
      pronam( 47)='Dielectric constant                               '
      pronam( 48)='Energy solvent restraints                   kJ/mol'
      pronam( 49)='Energy solute restraints                    kJ/mol'
      pronam( 50)='Energy solute umbrella potential            kJ/mol'
      lp(50)=nbias.gt.0
      pronam( 51)='Free energy solute umbrella potential       kJ/mol'
      pronam( 52)='Energy solvent intra electrostatic          kJ/mol'
      lp(52)=nwm.gt.0
      pronam( 53)='Energy solvent intra van der Waals          kJ/mol'
      lp(53)=nwm.gt.0
      pronam( 54)='Energy solvent inter electrostatic          kJ/mol'
      lp(54)=nwm.gt.0
      pronam( 55)='Energy solvent inter van der Waals          kJ/mol'
      lp(55)=nwm.gt.0
      pronam( 56)='Energy solvent bonds                        kJ/mol'
      lp(56)=nwm.gt.0
      pronam( 57)='Energy solvent angles                       kJ/mol'
      lp(58)=nwm.gt.0
      pronam( 58)='Energy solvent Urey-Bradley                 kJ/mol'
      lp(59)=nwm.gt.0
      pronam( 59)='Energy solvent proper dihedrals             kJ/mol'
      lp(60)=nwm.gt.0
      pronam( 60)='Energy solvent improper dihedrals           kJ/mol'
      lp(61)=nwm.gt.0
      pronam( 61)='Energy solvent polarization                 kJ/mol'
      pronam( 62)='Energy solvent external field               kJ/mol'
      pronam( 63)='Energy solvent reaction field               kJ/mol'
      pronam( 64)='Energy solvent kinetic translations         kJ/mol'
      pronam( 65)='Energy solvent kinetic                      kJ/mol'
      lp(65)=nwm.gt.0
      pronam( 66)='Energy solvent potential long range         kJ/mol'
      pronam( 67)='Energy solvent potential                    kJ/mol'
      pronam( 68)='Energy solvent total                        kJ/mol'
      pronam( 69)='Virial solvent                              kJ/mol'
      pronam( 70)='Mean square displacement solvent             nm**2'
      lp(70)=nwm.gt.0.and.ntype.gt.1.and.idifco.gt.0
      pronam( 71)='Diffusion constant solvent                  m**2/s'
      lp(71)=nwm.gt.0.and.ntype.gt.1.and.idifco.gt.0
      pronam( 72)='Dipole moment solvent permanent            e30 C m'
      pronam( 73)='Dipole moment solvent induced              e30 C m'
      pronam( 74)='Dipole moment solvent                      e30 C m'
      pronam( 75)='Energy restraint potentials                 kJ/mol'
      lp(75)=lnoe
      pronam( 76)='Energy umbrella potentials                  kJ/mol'
      lp(76)=lpmf
      if(ntype.ne.3) then
      pronam( 77)='Particle mesh Ewald reciprocal energy       kJ/mol'
      lp(77)=lpme
      pronam( 80)='Particle mesh Ewald correction energy       kJ/mol'
      lp(80)=lpme
      else
      pronam( 77)='Particle mesh Ewald reciprocal energy       kJ/mol'
      lp(77)=lpme
      pronam( 78)='PME reciprocal energy perturbation 2        kJ/mol'
      lp(78)=lpme
      pronam( 79)='PME reciprocal energy perturbation 3        kJ/mol'
      lp(79)=lpme
      pronam( 80)='Particle mesh Ewald correction energy       kJ/mol'
      lp(80)=lpme
      pronam( 81)='PME correction energy perturbation 2        kJ/mol'
      lp(81)=lpme
      pronam( 82)='PME correction energy perturbation 3        kJ/mol'
      lp(82)=lpme
      endif
      pronam( 83)='Energy distance restraints                  kJ/mol'
      pronam( 84)='Synchronization time                             s'
      lp( 84)=ntype.ge.1
      pronam( 85)='Efficiency load balancing                         '
      lp( 85)=ntype.ge.1
      pronam( 86)='Time per step                                    s'
      lp( 86)=ntype.ge.1
      pronam( 87)='Virial tensor element (1,1)                 kJ/mol'
      lp( 87)=.true.
      pronam( 88)='Virial tensor element (1,2)                 kJ/mol'
      lp( 88)=.true.
      pronam( 89)='Virial tensor element (1,3)                 kJ/mol'
      lp( 89)=.true.
      pronam( 90)='Virial tensor element (2,1)                 kJ/mol'
      lp( 90)=.true.
      pronam( 91)='Virial tensor element (2,2)                 kJ/mol'
      lp( 91)=.true.
      pronam( 92)='Virial tensor element (2,3)                 kJ/mol'
      lp( 92)=.true.
      pronam( 93)='Virial tensor element (3,1)                 kJ/mol'
      lp( 93)=.true.
      pronam( 94)='Virial tensor element (3,2)                 kJ/mol'
      lp( 94)=.true.
      pronam( 95)='Virial tensor element (3,3)                 kJ/mol'
      lp( 95)=.true.
      pronam( 96)='Pressure tensor element (1,1)                   Pa'
      lp( 96)=lpbc
      pronam( 97)='Pressure tensor element (1,2)                   Pa'
      lp( 97)=lpbc
      pronam( 98)='Pressure tensor element (1,3)                   Pa'
      lp( 98)=lpbc
      pronam( 99)='Pressure tensor element (2,1)                   Pa'
      lp( 99)=lpbc
      pronam(100)='Pressure tensor element (2,2)                   Pa'
      lp(100)=lpbc
      pronam(101)='Pressure tensor element (2,3)                   Pa'
      lp(101)=lpbc
      pronam(102)='Pressure tensor element (3,1)                   Pa'
      lp(102)=lpbc
      pronam(103)='Pressure tensor element (3,2)                   Pa'
      lp(103)=lpbc
      pronam(104)='Pressure tensor element (3,3)                   Pa'
      lp(104)=lpbc
      pronam(105)='Kinetic energy element (1,1)                kJ/mol'
      lp(105)=.true.
      pronam(106)='Kinetic energy element (1,2)                kJ/mol'
      lp(106)=.true.
      pronam(107)='Kinetic energy element (1,3)                kJ/mol'
      lp(107)=.true.
      pronam(108)='Kinetic energy element (2,1)                kJ/mol'
      lp(108)=.true.
      pronam(109)='Kinetic energy element (2,2)                kJ/mol'
      lp(109)=.true.
      pronam(110)='Kinetic energy element (2,3)                kJ/mol'
      lp(110)=.true.
      pronam(111)='Kinetic energy element (3,1)                kJ/mol'
      lp(111)=.true.
      pronam(112)='Kinetic energy element (3,2)                kJ/mol'
      lp(112)=.true.
      pronam(113)='Kinetic energy element (3,3)                kJ/mol'
      lp(113)=.true.
      pronam(114)='Box x-dimension                                 nm'
      lp(114)=lpbc
      pronam(115)='Box y-dimension                                 nm'
      lp(115)=lpbc
      pronam(116)='Box z-dimension                                 nm'
      lp(116)=lpbc
      pronam(117)='Free energy perturbation 1                  kJ/mol'
      pronam(118)='Free energy perturbation 2                  kJ/mol'
      pronam(119)='Free energy perturbation 1 including mass   kJ/mol'
      pronam(120)='Free energy perturbation 2 including mass   kJ/mol'
c
      isprop=120
      it=isprop
c
      if(nsf.gt.0) then
      it=0
      do 3 i=1,nsf
      it=isprop+(i-1)*30
      if(it+27.gt.maxpro) call md_abort('Increase dimension maxpro',0)
      pronam(it+1)='Energy solute     intra-polarization        kJ/mol'
      pronam(it+2)='Energy solute     intra-electrostatic       kJ/mol'
      lp(it+2)=.true.
      pronam(it+3)='Energy solute     intra-van der Waals       kJ/mol'
      lp(it+3)=.true.
      pronam(it+4)='Energy solute    -solvent polarization      kJ/mol'
      pronam(it+5)='Energy solute    -solvent electrostatic     kJ/mol'
      lp(it+5)=nwm.gt.0
      pronam(it+6)='Energy solute    -solvent van der Waals     kJ/mol'
      lp(it+6)=nwm.gt.0
      pronam(it+7)='Energy solute     bonds                     kJ/mol'
      lp(it+7)=.true.
      pronam(it+8)='Energy solute     angles                    kJ/mol'
      lp(it+8)=.true.
      pronam(it+9)='Energy solute     Urey-Bradley              kJ/mol'
      lp(it+9)=.true.
      pronam(it+10)='Energy solute     proper dihedrals          kJ/mol'
      lp(it+10)=.true.
      pronam(it+11)='Energy solute     improper dihedrals        kJ/mol'
      lp(it+11)=.true.
      pronam(it+12)='Energy solute     external field            kJ/mol'
      pronam(it+13)='Energy solute     reaction field            kJ/mol'
      pronam(it+14)='Energy solute     van der Waals 3rd         kJ/mol'
      lp(it+14)=.true.
      pronam(it+15)='Energy solute     electrostatic 3rd         kJ/mol'
      lp(it+15)=.true.
      pronam(it+16)='Energy solute     van der Waals nb          kJ/mol'
      lp(it+16)=.true.
      pronam(it+17)='Energy solute     electrostatic nb          kJ/mol'
      lp(it+17)=.true.
      pronam(it+18)='Energy solute     kinetic translations      kJ/mol'
      pronam(it+19)='Energy solute     kinetic                   kJ/mol'
      lp(it+19)=.true.
      pronam(it+20)='Energy solute     potential                 kJ/mol'
      lp(it+20)=.true.
      pronam(it+21)='Energy solute     total                     kJ/mol'
      pronam(it+22)='Virial solute                               kJ/mol'
      pronam(it+23)='Mean square displacement solute              nm**2'
      pronam(it+24)='Diffusion coefficient solute                m**2/s'
      pronam(it+25)='Dipole moment solute     permanent         e30 C m'
      pronam(it+26)='Dipole moment solute     induced           e30 C m'
      pronam(it+27)='Dipole moment solute                       e30 C m'
      pronam(it+28)='Moment of inertia solute                          '
      pronam(it+29)='Reaction field energy solute    -solvent    kJ/mol'
      pronam(it+30)='Reaction field energy solute                kJ/mol'
      do 4 j=1,21
      write(pronam(it+j)(15:17),'(i3)') i
    4 continue
      write(pronam(it+22)(33:35),'(i3)') i
      write(pronam(it+23)(30:32),'(i3)') i
      write(pronam(it+24)(22:24),'(i3)') i
      write(pronam(it+25)(22:24),'(i3)') i
      write(pronam(it+26)(22:24),'(i3)') i
      write(pronam(it+27)(26:28),'(i3)') i
      write(pronam(it+28)(30:32),'(i3)') i
      write(pronam(it+29)(30:32),'(i3)') i
    3 continue
      it=it+30
c
      if(nsf.gt.1) then
      if(it+nsf*(nsf-1)*2.gt.maxpro)
     + call md_abort('Increase dimension maxpro',0)
      do 5 i=1,nsf-1
      do 6 j=i+1,nsf
      it=it+1
      pronam(it)='Energy solute    -    polarization          kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      it=it+1
      pronam(it)='Energy solute    -    electrostatic         kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
      it=it+1
      pronam(it)='Energy solute    -    van der Waals         kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
      it=it+1
      pronam(it)='Energy solute    -    restraint             kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      it=it+1
      pronam(it)='Force  solute    -    x                  kJ/mol nm'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
      it=it+1
      pronam(it)='Force  solute    -    y                  kJ/mol nm'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
      it=it+1
      pronam(it)='Force  solute    -    z                  kJ/mol nm'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
      it=it+1
      pronam(it)='Force  solute    -                       kJ/mol nm'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      lp(it)=.true.
    6 continue
    5 continue
      endif
c
      if(msm.gt.0.and.npstat.gt.0) then
      if(it+2*nsm.gt.maxpro)
     + call md_abort('Increase dimension maxpro',0)
      do 7 i=1,nsm
      it=it+1
      pronam(it)='Moment of inertia solute molecule       kg nm2/mol'
      write(pronam(it)(34:38),'(i5)') i
      lp(it)=.true.
    7 continue
      do 8 i=1,nsm
      it=it+1
      pronam(it)='Radius of gyration solute molecule              nm'
      write(pronam(it)(34:38),'(i5)') i
      lp(it)=iradgy.gt.0
    8 continue
      do 88 i=1,nsm-1
      do 89 j=i+1,nsm
      it=it+1
      pronam(it)='Distance solute    -                            nm'
      write(pronam(it)(17:19),'(i3)') i
      write(pronam(it)(21:23),'(i3)') j
      lp(it)=.true.
   89 continue
   88 continue
      endif
c
      endif
c
      if(npmfi.gt.0) then
      if(it+npmfi.gt.maxpro)
     + call md_abort('Increase dimension maxpro',0)
      do 43 i=1,npmfi
      it=it+1
      pronam(it)='Energy pmf                                  kJ/mol'
      write(pronam(it)(12:16),'(i5)') i
      lp(it)=.true.
   43 continue
      endif
c
      if(iprof.eq.1) then
      if(it+42.gt.maxpro)
     + call md_abort('Increase dimension maxpro',0)
      it=it+1
      pronam(it)='Prof: Initialization velocities                  s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Periodic Boundary Conditions               s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Dynamic Load Balancing                     s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Atom Redistribution                        s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Center of Mass Coordinates                 s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Recording Trajectory                       s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: CAFE Initialization                        s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: QMD Forces Evaluation                      s'
      lp(it)=lqmd
      it=it+1
      pronam(it)='Prof: Synchronization                            s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: SPACE Initialization                       s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Synchronization                            s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Collapse Option Calculation                s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: External Field Contribution                s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Multiprocs Data Initialization             s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Multiprocs Data Communication              s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Multiprocs Data Forces Evaluation          s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Restraints Data Initialization             s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Restraints Data Communication              s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Restraints Data Forces Evaluation          s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: PMF Data Initialization                    s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: PMF Data Communication                     s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: PMF Data Forces Evaluation                 s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Induced Dipoles Evaluation                 s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: PME Initialization                         s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Charge Grid Evaluation                 s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Charge Grid Communication              s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Charge Grid Retrieval                  s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Reverse 3D-pFFT                        s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Energy Evaluation                      s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Synchronization                        s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Forward 3D-pFFT                        s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Synchronization                        s'
      lp(it)=lpme
      if(nbget.eq.0) then
      it=it+1
      pronam(it)='Prof: Local Coordinates Retrieval                s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Local Forces Evaluation                    s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Local Force Accumulation                   s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Local Coordinates Retrieval            s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Local Forces Evaluation                s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Local Force Accumulation               s'
      lp(it)=.true.
      else
      it=it+1
      pronam(it)='Prof: Prefetching coordinates                    s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Blocking Wait for Coordinates          s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Local Forces Evaluation                    s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Local Forces Evaluation                s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Forces Accumulation                        s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Non-Blocking Wait for Forces               s'
      lp(it)=.true.
      endif
      it=it+1
      pronam(it)='Prof: PME Wait                                   s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Forces Evaluation                      s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: PME Forces Communication                   s'
      lp(it)=lpme
      it=it+1
      pronam(it)='Prof: Node Synchronization                       s'
      lp(it)=.false.
      it=it+1
      pronam(it)='Prof: PME Flag                                   s'
      lp(it)=.false.
      it=it+1
      pronam(it)='Prof: Synchronization                            s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: SPACE Finalization                         s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: QM/MM Forces Evaluation                    s'
      lp(it)=lqmmm
      it=it+1
      pronam(it)='Prof: Forces Normalization                       s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Guided MD Forces                           s'
      lp(it)=lguide
      it=it+1
      pronam(it)='Prof: MD Time Step                               s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: SHAKE                                      s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: CAFE Finalization                          s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Center of Mass Motion Removal              s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: SPACE Update                               s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Property Evaluation                        s'
      lp(it)=.true.
      it=it+1
      pronam(it)='Prof: Restart File Update                        s'
      lp(it)=.true.
      endif
c
      nprop=0
      do 2 i=1,maxpro
      if(lp(i)) then
      nprop=nprop+1
      ixp(nprop)=i
      endif
    2 continue
c
      maxp=it
c
      call prp_init()
c
      mda=mdacq
      if(ntype.eq.3.and.me.eq.0) then
      if(.not.ma_push_get(mt_dbl,mda,'dfr',l_dfr,i_dfr))
     + call md_abort('Failed to allocate memory for dfr',0)
      if(.not.ma_push_get(mt_dbl,mda,'dfrm',l_dfrm,i_dfrm))
     + call md_abort('Failed to allocate memory for dfrm',0)
      endif
      lfreem=.false.
c
      if(nserie.gt.0) call prp_rdrest(lfn,fil)
c
      return
      end
      subroutine prp_init()
c
      implicit none
c
#include "prp_common.fh"
c
      integer i
c
      nsum=0
      nsumt=0
      nsump=0
      do 1 i=1,maxpro
      psum(i)=zero
      p2sum(i)=zero
      pslop(i)=zero
      psumt(i)=zero
      p2sumt(i)=zero
      pslopt(i)=zero
      psump(i)=zero
    1 continue
      tsum=zero
      t2sum=zero
      tsumt=zero
      t2sumt=zero
c
      nderiv=0
      dfree=zero
      do 2 i=1,24
      deriv(i)=zero
    2 continue
c
      lhdr=.false.
      lfhdr=.false.
      lauto=.false.
      lappr=.false.
      samrat=one
c
      return
      end
      subroutine prp_wtrest(lfn)
c
      implicit none
c
#include "prp_common.fh"
c
      integer lfn
c
      integer i
c
      if(me.ne.0) return
c
      write(lfn,1000)
 1000 format('restart properties')
      write(lfn,1001) maxp,nsum,nsumt
 1001 format(3i7)
      write(lfn,1002) tsum,t2sum,tsumt,t2sumt
 1002 format(4e20.12)
      write(lfn,1003) (psum(i),i=1,maxp)
      write(lfn,1003) (p2sum(i),i=1,maxp)
      write(lfn,1003) (pslop(i),i=1,maxp)
      write(lfn,1003) (psumt(i),i=1,maxp)
      write(lfn,1003) (p2sumt(i),i=1,maxp)
      write(lfn,1003) (pslopt(i),i=1,maxp)
 1003 format(4e20.12)
c
      return
      end
      subroutine prp_rdrest(lfn,fil)
c
      implicit none
c
#include "prp_common.fh"
c
      integer lfn
      character*255 fil
c
      character*18 string
      integer i
c
      if(me.ne.0) return
c
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      read(lfn,1000,end=9998) string
 1000 format(a18)
      if(string.ne.'restart properties') goto 1
c
      read(lfn,1001) maxp,nsum,nsumt
 1001 format(3i7)
      read(lfn,1002) tsum,t2sum,tsumt,t2sumt
 1002 format(4e20.12)
      read(lfn,1003) (psum(i),i=1,maxp)
      read(lfn,1003) (p2sum(i),i=1,maxp)
      read(lfn,1003) (pslop(i),i=1,maxp)
      read(lfn,1003) (psumt(i),i=1,maxp)
      read(lfn,1003) (p2sumt(i),i=1,maxp)
      read(lfn,1003) (pslopt(i),i=1,maxp)
 1003 format(4e20.12)
c
      return
 9998 continue
      call md_abort('Unable to read restart file in prp_rdrest ',me)
      return
 9999 continue
      call md_abort('Unable to open restart file in prp_rdrest ',me)
      return
      end
      subroutine prp_finish()
c
      implicit none
c
#include "prp_common.fh"
#include "mafdecls.fh"
c
      if(ntype.eq.3.and.me.eq.0) then
      if(.not.ma_pop_stack(l_dfrm))
     + call md_abort('Failed to deallocate dfrm',me)
      if(.not.ma_pop_stack(l_dfr))
     + call md_abort('Failed to deallocate dfr',me)
      endif
c
      return
      end
      subroutine prp_setup(wbx)
c
      implicit none
c
#include "prp_common.fh"
c
      real*8 wbx
c
      wbox=wbx
c
      return
      end
