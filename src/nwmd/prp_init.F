      subroutine prp_start(iout,iprp,nf1,nf2,nf3,
     + lcf,lsp,ltwn,lpol,lpme,
     + nwm,msfi,nsfi,mpei)
c
      implicit none
c
#include "prp_common.fh"
#include "global.fh"
c
      integer iout,iprp,nwm,msfi,nsfi,mpei,nf1,nf2,nf3
      logical lcf,lsp,ltwn,lpol,lpme
c
      integer i,j,it
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      lfnout=iout
      lfnprp=iprp
      lcafe=lcf
      lspac=lsp
      ltwin=ltwn
      msf=msfi
      nsf=nsfi
      mpe=mpei
      nfoutp=nf1
      nfstat=nf2
      nfprop=nf3
c
      do 1 i=1,250
      lp(i)=.false.
    1 continue
c
      pronam(  1)='Step                                              '
      pronam(  2)='Time                                            ps'
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
      lp(19)=.true.
      pronam( 20)='Temperature rotations                            K'
      lp(20)=.true.
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
      pronam( 36)='Energy perturbation 1                       kJ/mol'
      pronam( 37)='Energy perturbation 2                       kJ/mol'
      pronam( 38)='Free energy perturbation 1                  kJ/mol'
      pronam( 39)='Free energy perturbation 2                  kJ/mol'
      pronam( 40)='Derivative free energy wrt lambda           kJ/mol'
      pronam( 41)='H*dH/dL                                (kJ/mol)**2'
      pronam( 42)='Boltzmann probability factor                      '
      pronam( 43)='System dipole moment                       e30 C m'
      pronam( 44)='Displacement by pressure scaling                nm'
      pronam( 45)='Finite system g-factor                            '
      pronam( 46)='Infinite system g-factor                          '
      pronam( 47)='Dielectric constant                               '
      pronam( 48)='Energy solvent restraints                   kJ/mol'
      pronam( 49)='Energy solute restraints                    kJ/mol'
      pronam( 50)='Energy solute umbrella potential            kJ/mol'
      pronam( 51)='Free energy solute umbrella potential       kJ/mol'
      pronam( 52)='Energy solvent intra electrostatic          kJ/mol'
      pronam( 53)='Energy solvent intra van der Waals          kJ/mol'
      pronam( 54)='Energy solvent inter electrostatic          kJ/mol'
      lp(54)=nwm.gt.0
      pronam( 55)='Energy solvent inter van der Waals          kJ/mol'
      lp(55)=nwm.gt.0
      pronam( 56)='Energy solvent bonds                        kJ/mol'
      lp(56)=nwm.gt.0
      pronam( 57)='Energy solvent angles                       kJ/mol'
      lp(57)=nwm.gt.0
      pronam( 58)='Energy solvent proper dihedrals             kJ/mol'
      lp(58)=nwm.gt.0
      pronam( 59)='Energy solvent improper dihedrals           kJ/mol'
      lp(59)=nwm.gt.0
      pronam( 60)='Energy solvent polarization                 kJ/mol'
      pronam( 61)='Energy solvent external field               kJ/mol'
      pronam( 62)='Energy solvent reaction field               kJ/mol'
      pronam( 63)='Energy solvent kinetic translations         kJ/mol'
      pronam( 64)='Energy solvent kinetic                      kJ/mol'
      pronam( 65)='Energy solvent potential long range         kJ/mol'
      pronam( 66)='Energy solvent potential                    kJ/mol'
      pronam( 67)='Energy solvent total                        kJ/mol'
      pronam( 68)='Virial solvent                              kJ/mol'
      pronam( 69)='Mean square displacement solvent /6          nm**2'
      pronam( 70)='Diffusion constant solvent * e9             m**2/s'
      pronam( 71)='Dipole moment solvent permanent            e30 C m'
      pronam( 72)='Dipole moment solvent induced              e30 C m'
      pronam( 73)='Dipole moment solvent                      e30 C m'
      pronam( 74)='Energy bias potentials                      kJ/mol'
      pronam( 75)='Energy umbrella potentials                  kJ/mol'
      pronam( 76)='Particle mesh Ewald reciprocal energy       kJ/mol'
      pronam( 77)='Particle mesh Ewald correction energy       kJ/mol'
      pronam( 78)='Virial tensor element (1,1)                 kJ/mol'
      pronam( 79)='Virial tensor element (1,2)                 kJ/mol'
      pronam( 80)='Virial tensor element (1,3)                 kJ/mol'
      pronam( 81)='Virial tensor element (2,2)                 kJ/mol'
      pronam( 82)='Virial tensor element (2,3)                 kJ/mol'
      pronam( 83)='Virial tensor element (3,3)                 kJ/mol'
      pronam( 84)='Pressure tensor element (1,1)                   Pa'
      pronam( 85)='Pressure tensor element (1,2)                   Pa'
      pronam( 86)='Pressure tensor element (1,3)                   Pa'
      pronam( 87)='Pressure tensor element (2,2)                   Pa'
      pronam( 88)='Pressure tensor element (2,3)                   Pa'
      pronam( 89)='Pressure tensor element (3,3)                   Pa'
      pronam( 90)='Kinetic energy element (1,1)                kJ/mol'
      pronam( 91)='Kinetic energy element (1,2)                kJ/mol'
      pronam( 92)='Kinetic energy element (1,3)                kJ/mol'
      pronam( 93)='Kinetic energy element (2,2)                kJ/mol'
      pronam( 94)='Kinetic energy element (2,3)                kJ/mol'
      pronam( 95)='Kinetic energy element (3,3)                kJ/mol'
      pronam( 96)='Energy distance restraints                  kJ/mol'
c
      isprop=96
      if(nsf.gt.0) then
      do 3 i=1,nsf
      it=isprop+(i-1)*27
      pronam(it+1)='Energy solute     intra-polarization        kJ/mol'
      pronam(it+2)='Energy solute     intra-Coulombic           kJ/mol'
      pronam(it+3)='Energy solute     intra-Lennard-Jones       kJ/mol'
      pronam(it+4)='Energy solute    -solvent polarization      kJ/mol'
      pronam(it+5)='Energy solute    -solvent Coulombic         kJ/mol'
      lp(it+5)=nwm.gt.0
      pronam(it+6)='Energy solute    -solvent Lennard-Jones     kJ/mol'
      lp(it+6)=nwm.gt.0
      pronam(it+7)='Energy solute     bonds                     kJ/mol'
      lp(it+7)=.true.
      pronam(it+8)='Energy solute     angles                    kJ/mol'
      lp(it+8)=.true.
      pronam(it+9)='Energy solute     proper dihedrals          kJ/mol'
      lp(it+9)=.true.
      pronam(it+10)='Energy solute     improper dihedrals        kJ/mol'
      lp(it+10)=.true.
      pronam(it+11)='Energy solute     external field            kJ/mol'
      pronam(it+12)='Energy solute     reaction field            kJ/mol'
      pronam(it+13)='Energy solute     Coulombic                 kJ/mol'
      pronam(it+14)='Energy solute     Lennard-Jones             kJ/mol'
      pronam(it+15)='Energy solute     kinetic translations      kJ/mol'
      pronam(it+16)='Energy solute     kinetic                   kJ/mol'
      pronam(it+17)='Energy solute     potential                 kJ/mol'
      pronam(it+18)='Energy solute     total                     kJ/mol'
      pronam(it+19)='Virial solute                               kJ/mol'
      pronam(it+20)='Mean square displacement solute     /6       nm**2'
      pronam(it+21)='Diffusion coefficient solute *e9            m**2/s'
      pronam(it+22)='Dipole moment solute     permanent         e30 C m'
      pronam(it+23)='Dipole moment solute     induced           e30 C m'
      pronam(it+24)='Dipole moment solute                       e30 C m'
      pronam(it+25)='Moment of inertia solute                          '
      pronam(it+26)='Reaction field energy solute    -solvent    kJ/mol'
      pronam(it+27)='Reaction field energy solute                kJ/mol'
      do 4 j=1,19
      write(pronam(it+j)(15:17),'(i3)') i
    4 continue
      write(pronam(it+20)(33:35),'(i3)') i
      write(pronam(it+21)(30:32),'(i3)') i
      write(pronam(it+22)(22:24),'(i3)') i
      write(pronam(it+23)(22:24),'(i3)') i
      write(pronam(it+24)(22:24),'(i3)') i
      write(pronam(it+25)(26:28),'(i3)') i
      write(pronam(it+26)(30:32),'(i3)') i
      write(pronam(it+27)(30:32),'(i3)') i
    3 continue
      it=it+27
      if(nsf.gt.1) then
      do 5 i=1,nsf-1
      do 6 j=i+1,nsf
      it=it+1
      pronam(it)='Energy solute    -    polarization          kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      it=it+1
      pronam(it)='Energy solute    -    Coulombic             kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      it=it+1
      pronam(it)='Energy solute    -    Lennard-Jones         kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
      it=it+1
      pronam(it)='Energy solute    -    restraint             kJ/mol'
      write(pronam(it)(15:17),'(i3)') i
      write(pronam(it)(19:21),'(i3)') j
    6 continue
    5 continue
      endif
      endif
c
      nprop=0
      do 2 i=1,250
      if(lp(i)) then
      nprop=nprop+1
      ixp(nprop)=i
      endif
    2 continue
c
      nsum=0
      nsumt=0
      do 7 i=1,250
      psum(i)=zero
      p2sum(i)=zero
      pslop(i)=zero
      psumt(i)=zero
      p2sumt(i)=zero
      pslopt(i)=zero
    7 continue
      tsum=zero
      t2sum=zero
      tsumt=zero
      t2sumt=zero
c
      lhdr=.false.
c
      return
      end
