      subroutine dft_input(rtdb)
c     $Id: dft_input.F,v 1.4 1996-10-05 00:52:17 d3h449 Exp $
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
c     
c>>>  added for constant ntags_bsmx
c     
#include "bas_staticP.fh"
c     
c     read input for DFT calculation.  All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
      integer nopt,  nopt_quad, nopt_munit, nopt_grid,
     &     nopt_system
      parameter (nopt=57,  nopt_quad=2, nopt_munit=6,
     &     nopt_grid=5, nopt_system=4)
      character*255 test
      character*10 gridtype

c     
c>>>  declarations from NWdft
c     
      integer infinity ! for infinite loop
      parameter (infinity=100000) 
c     
c>>>  variables from various NWdft routines declared here (perhaps
c>>>  temporarily)
c     

      integer  n,
     &     iqrot, nagrid, nrmax,
     &     nrinc, idecomp, mult, 
     &     itrscf, ipol, noio,
     &     iguess,
     &     ndamp, iscfcon, igcon, idirect, nfock, ifield,
     &     numset, igeotype, numgrp,
     &     iAOacc,iXCacc,irqmax,itol2e,imull,iangquad,
     &     ncydp, ncysh, ncyds,
     &     nquad_task
      double precision rlshift, TOLL, radius, rm
      double precision dampon, dampoff, levlon, levloff, 
     &                 diison, diisoff
      integer nfield, iline
      logical pseudo
      logical nodamping, nolevelshifting, nodiis
      logical ltight, lmaxov, incore
      logical store_wght, ldelley
      integer nrad(ntags_bsmx)
      double precision rscale,  cdist(3), cang(3)
c     
c     XC functionals coefficients
c     
      character*50 xname,cname
      double precision cfac(10),xfac(10)
c--   >  defaults from rdcalc.F
c     
      ipol=1
C
C     noio is the default
C
      noio=1
      iguess=0
      cang(1)=0.d0
      cang(2)=0.d0
      cang(3)=0.d0
      cdist(1)=0.d0
      cdist(2)=0.d0
      cdist(3)=0.d0
      mult=1
C
C     Convergence aids defaults:
C
      nodamping = .false.
      nolevelshifting = .false.
      nodiis = .false.
      dampon = 1.d8
      dampoff = 1.d-1
      levlon = 1.d-1
      levloff = 1.d-4
      diison = 1.d-1
      diisoff = 0.d0
      ncydp = 0
      ncysh = 0
      ncyds = 0
      ndamp=70
      rlshift=0.5d0
c      icvgtyp=2
c
      ltight=.false.
      incore=.true.
      lmaxov=.false.      
      itrscf=30
      igcon=-1
      iscfcon=7
      pseudo=.false.
      nfock=10
      idirect=0
c do NOT turn off diis unless specified by user
c      idisoff=15
C
C
      iAOacc=20
      iXCacc=20
      irqmax=40
C
C     bi-eletronic tolerance of 12 used as a default
C
      itol2e=12
      imull=0
      iangquad=0
      rm = 2.0d0
      TOLL=1.D-15
      radius=16.d0
      nquad_task = 1
      ldelley = .false.
      store_wght = .false.
      do iline=1,10
        xfac(iline)=0.d0
        cfac(iline)=0.d0
      enddo
c     
c     default XC: vwn and slater
c     
      cfac(1)=1.d0
      xfac(2)=1.d0
      xname='Slater            Exchange    Functional'
      cname='Vosko-Wilk-Nusair Correlation Functional'
c     
c>>>  sym input defaults
c     
      numgrp = 1 ! sym group P1 or C1
      numset = 1 ! space group setting number
      igeotype = 0 ! no periodic dimensions
c     
c>>>  defaults from rdprop.F
c     
      idecomp = 0
c     
c     set default grid to  'MEDIUM'
c     
      gridtype='medium'
      nagrid = 8
      do n = 1,ntags_bsmx
        nrad(n) = 40
      enddo 
c     
c>>>  read inputs; slam them in the rtdb
c     
      do iline = 1,infinity !  infinite loop over input lines

        if (.not. inp_read()) call errquit
     $       ('dft_input: inp_read failed', 0)

        nfield = inp_n_field()
        do ifield = 1,nfield

          if (inp_a(test)) then
c           
            if ( inp_compare(.false., 'dft', test)) then
c             
c>>>          dft
c             
              ipol=1
            elseif ( inp_compare(.false., 'odft', test)) then
c             
c>>>          odft
c             
              ipol=2
            elseif ( inp_compare(.false., 'guess', test)) then
c             
c>>>          guess; next keyword gives type of initial guess
c             
              if (inp_a(test)) then
                if ( inp_compare(.false., 'atomic', test) ) then
                  iguess = 0
                elseif ( inp_compare(.false., 'hcore', test) ) then
                  iguess = 1
                else 
                  write(6,'(a,a,a)') ' dft_input: guess option ',test,
     &                 ' not recognized'
                endif 
              else 
                write (6,*)
     &               ' dft_input:  WARNING:',
     &               'guess option not found; ignored'
                goto 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'hcore', test)) then
c             
c>>>          atomic
c             
              iguess = 0
            elseif ( inp_compare(.false., 'atomic', test)) then
c             
c>>>          hcore
c             
              iguess = 1
            elseif ( inp_compare(.false., 'mult', test)) then
c             
c>>>          mult
c             
              if ( .not. inp_i(mult) ) then
                write(6,*)
     &               ' dft_input: spin multiplicity value not found'
                goto 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'charg', test)) then
c             
c>>>          charg
c             
c****         
c****         this keyword has been disabled
c****         
              call errquit(
     &             'dft_input: keyword charg no longer supported'
     &             ,0)

            elseif ( inp_compare(.false., 'itrscf', test)) then
c             
c>>>          itrscf
c             
              if ( .not. inp_i(itrscf) ) then
                write(6,*) ' dft_input: itrscf value not found'
                goto 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'direct', test)) then
c             
c>>>          direct
c             
              idirect = 1
            elseif ( inp_compare(.false., 'scfcon', test)) then
c             
c>>>          scfcon
c             
              if ( .not. inp_i(iscfcon) ) then
                write(6,*) ' dft_input: iscfcon value not found'
                goto 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'gcon', test)) then
c             
c>>>          gcon
c             
              if ( .not. inp_i(igcon) ) then
                write(6,*) ' dft_input: igcon value not found'
                goto 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'ecp', test)) then
c             
c>>>          ecp
c             
              pseudo = .true.
c
c>>>        Next are input parameters to define convergence aids.
c           
            elseif ( inp_compare(.false., 'dampon', test)) then
c             
c>>>          dampon
c             
              if ( .not. inp_f(dampon) ) then
                write(6,*)
     &               ' dft_input: dampon value not found; defaulting ',
     &               ' to 1.d8 Hartree'
                dampon = 1.d8
              endif            
            elseif ( inp_compare(.false., 'dampoff', test)) then
c             
c>>>          dampoff
c             
              if ( .not. inp_f(dampoff) ) then
                write(6,*)
     &               ' dft_input: dampoff value not found; defaulting ',
     &               ' to 0.0d0 Hartree'
                dampoff = 0.0d0
              endif            
            elseif ( inp_compare(.false., 'levlon', test)) then
c             
c>>>          levlon
c             
              if ( .not. inp_f(levlon) ) then
                write(6,*)
     &               ' dft_input: levlon value not found; defaulting ',
     &               ' to 1.d8 Hartree'
                levlon = 1.d8
              endif            
            elseif ( inp_compare(.false., 'levloff', test)) then
c             
c>>>          levloff
c             
              if ( .not. inp_f(levloff) ) then
                write(6,*)
     &               ' dft_input: levloff value not found; defaulting ',
     &               ' to 0.0d0 Hartree'
                levloff = 0.0d0
              endif            
            elseif ( inp_compare(.false., 'diison', test)) then
c             
c>>>          diison
c             
              if ( .not. inp_f(diison) ) then
                write(6,*)
     &               ' dft_input: diison value not found; defaulting ',
     &               ' to 1.d8 Hartree'
                diison = 1.d8
              endif            
            elseif ( inp_compare(.false., 'diisoff', test)) then
c             
c>>>          diisoff
c             
              if ( .not. inp_f(diisoff) ) then
                write(6,*)
     &               ' dft_input: diisoff value not found; defaulting ',
     &               ' to 0.0d0 Hartree'
                diisoff = 0.0d0
              endif            
            elseif ( inp_compare(.false., 'lshift', test)) then
c             
c>>>          lshift
c             
              if ( .not. inp_f(rlshift) ) then
                write(6,*)
     &               ' dft_input: lshift value not found; defaulting ',
     &               ' to  0.5 Hartree'
                rlshift = 0.5
              endif            
            elseif ( inp_compare(.false., 'ncysh', test)) then
c             
c>>>          ncysh (# cycle level shifter is active)
c             
              if ( .not. inp_i(ncysh) ) then
                write(6,*) ' dft_input: ncysh value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'damp', test)) then
c             
c>>>          damp
c             
              if ( .not. inp_i(ndamp) ) then
                write(6,*) ' dft_input: ndamp value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'ncydp', test)) then
c             
c>>>          ncydp (# cycle damping is active)
c             
              if ( .not. inp_i(ncydp) ) then
                write(6,*) ' dft_input: ncydp value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'ncyds', test)) then
c             
c>>>          ncyds (# cycle DIIS is active)
c             
              if ( .not. inp_i(ncyds) ) then
                write(6,*) ' dft_input: ncyds value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'nodamping', test)) then
c             
c>>>          nodamping
c             
              nodamping = .true.
            elseif ( inp_compare(.false., 'nolevelshifting', test)) then
c             
c>>>          nolevelshifting
c             
              nolevelshifting = .true.
            elseif ( inp_compare(.false., 'nodiis', test)) then
c             
c>>>          nodiis
c             
              nodiis = .true.
c            elseif ( inp_compare(.false., 'idiisoff', test)) then
c             
c>>>          idiisoff
c             
c              icvgtyp=2
c              if ( .not. inp_i(idisoff) ) then
c                write(6,*) ' dft_input: idisoff value not found'
c                goto 10000 !  die with error
c              endif            
            elseif ( inp_compare(.false., 'nfock', test)) then
c             
c>>>          nfock
c             
              if ( .not. inp_i(nfock) ) then
                write(6,*) ' dft_input: nfock value not found'
                goto 10000 !  die with error
              endif            
c            elseif ( inp_compare(.false., 'diis', test)) then
c             
c>>>          nodiis
c             
c              icvgtyp=0
c              nfock=0
            elseif ( inp_compare(.false., 'noio', test)) then
c             
c>>>          noio
c             
              noio = 1
            elseif ( inp_compare(.false., 'decomp', test)) then
c             
c>>>          decomp
c             
              idecomp = 1
            elseif ( inp_compare(.false., 'grid', test)) then
c             
c>>>          grid
c             
c             
c>>>          nrad is indexed by atom type
c             
              if (inp_a(test)) then
                if ( test .eq. 'coarse' ) then ! coarse
                  gridtype='coarse'
                  do n = 1,ntags_bsmx
                    nrad(n) = 30
                  enddo 
                  nagrid = 6 ! was 10 in rdquad
                elseif ( test .eq. 'medium' ) then ! medium
                  gridtype='medium'
                  do n = 1,ntags_bsmx
                    nrad(n) = 40
                  enddo 
                  nagrid = 8
                elseif ( test .eq. 'fine' ) then ! fine
                  gridtype='fine'
                  do n = 1,ntags_bsmx
                    nrad(n) = 60
                  enddo 
                  nagrid = 10
                elseif ( test .eq. 'xfine' ) then ! xfine
                  gridtype='xfine'
                  do n = 1,ntags_bsmx
                    nrad(n) = 64
                  enddo 
                  nagrid = 24
                elseif ( test .eq. 'user' ) then ! user
                  write(6,*)
     &                 ' dft_input: grid option user not implemented'
                else 
                  write(6,*) ' dft_input: grid type not recognized'
                endif 
              else 
                write(6,*) ' dft_input: grid type not found'
              endif 
            elseif ( inp_compare(.false., 'coarse', test)) then
c             
c>>>          coarse
c             
              gridtype=test
              do n = 1,ntags_bsmx
                nrad(n) = 30
              enddo 
              nagrid = 6
            elseif ( inp_compare(.false., 'medium', test)) then
c             
c>>>          medium
c             
              gridtype=test
              do n = 1,ntags_bsmx
                nrad(n) = 40
              enddo 
              nagrid = 8
            elseif ( inp_compare(.false., 'fine', test)) then
c             
c>>>          fine
c             
              gridtype=test
              do n = 1,ntags_bsmx
                nrad(n) = 60
              enddo 
              nagrid = 10
            elseif ( inp_compare(.false., 'xfine', test)) then
c             
c>>>          xfine
c             
              gridtype=test
              do n = 1,ntags_bsmx
                nrad(n) = 90
              enddo 
              nagrid = 10
            elseif ( inp_compare(.false., 'user', test)) then
c             
c>>>          user
c             
              write(6,*)
     &             ' dft_input: grid option user not implemented'

c             
c>>>          next are inputs for translational symmetry routines
c             
c             numgrp:    Space group number.    Default=1 (P1)
c             numset:    Setting number for space group.  Default=1
c             igeotype:  Defines number of periodic dimensions in system.
c             0=Molecule, 1=Polymer, 2=Surface, 3=Crystal
c             Default = -1.  If not changed by user -1 flags abort.
c             cang(3):   Angles for cell parameters.
c             cdist(3):  Lattice constants.
c             
c>>>          keywords are:  system molecule polymer surface crystal 
c>>>          group setting A B C alpha beta gamma
c             
            elseif ( inp_compare(.false., 'system', test)) then
c             
c>>>          system 
c             
              if (inp_a(test)) then
                if ( test .eq. 'molecule' ) then ! molecule
                  igeotype = 0
                elseif ( test .eq. 'polymer' ) then ! polymer
                  igeotype = 1
                elseif ( test .eq. 'surface' ) then ! surface
                  igeotype = 2
                elseif ( test .eq. 'crystal' ) then ! crystal
                  igeotype = 3
                endif 
              else 
                write(6,*) ' dft_input: system type not recognized'
                go to 10000 !  die with error
              endif 
            elseif ( inp_compare(.false., 'molecule', test)) then
c             
c>>>          molecule 
c             
              igeotype = 0
            elseif ( inp_compare(.false., 'polymer', test)) then
c             
c>>>          polymer 
c             
              igeotype = 1
            elseif ( inp_compare(.false., 'surface', test)) then
c             
c>>>          surface 
c             
              igeotype = 2
            elseif ( inp_compare(.false., 'crystal', test)) then
c             
c>>>          crystal 
c             
              igeotype = 3
            elseif ( inp_compare(.false., 'group', test)) then
c             
c>>>          group 
c             
              if ( .not. inp_i(numgrp) ) then
                write(6,*) ' dft_input: numgrp value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'setting', test)) then
c             
c>>>          setting 
c             
              if ( .not. inp_i(numset) ) then
                write(6,*) ' dft_input: numset value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'A', test)) then
c             
c>>>          A (first lattice constant)
c             
              if ( .not. inp_f(cdist(1)) ) then
                write(6,*) ' dft_input: A value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'B', test)) then
c             
c>>>          B 
c             
            elseif ( inp_compare(.false., 'C', test)) then
              if ( .not. inp_f(cdist(2)) ) then
                write(6,*) ' dft_input: B value not found'
                goto 10000 !  die with error
              endif            
c             
c>>>          C 
c             
              if ( .not. inp_f(cdist(3)) ) then
                write(6,*) ' dft_input: C value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'alpha', test)) then
c             
c>>>          alpha 
c             
              if ( .not. inp_f(cang(1)) ) then
                write(6,*) ' dft_input: alpha value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'beta', test)) then
c             
c>>>          beta 
c             
              if ( .not. inp_f(cang(2)) ) then
                write(6,*) ' dft_input: beta value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'gamma', test)) then
c             
c>>>          gamma
c             
              if ( .not. inp_f(cang(3)) ) then
                write(6,*) ' dft_input: gamma value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'slater', test)) then
c             
c>>>          Slater f for Xchang
c             
              xname='Slater            Exchange    Functional'
              xfac(2)=1.d0
            elseif ( inp_compare(.false., 'vwn', test)) then
c             
c>>>          VWN f for Corr
c             
              cname='Vosko-Wilk-Nusair Correlation Functional'
              cfac(1)=1.d0
            elseif ( inp_compare(.false., 'becke88', test)) then
c             
c>>>          Becke88 f for Xchang
c             
              xname='Becke 1988        Exchange    Functional'
              xfac(2)=0.d0
              xfac(3)=1.d0
            elseif ( inp_compare(.false., 'lyp', test)) then
c             
c>>>          LYP f for Corr
c             
              cname='Lee-Yang-Parr     Correlation Functional'
              cfac(1)=0.d0
              cfac(2)=1.d0
            elseif ( inp_compare(.false., 'HFexch', test)) then
c             
c>>>          Exact Exchange (HF) and no corr (could change with
c             a subsequent keyword for correlation)
c             
              xname='Hartree-Fock (Exact) Exchange'
              cname=' No Correlation Potential'
              xfac(2)=0.d0
              xfac(1)=1.d0
              cfac(1)=0.d0
              
            elseif ( inp_compare(.false., 'perdew86', test)) then
c             
c>>>          P86 f for Corr
c             
              cname='Perdew 1986       Correlation Functional'
              cfac(1)=0.d0
              cfac(3)=1.d0
              cfac(4)=1.d0
            elseif ( inp_compare(.false., 'perdew91', test)) then
c             
c>>>          P91 f for Corr
c             
              cname='Perdew 1991       Correlation Functional'
              cfac(1)=0.d0
              cfac(5)=1.d0
            elseif ( inp_compare(.false., 'pw91lda', test)) then
c             
c>>>          P91LDA f for Corr
c             
              cname='Perdew 1991 LDA    Correlation Functional'
              cfac(1)=0.d0
              cfac(6)=1.d0
            elseif ( inp_compare(.false., 'vwnrpa', test)) then
c             
c>>>          VWNRPA for Corr
c             
              cname='VWN RPA    Correlation Functional'
              cfac(1)=0.d0
              cfac(7)=1.d0
            elseif ( inp_compare(.false., 'perdew81', test)) then
c             
c>>>          P81 f for Corr
c             
              cname='Perdew 1981       Correlation Functional'
              cfac(1)=0.d0
              cfac(3)=1.0d0
            elseif ( inp_compare(.false., 'acm', test)) then
c             
c>>>          Adiabatic connection method
c             AD Becke JCP 98, 5648 (1993)
c             Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + C*Ec(P91) 
c             A=0.20, B=0.72, C=0.81               
c             
              xname = 'Adiabatic Connection Method for XC potential'
              cname = 
     &        'EXC=EX(HF)+EXC(Slater&VWN)+EX(Becke88)+EC(Perdew91)'
              cfac(1)=1.00d0
              cfac(5)=0.81d0
              xfac(1)=0.20d0
              xfac(2)=0.80d0
              xfac(3)=0.72d0
              
            elseif ( inp_compare(.false., 'b3lyp', test)) then
c             
c>>>          Becke3LYP
c             AD Becke JCP 98, 5648 (1993)
c             Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + (1-C)*Ec(VWN)
c             + C*Ec(VWN) 
c             A=0.20, B=0.72, C=0.81               
c             
              xname = 'B3LYP'
              cname = 
     &             'EXC=EX(HF)+EXC(Slater&VWN)+EX(Becke88)+EC(LYP)'
              cfac(1)=0.00d0
              cfac(7)=0.19d0
              cfac(2)=0.81d0
              xfac(1)=0.20d0
              xfac(2)=0.80d0
              xfac(3)=0.72d0
              
            elseif ( inp_compare(.false., 'beckehandh', test)) then
c             
c>>>          Becke Half and Half
c             AD Becke JCP 98, 1372
c             
c             0.5*Exc(LSD) + 0.5*Ex(HF) 
c             
              xname = 'Becke half and half Exchange potential'
              cname = 'EXC=0.5*EX(HF)+0.5*EX(Slater+Becke88)'
              xfac(1)=0.5d0
              xfac(2)=0.5d0
              xfac(3)=0.5d0
            elseif ( inp_compare(.false., 'accAOfunc', test)) then
c             
c>>>          for use in int_eval_gbsets (called by quadvxcs)
c             acc for evaluating AO gaussians
c             
              if ( .not. inp_i(iAOacc) ) then
                write(6,*) ' dft_input: accAOfunc value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'accXCfunc', test)) then
c             
c>>>          for use in int_eval_gbsets (called by quadvxc1)
c             acc for evaluating XC gaussians
c             
              if ( .not. inp_i(iXCacc) ) then
                write(6,*) ' dft_input: accXCfunc value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'accQrad', test)) then
c             
c>>>          for use in setquad
c             r_max for quadrature points
c             
              if ( .not. inp_i(irqmax) ) then
                write(6,*) ' dft_input: accQrad value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'accCoul', test)) then
c             
c>>>          for use in fitcd fitvcoul
c             for schwarz screening
c             
              if ( .not. inp_i(itol2e) ) then
                write(6,*) ' dft_input: accCoul value not found'
                goto 10000 !  die with error
              endif            
            elseif ( inp_compare(.false., 'tight', test)) then
c             
c>>>          Set screening threshholds tight
c             
              ltight = .true.

            elseif ( inp_compare(.false., 'max_ovl', test)) then
c             
c>>>          Move orbitals around in movecs to maintain max 
c             overlap between iterations
c             
              lmaxov = .true.

            elseif ( inp_compare(.false., 'incore', test)) then
c             
c>>>          Try to do an incore coulomb construction
c             
              incore = .true.

            elseif ( inp_compare(.false., 'mulliken', test)) then
c             
c>>>          do Mulliken pop analysis
c             
              imull=1
            elseif ( inp_compare(.false., 'gausleg', test)) then
c             
c>>>          use Gauss-Legendre for angular quadr (Handy)
c             
              iangquad=0
            elseif ( inp_compare(.false., 'lebedev', test)) then
c             
c>>>          use Lebedev angular quadrature
c             
              if ( .not. inp_i(iangquad) ) iangquad=4
c
c a few not very well know options which were buried in the XC subs
c
            elseif ( inp_compare(.false., 'store_wght', test)) then
c             
              store_wght = .true.
c
            elseif ( inp_compare(.false., 'delley', test)) then
c             
              ldelley = .true.
c
            elseif ( inp_compare(.false., 'nquad_task', test)) then
c             
              if ( .not. inp_i(nquad_task) ) then
                write(6,*)
     &               ' dft_input: nquad_task value not found; ',
     &               ' defaulting to 1 '
                nquad_task = 1
              endif            
c
            elseif ( inp_compare(.false., 'tol_rho', test)) then
c             
c>>>          TOLL
c             
              if ( .not. inp_f(TOLL) ) then
                write(6,*)
     &               ' dft_input: tol_rho value not found; defaulting ',
     &               ' to 1.d-15 '
                TOLL=1.D-15
              endif            
c
            elseif ( inp_compare(.false., 'rm', test)) then
c             
c>>>          rm
c             
              if ( .not. inp_f(rm) ) then
                write(6,*)
     &               ' dft_input: rm value not found; defaulting ',
     &               ' to 2.0d0 '
                rm = 2.0d0
              endif            
c
            elseif ( inp_compare(.false., 'radius', test)) then
c             
c>>>          radius
c             
              if ( .not. inp_f(radius) ) then
                write(6,*)
     &               ' dft_input: radius value not found; defaulting ',
     &               ' to 16.d0 '
                radius=16.d0
              endif            
c
            else if (inp_compare(.false.,'print',test)) then
              call util_print_input(rtdb, 'dft')
            else if (inp_compare(.false.,'noprint',test)) then
              call util_print_input(rtdb, 'dft')
            elseif ( inp_compare(.false., 'end', test)) then
c             
c>>>          end 
c             
              goto 99999
            else
              write(6,'(a,a,a)') ' dft_input: keyword ',test,
     &             ' not recognized'
              goto 10000 !  die with error
            endif 
          else           ! no further field found in line
            goto 88 ! next line
          endif          ! terminates block    if (inp_a(test)) then
        enddo            ! end of loop over input fields 
   88   continue 
      enddo              ! end of infinite loop over input lines
99999 continue 
c     
c      if(icvgtyp.ne.2) nfock=0
      nrinc = 1
      nrmax = 0
      do n = 1,ntags_bsmx
        nrmax = max(nrmax,nrad(n))
      enddo 
c     Set the switch to rotate the quadrature points. Eventually this switch
c     will control whether the user wants the grid rotationally invariant.
      iqrot = 0
c     Output these parameters.
      rscale = 0.0d0
c     
c>>>  place logical parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:nodamping', mt_log, 1, nodamping))
     $     call errquit('dft_input: rtdb_put failed', 1)
      if (.not. rtdb_put(rtdb, 'dft:nolevelshifting', mt_log, 1, 
     $  nolevelshifting))call errquit('dft_input: rtdb_put failed', 2)
      if (.not. rtdb_put(rtdb, 'dft:nodiis', mt_log, 1, nodiis))
     $     call errquit('dft_input: rtdb_put failed', 3)
      if (.not. rtdb_put(rtdb, 'dft:tight', mt_log, 1, ltight))
     $     call errquit('dft_input: rtdb_put failed', 4)
      if(.not.rtdb_put(rtdb,'dft:incore', mt_log, 1, incore)) 
     $     call errquit('dft_input: rtdb_put failed', 5)
      if(.not. rtdb_put(rtdb, 'dft:max_ovl', mt_log, 1,lmaxov ))
     $     call errquit('dft_input: rtdb_put failed', 6)
      if(.not. rtdb_put(rtdb, 'dft:store_wght', mt_log, 1,store_wght ))
     $     call errquit('dft_input: rtdb_put failed', 7)
      if(.not. rtdb_put(rtdb, 'dft:delley', mt_log, 1,ldelley ))
     $     call errquit('dft_input: rtdb_put failed', 8)
c     
c>>>  place integer parameters in rtdb
c     
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     $     call errquit('dft_input: rtdb_put failed', 9) 
      if (.not. rtdb_cput(rtdb, 'dft:cname', 1, cname))
     $     call errquit('dft_input: rtdb_put failed', 10) 
      if (.not. rtdb_cput(rtdb, 'dft:xname', 1, xname))
     $     call errquit('dft_input: rtdb_put failed', 11) 
      if (.not. rtdb_put(rtdb, 'dft:itrscf', mt_int, 1, itrscf))
     $     call errquit('dft_input: rtdb_put failed', 12) 
      if (.not. rtdb_put(rtdb, 'dft:iscfcon', mt_int, 1, iscfcon))
     $     call errquit('dft_input: rtdb_put failed', 13)
      if (.not. rtdb_put(rtdb, 'dft:iguess', mt_int, 1, iguess))
     $     call errquit('dft_input: rtdb_put failed', 14)
      if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
     $     call errquit('dft_input: rtdb_put failed', 15)
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     $     call errquit('dft_input: rtdb_put failed', 16)
      if (.not. rtdb_put(rtdb, 'dft:igcon', mt_int, 1, igcon))
     $     call errquit('dft_input: rtdb_put failed', 17)
c      if (.not. rtdb_put(rtdb, 'dft:icvgtyp', mt_int, 1, icvgtyp))
c     $     call errquit('dft_input: rtdb_put failed', 16)
      if (.not. rtdb_put(rtdb, 'dft:nfock', mt_int, 1, nfock))
     $     call errquit('dft_input: rtdb_put failed', 18)
      if (.not. rtdb_put(rtdb, 'dft:ncydp', mt_int, 1, ncydp))
     $     call errquit('dft_input: rtdb_put failed', 19)
      if (.not. rtdb_put(rtdb, 'dft:ncysh', mt_int, 1, ncysh))
     $     call errquit('dft_input: rtdb_put failed', 20)
      if (.not. rtdb_put(rtdb, 'dft:ncyds', mt_int, 1, ncyds))
     $     call errquit('dft_input: rtdb_put failed', 21)
      if (.not. rtdb_put(rtdb, 'dft:idirect', mt_int, 1, idirect))
     $     call errquit('dft_input: rtdb_put failed', 22)
c      if (.not. rtdb_put(rtdb, 'dft:idisoff', mt_int, 1, idisoff))
c     $     call errquit('dft_input: rtdb_put failed', 22)
      if (.not. rtdb_put(rtdb, 'dft:noio', mt_int, 1, noio))
     $     call errquit('dft_input: rtdb_put failed', 23)
      if (.not. rtdb_put(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     $     call errquit('dft_input: rtdb_put failed', 24)
      if (.not. rtdb_put(rtdb, 'dft:idecomp', mt_int, 1, idecomp))
     $     call errquit('dft_input: rtdb_put failed', 25)
      if (.not. rtdb_put(rtdb, 'dft:nrad', mt_int, ntags_bsmx, nrad))
     $     call errquit('dft_input: rtdb_put failed', 26)
      if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     $     call errquit('dft_input: rtdb_put failed', 27)
      if (.not. rtdb_put(rtdb, 'dft:iqrot', mt_int, 1, iqrot))
     $     call errquit('dft_input: rtdb_put failed', 28)
      if (.not. rtdb_put(rtdb, 'dft:nrinc', mt_int, 1, nrinc))
     $     call errquit('dft_input: rtdb_put failed', 29)
      if (.not. rtdb_put(rtdb, 'dft:nrmax', mt_int, 1, nrmax))
     $     call errquit('dft_input: rtdb_put failed', 30)
      if (.not. rtdb_put(rtdb, 'dft:igeotype', mt_int, 1, igeotype))
     $     call errquit('dft_input: rtdb_put failed', 31)
      if (.not. rtdb_put(rtdb, 'dft:numgrp', mt_int, 1, numgrp))
     $     call errquit('dft_input: rtdb_put failed', 32)
      if (.not. rtdb_put(rtdb, 'dft:numset', mt_int, 1, numset))
     $     call errquit('dft_input: rtdb_put failed', 33)
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1,nquad_task))
     $     call errquit('dft_input: rtdb_put failed', 34)
      if (.not. rtdb_put(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))
     $     call errquit('dft_input: rtdb_put failed', 35)
      if (.not. rtdb_put(rtdb, 'dft:iXCacc', mt_int, 1, iXCacc))
     $     call errquit('dft_input: rtdb_put failed', 36)
      if (.not. rtdb_put(rtdb, 'dft:irqmax', mt_int, 1, irqmax))
     $     call errquit('dft_input: rtdb_put failed', 37)
      if (.not. rtdb_put(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     $     call errquit('dft_input: rtdb_put failed', 38)
      if (.not. rtdb_put(rtdb, 'dft:imull', mt_int, 1, imull))
     $     call errquit('dft_input: rtdb_put failed', 39)
      if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     $     call errquit('dft_input: rtdb_put failed', 40)
c     
c>>>  store real parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:dampon', mt_dbl, 1, dampon))
     $     call errquit('dft_input: rtdb_put failed', 41)
      if (.not. rtdb_put(rtdb, 'dft:dampoff', mt_dbl, 1, dampoff))
     $     call errquit('dft_input: rtdb_put failed', 42)
      if (.not. rtdb_put(rtdb, 'dft:levlon', mt_dbl, 1, levlon))
     $     call errquit('dft_input: rtdb_put failed', 43)
      if (.not. rtdb_put(rtdb, 'dft:levloff', mt_dbl, 1, levloff))
     $     call errquit('dft_input: rtdb_put failed', 44)
      if (.not. rtdb_put(rtdb, 'dft:diison', mt_dbl, 1, diison))
     $     call errquit('dft_input: rtdb_put failed', 45)
      if (.not. rtdb_put(rtdb, 'dft:diisoff', mt_dbl, 1, diisoff))
     $     call errquit('dft_input: rtdb_put failed', 46)
      if (.not. rtdb_put(rtdb, 'dft:rscale', mt_dbl, 1, rscale))
     $     call errquit('dft_input: rtdb_put failed', 47)
      if (.not. rtdb_put(rtdb, 'dft:cang', mt_dbl, 3, cang))
     $     call errquit('dft_input: rtdb_put failed', 48)
      if (.not. rtdb_put(rtdb, 'dft:cdist', mt_dbl, 3, cdist))
     $     call errquit('dft_input: rtdb_put failed', 49)
      if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, 10, cfac))
     $     call errquit('dft_input: rtdb_put failed', 50)
      if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, 10, xfac))
     $     call errquit('dft_input: rtdb_put failed', 51)
      if (.not. rtdb_put(rtdb, 'dft:lshift', mt_dbl, 1,rlshift))
     $     call errquit('dft_input: rtdb_put failed', 52)
      if (.not. rtdb_put(rtdb, 'dft:radius', mt_dbl, 1,radius))
     $     call errquit('dft_input: rtdb_put failed', 53)
      if (.not. rtdb_put(rtdb, 'dft:tol_rho', mt_dbl, 1,TOLL))
     $     call errquit('dft_input: rtdb_put failed', 54)
      if (.not. rtdb_put(rtdb, 'dft:rm', mt_dbl, 1,rm))
     $     call errquit('dft_input: rtdb_put failed', 55)
c     if (.not. rtdb_put(rtdb, 'dft:iEx', mt_int, 1, iEx))
c     $     call errquit('dft_input: rtdb_put failed', 46)
c     if (.not. rtdb_put(rtdb, 'dft:iEc', mt_int, 1, iEc))
c     $     call errquit('dft_input: rtdb_put failed', 47)
c
      return
c       
10000 continue !  die with error
      write (6,*) 'dft_input: error in input file, last field: '
      write (6,'(a)') test
      call errquit('dft_input: error in input file', 0)
c       
      end

