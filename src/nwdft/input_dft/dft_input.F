      subroutine dft_input(rtdb)
c
c     $Id: dft_input.F,v 1.12 1997-02-13 01:07:35 d3g681 Exp $
c
c     The main purpose of this subroutine is to parse the input
c     and place any "found" parameters into the rtdb.
c
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
c
c     The following parameters must match those in
c     bas_staticP.fh e.g., dft_ntags_bsmx = ntags_bsmx.
c
      integer dft_ntags_bsmx
      parameter (dft_ntags_bsmx = 20)
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
     &     nagrid,
     &     idecomp, mult, 
     &     itrscf, ipol, noio,
     &     iguess,
     &     ndamp, iscfcon, igcon, idirect, nfock, ifield,
     &     iAOacc,iXCacc,irqmax,itol2e,imull,iangquad,
     &     ncydp, ncysh, ncyds,
     &     nquad_task
      double precision rlshift, TOLL, radius, rm, hl_tol
      double precision dampon, dampoff, levlon, levloff, 
     &                 diison, diisoff
      integer nfield, iline
      logical pseudo, status
      logical nodamping, nolevelshifting, nodiis
      logical levelshift, damp, diis, direct
      logical ltight, lmaxov, incore
      logical store_wght, ldelley
c
c     Define logical flags for local (l) and nonlocal (nl)
c     parts of exchange-correlation functionals.
c
      logical  lcfac(15),  lxfac(15)
      logical nlcfac(15), nlxfac(15)
      logical xccomb(15)
c
      integer nrad(dft_ntags_bsmx)
c     
c     XC functionals coefficients
c     
      double precision cfac(15),xfac(15)
c
c     Any DFT input causes calculation to be marked as not converged
c
      status = rtdb_delete(rtdb,'dft:converged')
c
c     Only these few arrays are initialized (not defaults).
c
      do iline=1,15
          xfac(iline)=0.d0
          cfac(iline)=0.d0
         lcfac(iline)=.false.
        nlcfac(iline)=.false.
         lxfac(iline)=.false.
        nlxfac(iline)=.false.
        xccomb(iline)=.false.
      enddo
c     
c>>>  Read inputs; slam them in the rtdb.
c     
      do iline = 1,infinity !  infinite loop over input lines

        if (.not. inp_read()) call errquit
     $       ('dft_input: inp_read failed', 0)

        nfield = inp_n_field()
        do ifield = 1,nfield
c
          if (inp_a(test)) then
c           
            if ( inp_compare(.false., 'dft', test)) then
c             
c>>>          dft
c             
              ipol=1
              if (.not. rtdb_put(rtdb, 'dft:ipol',
     &           mt_int, 1, ipol))
     &           call errquit('dft_input: rtdb_put failed', 105)
c
            elseif ( inp_compare(.false., 'odft', test)) then
c             
c>>>          odft
c             
              ipol=2
              if (.not. rtdb_put(rtdb, 'dft:ipol',
     &           mt_int, 1, ipol))
     &           call errquit('dft_input: rtdb_put failed', 105)
c
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
              if (.not. rtdb_put(rtdb, 'dft:iguess',
     &           mt_int, 1, iguess))
     &           call errquit('dft_input: rtdb_put failed', 103)
c
            elseif ( inp_compare(.false., 'mult', test)) then
c             
c>>>          mult
c             
              if ( .not. inp_i(mult) ) then
                write(6,*)
     &               ' dft_input: spin multiplicity value not found'
                goto 10000 !  die with error
              endif 
              if (.not. rtdb_put(rtdb, 'dft:mult', 
     &           mt_int, 1, mult))
     &           call errquit('dft_input: rtdb_put failed', 104)
c
            elseif ( inp_compare(.false., 'itrscf', test)) then
c             
c>>>          itrscf
c             
              if ( .not. inp_i(itrscf) ) then
                write(6,*) ' dft_input: itrscf value not found'
                goto 10000 !  die with error
              endif 
              if (.not. rtdb_put(rtdb, 'dft:itrscf',
     &           mt_int, 1, itrscf))
     &           call errquit('dft_input: rtdb_put failed', 101) 
c
            elseif ( inp_compare(.false., 'direct', test)) then
c             
c>>>          direct
c             
              direct = .true.
              if (.not. rtdb_put(rtdb, 'dft:direct',
     &           mt_log, 1, direct))
     &           call errquit('dft_input: rtdb_put failed', 3)
              idirect = 1
              if (.not. rtdb_put(rtdb, 'dft:idirect',
     &           mt_int, 1, idirect))
     &           call errquit('dft_input: rtdb_put failed', 111)
c
            elseif ( inp_compare(.false., 'scfcon', test)) then
c             
c>>>          scfcon
c             
              if ( .not. inp_i(iscfcon) ) then
                write(6,*) ' dft_input: iscfcon value not found'
                goto 10000 !  die with error
              endif 
              if (.not. rtdb_put(rtdb, 'dft:iscfcon',
     &           mt_int, 1, iscfcon))
     &           call errquit('dft_input: rtdb_put failed', 102)
c
            elseif ( inp_compare(.false., 'gcon', test)) then
c             
c>>>          gcon
c             
              if ( .not. inp_i(igcon) ) then
                write(6,*) ' dft_input: igcon value not found'
                goto 10000 !  die with error
              endif 
              if (.not. rtdb_put(rtdb, 'dft:igcon',
     &           mt_int, 1, igcon))
     &           call errquit('dft_input: rtdb_put failed', 106)
c
            elseif ( inp_compare(.false., 'ecp', test)) then
c             
c>>>          ecp
c             
              pseudo = .true.
c
c>>>        Next are input parameters to define convergence aids.
c           
            elseif ( inp_compare(.false., 'hl_tol', test)) then
c             
c>>>          hl_tol
c             
              if ( .not. inp_f(hl_tol) ) then
                write(6,*)
     &               ' dft_input: hl_tol value not found; defaulting ',
     &               ' to 0.1 Hartree'
                hl_tol = 0.1d0
              endif            
              if (.not. rtdb_put(rtdb, 'dft:hl_tol',
     &           mt_dbl, 1, hl_tol))
     &           call errquit('dft_input: rtdb_put failed', 201)
c
            elseif ( inp_compare(.false., 'dampon', test)) then
c             
c>>>          dampon
c             
              if ( .not. inp_f(dampon) ) then
                write(6,*)
     &               ' dft_input: dampon value not found; defaulting ',
     &               ' to 0.d0 Hartree'
                dampon = 0.0d0
              endif            
              if (.not. rtdb_put(rtdb, 'dft:dampon',
     &           mt_dbl, 1, dampon))
     &           call errquit('dft_input: rtdb_put failed', 201)
c
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
              if (.not. rtdb_put(rtdb, 'dft:dampoff',
     &           mt_dbl, 1, dampoff))
     &           call errquit('dft_input: rtdb_put failed', 202)
c
            elseif ( inp_compare(.false., 'levlon', test)) then
c             
c>>>          levlon
c             
              if ( .not. inp_f(levlon) ) then
                write(6,*)
     &               ' dft_input: levlon value not found; defaulting ',
     &               ' to 0.d0 Hartree'
                levlon = 0.0d0
              endif            
              if (.not. rtdb_put(rtdb, 'dft:levlon',
     &           mt_dbl, 1, levlon))
     &           call errquit('dft_input: rtdb_put failed', 203)
c
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
              if (.not. rtdb_put(rtdb, 'dft:levloff',
     &           mt_dbl, 1, levloff))
     &           call errquit('dft_input: rtdb_put failed', 204)
c
            elseif ( inp_compare(.false., 'diison', test)) then
c             
c>>>          diison
c             
              if ( .not. inp_f(diison) ) then
                write(6,*)
     &               ' dft_input: diison value not found; defaulting ',
     &               ' to 0.d0 Hartree'
                diison = 0.0d0
              endif            
              if (.not. rtdb_put(rtdb, 'dft:diison',
     &           mt_dbl, 1, diison))
     &           call errquit('dft_input: rtdb_put failed', 205)
c
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
              if (.not. rtdb_put(rtdb, 'dft:diisoff',
     &           mt_dbl, 1, diisoff))
     &           call errquit('dft_input: rtdb_put failed', 206)
c
            elseif ( inp_compare(.false., 'lshift', test)) then
c             
c>>>          lshift
c             
              levelshift = .true.
              if (.not. rtdb_put(rtdb, 'dft:levelshift',
     &           mt_log, 1, levelshift))
     &           call errquit('dft_input: rtdb_put failed', 3)
              if ( .not. inp_f(rlshift) ) then
                write(6,*)
     &               ' dft_input: lshift value not found; defaulting ',
     &               ' to  0.5 Hartree'
                rlshift = 0.5
              endif            
              if (.not. rtdb_put(rtdb, 'dft:rlshift',
     &           mt_dbl, 1,rlshift))
     &           call errquit('dft_input: rtdb_put failed', 212)
c
            elseif ( inp_compare(.false., 'ncysh', test)) then
c             
c>>>          ncysh (# cycle level shifter is active)
c             
              if ( .not. inp_i(ncysh) ) then
                write(6,*) ' dft_input: ncysh value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:ncysh',
     &           mt_int, 1, ncysh))
     &           call errquit('dft_input: rtdb_put failed', 109)
c
            elseif ( inp_compare(.false., 'damp', test)) then
c             
c>>>          damp
c
              damp = .true.             
              if (.not. rtdb_put(rtdb, 'dft:damp',
     &           mt_log, 1, damp))
     &           call errquit('dft_input: rtdb_put failed', 3)
              if ( .not. inp_i(ndamp) ) then
                write(6,*)
     &               ' dft_input: ndamp value not found; defaulting ',
     &               ' to 70%.'
                ndamp = 70
              endif            
              if (.not. rtdb_put(rtdb, 'dft:ndamp',
     &           mt_int, 1, ndamp))
     &           call errquit('dft_input: rtdb_put failed', 113)
c
            elseif ( inp_compare(.false., 'ncydp', test)) then
c             
c>>>          ncydp (# cycle damping is active)
c             
              if ( .not. inp_i(ncydp) ) then
                write(6,*) ' dft_input: ncydp value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:ncydp',
     &           mt_int, 1, ncydp))
     &           call errquit('dft_input: rtdb_put failed', 108)
c
            elseif ( inp_compare(.false., 'ncyds', test)) then
c             
c>>>          ncyds (# cycle DIIS is active)
c             
              if ( .not. inp_i(ncyds) ) then
                write(6,*) ' dft_input: ncyds value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:ncyds',
     &           mt_int, 1, ncyds))
     &           call errquit('dft_input: rtdb_put failed', 110)
c
            elseif ( inp_compare(.false., 'nodamping', test)) then
c             
c>>>          nodamping
c             
              nodamping = .true.

              if (.not. rtdb_put(rtdb, 'dft:nodamping',
     &           mt_log, 1, nodamping))
     &           call errquit('dft_input: rtdb_put failed', 1)
c
            elseif ( inp_compare(.false., 'nolevelshifting', test)) then
c             
c>>>          nolevelshifting
c             
              nolevelshifting = .true.
              if (.not. rtdb_put(rtdb, 'dft:nolevelshifting',
     &           mt_log, 1, nolevelshifting))
     &           call errquit('dft_input: rtdb_put failed', 2)
c
            elseif ( inp_compare(.false., 'diis', test)) then
c             
c>>>          diis
c             
              diis = .true.
              if (.not. rtdb_put(rtdb, 'dft:diis',
     &           mt_log, 1, diis))
     &           call errquit('dft_input: rtdb_put failed', 3)
c
            elseif ( inp_compare(.false., 'nodiis', test)) then
c             
c>>>          nodiis
c             
              nodiis = .true.
              if (.not. rtdb_put(rtdb, 'dft:nodiis',
     &           mt_log, 1, nodiis))
     &           call errquit('dft_input: rtdb_put failed', 3)
c
            elseif ( inp_compare(.false., 'nfock', test)) then
c             
c>>>          nfock
c             
              if ( .not. inp_i(nfock) ) then
                write(6,*) ' dft_input: nfock value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:nfock',
     &           mt_int, 1, nfock))
     &           call errquit('dft_input: rtdb_put failed', 107)
c
            elseif ( inp_compare(.false., 'noio', test)) then
c             
c>>>          noio
c             
              noio = 1
              if (.not. rtdb_put(rtdb, 'dft:noio',
     &           mt_int, 1, noio))
     &           call errquit('dft_input: rtdb_put failed', 112)
c
            elseif ( inp_compare(.false., 'decomp', test)) then
c             
c>>>          decomp
c             
              idecomp = 1
              if (.not. rtdb_put(rtdb, 'dft:idecomp',
     &           mt_int, 1, idecomp))
     &           call errquit('dft_input: rtdb_put failed', 114)
c
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
                  do n = 1,dft_ntags_bsmx
                    nrad(n) = 30
                  enddo 
                  nagrid = 6 ! was 10 in rdquad
                elseif ( test .eq. 'medium' ) then ! medium
                  gridtype='medium'
                  do n = 1,dft_ntags_bsmx
                    nrad(n) = 40
                  enddo 
                  nagrid = 8
                elseif ( test .eq. 'fine' ) then ! fine
                  gridtype='fine'
                  do n = 1,dft_ntags_bsmx
                    nrad(n) = 60
                  enddo 
                  nagrid = 10
                elseif ( test .eq. 'xfine' ) then ! xfine
                  gridtype='xfine'
                  do n = 1,dft_ntags_bsmx
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
              if (.not. rtdb_put(rtdb,'dft:nrad',
     &           mt_int, dft_ntags_bsmx, nrad))
     &           call errquit('dft_input: rtdb_put failed', 115)
              if (.not. rtdb_put(rtdb, 'dft:nagrid',
     &           mt_int, 1, nagrid))
     &           call errquit('dft_input: rtdb_put failed', 116)
              if (.not. rtdb_cput(rtdb, 'dft:gridtype', 
     &           1, gridtype))
     &           call errquit('dft_input: rtdb_put failed', 301) 
c
c  Input for correlation functionals read here.
c
            elseif ( inp_compare(.false., 'vwn', test)) then
c             
c>>>          VWN f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn local descriptor found.'
                  lcfac(1) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn nonlocal descriptor found.'
                  write(6,*)' vwn is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(1) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(1) = .true.
                  nlcfac(1) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(1) = .true.
                nlcfac(1) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(1)) ) then
c                write(6,*)
c     &          ' dft_input: vwn multiplicative factor not found.'
                cfac(1) = 1.0d0
              endif
            elseif ( inp_compare(.false., 'lyp', test)) then
c             
c>>>          LYP f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: lyp local descriptor found.'
                 write(6,*)' The lyp functional cannot be conveniently',
     &                      ' split into local and nonlocal pieces; '
                  write(6,*)' total lyp is assumed desired. '                   
                   lcfac(2) = .true.
                  nlcfac(2) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: lyp nonlocal descriptor found.'
                 write(6,*)' The lyp functional cannot be conveniently',
     &                      ' split into local and nonlocal pieces; '
                  write(6,*)' total lyp is assumed desired. '                   
                   lcfac(2) = .true.
                  nlcfac(2) = .true.
                else
c
c                 assume total wanted
c
                  lcfac(2) = .true.
                  nlcfac(2) = .true.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(2) = .true.
                nlcfac(2) = .true.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(2)) ) then
c                write(6,*)
c     &          ' dft_input: lyp multiplicative factor not found.'
                cfac(2) = 1.0d0
              endif            
            elseif ( inp_compare(.false., 'perdew81', test)) then
c             
c>>>          P81 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew81 local descriptor found.'
                  lcfac(3) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew81 nonlocal descriptor found.'
                  write(6,*)' perdew81 is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(3) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(3) = .true.
                  nlcfac(3) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(3) = .true.
                nlcfac(3) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(3)) ) then
c                write(6,*)
c     &          ' dft_input: perdew81 multiplicative factor not found.'
                cfac(3) = 1.0d0
              endif            
            elseif ( inp_compare(.false., 'perdew86', test)) then
c             
c>>>          P86 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew86 local descriptor found.'
                  write(6,*)' perdew86 is a nonlocal functional;',
     &                      ' adding perdew81 local functional. '                   
                  lcfac(3) = .true.
                  lcfac(4) = .false.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew86 nonlocal descriptor found.'
                  nlcfac(3) = .false.
                  nlcfac(4) = .true.
                else
c
c                 assume total wanted
c
                  write(6,*)' perdew86 is a nonlocal functional;',
     &                      ' adding perdew81 local functional. '                   
                   lcfac(3) = .true.
                  nlcfac(3) = .false.
                   lcfac(4) = .false.
                  nlcfac(4) = .true.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                  write(6,*)' perdew86 is a nonlocal functional;',
     &                      ' adding perdew81 local functional. '                   
                 lcfac(3) = .true.
                nlcfac(3) = .false.
                 lcfac(4) = .false.
                nlcfac(4) = .true.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(4)) ) then
c                write(6,*)
c     &          ' dft_input: perdew86 multiplicative factor not found.'
                cfac(4) = 1.0d0
                if(lcfac(3))cfac(3) = cfac(4)
              endif            
            elseif ( inp_compare(.false., 'perdew91', test)) then
c             
c>>>          P91 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew91 local descriptor found.'
                  write(6,*)' perdew91 is a nonlocal functional;',
     &                      ' adding pw91lda local functional. '
c
c                 This is also done by default inside xc_perdew91.F
c                 when lcfac(5) is .true.
c
                  lcfac(6) = .true.
                  lcfac(5) = .false.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: perdew91 nonlocal descriptor found.'
                  nlcfac(5) = .true.
                else
c
c                 assume total wanted
c
                   lcfac(5) = .true.
                  nlcfac(5) = .true.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                 lcfac(5) = .true.
                nlcfac(5) = .true.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(5)) ) then
c                write(6,*)
c     &          ' dft_input: perdew91 multiplicative factor not found.'
                cfac(5) = 1.0d0
                if(lcfac(6))cfac(6) = cfac(5)
              endif            
            elseif ( inp_compare(.false., 'pw91lda', test)) then
c             
c>>>          P91LDA f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: p91LDA local descriptor found.'
                  lcfac(6) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: p91LDA nonlocal descriptor found.'
                  write(6,*)' pw91lda is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(6) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(6) = .true.
                  nlcfac(6) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(6) = .true.
                nlcfac(6) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(6)) ) then
c                write(6,*)
c     &          ' dft_input: p91LDA multiplicative factor not found.'
                cfac(6) = 1.0d0
              endif            
            elseif ( inp_compare(.false., 'vwnrpa', test)) then
c             
c>>>          VWNRPA f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwnrpa local descriptor found.'
                  lcfac(7) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwnrpa nonlocal descriptor found.'
                  write(6,*)' vwnrpa is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(7) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(7) = .true.
                  nlcfac(7) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(7) = .true.
                nlcfac(7) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(7)) ) then
c                write(6,*)
c     &          ' dft_input: vwnrpa multiplicative factor not found.'
                cfac(7) = 1.0d0
              endif            
c
            elseif ( inp_compare(.false., 'vwn_1', test)) then
c             
c>>>          VWN_1 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_1 local descriptor found.'
                  lcfac(8) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_1 nonlocal descriptor found.'
                  write(6,*)' vwn_1 is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(8) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(8) = .true.
                  nlcfac(8) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(8) = .true.
                nlcfac(8) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(8)) ) then
c                write(6,*)
c     &          ' dft_input: vwnrpa multiplicative factor not found.'
                cfac(8) = 1.0d0
              endif            
c
            elseif ( inp_compare(.false., 'vwn_2', test)) then
c             
c>>>          VWN_2 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_2 local descriptor found.'
                  lcfac(9) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_2 nonlocal descriptor found.'
                  write(6,*)' vwn_2 is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(9) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(9) = .true.
                  nlcfac(9) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(9) = .true.
                nlcfac(9) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(9)) ) then
c                write(6,*)
c     &          ' dft_input: vwnrpa multiplicative factor not found.'
                cfac(9) = 1.0d0
              endif            
c
            elseif ( inp_compare(.false., 'vwn_3', test)) then
c             
c>>>          VWN_3 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_3 local descriptor found.'
                  lcfac(10) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_3 nonlocal descriptor found.'
                  write(6,*)' vwn_3 is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(10) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(10) = .true.
                  nlcfac(10) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(10) = .true.
                nlcfac(10) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(10)) ) then
c                write(6,*)
c     &          ' dft_input: vwnrpa multiplicative factor not found.'
                cfac(10) = 1.0d0
              endif    
c
            elseif ( inp_compare(.false., 'vwn_4', test)) then
c             
c>>>          VWN_4 f for Corr
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_4 local descriptor found.'
                  lcfac(11) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: vwn_4 nonlocal descriptor found.'
                  write(6,*)' vwn_4 is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                  nlcfac(11) = .false.
                else
c
c                 assume total wanted
c
                  lcfac(11) = .true.
                  nlcfac(11) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lcfac(11) = .true.
                nlcfac(11) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(cfac(11)) ) then
c                write(6,*)
c     &          ' dft_input: vwnrpa multiplicative factor not found.'
                cfac(11) = 1.0d0
              endif    
c
c  Input for exchange functionals read here.
c
            elseif ( inp_compare(.false., 'HFexch', test)) then
c             
c>>>          HF exact exchange
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: hfexch local descriptor found.'
                  write(6,*)' HF exchange is by nature nonlocal',
     &                      ' and cannot be conveniently split'
                  write(6,*)' into local and nonlocal pieces;',
     &                      ' total HF exchange is assumed desired. '                   
                  lxfac(1) = .true.
                  nlxfac(1) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: hfexch nonlocal descriptor found.'
                  write(6,*)' HF exchange is by nature nonlocal',
     &                      ' and cannot be conveniently split'
                  write(6,*)' into local and nonlocal pieces;',
     &                      ' total HF exchange is assumed desired. '                   
                  lxfac(1) = .true.
                  nlxfac(1) = .true.
                else
c
c                 assume total wanted
c
                  lxfac(1) = .true.
                  nlxfac(1) = .true.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lxfac(1) = .true.
                nlxfac(1) = .true.
                call inp_prev_field()
              endif
              if ( .not. inp_f(xfac(1)) ) then
                write(6,*)
     &          ' dft_input: hfexch multiplicative factor not found.'
                xfac(1) = 1.0d0
              endif            
            elseif ( inp_compare(.false., 'slater', test)) then
c             
c>>>          Slater f for exchange
c             
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: slater local descriptor found.'
                  lxfac(2) = .true.
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: slater nonlocal descriptor found.'
                  nlxfac(2) = .false.
                  write(6,*)' slater is a local functional;',
     &                      ' "nonlocal" descriptor is ignored. '                   
                else
c
c                 assume total wanted
c
                  lxfac(2) = .true.
                  nlxfac(2) = .false.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lxfac(2) = .true.
                nlxfac(2) = .false.
                call inp_prev_field()
              endif
              if ( .not. inp_f(xfac(2)) ) then
c                write(6,*)
c     &          ' dft_input: slater multiplicative factor not found.'
                xfac(2) = 1.0d0
              endif            
            elseif ( inp_compare(.false., 'becke88', test)) then
c               
c>>>          B88 f for exchange
c               
              if ( inp_a(test) ) then
                if ( inp_compare(.false.,'local', test) ) then
c                  write(6,*)
c     &            ' dft_input: becke88 local descriptor found.'
                  write(6,*)' becke88 "local" is equivalent to the',
     &                      ' slater exchange local functional. '
                  write(6,*)' Resubmit using slater exchange. '
                  call errquit(' dft_input - functional description ',0)
                elseif ( inp_compare(.false.,'nonlocal', test) ) then
c                  write(6,*)
c     &            ' dft_input: becke88 nonlocal descriptor found.'
                  nlxfac(3) = .true.
                else
c
c                 assume total wanted
c
                  lxfac(3) = .true.
                  nlxfac(3) = .true.
                  call inp_prev_field()
                endif
              else
c
c               assume total wanted
c
                lxfac(3) = .true.
                nlxfac(3) = .true.
                call inp_prev_field()
              endif
              if ( .not. inp_f(xfac(3)) ) then
c                write(6,*)
c     &          ' dft_input: becke88 multiplicative factor not found.'
                xfac(3) = 1.0d0
              endif            
c
c  Input for special exchange-correlation combination functionals read here.
c
            elseif ( inp_compare(.false., 'acm', test)) then
c             
c>>>          Adiabatic connection method
c             AD Becke JCP 98, 5648 (1993)
c             Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + C*Ec(P91) 
c             A=0.20, B=0.72, C=0.81               
c             
              xccomb(1) = .true.
                cfac(1) = 1.00d0
               lcfac(1) = .true.
              nlcfac(1) = .false.
                cfac(5) = 0.81d0
               lcfac(5) = .false.
              nlcfac(5) = .true.
                xfac(1) = 0.20d0
               lxfac(1) = .true.
              nlxfac(1) = .true.
                xfac(2) = 0.80d0
               lxfac(2) = .true.
              nlxfac(2) = .false.
                xfac(3) = 0.72d0
               lxfac(3) = .false.
              nlxfac(3) = .true.
              
            elseif ( inp_compare(.false., 'b3lyp', test)) then
c             
c>>>          Becke3LYP
c             AD Becke JCP 98, 5648 (1993)
c             Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + (1-C)*Ec(VWN)
c             + C*Ec(VWN) 
c             A=0.20, B=0.72, C=0.81               
c             
              xccomb(2) = .true.
                cfac(2) = 0.81d0
               lcfac(2) = .true.
              nlcfac(2) = .true.
                cfac(7) = 0.19d0
               lcfac(7) = .true.
              nlcfac(7) = .false.
                xfac(1) = 0.20d0
               lxfac(1) = .true.
              nlxfac(1) = .true.
                xfac(2) = 0.80d0
               lxfac(2) = .true.
              nlxfac(2) = .false.
                xfac(3) = 0.72d0
               lxfac(3) = .false.
              nlxfac(3) = .true.
c              
            elseif ( inp_compare(.false., 'beckehandh', test)) then
c             
c>>>          Becke Half and Half
c             AD Becke JCP 98, 1372
c             
c             0.5*Exc(LSD) + 0.5*Ex(HF) 
c             
              xccomb(3) = .true.
                xfac(1) = 0.5d0
               lxfac(1) = .true.
              nlxfac(1) = .true.
                xfac(2) = 0.5d0
               lxfac(2) = .true.
              nlxfac(2) = .false.
                cfac(6) = 0.5d0
               lcfac(6) = .true.
              nlcfac(6) = .false.
c
            elseif ( inp_compare(.false., 'accAOfunc', test)) then
c             
c>>>          for use in int_eval_gbsets (called by quadvxcs)
c             acc for evaluating AO gaussians
c             
              if ( .not. inp_i(iAOacc) ) then
                write(6,*) ' dft_input: accAOfunc value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:iAOacc',
     &           mt_int, 1, iAOacc))
     &           call errquit('dft_input: rtdb_put failed', 124)
c
            elseif ( inp_compare(.false., 'accXCfunc', test)) then
c             
c>>>          for use in int_eval_gbsets (called by quadvxc1)
c             acc for evaluating XC gaussians
c             
              if ( .not. inp_i(iXCacc) ) then
                write(6,*) ' dft_input: accXCfunc value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:iXCacc',
     &           mt_int, 1, iXCacc))
     &           call errquit('dft_input: rtdb_put failed', 125)
c
            elseif ( inp_compare(.false., 'accQrad', test)) then
c             
c>>>          for use in setquad
c             r_max for quadrature points
c             
              if ( .not. inp_i(irqmax) ) then
                write(6,*) ' dft_input: accQrad value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:irqmax',
     &           mt_int, 1, irqmax))
     &           call errquit('dft_input: rtdb_put failed', 126)
c
            elseif ( inp_compare(.false., 'accCoul', test)) then
c             
c>>>          for use in fitcd fitvcoul
c             for schwarz screening
c             
              if ( .not. inp_i(itol2e) ) then
                write(6,*) ' dft_input: accCoul value not found'
                goto 10000 !  die with error
              endif            
              if (.not. rtdb_put(rtdb, 'dft:itol2e',
     &           mt_int, 1, itol2e))
     &           call errquit('dft_input: rtdb_put failed', 127)
c
            elseif ( inp_compare(.false., 'tight', test)) then
c             
c>>>          Set screening threshholds tight
c             
              ltight = .true.
              if (.not. rtdb_put(rtdb, 'dft:ltight',
     &           mt_log, 1, ltight))
     &           call errquit('dft_input: rtdb_put failed', 4)
c
            elseif ( inp_compare(.false., 'max_ovl', test)) then
c             
c>>>          Move orbitals around in movecs to maintain max 
c             overlap between iterations
c             
              lmaxov = .true.
              if(.not. rtdb_put(rtdb, 'dft:max_ovl',
     &          mt_log, 1, lmaxov))
     &          call errquit('dft_input: rtdb_put failed', 6)
c
            elseif ( inp_compare(.false., 'incore', test)) then
c             
c>>>          Try to do an incore coulomb construction
c             
              incore = .true.
              if(.not.rtdb_put(rtdb,'dft:incore',
     &          mt_log, 1, incore)) 
     &          call errquit('dft_input: rtdb_put failed', 5)
c
            elseif ( inp_compare(.false., 'mulliken', test)) then
c             
c>>>          do Mulliken pop analysis
c             
              imull=1
              if (.not. rtdb_put(rtdb, 'dft:imull',
     &           mt_int, 1, imull))
     &           call errquit('dft_input: rtdb_put failed', 128)
c
            elseif ( inp_compare(.false., 'gausleg', test)) then
c             
c>>>          use Gauss-Legendre for angular quadr (Handy)
c             
              iangquad=0
              if (.not. rtdb_put(rtdb, 'dft:iangquad',
     &           mt_int, 1, iangquad))
     &           call errquit('dft_input: rtdb_put failed', 129)
c
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
              if(.not. rtdb_put(rtdb, 'dft:store_wght',
     &          mt_log, 1, store_wght ))
     &          call errquit('dft_input: rtdb_put failed', 7)
c
            elseif ( inp_compare(.false., 'delley', test)) then
c             
              ldelley = .true.
              if(.not. rtdb_put(rtdb, 'dft:delley',
     &          mt_log, 1, ldelley ))
     &          call errquit('dft_input: rtdb_put failed', 8)
c
            elseif ( inp_compare(.false., 'nquad_task', test)) then
c             
              if ( .not. inp_i(nquad_task) ) then
                write(6,*)
     &               ' dft_input: nquad_task value not found; ',
     &               ' defaulting to 1 '
                nquad_task = 1
              endif            
              if (.not. rtdb_put(rtdb, 'dft:nquad_task',
     &           mt_int, 1,nquad_task))
     &           call errquit('dft_input: rtdb_put failed', 123)
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
              if (.not. rtdb_put(rtdb, 'dft:tol_rho',
     &           mt_dbl, 1,TOLL))
     &           call errquit('dft_input: rtdb_put failed', 214)
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
              if (.not. rtdb_put(rtdb, 'dft:rm',
     &           mt_dbl, 1,rm))
     &           call errquit('dft_input: rtdb_put failed', 215)
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
              if (.not. rtdb_put(rtdb, 'dft:radius',
     &           mt_dbl, 1,radius))
     &           call errquit('dft_input: rtdb_put failed', 213)
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
c     No alternative but to dump all array's associated
c     with type of XC into rtdb (whether specified or not).
c     
c>>>  place logical XC parameters in rtdb
c     
      if(.not. rtdb_put(rtdb, 'dft:lcfac',
     &   mt_log, 15, lcfac ))
     &   call errquit('dft_input: rtdb_put failed', 9)
      if(.not. rtdb_put(rtdb, 'dft:nlcfac',
     &   mt_log, 15, nlcfac ))
     &   call errquit('dft_input: rtdb_put failed', 10)
      if(.not. rtdb_put(rtdb, 'dft:lxfac',
     &   mt_log, 15, lxfac ))
     &   call errquit('dft_input: rtdb_put failed', 11)
      if(.not. rtdb_put(rtdb, 'dft:nlxfac',
     &   mt_log, 15, nlxfac ))
     &   call errquit('dft_input: rtdb_put failed', 12)
      if(.not. rtdb_put(rtdb, 'dft:xccomb',
     &   mt_log, 15, xccomb ))
     &   call errquit('dft_input: rtdb_put failed', 13)
c     
c>>>  place real XC parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:cfac',
     &   mt_dbl, 15, cfac))
     &   call errquit('dft_input: rtdb_put failed', 210)
      if (.not. rtdb_put(rtdb, 'dft:xfac',
     &   mt_dbl, 15, xfac))
     &   call errquit('dft_input: rtdb_put failed', 211)
      return
c       
10000 continue !  die with error
      write (6,*) 'dft_input: error in input file, last field: '
      write (6,'(a)') test
      call errquit('dft_input: error in input file', 0)
c       
      end

