      subroutine dft_inpana(rtdb)
c
c     $Id: dft_inpana.F,v 1.25 1997-03-31 23:46:43 d3h449 Exp $
c
c*********************************************************************
c     
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c     Write pertinent information to user output.
c     
c*********************************************************************
      implicit none
c
      logical even, cksetd
      integer noc(2)
      Integer rtdb              !  runtime database handle
      integer me,ichg,nel
      integer noc1
      integer n, na, nb
      integer nshells_cd, nshells_xc, nshells
      double precision anucl_charg, ckfac, t1, r1
      double precision EPSILON
      character*9 local_c, nonlocal_c, nineb_c
      character*3 on_off_1, on_off_2
      Parameter ( EPSILON=1.d-13 )
      logical oprint_general, oprint_grid, oprint_xc, 
     &     oprint_convergence, oprint_tolerances
      character*10 start_10c, NA_10c, asap_10c
      character*10 strng1, strng2, strng3, strng4, strng5, strng6
c     
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "geom.fh"
#include "bas.fh"
#include "cdft.fh"
#include "util.fh"
c     
      me=ga_nodeid()
c     
      oprint_general = util_print('general information',print_default)
      oprint_grid = util_print('grid information',print_default)
      oprint_xc = util_print('xc information',print_default)
      oprint_convergence = util_print('convergence information',
     &                                print_default)
      oprint_tolerances = util_print('screening tolerance information',
     &                               print_default)
c     
c     Determine number of electrons in system.
c     
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &     call errquit('scf: geom_nuc_charge failed', 0)
c     
      ichg = idint(1.000001d0 * rcharge) ! avoid bad truncation
      nel = int(anucl_charg+EPSILON) - ichg
c     
c     Check to see if calculation type is allowed.
c     
c     Even number of electrons required for RHF.
c     
      even=mod(nel,2).eq.0
c     
c     odd # of electrons or not a singlet state --> LSD
c     
      if ((.not.even).or.(mult.ne.1)) ipol=2
      
      noc(2)=0
c     
c     Calculate number of occupied orbitals.
c     
      if (ipol.eq.1)then
         noc1 = nel/2
         noc(2)= 0
         noc(1)= noc1
      else
c     
c        check consistency of no. elec and multiplicity
c     
         even=mod((nel+mult-1),2).eq.0
         if (.not.even) call errquit(
     &        ' no. of electrons and multiplicity not compatible',nel)
         noc(2) = (nel - mult + 1)/2
         noc(1) = nel - noc(2)
         noc1  = noc(1) + noc(2)
      endif
c     
c     Check to see if there are enough electrons for this
c     value of the multiplicity.
c     
      if (noc(2).lt.0)then
         call errquit('dft: #electrons not valid for multiplicity',mult)
      endif
c     
c     write noc (consistent with definition in ddscf) to rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('inpana: rtdb_put of noc failed', 0) 
c     
c     Write new data to checkpoint file.
c     
c     Analyze any user specified XC functionals ... set if none
c     
      cksetd = .true.
c     
c     Check if user has specified some type of functional,
c     if so, do not set defaults.
c     
      do n = 1, 15
         if (lcfac(n).or.nlcfac(n).or.lxfac(n).or.nlxfac(n))
     &        cksetd = .false.
      enddo
      if (cksetd)then
c     
c        Set functional defaults.
c     
         cfac(1) = 1.0d0
         lcfac(1) = .true.
         xfac(2) = 1.0d0
         lxfac(2) = .true.
c     
c        Update rtdb.
c     
         if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, 15, cfac))
     $        call errquit('dft_input: rtdb_put failed', 210)
         if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, 15, xfac))
     $        call errquit('dft_input: rtdb_put failed', 211)
         if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, 15, lcfac ))
     $        call errquit('dft_input: rtdb_put failed', 9)
         if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, 15, lxfac ))
     $        call errquit('dft_input: rtdb_put failed', 11)
      endif
c
c     Check/set defaults for convergence schemes
c     
c     Three types of convergence speedup schemes:
c     1) based on number of cycles performed,
c     2) based on differences of total energy less than some threshold,
c     3) DIIS plus levelshifting if homo-lumo gap is small.
c 
c     Examples of these might be:
c
c     1) based on number of cycles performed,
c     ncydp = 3
c     ndamp = 40
c     ncysh = iterations
c     rlshift = 0.5
c     ncyds = iterations        
c
c     2) based on differences of total energy being less than some threshold,
c     dampon = 1.d8
c     dampoff = 1.0d-1
c     levlon = 1.d-1
c     levloff = 1.0d-3
c     rlshift = 0.5
c     diison = 1.d-1
c     diisoff = 0.0d0
c
c     3) DIIS plus levelshifting if homo-lumo gap is small (default)
c     nodamping = .true.
c     ncydp = 0
c     ndamp = 0
c     ncysh = iterations
c     rlshift = 0.5
c     ncyds = iterations        
c     
      if(nodamping)then
         damp = .false.
         ncydp = 0
         ndamp = 0
      endif
      if(nolevelshifting)then
         levelshift = .false.
         ncysh = 0
         rlshift = 0.0
      endif
      if(nodiis)then
         diis = .false.
         ncyds = 0
      endif
      if (damp)then
c
c        check to make sure either number of damping iterations or
c        energy criterion has been specified.
c
         if (ncydp.eq.0.and.dampon.eq.0.0d0)then
            ncydp = iterations
         endif
      endif
      if (levelshift)then
c
c        check to make sure either number of levelshifting iterations or
c        energy criterion has been specified.
c
         if (ncysh.eq.0.and.levlon.eq.0.0d0)then
            ncysh = iterations
         endif
      endif
      if (diis)then
c
c        check to make sure either number of diis iterations or
c        energy criterion has been specified.
c
         if (ncyds.eq.0.and.diison.eq.0.0d0)then
            ncyds = iterations
         endif
      endif
c     
c     If convergence input based upon #cycles then turn off energy constraints.
c     
      if (ncydp.ne.0)then
         dampon  = -999.9
         dampoff = -999.9
      endif
      if (ncysh.ne.0)then
         levlon  = -999.9
         levloff = -999.9
      endif
      if (ncyds.ne.0)then
         diison  = -999.9
         diisoff = -999.9
      endif
c
c     check special case with damping - change
c     default of 2 to "iterations" if no other
c     convergence control specified
c
      if (ncysh.eq.0 .and. ncyds.eq.0 .and. ncydp.eq.2) 
     &   ncydp = iterations
c
      if (me.eq.0) then
c     
c        Write to output.
c
         if (oprint_general)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'General Information',20,.true.)
            write(LuOut,9020)
            if (ipol.eq.1)then
               write(LuOut,9050)
            elseif (ipol.eq.2)then
               write(LuOut,9055)
            endif
            write(LuOut,8150)ncenters
            na = noc(1)
            nb = noc(2)
            if (ipol.eq.1)nb = noc(1)
            write(LuOut,8200)nel,na,nb,ichg,mult
            on_off_1 = 'off'
            if (oskel)on_off_1 = 'on'
            on_off_2 = 'off'
            if (oadapt)on_off_2 = 'on'
            write(LuOut,9056)on_off_1, on_off_2
            write(LuOut,9030)iterations
            if (idirect.eq.1)write(LuOut,9110)
            if (.not. bas_numcont(AO_bas_han, nshells))
     &         call errquit('rdinput:rdinput:',86)
            write(LuOut,4001)nbf, nshells
            if (nbf_cd.gt.0)then
               write(LuOut,9130)
              if (.not. bas_numcont(CD_bas_han, nshells_cd))
     &           call errquit('rdinput:rdinput:',87)
               write(LuOut,4002)nbf_cd, nshells_cd
            endif
            if (nbf_xc.gt.0)then
               write(LuOut,9120)
               if (.not. bas_numcont(XC_bas_han, nshells_xc))
     &            call errquit('rdinput:rdinput:',88)
               write(LuOut,4003)nbf_xc, nshells_xc
            endif
            write(LuOut,9035)e_conv
            if (d_conv.gt.0)then
               write(LuOut,9040)d_conv
            endif
            if (g_conv.gt.0)then
               write(LuOut,9045)g_conv
            endif
         endif
         if (oprint_xc)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'XC Information',20,.true.)
c     
c           Write out XC info. Combo info first, than X components,
c           than C components.
c     
            local_c = 'local    '
            nonlocal_c = 'non-local'
            nineb_c = '         '
            do n = 1, 15
               if (xccomb(n))write(LuOut,9223) xcname(n)
            enddo
c     
c           Do exact exchange differently.
c     
            if (lxfac(1).or.nlxfac(1))
     &           write(LuOut,9224) xname(1), xfac(1), nineb_c
            do n = 2, 15
               if (lxfac(n).and.nlxfac(n))then
                  write(LuOut,9224) xname(n), xfac(n), nineb_c
               elseif (lxfac(n).and.(.not.nlxfac(n)))then
                  write(LuOut,9224) xname(n), xfac(n), local_c
               elseif ((.not.lxfac(n)).and.nlxfac(n))then
                  write(LuOut,9224) xname(n), xfac(n), nonlocal_c
               endif
            enddo
            do n = 1, 15
               if (lcfac(n).and.nlcfac(n))then
                  write(LuOut,9224) cname(n), cfac(n), nineb_c
               elseif (lcfac(n).and.(.not.nlcfac(n)))then
                  write(LuOut,9224) cname(n), cfac(n), local_c
               elseif ((.not.lcfac(n)).and.nlcfac(n))then
                  write(LuOut,9224) cname(n), cfac(n), nonlocal_c
               endif
            enddo
c     
c           Check XC coefficients to make sure appropriate components
c           sum to 1.0
c     
            ckfac = 0.0d0
            do n = 1, 15
               if (lcfac(n))ckfac = ckfac + cfac(n)
            enddo
            if (abs(ckfac-1.0d0).gt.1.d-8)then
               write(LuOut,*)
     &            ' WARNING: Sum of local correlation is ',ckfac
               write(LuOut,*)' Sum of components do not equal unity. '
            endif
            ckfac = 0.0d0
            do n = 1, 15
               if (nlcfac(n))ckfac = ckfac + cfac(n)
            enddo
            if (abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
               write(LuOut,*)
     &            ' WARNING: Sum of nonlocal correlation is ',ckfac
               write(LuOut,*)
     &            ' Sum of components do not equal unity or 0. '
            endif
            ckfac = 0.0d0
            do n = 1, 15
               if (lxfac(n))ckfac = ckfac + xfac(n)
            enddo
            if (abs(ckfac-1.0d0).gt.1.d-8)then
               write(LuOut,*)
     &            ' WARNING: Sum of local exchange is ',ckfac
               write(LuOut,*)' Sum of components do not equal unity. '
            endif
            ckfac = 0.0d0
            do n = 1, 15
               if (nlxfac(n))ckfac = ckfac + xfac(n)
            enddo
            if (abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
               write(LuOut,*)
     &            ' WARNING: Sum of nonlocal exchange is ',ckfac
               write(LuOut,*)
     &            ' Sum of components do not equal unity or 0. '
            endif
         endif
         if (oprint_grid)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'Grid Information',20,.true.)
            write(LuOut,9135)gridtype
            write(LuOut,9136)
            if (iangquad.eq.0)then
              write(LuOut,9137)nradpts,nagrid,2*nagrid
            else
              write(LuOut,9138)nradpts,nang_leb_pts
            endif
            write(LuOut,4004)nqshells
            if (ldelley)then
               write(LuOut,9139)
            else
               write(LuOut,9140)
            endif
            if (nquad_task.ne.1)then
               write(LuOut,9141)
            endif
         endif
         if (oprint_convergence)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'Convergence Information',20,.true.)
            write(LuOut,3231)hl_tol, nfock
            write(LuOut,3232)ndamp, rlshift
            asap_10c  = '  ASAP    '
            start_10c = '  start   '
            NA_10c = '   N/A    '
            if(ncydp.ne.0)then
               strng1 = start_10c
               write(strng4,'(i3,7h iters )')ncydp
            elseif(nodamping)then
               strng1 = NA_10c
               strng4 = NA_10c
            else
               write(strng1,'(d10.2)')dampon
               write(strng4,'(d10.2)')dampoff
            endif
c     
            if(ncysh.ne.0)then
               strng2 = asap_10c
               write(strng5,'(i3,7h iters )')ncysh
            elseif(nolevelshifting)then
               strng2 = NA_10c
               strng5 = NA_10c
            else
               write(strng2,'(d10.2)')levlon
               write(strng5,'(d10.2)')levloff
            endif
c     
            if(ncyds.ne.0)then
               strng3 = start_10c
               write(strng6,'(i3,7h iters )')ncyds
            elseif(nodiis)then
               strng3 = NA_10c
               strng6 = NA_10c
            else
               write(strng3,'(d10.2)')diison
               write(strng6,'(d10.2)')diisoff
            endif
            write(LuOut,3233)strng1,strng2,strng3,strng4,strng5,strng6
            call flush(LuOut)
         endif
         if (oprint_tolerances)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'Screening Tolerance Information',20,.true.)
            t1 = 10.d0**(-itol2e)
            r1 = dble(irqmax)
            write(LuOut,9372)tol_rho, iaoacc, itol2e, t1, ixcacc, 
     &                    r1, radius
         endif
      endif
      return
 3231 format(10x,'Convergence aids based upon iterative change in ',/,
     &       10x,'total energy or number of iterations. ',/,
     &       10x,'Levelshifting, if invoked, occurs when the ',/,
     &       10x,'HOMO/LUMO gap drops below (HL_TOL): ',d8.2,/,
     &       10x,'DIIS, if invoked, will attempt to extrapolate ',/,
     &       10x,'using up to (NFOCK): ',i2,' stored Fock matrices.',/)
 3232 format(10x,
     &       10x,'Damping(',i2,'%)  Levelshifting(',f3.1,
     &           ')       DIIS',/,
     &       10x,8x,15('-'),1x,19('-'),1x,15('-'))
 3233 format(10x,'dE  on:',2x,a10,7x,a10,10x,a10,/,
     &       10x,'dE off:',2x,a10,7x,a10,10x,a10,/)
 4001 format(10x,'AO basis - number of functions: ',i5,/,
     &       10x,'           number of shells: ',i5)
 4002 format(10x,'CD basis - number of functions: ',i5,/,
     &       10x,'           number of shells: ',i5)
 4003 format(10x,'XC basis - number of functions: ',i5,/,
     &       10x,'           number of shells: ',i5)
 4004 format(10x,'Number of quadrature shells: ',i5)
 9020 format(10x,'SCF calculation type: DFT')
 9030 format(10x,'Maximum number of iterations: ',I3)
 9035 format(10x,'Convergence on energy requested: ',d8.2)
 9040 format(10x,'Convergence on density requested: ',d8.2)
 9045 format(10x,'Convergence on gradient requested: ',d8.2)
 9050 format(10x,'Wavefunction type:  closed shell.')
 9055 format(10x,'Wavefunction type:  spin polarized.')
 9056 format(10x,'Use of symmetry is: ',a3,
     &           '; symmetry adaption is: ',a3)
 9110 format(10x,'This is a Direct SCF calculation.')
 9120 format(10x,'An Exch-Corr fitting basis will be used.')
 9130 format(10x,'A Charge density fitting basis will be used.')
 9135 format(10x,'Grid used for XC integration:  ',a)
 9136 format(10x,'Radial quadrature: Euler-MacLaurin, ')
 9137 format(10x,'Angular quadrature: Gauss-Legendre, with ',/,
     &       10x,i4,' radial points, ',i3,' theta points, and ',
     &           i3,' phi points.')
 9138 format(10x,'Angular quadrature: Lebedev, with ',/,
     &       10x,i4,' radial points and ',
     &           i3,' angular points.')
 9139 format(10x,'Spatial weights used: Delley. ')
 9140 format(10x,'Spatial weights used: Becke. ')
 9141 format(10x,'Parallel task size associated with evaluation of ',/,
     &           'grid based components has been modified to: ',i2)
 9223 format(10x,a40)
 9224 format(10x,a40,1x,f6.2,1x,a9)
 8150 format(10x,'No. of atoms     :',2x,i4)
 8200 format(10x,'No. of electrons :',2x,i4,/,
     &       10x,' Alpha electrons :',2x,i4,/,
     &       10x,'  Beta electrons :',2x,i4,/,
     &       10x,'Charge           :',2x,i4,/,
     &       10x,'Spin multiplicity:',2x,i4)
 9372 format(10x,'Density screening/tol_rho: ',d8.2,/,
     &       10x,'AO Gaussian exp screening on grid/accAOfunc: ',i3,/,
     &       10x,'CD Gaussian exp screening on grid/accCoul: ',i3,/,
     &       10x,'Schwarz screening/accCoul: ',d8.2,/,
     &       10x,'XC Gaussian exp screening on grid/accXCfunc: ',i3,/,
     &       10x,'Max radial quadrature (shell size)/accQrad: ',d8.2,/,
     &       10x,'Spatial weight screening/radius: ',d8.2,/)
      return
      end
