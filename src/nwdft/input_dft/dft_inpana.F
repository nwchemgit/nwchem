      subroutine dft_inpana(rtdb)
c
c     $Id: dft_inpana.F,v 1.71 2008-12-06 23:09:50 niri Exp $
c
c*********************************************************************
c     
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c     Write pertinent information to user output.
c     
c*********************************************************************
      implicit none
#include "errquit.fh"
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
      double precision AUTOANG, ANGTOAU
      Parameter (AUTOANG = 0.529177249d0, ANGTOAU = 1.d0/AUTOANG)
      logical even, cksetd, no_prune
      integer noc(2)
      Integer rtdb              !  runtime database handle
      integer me,ichg,nel
      integer noc1
      integer n, na, nb, test_sic
      integer nshells
      integer ictr
      double precision bsr, radang
      double precision anucl_charg, ckfac, t1
      double precision ictr_coord(3), ictr_chg,smear_sigma
      character*3 on_off_1, on_off_2
      logical oprint_general, oprint_grid, oprint_xc, 
     &     oprint_convergence, oprint_tolerances, oprint_sic,
     ,     lsigma
      character*9 local_c, nonlocal_c, nineb_c
      character*10 start_10c, NA_10c, asap_10c
      character*10 strng1, strng2, strng3, strng4, strng5, strng6
      character*16 tag, theory
      character*20 rgridnames(8)
      logical xc_gotxc
      external xc_gotxc
      data rgridnames /'Euler-MacLaurin','Mura-Knowles',
     .     'Treutler-Ahlrichs','Gauss-Legendre','G-C-interv',
     ,     'Lindh','Chebyshev','Legendre'/
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
      oprint_sic = util_print('sic information',print_default)
c     
c     Figure out the number of electrons from the required total
c     charge and the sum of nuclear charges
c     
      if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $        call errquit('dft_inpana: theory not specified',0,
     &       RTDB_ERR)
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &     call errquit('dft_inpana: geom_nuc_charge failed', 0,
     &       GEOM_ERR)
      nel = nint(anucl_charg - rcharge)
      if (nel .le. 0) call errquit
     $     ('dft_inpana: negative no. of electrons ?', nel, INPUT_ERR)
      if (abs(anucl_charg - rcharge - dble(nel)) .gt. 1d-8)
     $     call errquit('dft_inpana: non-integral # of electrons ?', 0,
     &       INPUT_ERR)
c     
c     Check to see if calculation type is allowed.
c     
c     Even number of electrons required for RHF.
c     
      even=mod(nel,2).eq.0
c     
c     odd # of electrons or not a singlet state --> LSD
c     
      if ((.not.even).or.(mult.ne.1).or.(theory.eq.'sodft')) ipol=2
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
         if (.not.even) then
           write(LuOut,*)' number of electrons :',nel
           write(LuOut,*)' multiplicity        :',mult
           call errquit(
     &         ' no. of electrons and multiplicity not compatible',nel,
     &       INPUT_ERR)
         endif
         if(mult.gt.0) then
            noc(2) = (nel - mult + 1)/2
            noc(1) = nel - noc(2)
            noc1  = noc(1) + noc(2)
         else
            noc(1) = (nel + mult + 1)/2
            noc(2) = nel - noc(1)
            noc1  = noc(1) + noc(2)
         endif
      endif
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     $     call errquit('inpana: dft:ipol put failed', 0, RTDB_ERR)
c     
c     Check to see if there are enough electrons for this
c     value of the multiplicity.
c     
      if (noc(2).lt.0)then
         call errquit('dft: #electrons not valid for multiplicity',mult,
     &       INPUT_ERR)
      endif
c     
c     write noc (consistent with definition in ddscf) to rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('inpana: rtdb_put of noc failed', 0, RTDB_ERR)
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
      do n = 1, numfunc
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
         if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, numfunc, cfac))
     $        call errquit('dft_input: rtdb_put failed', 210, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('dft_input: rtdb_put failed', 211, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, numfunc, lcfac ))
     $        call errquit('dft_input: rtdb_put failed', 9, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, numfunc, lxfac ))
     $        call errquit('dft_input: rtdb_put failed', 11, RTDB_ERR)
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
c     check on 2-e integral and XC grid tolerances
c
      call dft_inpanae(rtdb)
c
c     Check for no pruning; no_prune
c
      if (.not. rtdb_get(rtdb, 'dft:no_prune', mt_log, 1,
     &         no_prune))no_prune = .false.
      if (.not. rtdb_get(rtdb, 'dft:test_sic', mt_int, 1,
     &     test_sic))test_sic = 0
      lsigma=rtdb_get(rtdb, 'dft:smear_sigma',mt_dbl,1,smear_sigma)
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
            ichg = nint(rcharge)
            write(LuOut,8200)nel,na,nb,ichg,mult
            on_off_1 = 'off'
            if (oskel)on_off_1 = 'on'
            on_off_2 = 'off'
            if (oadapt)on_off_2 = 'on'
            write(LuOut,9056)on_off_1, on_off_2
            write(LuOut,9030)iterations
            if (direct)write(LuOut,9110)
            if (.not. bas_numcont(AO_bas_han, nshells))
     &         call errquit('rdinput:rdinput:',86, BASIS_ERR)
            write(LuOut,4001)nbf, nshells
            if (nbf_cd.gt.0)then
               write(LuOut,9130)
              if (.not. bas_numcont(CD_bas_han, nshells_cd))
     &           call errquit('rdinput:rdinput:',87, BASIS_ERR)
               write(LuOut,4002)nbf_cd, nshells_cd
            endif
            if (nbf_xc.gt.0)then
               write(LuOut,9120)
               if (.not. bas_numcont(XC_bas_han, nshells_xc))
     &            call errquit('rdinput:rdinput:',88, BASIS_ERR)
               write(LuOut,4003)nbf_xc, nshells_xc
            endif
            write(LuOut,9035)e_conv
            if (d_conv.gt.0)then
               write(LuOut,9040)d_conv
            endif
            if (g_conv.gt.0)then
               write(LuOut,9045)g_conv
            endif
            call util_flush(LuOut)
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
            do n = 1, numfunc
               if (xccomb(n))write(LuOut,9223) xcname(n)
            enddo
c     
c           Do exact exchange differently.
c     
            if (lxfac(1).or.nlxfac(1))
     &           write(LuOut,9224) xname(1), xfac(1), nineb_c
            do n = 2, numfunc
               if (lxfac(n).and.nlxfac(n))then
                  write(LuOut,9224) xname(n), xfac(n), nineb_c
               elseif (lxfac(n).and.(.not.nlxfac(n)))then
                  write(LuOut,9224) xname(n), xfac(n), local_c
               elseif ((.not.lxfac(n)).and.nlxfac(n))then
                  write(LuOut,9224) xname(n), xfac(n), nonlocal_c
               endif
            enddo
            do n = 1, numfunc
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
c            ckfac = 0.0d0
c            do n = 1, numfunc
c               if (lcfac(n))ckfac = ckfac + cfac(n)
c            enddo
c            if (abs(ckfac-1.0d0).gt.1.d-8)then
c               write(LuOut,*)
c     &            ' WARNING: Sum of local correlation is ',ckfac
c               write(LuOut,*)' Sum of components do not equal unity. '
c            endif
c            ckfac = 0.0d0
c            do n = 1, numfunc
c               if (nlcfac(n))ckfac = ckfac + cfac(n)
c            enddo
c            if (abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
c               write(LuOut,*)
c     &            ' WARNING: Sum of nonlocal correlation is ',ckfac
c               write(LuOut,*)
c     &            ' Sum of components do not equal unity or 0. '
c            endif
c            ckfac = 0.0d0
c            do n = 1, numfunc
c               if (lxfac(n))ckfac = ckfac + xfac(n)
c            enddo
c            if (abs(ckfac-1.0d0).gt.1.d-8)then
c               write(LuOut,*)
c     &            ' WARNING: Sum of local exchange is ',ckfac
c               write(LuOut,*)' Sum of components do not equal unity. '
c            endif
c            ckfac = 0.0d0
c            do n = 1, numfunc
c               if (nlxfac(n))ckfac = ckfac + xfac(n)
c            enddo
c            if (abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
c               write(LuOut,*)
c     &            ' WARNING: Sum of nonlocal exchange is ',ckfac
c               write(LuOut,*)
c     &            ' Sum of components do not equal unity or 0. '
c            endif
c
c           Check if asymptotic correction will be added to potential
c           If both LB94 and CS00 are .true., it is assumed that the 
c           user meant to use CS00 (since CS00 uses LB94)
c
            if (cs00) then
               if (delta_ac.gt.1.0d90) then
               write(LuOut,*)
               write(LuOut,9226) 
     &         '  CS with a Zhan-Nichols-Dixon shift        '
               else
               write(LuOut,*)
               write(LuOut,9227) 
     &         '  Casida-Salahub correction with a shift    ',
     &         delta_ac,'au'
               endif
            else if (lb94) then
               write(LuOut,*)
               write(LuOut,9226) 
     &         '         van Leeuwen-Baerends correction    '
            endif
            call util_flush(LuOut)
         endif
         if (oprint_sic) then
           if (test_sic.eq.1) then
             write(LuOut,'(/,14x,"SIC perturbative approximation")')
           else
             if (test_sic.eq.2) then
               write(LuOut,
     .               '(/14x,"SIC/OEP without localized orbitals")')
             else
               if (test_sic.eq.4) then
                 write(LuOut,
     .                 '(/14x,"SIC/OEP with localized orbitals")')
               end if
             end if
           end if
         end if
         if (oprint_grid.and.xc_gotxc())then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'Grid Information',20,.true.)
            write(LuOut,9135)gridtype
            write(LuOut,9136) rgridnames(wradgrid)
            if (.not.leb)then
               write(LuOut,9142)
            else
               write(LuOut,9143)
            endif
            write(LuOut,9144)
            do n = 1, ntypes
c
c              Find an atom of this kind in the complete list.
c
               do ictr = 1, ncenters
                  if (iatype(ictr).eq.n) then
                     if (.not. geom_cent_get(geom, ictr, tag, 
     &                  ictr_coord, ictr_chg))call errquit
     &                  ('dft_inpana: geom_cent_get failed', 0,
     &       GEOM_ERR)
                     goto 40
                  endif
               enddo
   40          continue
               bsr = bsrad_atom_type(n)*AUTOANG
               radang = dble(nint(autoang*dble(rad_cutoff(1,n))))
               if (leb)then
                  write(LuOut,9138)tag,bsr,nrad(n),radang,
     &                nang(n)
               else
                  write(LuOut,9137)tag,bsr,nrad(n),radang,
     &                nang(n),2*nang(n)
               endif
            enddo
            if (no_prune)then
               on_off_1 = 'off'
            else                
               on_off_1 = 'on'
            endif
            write(LuOut,4005)on_off_1
            write(LuOut,4004)nqshells
c
            if (ldelley)then
               write(LuOut,9140) 'Delley'
            elseif(lssw) then
               if(whichssw.eq.'ssf ') then
                  write(LuOut,9140) ' Straatman-Scuseria-Fritsch'
               elseif(whichssw.eq.'erf1') then
                  write(LuOut,9140) ' Erf1'
               elseif(whichssw.eq.'erf2') then
                  write(LuOut,9140) ' Erf2'
               else
                  write(LuOut,9140) whichssw
               endif
            else
               write(LuOut,9140) 'Becke'
            endif
            if (nquad_task.ne.1)then
               write(LuOut,9141)nquad_task
            endif
            if (nq_chunk.ne.0)then
               write(LuOut,9145)nq_chunk
            endif
            call util_flush(LuOut)
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
            call util_flush(LuOut)
         endif
         if(lsigma .and. oprint_general) then
            write(luout,
     .           "(10x,'Smearing applied: ',d8.2,' (hartree)')"
     .           ) smear_sigma
         endif
         

         if (oprint_tolerances)then     
            write(LuOut,*)
            call util_print_centered
     &         (LuOut,'Screening Tolerance Information',20,.true.)
            t1 = 10.d0**(-itol2e)
c            write(LuOut,9372)tol_rho, iaoacc, icdacc, ixcacc, t1,  
c     &                    r1
            write(LuOut,9372)tol_rho, iaoacc, icdacc, ixcacc, t1
            call util_flush(LuOut)
         endif
      endif
c
      return
 3231 format(10x,'Convergence aids based upon iterative change in ',/,
     &       10x,'total energy or number of iterations. ',/,
     &       10x,'Levelshifting, if invoked, occurs when the ',/,
     &       10x,'HOMO/LUMO gap drops below (HL_TOL): ',1Pd8.2,/,
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
 4005 format(10x,'Grid pruning is: ',a3)
 9020 format(10x,'SCF calculation type: DFT')
 9030 format(10x,'Maximum number of iterations: ',I3)
 9035 format(10x,'Convergence on energy requested: ',1Pd8.2)
 9040 format(10x,'Convergence on density requested: ',1Pd8.2)
 9045 format(10x,'Convergence on gradient requested: ',1Pd8.2)
 9050 format(10x,'Wavefunction type:  closed shell.')
 9055 format(10x,'Wavefunction type:  spin polarized.')
 9056 format(10x,'Use of symmetry is: ',a3,
     &           '; symmetry adaption is: ',a3)
 9110 format(10x,'This is a Direct SCF calculation.')
 9120 format(10x,'An Exch-Corr fitting basis will be used.')
 9130 format(10x,'A Charge density fitting basis will be used.')
 9135 format(10x,'Grid used for XC integration:  ',a)
cedo 9136 format(10x,'Radial quadrature: Euler-MacLaurin. ')
 9136 format(10x,'Radial quadrature: ',A)
 9137 format(10x,a16,2x,f6.2,6x,i3,6x,f6.1,3x,i2,1x,'*',1x,i2)
 9138 format(10x,a16,2x,f6.2,6x,i3,8x,f6.1,3x,2x,i5)
c 9139 format(10x,'Spatial weights used: Delley. ')
 9140 format(10x,'Spatial weights used: ',A)
 9141 format(10x,'Parallel task size associated with evaluation of ',/,
     &       10x,'grid based components has been modified to: ',i2)
 9145 format(10x,'Chunking of the angular grid is being used; ',
     &           'nq_chunk = ',i4)
 9142 format(10x,'Angular quadrature: Gauss-Legendre. ')
 9143 format(10x,'Angular quadrature: Lebedev. ')
 9144 format(10x,'Tag',14x,'B.-S. Rad.',1x,'Rad. Pts.',1x,'Rad. Cut.',
     &        1x,'Ang. Pts.',/,
     &       10x,'---',14x,'----------',1x,'---------',1x,'---------',
     &        1x,'---------')
 9223 format(10x,a40)
 9224 format(10x,a40,1x,f6.3,1x,a9)
 9225 format(10x,a40,1x,f8.4)
 9226 format(10x,a44)
 9227 format(10x,a44,1x,f10.6,1x,a2)
 8150 format(10x,'No. of atoms     :',2x,i4)
 8200 format(10x,'No. of electrons :',2x,i4,/,
     &       10x,' Alpha electrons :',2x,i4,/,
     &       10x,'  Beta electrons :',2x,i4,/,
     &       10x,'Charge           :',2x,i4,/,
     &       10x,'Spin multiplicity:',2x,i4)
 9372 format(10x,'Density screening/tol_rho: ',1Pd8.2,/,
     &       10x,'AO Gaussian exp screening on grid/accAOfunc: ',i3,/,
     &       10x,'CD Gaussian exp screening on grid/accCDfunc: ',i3,/,
     &       10x,'XC Gaussian exp screening on grid/accXCfunc: ',i3,/,
     &       10x,'Schwarz screening/accCoul: ',1Pd8.2,/)
      end
      subroutine dft_inpanae(rtdb)
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb
c
      logical oprint_general,oprint_tolgr
      double precision tol2e
c
      oprint_general = util_print('general information',print_default)
      oprint_tolgr = util_print('grid_tol_info', print_high)
c
c     get e_conv from rtdb because hess_init
c
      if (.not.rtdb_get(rtdb,'dft:e_conv',mt_dbl,1,e_conv))
     .     call errquit('dftinpanae: rtdbget econv failed',0, RTDB_ERR)
c
c     check integral tolerances to make sure they match 
c     requested convergence tolerances.
c
c     make sure itol2e is less than 0.1*e_conv or (0.01*g_conv**2)
c
      if (10.d0**(-itol2e).gt.(1.0d-1*e_conv))then
         itol2e = -nint(log10(1.0d-1*e_conv))
         if (.not. rtdb_put(rtdb, 'dft:itol2e',
     &      mt_int, 1, itol2e))
     &      call errquit('dft_inpanae: rtdb_put failed', 127, RTDB_ERR)
         if (ga_nodeid().eq.0.and.oprint_general)then
            write(LuOut,*)' itol2e modified to match energy'
            write(LuOut,*)' convergence criterion.'
         endif
      endif
c
c     check density tolerance to make sure it matches
c     requested convergence tolerances.
c
c     make sure tol_rho is less than 0.01*e_conv or (0.01*g_conv**2)
c
      if (tol_rho.gt.(1.0d-3*e_conv))then
         tol_rho = 1.0d-3*e_conv
         if (.not. rtdb_put(rtdb, 'dft:tol_rho',
     &      mt_dbl, 1, tol_rho))
     &      call errquit('dft_inpanae: rtdb_put failed', 127, RTDB_ERR)
         if (ga_nodeid().eq.0.and.oprint_general)then
            write(LuOut,*)' tol_rho modified to match energy'
            write(LuOut,*)' convergence criterion.'
         endif
      endif
      if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1,
     &   iAOacc))then
         iAOacc=-nint(log(e_conv))
      else
        iAOacc=max(iAOacc,-nint(log(e_conv)))
      endif	
       if (.not. rtdb_put(rtdb, 'dft:iAOacc',
     &      mt_int, 1, iAOacc))
     &      call errquit('dft_inpanae: rtdb_put failed', 124, RTDB_ERR)
      tol2e = 10.d0**(-itol2e)
      if (.not. rtdb_put(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) 
     .     call errquit('dftinpanae:rtdbput failed',0, RTDB_ERR)
      if(oprint_tolgr.and.ga_nodeid().eq.0) then
        write(luout,*) ' dftinpanae: itol2e ',itol2e,
     ,     ' iaoacc ',iaoacc,
     ,    ' tol_rho ',tol_rho
       call util_flush(luout)
      endif
      return
      end
