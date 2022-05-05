      subroutine prop_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "context.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
c>>>  Read input for properties calculations.    
c     
c     All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
      integer iline, ifield, nfield, ii, jj, icount
      
      character*255 movecs_in
      character*255 movecs_out
      character*255 test
      integer gshiftfile,gshiftopt
      integer shldfile,shldopt
      integer hypfile, hypopt
      integer efgfile, efgopt
      integer nbofile, nboopt
      integer aimfile
      integer moldenfile
      integer polfromsos
      character*32 norm_style
      integer dipole
      integer quadrupole
      integer octupole
      integer mulliken
      integer spinpopulation
      integer esp
      integer espiso, iso, tol
      integer efield
      integer efieldmap
      integer efieldgrad
      integer efieldgradZ4 ! Added by FA
      integer gshift
      integer efieldgradmap
      integer electrondensity
      integer electrondensitymap
      integer spindensity
      integer spindensitymap
      integer stoneanalysis
      integer localization
      integer hyperfine
      integer spinspin, spnspnpairs
      integer giao, giaoatoms
      integer hfineatoms
      integer efgz4atoms ! Added by FA
      integer center
      integer cphf3write,checkorth
      integer l_atefg,k_atefg,
     &        l_atshi,k_atshi,
     &        l_athyp,k_athyp, 
     &        l_pr, k_pr
      double precision point(3)
      logical status

      integer atomnr(100) ! FA-04-22-10

c ... jochen:
      integer aoresponse, nfreq, l_fr, k_fr, 
     &        velocity, magneticpert, aor_giao, aor_damp, aor_conv,
     &        aor_orbeta, aor_gprime, aor_quad, aor_analyze,
     &        loc_boys, loc_pm, loc_ibo
      double precision gamwidth, convergence
c     
c>>>  declarations from NWdft
c     
      integer infinity ! for infinite loop
      parameter (infinity=100000) 
c
      movecs_in     = ' '
      movecs_out    = ' '
c
c>>>  Default property settings.
c
      aimfile           = 0
      moldenfile        = 0
      polfromsos        = 0
      nbofile           = 0
      nboopt            = 0
      efgfile           = 0 
      efgopt            = 0 
      hypfile           = 0 
      hypopt            = 0 
      gshiftfile        = 0 
      gshiftopt         = 0 
      shldfile          = 0 
      shldopt           = 0 
c
      dipole            = 1
      quadrupole        = 1
      octupole          = 1
      mulliken          = 1
      spinpopulation    = 1
      esp               = 1
      espiso            = 1
      efield            = 1
      efieldmap         = 1
      efieldgrad        = 1
      efieldgradZ4      = 1
      efieldgradmap     = 1
      electrondensity   = 1
      electrondensitymap= 1
      spindensity       = 1
      spindensitymap    = 1
      stoneanalysis     = 1
      hyperfine         = 1
      spinspin          = 1
      giao              = 1
      gshift            = 1 
      center            = 1 ! default of center of charge
      point(1)          = 0.0d+0
      point(2)          = 0.0d+0
      point(3)          = 0.0d+0
      giaoatoms         = 0
      hfineatoms        = 0
      efgz4atoms        = 0 ! FA
      spnspnpairs       = 0

c ... jochen:
      aoresponse        = 1 
      velocity          = 1
      magneticpert      = 1
      aor_giao          = 1
      aor_damp          = 1
      aor_conv          = 1
      aor_orbeta        = 1
      aor_gprime        = 1
      aor_quad          = 1
      aor_analyze       = 1
      nfreq             = 0
      cphf3write        = 0 ! FA
      checkorth         = 0 ! FA
      gamwidth          = 0d0
      localization      = 1
      loc_boys          = 1
      loc_pm            = 1
      loc_ibo           = 1
c     
c     read inputs; store in rtdb
c     
      do iline = 1,infinity     !  infinite loop over input lines
        
      if (.not. inp_read()) call errquit
     $   ('prop_input: inp_read failed', 0, INPUT_ERR)
        
      nfield = inp_n_field()
      do ifield = 1,nfield
        
      if (inp_a(test)) then
c       
        if ( inp_compare(.false., 'hyperfine', test)) then
          hyperfine = 0
          if (inp_i(hfineatoms)) then
            if (.not. ma_push_get(mt_int,hfineatoms,'atoms',
     &                            l_athyp,k_athyp))
     &         call errquit('prp_input ma_push_get fail k_pr',0,MA_ERR)
            icount = 0
            do ii = 1, hfineatoms
              status = inp_i(int_mb(k_athyp+ii-1))
              if (status) then
                icount = icount+1
              else
                write(6,'(a,a,i4,a)') 
     &             'NMR hyperfine property input found',
     &             ' fewer atoms than expected, only',icount,
     &             ' will be used'
                hfineatoms = icount
                if (icount.eq.0) then 
                  if (.not.ma_pop_stack(l_athyp)) call
     &        errquit('prop_input: ma_pop_stack failed k_athyp',
     &                0,MA_ERR)
                endif
                goto 14
              endif
            enddo
          endif
   14     continue      
        elseif ( inp_compare(.false., 'spinspin', test)) then
          spinspin = 0
          if (inp_i(spnspnpairs)) then
            if (.not. ma_push_get(mt_int,spnspnpairs*2,'prs',l_pr,k_pr))
     &         call errquit('prp_input ma_push_get fail k_pr',0,MA_ERR)
            icount = 0
            do ii = 1, spnspnpairs
              status = inp_i(int_mb(k_pr+ii-1))
              status = inp_i(int_mb(k_pr+ii+spnspnpairs-1))
              if (status) then
                icount = icount+1
              else
                write(6,'(a,a,i4,a)') 
     &             'Spin-spin property input found',
     &             ' fewer pairs than expected, only',icount,
     &             ' will be used'
                if (icount.gt.0) then
                  do jj = 0, icount-1
                    int_mb(k_pr+jj+icount) = int_mb(k_pr+jj+spnspnpairs)
                  enddo
                endif
                spnspnpairs = icount
                if (icount.eq.0) then 
                  if (.not.ma_pop_stack(l_pr)) call
     &         errquit('prop_input: ma_pop_stack failed k_pr',0,MA_ERR)
                endif
                goto 10
              endif
            enddo
          endif
   10     continue                
        elseif ( inp_compare(.false., 'shielding', test)) then
          giao = 0
          if (inp_i(giaoatoms)) then
             if (.not.ma_alloc_get(
     &          mt_int,giaoatoms,'atoms',l_atshi,k_atshi))
     &        call errquit('prp_input ma_alloc_get failed k_ashi',
     &                     0,MA_ERR)
            icount = 0
            do ii = 1, giaoatoms
              status = inp_i(int_mb(k_atshi+ii-1))
              if (status) then
                icount = icount+1
              else
                write(6,'(a,a,i4,a)') 
     &             'NMR shielding property input found',
     &             ' fewer atoms than expected, only',icount,
     &             ' will be used'
                giaoatoms = icount
                if (icount.eq.0) then 
                  if (.not.ma_free_heap(l_atshi)) call
     &              errquit('prop_input: ma_free_heap l_atshi',
     &                      0, MA_ERR)
                endif
                goto 11
              endif
            enddo
          endif
   11     continue                
        elseif ( inp_compare(.false., 'gshift', test)) then
          gshift = 0   
c  ... jochen: parse input for response/aoresponse key
      elseif ( 
     &    inp_compare(.false.,   'response', test) .or.
     &    inp_compare(.false., 'aoresponse', test)) then
        aoresponse = 0
        if (inp_i(nfreq)) then
          if (.not. ma_push_get(mt_dbl,nfreq,'freqs',l_fr,k_fr))
     &       call errquit('freq_input ma_push_get fail l_fr',0,MA_ERR)
          icount = 0
          do ii = 1, nfreq
            status = inp_f(dbl_mb(k_fr+ii-1))
            if (status) then
              icount = icount+1
            else
              write(6,'(a,a,i4,a)') 
     &           'Response property input found',
     &           ' fewer frequencies than expected, only',icount,
     &           ' will be used'
              nfreq = icount
              if (icount.eq.0) then 
                if (.not.ma_pop_stack(l_fr)) call
     &          errquit('prop_input: ma_pop_stack failed k_fr',0,MA_ERR)
              endif
              goto 12
            endif
          enddo
        endif
   12   continue  
      elseif ( inp_compare(.false., 'velocity', test)) then
         velocity = 0
      elseif ( inp_compare(.false., 'magneticpert', test)) then
         magneticpert = 0
      elseif ( inp_compare(.false., 'giao', test)) then
         aor_giao = 0
      elseif ( inp_compare(.false., 'orbeta', test)) then
        aor_orbeta = 0
      elseif ( inp_compare(.false., 'gprime', test)) then
         aor_gprime = 0
      elseif ( inp_compare(.false., 'bdtensor', test)) then
         aor_quad = 0
      elseif ( inp_compare(.false., 'analysis', test)) then
         aor_analyze = 0
      elseif ( inp_compare(.false., 'gshift', test)) then
         aor_giao = 0
      elseif ( inp_compare(.false., 'damping', test)) then
         aor_damp = 0
         if (.not.inp_f(gamwidth)) then
           call errquit('Response: error reading damping constant',
     &        0, INPUT_ERR)
         endif
c cphf3write=1 ==> write (z1,Az1) to file
c z1 , history matrix of solutions 
c Az1, history matrix of products Az1       
        status=inp_i(cphf3write)
c        write(*,*) 'Reading cphf3write=',cphf3write
        status=inp_i(checkorth)
c        write(*,*) 'Reading checkorth=',checkorth
      elseif ( inp_compare(.false., 'convergence', test)) then
         aor_conv = 0
         if (.not.inp_f(convergence)) then
           call errquit('Response: error reading convergence crit.',
     &        0, INPUT_ERR)
         endif
c        ... jochen end
c ... jochen: localization input         
      elseif ( inp_compare(.false., 'localization', test)) then
        localization = 0
        if (inp_a(test)) then
          if ( inp_compare(.false., 'boys', test)) then
            loc_boys = 0
          elseif ( inp_compare(.false., 'ibo', test)) then
            loc_ibo = 0
          elseif ( inp_compare(.false., 'pm', test)) then
            loc_pm = 0
          else ! default is PM
            loc_pm = 0
          end if
        else
          loc_pm = 0
        end if
c ... jochen end        
      elseif ( inp_compare(.false., 'stoneanalysis', test)) then
         stoneanalysis  = 0
         stoneanalysis  = 1
      elseif ( inp_compare(.false., 'spindensitymap', test)) then
         spindensitymap = 0
         spindensity    = 0
      elseif ( inp_compare(.false., 'spindensity', test)) then
         spindensity = 0
      elseif ( inp_compare(.false., 'electrondensitymap', test)) then
         electrondensitymap = 0
         electrondensity    = 0
      elseif ( inp_compare(.false., 'electrondensity', test)) then
         electrondensity = 0
      elseif ( inp_compare(.false., 'efieldgradmap', test)) then
         efieldgradmap   = 0
         efieldgrad      = 0
      elseif ( inp_compare(.false., 'efieldgrad', test)) then
         efieldgrad      = 0
      elseif ( inp_compare(.false., 'efieldgradZ4', test)) then
         efieldgradZ4    = 0
c
c        calculate the efg for requested atoms
c
         if (inp_i(efgz4atoms)) then
            if (.not. ma_push_get(mt_int,efgz4atoms,'atoms',
     &                            l_atefg,k_atefg))
     &         call errquit('prp_input ma_push_get fail k_pr',0,MA_ERR)
            icount = 0
            do ii = 1, efgz4atoms
              status = inp_i(int_mb(k_atefg+ii-1))
              if (status) then
                icount = icount+1
              else
                write(6,'(a,a,i4,a)') 
     &             'EFGZ4 property input found',
     &             ' fewer atoms than expected, only',icount,
     &             ' will be used'
                efgz4atoms = icount
                if (icount.eq.0) then 
                  if (.not.ma_pop_stack(l_atefg)) call
     &        errquit('prop_input: ma_pop_stack failed k_atefg',
     &                0,MA_ERR)
                endif
                goto 13
              endif
            enddo
         endif
   13    continue      
c --------- storing selected atmnr in rtdb ------- START
          if (efgz4atoms.gt.0) then
           if (.not. rtdb_put(rtdb, 'efgz4:natoms',mt_int,
     $                        1,efgz4atoms))
     $      call errquit('prop_input-EFGZ4: rtdb_put failed',
     $                   555, RTDB_ERR)
           if (.not. rtdb_put(rtdb, 'efgz4:atom list',mt_int,
     $                        efgz4atoms,int_mb(k_atefg)))
     $      call errquit('prop_input-EFGZ4: rtdb_put failed',
     $                   555, RTDB_ERR)
          if (.not.ma_pop_stack(l_atefg)) call
     &        errquit('prop_input: ma_pop_stack failed k_atefg',
     &                0,MA_ERR)
          endif
c --------- storing selected atmnr in rtdb ------- END
c ==========================================
      elseif ( inp_compare(.false., 'efieldmap', test)) then
         efieldmap       = 0
         efield          = 0
      elseif ( inp_compare(.false., 'efield', test)) then
         efield          = 0
      elseif ( inp_compare(.false., 'espiso', test)) then
         espiso          = 0

         if (inp_a(test)) then
           if( inp_compare(.false., 'iso', test))
     >          write(*,*) "iso found"
                if (inp_f(iso)) then
                   continue
                end if
         end if
         if (inp_a(test)) then
           if( inp_compare(.false., 'tol', test))
     >          write(*,*) "tol found"
                if (inp_f(tol)) then
                   continue
                end if
         end if
      elseif ( inp_compare(.false., 'esp', test)) then
         esp             = 0
         if (inp_a(test)) then
           if( inp_compare(.false., 'grid', test))
     >          call prop_grid_input(rtdb)
         end if
      elseif ( inp_compare(.false., 'grid', test)) then
         call prop_grid_input(rtdb)
      elseif ( inp_compare(.false., 'spinpopulation', test)) then
         spinpopulation  = 0
         spinpopulation  = 1
      elseif ( inp_compare(.false., 'mulliken', test)) then
         mulliken  = 0
      elseif ( inp_compare(.false., 'octupole', test)) then
         octupole  = 0
      elseif ( inp_compare(.false., 'quadrupole', test)) then
         quadrupole  = 0
      elseif ( inp_compare(.false., 'vectors', test)) then
         if(.not.inp_a(movecs_in)) goto 10000
      elseif ( inp_compare(.false., 'dipole', test)) then
         dipole  = 0
      elseif ( inp_compare(.false., 'nbofile', test)) then
         nbofile = 1
c ... jochen: read an option (integer value) to define 
c     what NBO options will appear    
          if (inp_i(nboopt)) then 
             continue           ! just read the number
          end if
      elseif ( inp_compare(.false., 'aimfile', test)) then
         aimfile = 1
      elseif ( inp_compare(.false., 'moldenfile', test)) then
         moldenfile = 1
      elseif ( inp_compare(.false., 'polfromsos', test)) then
         polfromsos = 1
      elseif ( inp_compare(.false., 'molden_norm', test)) then
         if(.not. inp_a(norm_style))
     $     call errquit(' prop_input: missing val for norm_stl ',0,
     I           INPUT_ERR)
         
      elseif ( inp_compare(.false., 'efgfile', test)) then
         efgfile = 1
c ... Fredy: read an option (integer value) to define 
c     what EFG options will appear    
          if (inp_i(efgopt)) then 
             continue           ! just read the number
          end if
      elseif ( inp_compare(.false., 'hypfile', test)) then
         hypfile = 1
          if (inp_i(hypopt)) then 
             continue           ! just read the number
          end if
      elseif ( inp_compare(.false., 'gshiftfile', test)) then
         gshiftfile = 1
          if (inp_i(gshiftopt)) then
             continue           ! just read the number
          end if
      elseif ( inp_compare(.false., 'shldfile', test)) then
         shldfile = 1
          if (inp_i(shldopt)) then
             continue           ! just read the number
          end if
      elseif ( inp_compare(.false., 'all', test)) then
         giao   = 0
         gshift = 0
         localization = 0
         loc_pm = 0
         stoneanalysis  = 0
         stoneanalysis  = 1
         spindensitymap = 0
         spindensity    = 0
         electrondensitymap = 0
         electrondensity    = 0
         efieldgradmap   = 0
         efieldgrad      = 0 ! regular efg
         efieldgradZ4    = 0 ! picture change efg with zora/so-zora
         efieldmap       = 0
         efield          = 0
         esp             = 0
         hyperfine       = 0
         spinspin        = 0
         spinpopulation  = 0
         spinpopulation  = 1
         mulliken  = 0
         octupole  = 0
         quadrupole  = 0
         dipole    = 0
         nbofile   = 1
         aimfile = 0
         efgfile   = 1
         hypfile   = 1
         gshiftfile= 1
         hypfile   = 1
         gshiftfile= 1
         shldfile  = 1
      elseif ( inp_compare(.false., 'center', test)) then
         if (inp_a(test)) then
           if (inp_compare(.false., 'coc', test)) then
             center = 1  ! default of center of charge
           elseif (inp_compare(.false., 'com', test)) then
             center = 2  ! center of mass
           elseif (inp_compare(.false., 'origin', test)) then
             center = 3  ! origin
           elseif (inp_compare(.false., 'arb', test)) then
             center = 4  ! arbitrary point
             status = .true.
             status = status .and. inp_f(point(1))
             status = status .and. inp_f(point(2))
             status = status .and. inp_f(point(3))
             if (.not.status)
     $         call errquit 
     $          ('prop_input: invalid arbitrary point',555, INPUT_ERR)
           else
             call errquit ('prop_input: not a valid center type',555,
     &       INPUT_ERR)
           endif
         else
            call errquit ('prop_input: no center information',555,
     &       INPUT_ERR)
         endif
c             
      else if (inp_compare(.false.,'print', test) .or. 
     $        inp_compare(.false.,'noprint', test)) then
c          
         call util_print_input(rtdb, 'prop')
c
      elseif ( inp_compare(.false., 'end', test)) then
         goto 99999
      else
         write(6,'(a,a,a)') ' prop_input: keyword ',test,
     &             ' not recognized'
         goto 10000 !  die with error
      endif 
          else           ! no further field found in line
            goto 88      ! next line
      endif              ! terminates block    if (inp_a(test)) then
      enddo              ! end of loop over input fields 
   88   continue 
      enddo              ! end of infinite loop over input lines
99999 continue 
c
c    place integer parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'prop:nbofile', mt_int, 1,
     $                               nbofile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... jochen: also write NBO option key to RTDB:
      if (.not. rtdb_put(rtdb, 'prop:nboopt', mt_int, 1,
     $                               nboopt           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:aimfile', mt_int, 1,
     $                               aimfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:moldenfile', mt_int, 1,
     $                               moldenfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c         endif
      if(moldenfile.ne.0) then
       if(norm_style(1:inp_strlen(norm_style)).eq.'janpa') then
      if (.not. rtdb_put(rtdb, 'molden:do_norm_janpa', mt_log, 1,
     $                               .true.          ))
     $     call errquit('prop_input: rtdb_put failed', 1, RTDB_ERR)
       elseif (norm_style(1:inp_strlen(norm_style)).eq.'nwchem') then
      if (.not. rtdb_put(rtdb, 'molden:do_norm_nwchem', mt_log, 1,
     $                               .true.          ))
     $     call errquit('prop_input: rtdb_put failed', 2, RTDB_ERR)
       elseif (norm_style(1:inp_strlen(norm_style)).eq.'none') then
      if (.not. rtdb_put(rtdb, 'molden:do_norm_nwchem', mt_log, 1,
     $                               .false.          ))
     $     call errquit('prop_input: rtdb_put failed', 2, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'molden:do_norm_janpa', mt_log, 1,
     $                               .false.          ))
     $     call errquit('prop_input: rtdb_put failed', 1, RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'prop:do_fockbuild',mt_log,1,.false.))
c     $     call errquit('prop_input: rtdb_put failed',9, RTDB_ERR)
       else
          call errquit('prop_input: wrong molden_norm value',2,RTDB_ERR)
       endif
       endif
      if (.not. rtdb_put(rtdb, 'prop:polfromsos', mt_int, 1,
     $                               polfromsos           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c .... Fredy: write efg output for NBO program 04-24-10
      if (.not. rtdb_put(rtdb, 'prop:efgfile', mt_int, 1,
     $                               efgfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... Fredy: also write EFG option key to RTDB:
      if (.not. rtdb_put(rtdb, 'prop:efgopt', mt_int, 1,
     $                               efgopt           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:hypfile', mt_int, 1,
     $                               hypfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:hypopt', mt_int, 1,
     $                               hypopt           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:gshiftfile', mt_int, 1,
     $                               gshiftfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:gshiftopt', mt_int, 1,
     $                               gshiftopt           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:shldfile', mt_int, 1,
     $                               shldfile           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:shldopt', mt_int, 1,
     $                               shldopt           ))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:dipole', mt_int, 1,
     $                               dipole            ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:quadrupole', mt_int, 1,
     $                               quadrupole        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:octupole', mt_int, 1, 
     $                               octupole          ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:mulliken', mt_int, 1,
     $                               mulliken          ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:spinpopulation', mt_int, 1,
     $                               spinpopulation    ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:esp', mt_int, 1,
     $                               esp               ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:espiso', mt_int, 1,
     $                               espiso             ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:espiso:iso', mt_dbl, 1,
     $                               iso             ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:espiso:tol', mt_dbl, 1,
     $                               tol             ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efield', mt_int, 1,
     $                               efield            ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efieldmap', mt_int, 1,
     $                               efieldmap         ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efieldgrad', mt_int, 1,
     $                               efieldgrad        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efieldgradZ4', mt_int, 1,! Added FA
     $                               efieldgradZ4        ))   ! Added FA  
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR) !FA
      if (.not. rtdb_put(rtdb, 'prop:efieldgradmap', mt_int, 1,
     $                               efieldgradmap     ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:electrondensity', mt_int, 1,
     $                               electrondensity   ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:electrondensitymap', mt_int, 1,
     $                               electrondensitymap))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:spindensity', mt_int, 1,
     $                               spindensity       ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:spindensitymap', mt_int, 1,
     $                               spindensitymap    ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:stoneanalysis', mt_int, 1,
     $                               stoneanalysis     ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'prop:loc_boys', mt_int, 1,
     $                               loc_boys  ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'prop:loc_pm', mt_int, 1,
     $                               loc_pm  ))
     $  call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c      
      if (.not. rtdb_put(rtdb, 'prop:loc_ibo', mt_int, 1,
     $                               loc_ibo  ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)      
      if (.not. rtdb_put(rtdb, 'prop:spinspin', mt_int, 1,
     $                               spinspin         ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (spnspnpairs.gt.0) then
          if (.not.rtdb_put(rtdb,'spinspin:pairs',mt_int,1,spnspnpairs))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
          if (.not. rtdb_put(rtdb,'spinspin:pair list',mt_int,
     &                       2*spnspnpairs,int_mb(k_pr)))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
          if (.not.ma_pop_stack(l_pr)) call
     &        errquit('prop_input: ma_pop_stack failed k_pr',0,MA_ERR)
      endif
      if (.not. rtdb_put(rtdb, 'prop:giao', mt_int, 1,
     $                               giao              ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (giaoatoms.gt.0) then
          if (.not. rtdb_put(rtdb, 'giao:natoms', mt_int, 1,giaoatoms))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
          if (.not. rtdb_put(rtdb, 'giao:atom list', mt_int, giaoatoms,
     &                       int_mb(k_atshi)))
     $    call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
          if (.not.ma_free_heap(l_atshi)) 
     &    call errquit('prop_input: ma_free_heap l_atshi-1',
     &                      0, MA_ERR)
      endif
      if (.not. rtdb_put(rtdb, 'prop:hyperfine', mt_int, 1,
     $                               hyperfine         ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (hfineatoms.gt.0) then
         if (.not. rtdb_put(rtdb, 'hfine:natoms', mt_int, 1,hfineatoms))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'hfine:atom list', mt_int, hfineatoms,
     &                      int_mb(k_athyp)))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
         if (.not.ma_pop_stack(l_athyp)) call
     &       errquit('prop_input: ma_pop_stack failed k_athyp',
     &               0,MA_ERR)
      endif

      if (.not. rtdb_put(rtdb, 'prop:gshift', mt_int, 1,
     $                               gshift             ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... jochen: add aoresponse key, frequencies, and options to RTDB
      if (.not. rtdb_put(rtdb, 'prop:aoresponse', mt_int, 1,
     $                               aoresponse        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... jochen: I tested response with DFT where the convergence
c             limits made a lot of trouble. Therefore we set the
c             "nocheck" option which is used later in hnd_property.F
      if (aoresponse.eq.0) then
       if (.not. rtdb_put(rtdb, 'prop:nocheck', mt_log, 1, .true.))
     $      call errquit('prop_input: rtdb_put nocheck failed', 0, 
     $      RTDB_ERR)
      endif
      if (.not. rtdb_put(rtdb, 'aoresponse:velocity', mt_int, 1,
     $   velocity))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:magneticpert', mt_int, 1,
     $   magneticpert))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:giao', mt_int, 1,
     $   aor_giao ))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:orbeta', mt_int, 1,
     $   aor_orbeta))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:gprime', mt_int, 1,
     $   aor_gprime))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)      
      if (.not. rtdb_put(rtdb, 'aoresponse:quad', mt_int, 1,
     $   aor_quad))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:analyze', mt_int, 1,
     $   aor_analyze))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:damping', mt_int, 1,
     $   aor_damp))
     $   call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... jochen: write cphf iteration threshold if set by input
      if (aor_conv.eq.0) then
      if (.not. rtdb_put(rtdb, 'cphf:thresh',    mt_dbl, 1,
     $                               convergence     ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      endif ! aor_conv
c ----- FA-11-30-12: Adding new parameter: aor_damp ---- START
      if (aor_damp.eq.0) then ! FA-06-06-12
        if (.not. rtdb_put(rtdb, 'cphf:cphf3write',    mt_int, 1,
     $                               cphf3write))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'cphf:checkorth',    mt_int, 1,
     $                               checkorth))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      endif ! aor_damp
c ----- FA-11-30-12: Adding new parameter: aor_damp ---- END
      if (.not. rtdb_put(rtdb, 'aoresponse:gamwidth', mt_dbl, 1,
     $                               gamwidth        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)

      if (nfreq.gt.0) then
        if (.not. rtdb_put(rtdb, 'aoresponse:nfreq', mt_int, 1,nfreq))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'aoresponse:freqlist', mt_dbl, nfreq,
     &     dbl_mb(k_fr)))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
        if (.not.ma_pop_stack(l_fr)) call
     &     errquit('prop_input: ma_pop_stack failed k_fr',0,MA_ERR)
      else
        if (.not. rtdb_put(rtdb, 'aoresponse:nfreq', mt_int, 1,0))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      end if
c ... jochen: end response section
      if (.not. rtdb_put(rtdb, 'prop:center',mt_int, 1, center))
     $      call errquit('prop_input: rtdb_put failed', 555, RTDB_ERR)
      if (center.eq.4) then
        if (.not. rtdb_put(rtdb, 'prop:center_val',mt_dbl, 3, point))
     $      call errquit('prop_input: rtdb_put failed', 555, RTDB_ERR)
      endif
c       
      if(movecs_in.ne.' ') then
         if (.not. rtdb_cput(rtdb, 'prop:vectors',1, movecs_in))
     $      call errquit('prop_input: rtdb_cput failed', 555, RTDB_ERR)
      end if

      return
c       
10000 continue !  die with error
      write (6,*) 'prop_input: error in input file, last field: '
      write (6,'(a)') test
      call errquit('prop_input: error in input file', 0, INPUT_ERR)
c       
      end

      subroutine prop_grid_input(rtdb)
c
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "util_params.fh"
c
      integer rtdb
c
      character*32 tag
      character*32 pname
      character*255 token
      character*255 buffer
      character*255 token0
      integer i,ic
      double precision f3(3),f
      integer itok,itok3(3)
c
      pname = "prop_input_grid"
c
      tag = 'prop:cubefile'
      if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >    goto 999
c
600   continue
      if(.not.inp_a(token)) goto 20
      if (inp_compare(.false.,'step',token)) then
        ic = 0
        do i=1,1
          if(.not.inp_f(f3(i))) exit
          ic = ic + 1
        end do
        if(ic.eq.1) then
          f3 = f3(1)
        else if(ic.ne.3) then
           call errquit(pname//'need 1 float '//token,
     >                  0, INPUT_ERR)
        end if
        tag="prop:grid:step"
        if(.not.rtdb_put(rtdb,tag,mt_dbl,1,f3))
     >       goto 999

      else if (inp_compare(.false.,'pad',token)) then
        ic = 0
        do i=1,3
          if(.not.inp_f(f3(i))) exit
          ic = ic + 1
        end do
        if(ic.eq.1) then
          f3 = f3(1)
        else if(ic.ne.3) then
           call errquit(pname//'need 1 or 3 floats '//token,
     >                  0, INPUT_ERR)
        end if
        tag="prop:grid:pad"
        f3 = f3*cang2au
        if(.not.rtdb_put(rtdb,tag,mt_dbl,3,f3)) 
     >       goto 999
      else if (inp_compare(.false.,'rmin',token)) then
        ic = 0
        do i=1,3
          if(.not.inp_f(f3(i))) exit
          ic = ic + 1
        end do
        if(ic.ne.3) then
           call errquit(pname//'need 3 floats '//token,0, INPUT_ERR)
        end if
        tag="prop:grid:rmin"
        f3 = f3*cang2au
        if(.not.rtdb_put(rtdb,tag,mt_dbl,3,f3)) 
     >       goto 999
      else if (inp_compare(.false.,'rmax',token)) then
        ic = 0
        do i=1,3
          if(.not.inp_f(f3(i))) exit
          ic = ic + 1
        end do
        if(ic.ne.3) then
           call errquit(pname//'need 3 floats '//token,0, INPUT_ERR)
        end if
        tag="prop:grid:rmax"
        f3 = f3*cang2au
        if(.not.rtdb_put(rtdb,tag,mt_dbl,3,f3)) 
     >       goto 999
       else if (inp_compare(.false.,'ngrid',token)) then
        ic = 0
        do i=1,3
          if(.not.inp_i(itok3(i))) exit
          ic = ic + 1
        end do
        if(ic.eq.1) then
          itok3 = itok3(1)
        else if(ic.ne.3) then
           call errquit(pname//'need 1/3 integers '//token,0, INPUT_ERR)
        end if
        tag="prop:grid:ngrid"
        if(.not.rtdb_put(rtdb,tag,mt_int,3,itok3)) 
     >       goto 999
c      disabling filename for now
c      --------------------------
      else if (inp_compare(.false.,'output',token)) then
        if(.not.inp_a(buffer)) 
     >       goto 997
        tag="prop:grid:output"
        if (.not. rtdb_cput(rtdb, tag,1, trim(buffer)))
     >       goto 999
      else
        goto 997
      end if
      goto 600
20    continue
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)
997   call errquit(pname//'invalid token '//token,0,
     >        RTDB_ERR)

      end
c $Id$
