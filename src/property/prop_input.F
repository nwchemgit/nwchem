      subroutine prop_input(rtdb)
*
* $Id: prop_input.F,v 1.16 2006-11-16 18:23:01 jochen Exp $
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
      integer nbofile
      integer dipole
      integer quadrupole
      integer octupole
      integer mulliken
      integer spinpopulation
      integer esp
      integer efield
      integer efieldmap
      integer efieldgrad
      integer efieldgradmap
      integer electrondensity
      integer electrondensitymap
      integer spindensity
      integer spindensitymap
      integer stoneanalysis
      integer boyslocalization
      integer hyperfine
      integer spinspin, spnspnpairs
      integer giao, giaoatoms
      integer center
      integer l_at, k_at, l_pr, k_pr
      double precision point(3)
      logical status

c ... jochen:
      integer aoresponse, nfreq, l_fr, k_fr, 
     &        velocity, magneticpert, aor_giao, aor_damp, aor_conv
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
      nbofile           = 0
c
      dipole            = 1
      quadrupole        = 1
      octupole          = 1
      mulliken          = 1
      spinpopulation    = 1
      esp               = 1
      efield            = 1
      efieldmap         = 1
      efieldgrad        = 1
      efieldgradmap     = 1
      electrondensity   = 1
      electrondensitymap= 1
      spindensity       = 1
      spindensitymap    = 1
      stoneanalysis     = 1
      boyslocalization  = 1
      hyperfine         = 1
      spinspin          = 1
      giao              = 1
      
      center            = 1   ! default of center of charge
      point(1)          = 0.0d+0
      point(2)          = 0.0d+0
      point(3)          = 0.0d+0
      giaoatoms         = 0
      spnspnpairs       = 0

c ... jochen:
      aoresponse        = 1 
      velocity          = 1
      magneticpert      = 1
      aor_giao          = 1
      aor_damp          = 1
      aor_conv          = 1
      nfreq             = 0

c     
c>>>  read inputs; store in rtdb
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
     &             ' viewer pairs than expected, only',icount,
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
            if (.not. ma_push_get(mt_int,giaoatoms,'atoms',l_at,k_at))
     &         call errquit('prp_input ma_push_get fail k_pr',0,MA_ERR)
            icount = 0
            do ii = 1, giaoatoms
              status = inp_i(int_mb(k_at+ii-1))
              if (status) then
                icount = icount+1
              else
                write(6,'(a,a,i4,a)') 
     &             'NMR shielding property input found',
     &             ' viewer atoms than expected, only',icount,
     &             ' will be used'
                giaoatoms = icount
                if (icount.eq.0) then 
                  if (.not.ma_pop_stack(l_at)) call
     &        errquit('prop_input: ma_pop_stack failed k_at',0,MA_ERR)
                endif
                goto 11
              endif
            enddo
          endif
   11     continue                
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
      elseif ( inp_compare(.false., 'damping', test)) then
         aor_damp = 0
         if (.not.inp_f(gamwidth)) then
           call errquit('Response: error reading damping constant',
     &        0, INPUT_ERR)
         endif
      elseif ( inp_compare(.false., 'convergence', test)) then
         aor_conv = 0
         if (.not.inp_f(convergence)) then
           call errquit('Response: error reading convergence crit.',
     &        0, INPUT_ERR)
         endif
c ... jochen end
      elseif ( inp_compare(.false., 'boyslocalization', test)) then
         boyslocalization = 0
         boyslocalization = 1
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
      elseif ( inp_compare(.false., 'efieldmap', test)) then
         efieldmap       = 0
         efield          = 0
      elseif ( inp_compare(.false., 'efield', test)) then
         efield          = 0
      elseif ( inp_compare(.false., 'esp', test)) then
         esp             = 0
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
      elseif ( inp_compare(.false., 'all', test)) then
         giao = 0
         boyslocalization = 0
         boyslocalization = 1
         stoneanalysis  = 0
         stoneanalysis  = 1
         spindensitymap = 0
         spindensity    = 0
         electrondensitymap = 0
         electrondensity    = 0
         efieldgradmap   = 0
         efieldgrad      = 0
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
         dipole  = 0
         nbofile = 1
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
c>>>    place integer parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'prop:nbofile', mt_int, 1,
     $                               nbofile           ))
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
      if (.not. rtdb_put(rtdb, 'prop:efield', mt_int, 1,
     $                               efield            ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efieldmap', mt_int, 1,
     $                               efieldmap         ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:efieldgrad', mt_int, 1,
     $                               efieldgrad        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
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
      if (.not. rtdb_put(rtdb, 'prop:boyslocalization', mt_int, 1,
     $                               boyslocalization  ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'prop:hyperfine', mt_int, 1,
     $                               hyperfine         ))
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
     &                       int_mb(k_at)))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
          if (.not.ma_pop_stack(l_at)) call
     &        errquit('prop_input: ma_pop_stack failed k_at',0,MA_ERR)
      endif
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
     $                               velocity        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:magneticpert', mt_int, 1,
     $                               magneticpert        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:giao', mt_int, 1,
     $                               aor_giao        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'aoresponse:damping', mt_int, 1,
     $                               aor_damp        ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
c ... jochen: write cphf iteration threshold if set by input
      if (aor_conv.eq.0) then
      if (.not. rtdb_put(rtdb, 'cphf:thresh',    mt_dbl, 1,
     $                               convergence     ))
     $      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      endif ! aor_conv
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
c>>>    store real parameters in rtdb
c       
cc        if (.not. rtdb_put(rtdb, 'inp:xxx', mt_dbl, 1, xxx))
cc     $       call errquit('dft_input: rtdb_put failed', 0)
c
       if(movecs_in.ne.' ') then
         if (.not. rtdb_cput(rtdb, 'prop:vectors',1, movecs_in))
     $      call errquit('prop_input: rtdb_cput failed', 555, RTDB_ERR)
       end if

        return
c       
10000   continue !  die with error
        write (6,*) 'prop_input: error in input file, last field: '
        write (6,'(a)') test
        call errquit('prop_input: error in input file', 0, INPUT_ERR)
c       
        end
