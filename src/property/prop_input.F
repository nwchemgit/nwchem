      subroutine prop_input(rtdb)
*
* $Id: prop_input.F,v 1.7 1999-12-09 17:15:40 windus Exp $
*
      implicit none
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
      integer iline, ifield, nfield
      
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
      integer iglo
      integer giao
      integer center
      double precision point(3)
      logical status

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
      iglo              = 1
      giao              = 1
      center            = 1   ! default of center of charge
      point(1)          = 0.0d+0
      point(2)          = 0.0d+0
      point(3)          = 0.0d+0
c     
c>>>  read inputs; store in rtdb
c     
      do iline = 1,infinity !  infinite loop over input lines

        if (.not. inp_read()) call errquit
     $       ('prop_input: inp_read failed', 0)

        nfield = inp_n_field()
        do ifield = 1,nfield

      if (inp_a(test)) then
c             
      if ( inp_compare(.false., 'vectors', test)) then
           call vectors_input(rtdb,'prop')
c             
      elseif ( inp_compare(.false., 'iglo', test)) then
         iglo = 0
         iglo = 1
      elseif ( inp_compare(.false., 'giao', test)) then
         giao = 0
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
      elseif ( inp_compare(.false., 'dipole', test)) then
         dipole  = 0
      elseif ( inp_compare(.false., 'nbofile', test)) then
         nbofile = 1
      elseif ( inp_compare(.false., 'all', test)) then
         iglo = 0
         iglo = 1
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
     $          ('prop_input: invalid arbitrary point',555)
           else
             call errquit ('prop_input: not a valid center type',555)
           endif
           write(6,*) 'center is ',center
           write(6,*) 'point is ',point(1), point(2), point(3)
         else
            call errquit ('prop_input: no center information',555)
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
     $     call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:dipole', mt_int, 1,
     $                               dipole            ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:quadrupole', mt_int, 1,
     $                               quadrupole        ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:octupole', mt_int, 1, 
     $                               octupole          ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:mulliken', mt_int, 1,
     $                               mulliken          ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:spinpopulation', mt_int, 1,
     $                               spinpopulation    ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:esp', mt_int, 1,
     $                               esp               ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:efield', mt_int, 1,
     $                               efield            ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:efieldmap', mt_int, 1,
     $                               efieldmap         ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:efieldgrad', mt_int, 1,
     $                               efieldgrad        ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:efieldgradmap', mt_int, 1,
     $                               efieldgradmap     ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:electrondensity', mt_int, 1,
     $                               electrondensity   ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:electrondensitymap', mt_int, 1,
     $                               electrondensitymap))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:spindensity', mt_int, 1,
     $                               spindensity       ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:spindensitymap', mt_int, 1,
     $                               spindensitymap    ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:stoneanalysis', mt_int, 1,
     $                               stoneanalysis     ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:boyslocalization', mt_int, 1,
     $                               boyslocalization  ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:iglo', mt_int, 1,
     $                               iglo              ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:giao', mt_int, 1,
     $                               giao              ))
     $      call errquit('prop_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'prop:center',mt_int, 1, center))
     $      call errquit('prop_input: rtdb_put failed', 555)
      if (center.eq.4) then
        if (.not. rtdb_put(rtdb, 'prop:center_val',mt_dbl, 3, point))
     $      call errquit('prop_input: rtdb_put failed', 555)
      endif
c       
c>>>    store real parameters in rtdb
c       
cc        if (.not. rtdb_put(rtdb, 'inp:xxx', mt_dbl, 1, xxx))
cc     $       call errquit('dft_input: rtdb_put failed', 0)
        return
c       
10000   continue !  die with error
        write (6,*) 'prop_input: error in input file, last field: '
        write (6,'(a)') test
        call errquit('prop_input: error in input file', 0)
c       
        end
