      subroutine prop_input(rtdb)
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
      
      character*255 test
      integer nbofile

c     
c>>>  declarations from NWdft
c     
      integer infinity ! for infinite loop
      parameter (infinity=100000) 
c
c>>>  Default property settings.
c
      nbofile = 0
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
            if ( inp_compare(.false., 'nbofile', test)) then
c             
c>>>          print out nbo file for standalone nbo program.
c             
              nbofile = 1
            elseif ( inp_compare(.false., 'end', test)) then
c             
c>>>          end 
c             
              goto 99999
            else
              write(6,'(a,a,a)') ' prop_input: keyword ',test,
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
c>>>    place integer parameters in rtdb
c       
        if (.not. rtdb_put(rtdb, 'prop:nbofile', mt_int, 1, nbofile))
     $       call errquit('prop_input: rtdb_put failed', 0)
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

