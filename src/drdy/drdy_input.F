      subroutine drdy_input(rtdb)
      implicit none
*: typical drdy input deck
*:==============================================================================
*:drdy
*:------------------------------------------------------------------------------
*: *GENERAL
*:   TITLE
*:     Test run: H+H2 reaction, Page-McIver CUBE algorithm, no restart
*:   
*:   ATOMS
*:      1    H
*:      2    H
*:      3    H
*:   END
*: 
*: *REACT1
*:    GEOM
*:      1  0.0   0.0   0.0
*:      2  0.0   0.0   1.3886144
*:    END
*: 
*:    SPECIES LINRP
*: 
*: *REACT2
*:   GEOM
*:     3    0.0   0.0 190.3612132
*:   END
*: 
*:   SPECIES  ATOMIC
*: 
*: *PROD2
*:   GEOM
*:     1   0.0   0.0 190.3612132
*:   END
*:  
*:   SPECIES   ATOMIC
*: 
*: *PROD1
*: 
*:   GEOM
*:     2  0.0   0.0   1.3886144
*:     3   0.0   0.0   0.0
*:   END
*: 
*:   SPECIES  LINRP
*: 
*: *START
*: 
*:   GEOM
*:     1  0.0   0.0  -1.76531973
*:     2  0.0   0.0   0.0
*:     3  0.0   0.0   1.76531973
*:   END
*: 
*:   SPECIES  LINTS
*: 
*: *PATH
*:   SSTEP  0.05 
*:   SSAVE  0.05 
*:   SLP    0.50 
*:   SLM   -0.50
*:   SCALEMASS 0.6718993
*:   INTEGRA   CUBE
*:------------------------------------------------------------------------------
*:end:drdy
*:==============================================================================
*: 
*::includes::
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
*::passed:      
      integer rtdb ! [input] RTDB handle
*::locals:
      integer count ! number of lines read in block
      character*255 rtdbstring
      character*255 field
      logical block_end_found
*
      if (ga_nodeid().ne.0) return
      call inp_set_field(0) ! goto start of the line
      if (.not.inp_a(field)) call errquit
     &      ('drdy_input: no input present',911)
      if (.not.(inp_compare(.false.,'drdy',field))) call errquit
     &      ('drdy_input:fatal error: not the DRDY input block',911)
      count = 0
00100 continue
      if (inp_read()) then  ! read first line
        if (inp_a(field)) then ! read first line token
          if(inp_compare(.false.,'end',     field(1:3))) then
            goto 00200 ! done get out 
          elseif(inp_compare(.false.,'theory',field(1:6))) then
            if (.not.inp_a(field)) call errquit(
     &            'drdy_input: no theory on line',911)
            if (.not.rtdb_cput(rtdb,'drdy:theory',1,field))
     &            call errquit(
     &            'drdy_input: cput of theory failed',911)
            goto 00100
          elseif(inp_compare(.false.,'*',field(1:1))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'species',field(1:7))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'re-start',field(1:8))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'scalemass',field(1:9))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'sstep',field(1:5))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'ssave',field(1:5))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'slp',field(1:3))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'slm',field(1:3))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'saddle',field(1:6))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'nosaddle',field(1:8))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'sign',field(1:4))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'integra',field(1:8))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            goto 00100
          elseif(inp_compare(.false.,'title',field(1:5))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            if (inp_read()) then
              call drdy_input_add_line(rtdb,count) ! count incremented
            else
              call errquit('drdy_input: could not read title',911)
            endif
            goto 00100
          elseif(inp_compare(.false.,'geom',field(1:4))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            block_end_found = .false.
            do while (.not.block_end_found)
              if (inp_read()) then
                if (inp_a(field)) then ! read first line token              
                  block_end_found =
     &                  inp_compare(.false.,'end',field(1:3))
                endif
                call drdy_input_add_line(rtdb,count) ! count incremented
              else
                call errquit(
     &                'drdy_input: could not read geom block',911)
              endif
            enddo
            goto 00100
          elseif(inp_compare(.false.,'atoms',field(1:5))) then
            call drdy_input_add_line(rtdb,count)  ! count incremented
            block_end_found = .false.
            do while (.not.block_end_found)
              if (inp_read()) then
                if (inp_a(field)) then ! read first line token              
                  block_end_found =
     &                  inp_compare(.false.,'end',field(1:3))
                endif
                call drdy_input_add_line(rtdb,count) ! count incremented
              else
                call errquit(
     &                'drdy_input: could not read geom block',911)
              endif
            enddo
            goto 00100
          else
            write(luout,*)' line not recognized as part of drdy input'
            call errquit('drdy_input:fatal error: unexpected input',911)
          endif
        endif
      else
        call errquit('drdy_input:fatal error: unexpected eof',911)
      endif
*
00200 continue
*
      if (count.gt.0) then
        rtdbstring = 'drdy:nwchem:line:count'
        if (.not.rtdb_put(rtdb,rtdbstring,mt_int,1,count))
     &        call errquit
     &        ('drdy_input:fatal error: rtdb_put:count:failed',911)
      endif
*
      end
      subroutine drdy_input_add_line(rtdb,count)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
*
      integer rtdb
      integer count
*
      integer len_drdy
      character*255 rtdbstring
      character*255 drdyline
*
      call inp_set_field(0)
      drdyline = ' '
      if (.not.inp_line(drdyline)) call errquit
     &      ('drdy_input_add_line: inp_line failed',911)
      len_drdy = inp_strlen(drdyline)
      if (len_drdy.gt.80) call errquit
     &      ('drdy_input:fatal error: line gt 80 chars ',911)
      rtdbstring = ' '
      count = count + 1
      write(rtdbstring,10000)count
*      write(6,*)rtdbstring(1:inp_strlen(rtdbstring)),
*     &      '|||',drdyline(1:inp_strlen(drdyline))
      if (.not.rtdb_cput(rtdb,rtdbstring,1,drdyline))
     &      call errquit(
     &      'drdy_input_add_line:rtdb_cput failed on drdyline',911)
10000 format('drdy:',i10)
      end
