C> \ingroup nwdft
C> @{
C>
C> \file occup_input.F
C> Read the occupation numbers
C>
C> @}
C>
C> \ingroup nwdft
C> @{
C>
C> \brief Orbital occupations input reader
C>
C> Read the orbital occupation numbers either from the input or from
C> a file. An example of the input block is
C>
C> \code
C>   occup
C>     5   3
C>     1.0 1.0
C>     1.0 1.0
C>     1.0 1.0
C>     1.0
C>     1.0
C>   end
C> \endcode
C>
C> The first line with two integers specifies how many occupation
C> numbers there are to read for each spin channel. Next there are
C> a number of lines specifying the orbital occupations.
C>
C> Similary for reading the occupation numbers from a file
C>
C> \code
C>   occup
C>     5   3
C>     load file.occup
C>   end
C> \endcode
C>
C> After reading the occupation numbers they are stored on the 
C> runtime database in the fields:
C>
C> - `focc:occup_switch` -- there was an "occup" input block
C>
C> - `focc:occup` -- the number of alpha- and beta-occupation numbers
C>
C> - `focc:occup_list` -- the list of occupation numbers
C>
C> The occupation number list is stored essentially as read. E.g. for
C> the input block given above this list is stored as
C>
C> \code
C>   alpha(1) beta(1) alpha(2) beta(2) alpha(3) beta(3) alpha(4) alpha(5)
C> \endcode
C>
c     orbital occupations handler
c
      subroutine occup_input(rtdb)
c
C     $Id$

c     occup 
c      5 3
c      1.0 1.0
c      1.0 1.0
c      1.0 1.0
c      1.0
c      1.0
c      or
c      load file.occup
c     end
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer rtdb              !< [Input] The RTDB handle
      character*255 field       ! for character input
      integer offset,nocc1,nocc2    ! counts no. of centers as input
      logical status            ! scratch for return codes
      integer i,imin,imaj
      integer nfocc_in(2),nfocc ! input of nfocc
      integer l_occ,k_occ
      integer switch_focc
c
      character*255 fileoccup
      character*32 pname
      logical lfromfile
      pname = 'occup_input: '
      lfromfile = .false.
c     
c     read occup from the input deck and persist it in the rtdb
      if (ga_nodeid() .ne. 0) return
c     
c     Check that this is indeed a occupations line
c     
      call inp_set_field(0)     ! goto start of line
      if (.not. inp_a(field))
     $ call errquit(pname//'occup input error', 0, INPUT_ERR)
      if (.not. inp_compare(.false., 'occup', field))
     $ call errquit(pname//'occup input error', 0,INPUT_ERR)   
c
      nfocc_in(1) = 0  ! spin-1
      nfocc_in(2) = 0  ! spin-2
      status = .true.
      if (inp_read()) then
         if (.not. inp_i(nfocc_in(1))) nfocc_in(1) = 0
         if (.not. inp_i(nfocc_in(2))) nfocc_in(2) = 0
      end if ! inp_read()
c
      nfocc=nfocc_in(1)+nfocc_in(2) ! total fractional occupation
      if (nfocc.eq.0) return
      imin=2
      imaj=1   ! spin-1 is dominant
      if (nfocc_in(1).lt.nfocc_in(2)) then
        imin=1
        imaj=2  ! spin-2 is dominant
      endif
c
c     set the memory based on the number of fractional occupied states
      if (.not. ma_push_get(mt_dbl,nfocc,'occup', l_occ,k_occ))
     &    call errquit(pname//'ma_push_get fail k_occ',0,MA_ERR)  

c     read occupations <occ-A> <occ-B> sets
      nocc1=1
      nocc2=1
      offset = 0
 20   if (inp_read()) then
        status = inp_a(field)
        if (inp_compare(.false., 'end', field)) then
          goto 30
        else if (inp_compare(.false., 'load', field)) then
         if (.not. inp_a(fileoccup)) 
     &    call errquit(pname//'missing occups file',offset,INPUT_ERR)
          call occ_from_file(fileoccup,imin,imaj,nfocc,
     &              nfocc_in,dbl_mb(k_occ))
          lfromfile = .true.
        else
         if (nocc1 .gt. nfocc_in(imaj)) 
     &    call errquit(pname//'too many occupations?',offset,INPUT_ERR)
c        
c         input format  <occ-A>  <occ-B> 
          call inp_set_field(0)
          status = .true.
c
          if (nocc1 .le. nfocc_in(imin)) then
c
c           Reading occ-1 occ-2
            if ((.not.inp_f(dbl_mb(k_occ+offset))) .or. 
     &               (.not. inp_f(dbl_mb(k_occ+offset+1))))
     &          call errquit(pname//'<occ>',offset+1,INPUT_ERR)
            offset = offset + 2
            nocc1=nocc1+1
            nocc2=nocc2+1
          else if (offset/2+1 .gt. nfocc_in(imin)) then
c
c           Reading occ-1 (unpaired set)
            if (.not.inp_f(dbl_mb(k_occ+offset))) 
     &          call errquit(pname//'<occ>',offset+1,INPUT_ERR)
            if (imaj.eq.1) nocc1=nocc1+1
            if (imaj.eq.2) nocc2=nocc2+1
            offset = offset + 1
          end if  ! nocc1 check
        end if
        goto 20
       else
        call errquit(pname//'missing occupations',0,INPUT_ERR) 
       end if
c
 30   continue    
c
c     consistency check reading from input file
      if (.not.lfromfile) then
        if ((nocc1-1.ne.nfocc_in(1)) .or. (nocc2-1.ne.nfocc_in(2))) then
         write(luout,*) "nocc1-1: ",nocc1-1
         write(luout,*) "nfocc_in(1): ",nfocc_in(1)
         write(luout,*) "nocc2-1: ",nocc2-1
         write(luout,*) "nfocc_in(2): ",nfocc_in(2)
         call errquit(pname//'occupancy mismatch',0,INPUT_ERR)
        end if ! consistency check
      end if
c
c     store in rtdb 
      switch_focc=1 ! 1 means using occup keyword
c
      if (.not. rtdb_put(rtdb,'focc:occup_switch',mt_int,1,switch_focc))
     &  call errquit(pname//'rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup',mt_int,2,nfocc_in))
     &  call errquit(pname//'rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup_list',
     &    mt_dbl,nfocc,dbl_mb(k_occ)))
     &  call errquit(pname//'rtdb_put failed', 1201, RTDB_ERR)
      if (.not.ma_pop_stack(l_occ)) ! Free memory
     &  call errquit(pname//'ma_pop_stack failed k_occ',0,MA_ERR)
c
      end
c
c     load occupations from file
C>
C> \brief Read the orbital occupations from an external file
C>
C> Reads the occupation numbers from a file in exactly the same
C> format as read from the input file.
C>
      subroutine occ_from_file(fileoccup,imin,imaj,nfocc,nfocc_in,occ)

      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      character*255 fileoccup !< [Input] The filename
      integer imin  !< [Input] The minority spin channel
                    !< (1=alpha, 2=beta)
      integer imaj  !< [Input] The majority spin channel
      integer nfocc !< [Input] The total number of occupation numbers
      integer nfocc_in(2) !< [Input] The number of occupation numbers
                          !< per spin channel
      double precision occ(0:nfocc-1) !< [Output] The occupation numbers
c      
      character*32 pname
      integer unitno
c
      integer i,nocc1,nocc2,offset
c
      pname = 'occup_from_file: '
      unitno = 77
c
c     open occupations file
      open(unitno,file=fileoccup,err=1000)
c      write(6,*) "opening file"
c      write(6,*) "imaj: ",imaj
c      write(6,*) "imin: ",imin
c      write(6,*) imin,nfocc_in(imin)
c      write(6,*) imaj,nfocc_in(imaj)

c     read occupations <occ-A> <occ-B> sets
      nocc1=1
      nocc2=1
      offset = 0
      do i = 1,nfocc_in(imaj)
       if (nocc1 .le. nfocc_in(imin)) then
         read(unitno,*,err=1001,end=1001) occ(offset),occ(offset+1)
c         write(6,*) offset,occ(offset),occ(offset+1)
         offset = offset+2
         nocc1=nocc1+1
         nocc2=nocc2+1
       else if (offset/2+1 .gt. nfocc_in(imin)) then
         read(unitno,*,err=1001,end=1001) occ(offset)
c         write(6,*) offset,occ(offset)
         if (imaj.eq.1) nocc1=nocc1+1
         if (imaj.eq.2) nocc2=nocc2+1
         offset = offset+1
       end if
      end do
c
c     consistency check
      if ((nocc1-1.ne.nfocc_in(1)) .or. (nocc2-1.ne.nfocc_in(2))) then
         write(luout,*) "nocc1-1: ",nocc1-1
         write(luout,*) "nfocc_in(1): ",nfocc_in(1)
         write(luout,*) "nocc2-1: ",nocc2-1
         write(luout,*) "nfocc_in(2): ",nfocc_in(2)
         call errquit(pname//'occupancy mismatch',0,INPUT_ERR)
      end if ! consistency check
c
      close(unitno,err=1002)
c
      return
c
 1000 call errquit(pname//'failed to open file',0,-1)
 1001 call errquit(pname//'failed to read file',0,-1)
 1002 call errquit(pname//'failed to close file',0,-1)
c
      end
C> @}
