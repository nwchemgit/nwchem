      subroutine int_init(nbas, bases)
c $Id: int_init.F,v 1.6 1994-06-02 18:47:25 d3e129 Exp $
c
c initializes integral code to data structers for a integral computation
c
      implicit none
#include "bas.fh"
#include "apiP.fh"
c::passed
      integer nbas        ! number of basis sets to be used
      integer bases(nbas) ! basis set handles
c::local
      integer ibas, ang2use, angm
      logical status
c
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'int_init')
00100 continue

      if (.not.status) then
        write(6,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(6,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('int_init: basis handles hosed ',nbas)
      endif
      write(6,*)' int_init: basis set handles valid '
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('int_init: angm error',angm)
        ang2use = max(ang2use,angm)
00300 continue
c.      call defNxyz(ang2use)
c.. for derivatives
      call defNxyz(ang2use+1)
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(bases,nbas)
      init_int = 1
      end
      block data api_data
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
      implicit none
#include "apiP.fh"
c     
      data mem_2e4c  /0/
      data init_int  /0/
      data mem_1e    /0/
      data isz_2e4c  /0/
      data isz_1e    /0/
      data mem_1e_h1 /0/
      data isz_2e3c  /0/
      data isz_2e2c  /0/
      data mem_2e3c  /0/
      data mem_2e2c   /0/
      data isz_3ov   /0/
      data mem_3ov   /0/
      data DCexp     /0.0000000000D00/
      data DCcoeff   /1.0000000000D00/
      end
c
