      subroutine int_init(nbas, bases)
c $Id: int_init.F,v 1.2 1994-04-04 20:27:36 d3e129 Exp $
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
      integer ibas
      logical status
c
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'int_init')
00100 continue

      if (.not.status) then
        write(6,*)' one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(6,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        return
      endif
      write(6,*)' int_init: basis set handles valid '
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(bases,nbas)
      init_int = 0
      end
      block data api_data
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
      implicit none
#include "apiP.fh"
c     
      data init_int /0/
      data mem_2e4c /0/
      data mem_1e   /0/
      data isz_2e4c /0/
      data isz_1e   /0/

c     
      end
      
