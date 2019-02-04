      subroutine errquit(string, icode, errcode)
      implicit none
C$Id$
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      character*(*) string
      character*255 man_sec
      integer icode
      integer errcode
c     
c     error termination 
c     
c     string = error message printed to stdout
c     icode  = informative value printed to stdout
c
c     try to minimize the amount of output if all nodes 
c     detect an error
c     
      call ecce_print_module_entry('errquit')
      call ecce_print1_char('errquit message',string,1)
      call ecce_print1('errquit code', mt_int,icode,1)
      call ecce_print_module_exit('errquit', 'failed')
c
      man_sec = ' '
      if (ga_nodeid() .eq. 0) write(LuOut,1)
    1 format(1x,72('-'))
      write(LuOut, 2) string, icode
    2 format(1x,a,i20)
    3 format(/)
      write(LuOut,1)
      call inp_errout
      write(LuOut,1)
      write(LuOut,*) ' current input line : '
      call inp_outrec
      write(LuOut,1)
c
c       categorization of errors
c     
      write(LuOut,1) 
      if (errcode .eq. UNKNOWN_ERR) then
         write(LuOut,*) 'This error has not yet been assigned to a ',
     &        'category'
         man_sec = 'No section for this category'
      endif
c     
#if 0
      if (errcode .eq. MEM_ERR) then
         write(LuOut,*) 'You need to assign at least ', icode, 
     &        'bytes more memory'
         write(LuOut,3)
         write(LuOut,3)
         write(LuOut,*) 'Add a memory card of the following format: '
         write(LuOut,3)
         write(LuOut,*) 'memory total <default + ', icode, '>'
         man_sec = 'NWChem -> Top Level Directives -> Memory'
      endif
c
      if (errcode .eq. STACK_ERR) then
         write(LuOut,*) 'You need to assign at least ', icode,
     &        'bytes more stack memory'
         write(LuOut,3)
         write(LuOut,3)
         write(LuOut,*) 'Add a memory card of the following format: '
         write(LuOut,3)
         write(LuOut,*) 'memory stack <default + ', icode, '>'
         man_sec = 'NWChem -> Top Level Directives -> Memory'
      endif
c     
      if (errcode .eq. HEAP_ERR) then
         write(LuOut,*) 'You need to assign at least ', icode,
     &        'bytes more heap memory'
         write(LuOut,3)
         write(LuOut,3)
         write(LuOut,*) 'Add a memory card of the following format: '
         write(LuOut,3)
         write(LuOut,*) 'memory heap <default + ', icode, '>'
         man_sec = 'NWChem -> Top Level Directives -> Memory'
      endif
#endif
c     
      if (errcode .eq. RTDB_ERR) then
         write(LuOut,*) 'An error occured in the Runtime Database'
      endif
c     
      if (errcode .eq. BASIS_ERR) then
         write(LuOut,*) 'There is an error in the specified basis set'
      endif
c     
      if (errcode .eq. GEOM_ERR) then
         write(LuOut,*) 'There is an error related to the specified ',
     &        'geometry'
      endif
c
      if (errcode .eq. INPUT_ERR) then
         write(LuOut,*) 'There is an error in the input file'
      endif
c
      if (errcode .eq. CAPMIS_ERR) then
         write(LuOut,*) 'A feature requested has not yet been ',
     &        'implemented'
         write(LuOut,*) 
      endif
c     
#if 0
      if (errcode .eq. GA_ERR) then
         write(LuOut,*) 'You need to assign at least ', icode,
     &        ' bytes more global memory'
         man_sec = 'NWChem -> Top Level Directives -> Memory'
      endif
c     
      if (errcode .eq. MA_ERR) then
         write(LuOut,*) 'You need to assign at least ', icode,
     &        ' bytes more local memory'
         man_sec = 'NWChem -> Top Level Directives -> Memory'
      endif
#endif
c
      if (errcode .eq. INT_ERR) then
         write(LuOut,*) 'An error occured while computing integrals'
      endif
c
      if (errcode .eq. DISK_ERR) then
         write(LuOut,*) 'An error occured while trying to read or ',
     &        'write to disk space'
      endif
c
      if (errcode .eq. CALC_ERR) then
         write(LuOut,*) 'This type of error is most commonly ',
     &        'associated with calculations not reaching ',
     &        'convergence criteria'
      endif 
c
      write(LuOut,1) 
      write(LuOut,*) 'For more information see the NWChem manual at ',
     &     'https://github.com/nwchemgit/nwchem/wiki'
      write(LuOut,3)
      write(LuOut,*) 'For further details see manual section: ',
     &     man_sec

      call util_flush(LuOut)
      call ga_error(string, icode)
c
      end

      subroutine util_warning(string, icode, errcode)
      implicit none
C$Id$
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      character*(*) string
      character*255 man_sec
      integer icode
      integer errcode

      if (ga_nodeid() .eq. 0) 
     >    write(LuOut,*) string
      call util_flush(LuOut)
c
      end
