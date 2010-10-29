      subroutine dia_getref
c
c $Id$
c
c     read the reference coordinates
c     ------------------------------
c
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      filref=card(8:80)
c
c     get the size from the reference file
c
      call dia_sizref(nsai,msai,nwmi,mwmi,nwai,mwai)
c
c     allocate memory and initialize
c
      call dia_alloc(.true.)
c
c     read the reference coordinates
c
      call dia_rdref(dbl_mb(i_xref))
c
      call dia_all(int_mb(i_ssel),1,dbl_mb(i_swt),nsel)
      lsel=.false.
      call dia_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call dia_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
c
      if(me.eq.0) then
      write(*,'(a,a)') ' Reference coordinates read from ',
     + filref(1:index(filref,' ')-1)
      write(*,'(/,a,i5,/)') ' Number of atoms is ',nsai
      endif
      if(.not.ma_verify_allocator_stuff()) print*,'Oopsi 2'
c
      return
      end
