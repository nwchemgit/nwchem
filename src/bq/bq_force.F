c $Id$


      subroutine bq_force(rtdb)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"
#include "bq_params.fh"
#include "stdio.fh"
c
      integer rtdb
c
      character*32 pname
      integer bq_handle
      integer bq_ncent
      logical ostatus
c
      integer i,j
      integer i_cbq,i_qbq
      integer h_f
      integer i_f
      character*255 fname
c
      logical  bq_ncenter
      external bq_ncenter
c
      logical  bq_index_coord
      external bq_index_coord
c
      logical  bq_index_charge
      external bq_index_charge
c
      logical  bq_destroy
      external bq_destroy
c
      logical  bq_get_active
      external bq_get_active
c   
      pname = "bq_force"
c
      call bq_force_status(rtdb,ostatus)
      if(.not.ostatus) return
c      
      if(.not.rtdb_cget(rtdb, "bq:force_filename",1,fname))
     >     call util_file_prefix("bqforce.dat",fname)
c
c     get bq information
c     ------------------
      if(.not.bq_get_active(bq_handle)) then
         return
      end if
      if(.not.bq_ncenter(bq_handle,bq_ncent))
     >   call errquit(pname//':no bq centers',0,0)
      if(.not.bq_index_coord(bq_handle,i_cbq))
     >   call errquit(pname//':no bq coords',0,0)
      if(.not.bq_index_charge(bq_handle,i_qbq))
     >   call errquit(pname//':no bq coords',0,0)
c
      if(.not.ma_push_get(mt_dbl,3*bq_ncent,'bq f',h_f,i_f))
     + call errquit(pname//': Failed to allocate memory for c',
     + 3*bq_ncent, MA_ERR)
c
      call qmmm_bq_forces(rtdb,bq_ncent,
     >                dbl_mb(i_cbq),
     >                dbl_mb(i_qbq),
     >                dbl_mb(i_f))
c
      call bq_force_print(bq_ncent,dbl_mb(i_f),fname)
      if(.not.ma_pop_stack(h_f))
     & call errquit(pname//': Failed to deallocate stack ',0,
     &       MA_ERR)
      end 

      subroutine bq_force_print(n,f,fname)
c
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
c
      integer n
      double precision f(3,n)
      character*(*) fname
c
      integer i,k

      integer un
      character*30 pname
      logical util_get_io_unit
      external util_get_io_unit

      pname = "bq_force_print"
      if(.not.util_get_io_unit(un))
     >    call errquit(pname//"no free file unit",0,0)

      open(unit=un,file=fname,
     >       form='formatted',status='unknown')


      write(un,'(A)')"#Forces on Bqs fx fy fz"
9498  format(3F12.6)
      do i=1,n
        write(un,9498) (f(k,i),k=1,3)
      end do
      close(un)
      end subroutine
