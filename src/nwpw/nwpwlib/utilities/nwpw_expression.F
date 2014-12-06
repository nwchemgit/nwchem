
*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_compile                  *
*     *                                                        *
*     **********************************************************
*  
*   This function compiles the eqnstring into emachine code and stores it
* on the rtdb
*
*     Entry - rtdb: rtdb id
*             eqnstring: string that constains expression
*     Exit - returns an Id identifying the emachine
*
      integer function nwpw_expression_compile(rtdb,eqnstring)
      implicit none
      integer rtdb
      character*(*) eqnstring

#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      logical value,dprint
      real*8  fconst(50)
      integer code(5,200),i,nf,nc,print_level,nemach,n1
      character*80 rtdb_name

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      call util_print_get_level(print_level)
      dprint = (print_debug.le.print_level)

*     **** parse the eqnstring ****
      n1 = inp_strlen(eqnstring)
      call nwpw_emachine_parse(eqnstring,n1,nc,code,nf,fconst)

      if (dprint) then
         write(luout,'(A,A)') "eqnstring = ",trim(eqnstring)
         write(luout,*)
         do i=1,nc
            write(luout,'(A,I3,A,5I6)') "emachine: ln=",i-1," code=",
     >                             code(1,i),code(2,i),
     >                             code(3,i),code(4,i),code(5,i)
         end do
         write(luout,*)
         do i=1,nf
            write(luout,'("emachine: fconst(",I3,")=",E15.6)')
     >          i-1,fconst(i)
         end do
         write(luout,*)
         write(luout,*)
      end if

*     ***** read number of emachines on rtdb ****
      rtdb_name = 'nwpw:emachine:nemach'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nemach)) nemach = 0

*     ***** write emachine data to rtdb *****
      nemach = nemach + 1
      value = btdb_put(rtdb,rtdb_name,mt_int,1,nemach)

      rtdb_name = 'nwpw:emachine:eqnstring'//c_index_name(nemach)
      value = value.and.btdb_cput(rtdb,rtdb_name,1,eqnstring)

      rtdb_name = 'nwpw:emachine:nc'//c_index_name(nemach)
      value = value.and.btdb_put(rtdb,rtdb_name,mt_int,1,nc)

      rtdb_name = 'nwpw:emachine:nf'//c_index_name(nemach)
      value = value.and.btdb_put(rtdb,rtdb_name,mt_int,1,nf)

      rtdb_name = 'nwpw:emachine:code'//c_index_name(nemach)
      value = value.and.btdb_put(rtdb,rtdb_name,mt_int,5*nc,code)

      rtdb_name = 'nwpw:emachine:fconst'//c_index_name(nemach)
      value = value.and.btdb_put(rtdb,rtdb_name,mt_dbl,nf,fconst)
      if (.not.value)
     >   call errquit(
     >   'nwpw_expression_compile:error writing emachine',0,RTDB_ERR)

      nwpw_expression_compile = nemach
      return
      end

*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_reset                    *
*     *                                                        *
*     **********************************************************
      subroutine nwpw_expression_reset(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "errquit.fh"

      logical value
      character*80 rtdb_name

      rtdb_name = 'nwpw:emachine:nemach'
      value = rtdb_delete(rtdb,rtdb_name)

      return
      end



*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_start                    *
*     *                                                        *
*     **********************************************************
      subroutine nwpw_expression_start(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"

      logical started
      integer nc(2),nf(2),code(2),fconst(2)
      common /nwpw_expression/ nc,nf,code,fconst,started

*     **** local variables ****
      logical value
      integer i,nemach
      character*80 rtdb_name

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name
 
      rtdb_name = 'nwpw:emachine:started'
      if (.not.btdb_get(rtdb,rtdb_name,mt_log,1,started))
     >   started=.false.

*     ***** read number of emachines on rtdb ****
      rtdb_name = 'nwpw:emachine:nemach'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nemach)) nemach = 0

      if ((.not.started).and.(nemach.gt.0)) then
         started = .true.


*        ***** allocate and read pmeta,ameta,bmeta,parammeta, and ymeta ****
         value = BA_alloc_get(mt_int,nemach,'emachine:nc',nc(2),nc(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nemach,'emachine:nf',nf(2),nf(1))
         value = value.and.
     >           BA_alloc_get(mt_int,5*200*nemach,'emachine:code',
     >                        code(2),code(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,50*nemach,'emachine:fconst',
     >                        fconst(2),fconst(1))
         if (.not.value)
     >      call errquit(
     >      'nwpw_expression_start:out of heap',0,MA_ERR)

         do i=1,nemach
            rtdb_name = 'nwpw:emachine:nc'//c_index_name(i)
            value = value.and.
     >              btdb_get(rtdb,rtdb_name,mt_int,1,int_mb(nc(1)+i-1))

            rtdb_name = 'nwpw:emachine:nf'//c_index_name(i)
            value = value.and.
     >              btdb_get(rtdb,rtdb_name,mt_int,1,int_mb(nf(1)+i-1))

            rtdb_name = 'nwpw:emachine:code'//c_index_name(i)
            value = value.and.
     >              btdb_get(rtdb,rtdb_name,mt_int,
     >                       5*int_mb(nc(1)+i-1),
     >                       int_mb(code(1)+(i-1)*5*200))

            rtdb_name = 'nwpw:emachine:fconst'//c_index_name(i)
            value = value.and.
     >              btdb_get(rtdb,rtdb_name,mt_dbl,int_mb(nf(1)+i-1),
     >                       dbl_mb(fconst(1)+(i-1)*50))
         end do

*        **** write started ***
         rtdb_name = 'nwpw:emachine:started'
         value = value.and.btdb_put(rtdb,rtdb_name,mt_log,1,started)
         if (.not.value)
     >      call errquit(
     >      'nwpw_expression_start:starting emachine',0,RTDB_ERR)

      end if

      return
      end


*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_end                      *
*     *                                                        *
*     **********************************************************
      subroutine nwpw_expression_end(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

      logical started
      integer nc(2),nf(2),code(2),fconst(2)
      common /nwpw_expression/ nc,nf,code,fconst,started

*     **** local variables ****
      logical      value
      character*80 rtdb_name

      rtdb_name = 'nwpw:emachine:started'
      value = rtdb_delete(rtdb,rtdb_name)

      if (started) then
         value =           BA_free_heap(nc(2))
         value = value.and.BA_free_heap(nf(2))
         value = value.and.BA_free_heap(code(2))
         value = value.and.BA_free_heap(fconst(2))
         if (.not.value)
     >   call errquit(
     >   'nwpw_expression_end:error calling rtdb_delete',0,RTDB_ERR)
      end if

      return
      end

*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_eqnstring                *
*     *                                                        *
*     **********************************************************
*
*     returns eqnstring. 
*     Warning - If called from serial than btdb_parallel must me turned off
*
      subroutine nwpw_expression_eqnstring(rtdb,i,eqnstring)
      implicit none
      integer rtdb
      integer i
      character*(*) eqnstring

#include "btdb.fh"
#include "errquit.fh"

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:emachine:eqnstring'//c_index_name(i),
     >                    1,eqnstring))
     >   call errquit(
     >   'nwpw_expression_eqnstring:error calling btdb_cget',0,RTDB_ERR)

      return
      end

*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_f                        *
*     *                                                        *
*     **********************************************************
*
      real*8 function nwpw_expression_f(i,nion,rion)
      implicit none
      integer i
      integer nion
      real*8 rion(3,*)

#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"

      logical started
      integer nc(2),nf(2),code(2),fconst(2)
      common /nwpw_expression/ nc,nf,code,fconst,started

*     **** local variables ****
      real*8 f

*     **** external functions ****
      real*8   nwpw_emachine_f
      external nwpw_emachine_f

      f = 0.0
      if (started) then
         f = nwpw_emachine_f(int_mb(nc(1)+i-1),
     >                       int_mb(code(1)+(i-1)*5*200),
     >                       int_mb(nf(1)+i-1),
     >                       dbl_mb(fconst(1)+(i-1)*50),
     >                       nion,rion)
      end if

      nwpw_expression_f = f
      return
      end


*     **********************************************************
*     *                                                        *
*     *               nwpw_expression_fion                     *
*     *                                                        *
*     **********************************************************
*
      subroutine nwpw_expression_fion(alpha,i,nion,rion,fion)
      implicit none
      real*8 alpha
      integer i
      integer nion
      real*8 rion(3,*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"

      logical started
      integer nc(2),nf(2),code(2),fconst(2)
      common /nwpw_expression/ nc,nf,code,fconst,started

*     **** local variables ****
      integer i0,ii0
      real*8 df

*     **** external functions ****
      real*8   nwpw_emachine_df
      external nwpw_emachine_df

      if (started) then
         do ii0 = 1,nion
         do i0=1,3
            df = nwpw_emachine_df(i0,ii0,
     >                            int_mb(nc(1)+i-1),
     >                            int_mb(code(1)+(i-1)*5*200),
     >                            int_mb(nf(1)+i-1),
     >                            dbl_mb(fconst(1)+(i-1)*50),
     >                            nion,rion)
            fion(i0,ii0) = fion(i0,ii0) - alpha*df
         end do
         end do
      end if

      return
      end

