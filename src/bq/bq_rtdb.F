c $Id: bq_rtdb.F,v 1.3 2005-04-06 19:53:06 edo Exp $

      subroutine bq_tag_ncent(namespace,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 namespace
c
      n=inp_strlen(namespace)
      tag="bq:"//namespace(1:n)//":ncenter"

      end

      subroutine bq_tag_coord(namespace,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 namespace
c
      n=inp_strlen(namespace)
      tag="bq:"//namespace(1:n)//":coord"

      end

      subroutine bq_tag_charge(namespace,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 namespace
c
      n=inp_strlen(namespace)
      tag="bq:"//namespace(1:n)//":charge"

      end

      function bq_rtdb_store(irtdb,handle)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer irtdb
       integer handle
       logical bq_rtdb_store
c      local variables
       integer ncent
       integer i_c
       integer i_q
       character*(32) pname
       character*(255) tag
       character*(255) namespace

       logical  bq_ncenter
       external bq_ncenter

       logical  bq_check_handle
       external bq_check_handle

       logical  bq_namespace
       external bq_namespace

       logical  bq_index_coord
       external bq_index_coord

       logical  bq_index_charge
       external bq_index_charge

       pname = "bq_rtdb_store"

       bq_rtdb_store= .true.
       if(.not.bq_check_handle(handle)) then
         bq_rtdb_store = .false.
         return
       end if
         
       if(.not.bq_ncenter(handle,ncent))
     &       call errquit('bq_ncenter failed',0,0)
       if(.not.bq_namespace(handle,namespace))
     &       call errquit('bq_namespace failed',0,0)
       if(.not.bq_index_coord(handle,i_c))
     &       call errquit('bq_index_coord failed',0,0)
       if(.not.bq_index_charge(handle,i_q))
     &       call errquit('bq_index_charge failed',0,0)

       call bq_tag_ncent(namespace,tag)
       if(.not. rtdb_put(irtdb,tag,mt_int,1,ncent)) 
     &      call errquit(
     &      pname//' unable to store ncenter',
     &      0, RTDB_ERR)

       
       call bq_tag_coord(namespace,tag)
       if(.not. rtdb_put(irtdb,tag,mt_dbl,3*ncent,dbl_mb(i_c))) 
     &      call errquit(
     &      pname//' unable to store coord',
     &      0, RTDB_ERR)

       call bq_tag_charge(namespace,tag)
       if(.not. rtdb_put(irtdb,tag,mt_dbl,ncent,dbl_mb(i_q))) 
     &      call errquit(
     &      pname//' unable to store charge',
     &      0, RTDB_ERR)

       return
       end

      function bq_rtdb_load(irtdb,handle)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer irtdb
       integer handle
       logical bq_rtdb_load
c      local variables
       integer ncent
       integer h_c
       integer h_q
       character*(32) pname
       character*(255) tag
       character*(255) namespace
 
       integer ma_type,n

       logical  bq_ncenter
       external bq_ncenter

       logical  bq_check_handle
       external bq_check_handle

       logical  bq_namespace
       external bq_namespace

       logical  bq_pset
       external bq_pset

       pname = "bq_rtdb_load"

       bq_rtdb_load= .true.
       if(.not.bq_check_handle(handle)) then
         bq_rtdb_load = .false.
         return
       end if
         
       if(.not.bq_ncenter(handle,ncent))
     &       call errquit('bq_ncenter failed',0,0)
       if(ncent.ne.0) then
         bq_rtdb_load = .false.
         write(*,*) pname//"empty bq set first"
         return
       end if
        if(.not.bq_namespace(handle,namespace))
     &       call errquit('bq_namespace failed',0,0)

       call bq_tag_ncent(namespace,tag)
       if(.not. rtdb_get(irtdb,tag,mt_int,1,ncent)) 
     &      call errquit(
     &      pname//' unable to get ncenter',
     &      0, RTDB_ERR)

       call bq_tag_coord(namespace,tag)
       write(*,*) tag
       if(.not. rtdb_ma_get(irtdb,tag,ma_type,n,h_c))
     &      call errquit(
     &      pname//' uuunable to get coord',
     &      0, RTDB_ERR)


       call bq_tag_charge(namespace,tag)
       if(.not. rtdb_ma_get(irtdb,tag,ma_type,n,h_q))
     &      call errquit(
     &      pname//' unable to get charge',
     &      0, RTDB_ERR)

       if(.not. bq_pset(handle,ncent,h_q,h_c))
     &      call errquit(
     &      pname//' unable to pset bq',
     &      0, 0)

       
       return
       end

