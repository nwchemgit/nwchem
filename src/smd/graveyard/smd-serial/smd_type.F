      subroutine smd_type_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_type,sp_param,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_type_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_type_init(sp_type)
c
      return
      end

      subroutine smd_type_init(sp_type)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_type
c
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_ta,i_tr,i_it
      integer i_tta,i_ttr
      logical result
c
      pname = "smd_type_init"
c
c      write(*,*) "in "//pname
c
      tag = "atom:name"
      call smd_get_ind_size(tag,i_ta,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)

      tag = "atom:resname"
      call smd_get_ind(tag,i_tr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:atomname"
      call smd_get_ind_size(tag,i_tta,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      if(nt.le.0)
     >  call errquit(
     >       pname//'no types ',nt, RTDB_ERR)

      tag = "param:resname"
      call smd_get_ind(tag,i_ttr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

c
c     create type structure
c     ---------------------
c
      call smd_namespace_create(sp_type)
      tag = "type:id"
      call smd_data_create_get(sp_type,tag,na,MT_INT,i_it)

      call smd_type_set(nt,
     +                  byte_mb(i_tta),
     +                  byte_mb(i_ttr),
     +                  na,
     +                  byte_mb(i_ta),
     +                  byte_mb(i_tr),
     +                  int_mb(i_it))


      return
      end

      subroutine smd_type_set(nt,tta,ttr,na,ta,tr,it)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "smd_const_data.fh"
      integer na,nt
      character*(smd_string_size) ta(na)
      character*(smd_string_size) tr(na)
      character*(smd_string_size) tta(nt)
      character*(smd_string_size) ttr(nt)
      integer it(na)
c
      integer i,j
      character*(smd_string_size) atag,rtag
      character*30 pname
      logical lr,la,match
 
      pname = "smd_type_set"

      do i=1,na
       match = .false.
       atag = ta(i)
       rtag = tr(i)
       do j=1,nt
        lr=rtag.eq.ttr(j)
        la=atag.eq.tta(j)
c        write(*,*) "check",atag,tta(j),rtag,ttr(j)
c        write(*,*) "check",la,lr
        if(la.and.lr) then
          match = .true.
          it(i)=j
c          write(*,*) "found match"
          goto 110
        end if
       end do
110    continue
       if(.not.match) then
        call util_warning
     >         ("no type defined for "//atag//"in "//rtag,0,0)
       end if
      end do

      return
      end

c $Id$
