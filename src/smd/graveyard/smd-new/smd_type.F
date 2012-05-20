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


      call smd_type_init(sp_type,sp_atom,sp_param)
c
      return
      end

      subroutine smd_type_init(sp_type,sp_atom,sp_param)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_atom
      character*(*) sp_param
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
c     get total number of atoms and types
c     ---------------------------------------
      call smd_atom_ntot(sp_atom,na)
      call smd_param_ntot(sp_param,nt)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
      if(nt.le.0)
     >  call errquit(
     >       pname//'no types ',0, RTDB_ERR)
c
c     create type structure
c     ---------------------
c
      call smd_namespace_create(sp_type)
      tag = "type:id"
      call smd_data_create(sp_type,tag,na,MT_INT)
      call smd_data_get_index(sp_type,tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "atom:name"
      call smd_data_get_index(sp_atom,tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "atom:resname"
      call smd_data_get_index(sp_atom,tag,i_tr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:atomname"
      call smd_data_get_index(sp_param,tag,i_tta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:resname"
      call smd_data_get_index(sp_param,tag,i_ttr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)



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
      integer na,nt
      character*16 ta(na)
      character*16 tr(na)
      character*16 tta(nt)
      character*16 ttr(nt)
      integer it(na)
c
      integer i,j
      character*16 atag,rtag
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
