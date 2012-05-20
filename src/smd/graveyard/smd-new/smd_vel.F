      subroutine smd_vel_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_vel,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_vel_init_system"




      write(*,*) "in "//pname
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_vel_init(sp_vel,sp_atom)
c      call smd_vel_set(sp_vel,sp_type,sp_param)
c
      return
      end

      subroutine smd_vel_init(sp_vel,sp_atom)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_vel
      character*(*) sp_atom
c
      character*32 pname
      integer na
c
      pname = "smd_vel_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(sp_atom,na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create vel data structures
c     ---------------------------
      call smd_namespace_create(sp_vel)
      call smd_data_create(sp_vel,"vel",3*na,MT_DBL)
      return
      end
c
      subroutine smd_vel_random_guess()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*72 sp_vel
      character*72 sp_param
      character*72 sp_type
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_v,i_it
      integer i_m
      logical result
      double precision targetke,ke
c
      pname = "smd_type_init"
c
      write(*,*) "in "//pname
c
c     get components
c     --------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_data_get_index(sp_vel,tag,i_v,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      call  smd_data_get_size(sp_vel,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

c
c     get mass array
c     ------------------
      tag = "param:mass"
      call smd_data_get_index(sp_param,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get type array
c     --------------
      tag = "type:id"
      call smd_data_get_index(sp_type,tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      call smd_temper_get_ke_target(targetke)
      if(targetke.eq.0.0d0) 
     >   call smd_temper_ke_compute(298.15,targetke)
      call smd_vel_random0(na,
     >                     targetke,
     >                     dbl_mb(i_v),
     >                     dbl_mb(i_m),
     >                     int_mb(i_it))

      call smd_vel_ke_compute(ke)
c      write(*,*) "current ke compared",targetke,ke
      return
      end

      subroutine smd_vel_write(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      integer un
c     
      character*72 sp_vel
      character*72 sp_param
      character*72 sp_type
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_v,i_it
      integer i_m
      logical result
      double precision targetke,ke
c
      pname = "smd_type_init"
c
c      write(*,*) "in "//pname
c
c     get components
c     --------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_data_get_index(sp_vel,tag,i_v,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      call  smd_data_get_size(sp_vel,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

      call smd_vel_write0(un,
     >                   na,
     >                   dbl_mb(i_v))

      return
      end

      subroutine smd_vel_read(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      integer un
c     
      character*72 sp_vel
      character*72 sp_param
      character*72 sp_type
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_v,i_it
      integer i_m
      logical result
      double precision targetke,ke
c
      pname = "smd_type_init"
c
c      write(*,*) "in "//pname
c
c     get components
c     --------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_data_get_index(sp_vel,tag,i_v,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      call  smd_data_get_size(sp_vel,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

      call smd_vel_read0(un,
     >                   na,
     >                   dbl_mb(i_v))

      return
      end

      subroutine smd_vel_ke_compute(targetke)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*72 sp_vel
      character*72 sp_param
      character*72 sp_type
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_v,i_it
      integer i_m
      logical result
      double precision targetke
c
      pname = "smd_type_init"
c
c      write(*,*) "in "//pname
c
c     get components
c     --------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_data_get_index(sp_vel,tag,i_v,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      call  smd_data_get_size(sp_vel,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

c
c     get mass array
c     ------------------
      tag = "param:mass"
      call smd_data_get_index(sp_param,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get type array
c     --------------
      tag = "type:id"
      call smd_data_get_index(sp_type,tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      call smd_vel_ke0(na,
     >                     targetke,
     >                     dbl_mb(i_v),
     >                     dbl_mb(i_m),
     >                     int_mb(i_it))

c      write(*,*) "current ke 1",targetke
      return
      end

      SUBROUTINE smd_vel_random0(natms,targetke,vvv,typmass,atmtype)

      implicit none

      integer natms
      double precision targetke
      double precision vvv(natms,3)
      double precision typmass(*)
      integer atmtype(*)
c
      integer i,iatm,iseed

      double precision   x
      double precision  commass,comxvv,comyvv,comzvv
      double precision  instanke,xscale

      iseed=620419483
      comxvv=0.d0
      comyvv=0.d0
      comzvv=0.d0
      commass = 0.0d0

      do i=1,natms

       iatm=atmtype(i)
       call tool_randm(iseed,x)
       vvv(i,1)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(i,2)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(i,3)=(x-0.5)/sqrt(typmass(iatm))

       comxvv=comxvv+vvv(i,1)*typmass(iatm)
       comyvv=comyvv+vvv(i,2)*typmass(iatm)
       comzvv=comzvv+vvv(i,3)*typmass(iatm)
       commass=commass+typmass(iatm)

      enddo

      comxvv=comxvv/commass
      comyvv=comyvv/commass
      comzvv=comzvv/commass

      do i=1,natms

       vvv(i,1)=vvv(i,1)-comxvv
       vvv(i,2)=vvv(i,2)-comyvv
       vvv(i,3)=vvv(i,3)-comzvv

      enddo

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $         +typmass(iatm)*((vvv(i,1)**2+vvv(i,2)**2+vvv(i,3)**2))

      enddo

      instanke=0.5*instanke

      xscale=sqrt(targetke/instanke)

      do i=1,natms

       vvv(i,1)=xscale*vvv(i,1)
       vvv(i,2)=xscale*vvv(i,2)
       vvv(i,3)=xscale*vvv(i,3)

      enddo

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $         +typmass(iatm)*((vvv(i,1)**2+vvv(i,2)**2+vvv(i,3)**2))

      enddo

      instanke=0.5*instanke

c      write(*,*) "current ke comp1",targetke,instanke

      return

      END

      SUBROUTINE smd_vel_read0(un,natms,vvv)

      implicit none

      integer un,natms
      double precision vvv(natms,3)
c
      integer i

      do i=1,natms

       read(un,*) vvv(i,1),vvv(i,2),vvv(i,3)

      enddo

      return

      END

      SUBROUTINE smd_vel_write0(un,natms,vvv)

      implicit none

      integer un,natms
      double precision vvv(natms,3)
c
      integer i

      do i=1,natms

       write(un,*) vvv(i,1),vvv(i,2),vvv(i,3)

      enddo

      return

      END

      SUBROUTINE smd_vel_ke0(natms,ke,vvv,typmass,atmtype)

      implicit none

      integer natms
      double precision targetke
      double precision vvv(natms,3)
      double precision typmass(*)
      integer atmtype(*)
c
      integer i,iatm

      double precision ke
c
      ke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       ke=ke
     $          +typmass(iatm)*((vvv(i,1)**2+vvv(i,2)**2+vvv(i,3)**2))

      enddo

      ke=0.5*ke

      return

      END


c $Id$
