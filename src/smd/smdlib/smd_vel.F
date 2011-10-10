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

      call smd_vel_init(sp_vel)
c
      return
      end

      subroutine smd_vel_init(sp_vel)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_vel
c
      character*32 pname
      integer na
      character*255 filename
      logical result
c
      pname = "smd_vel_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create vel data structures
c     ---------------------------
      call smd_namespace_create(sp_vel)
      call smd_data_create(sp_vel,"vel",3*na,MT_DBL)
c
      call smd_velocfile_input(filename,result)
c      call smd_rtdb_get_veloc_input(filename,result)
      if((.not.result).or.(filename.eq."random")) then
        call smd_vel_random_guess()
      else
        call smd_vel_read(filename)
      end if
c
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
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_get_ind_dim(tag,i_v,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

c
c     get mass array
c     ------------------
      tag = "param:mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get type array
c     --------------
      tag = "type:id"
      call smd_get_ind(tag,i_it,result)
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

      subroutine smd_vel_write(filename)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      character*(*) filename
c
      integer un
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
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="unknown",form="formatted",file=filename)
c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_get_ind_dim(sp_vel,tag,i_v,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

      call smd_vel_write0(un,
     >                   na,
     >                   dbl_mb(i_v))

      close(un)
      return
      end

      subroutine smd_vel_read(filename)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      character*(*) filename
c
      integer un
      character*72 sp_vel
      character*72 sp_param
      character*72 sp_type
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_v,i_it,n
      integer i_m
      logical result
c
      pname = "smd_vel_read"
c
c      write(*,*) "in "//pname
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
c
      open(unit=un,status="old",form="formatted",file=filename)

      call smd_get_ind_dim("vel",i_v,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na=na/3
      call smd_vel_read0(un,
     >                   na,
     >                   dbl_mb(i_v))

      close(un)
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

c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_get_ind_dim(tag,i_v,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

c
c     get mass array
c     ------------------
      tag = "param:mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get type array
c     --------------
      tag = "type:id"
      call smd_get_ind(tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      call smd_vel_ke0(na,
     >                     targetke,
     >                     dbl_mb(i_v),
     >                     dbl_mb(i_m),
     >                     int_mb(i_it))

      call smd_energy_set_component("kinetic",targetke)
      write(*,*) "current ke 1",targetke
      return
      end

      SUBROUTINE smd_vel_random0(natms,targetke,vvv,typmass,atmtype)

      implicit none

      integer natms
      double precision targetke
      double precision vvv(3,natms)
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
       vvv(1,i)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(2,i)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(3,i)=(x-0.5)/sqrt(typmass(iatm))

       comxvv=comxvv+vvv(1,i)*typmass(iatm)
       comyvv=comyvv+vvv(2,i)*typmass(iatm)
       comzvv=comzvv+vvv(3,i)*typmass(iatm)
       commass=commass+typmass(iatm)

      enddo

      comxvv=comxvv/commass
      comyvv=comyvv/commass
      comzvv=comzvv/commass

      do i=1,natms

       vvv(1,i)=vvv(1,i)-comxvv
       vvv(2,i)=vvv(2,i)-comyvv
       vvv(3,i)=vvv(3,i)-comzvv

      enddo

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $         +typmass(iatm)*((vvv(1,i)**2+vvv(2,i)**2+vvv(3,i)**2))

      enddo

      instanke=0.5*instanke

      xscale=sqrt(targetke/instanke)

      do i=1,natms

       vvv(1,i)=xscale*vvv(1,i)
       vvv(2,i)=xscale*vvv(2,i)
       vvv(3,i)=xscale*vvv(3,i)

      enddo

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $         +typmass(iatm)*((vvv(1,i)**2+vvv(2,i)**2+vvv(3,i)**2))

      enddo

      instanke=0.5*instanke

c      write(*,*) "current ke comp1",targetke,instanke

      return

      END

      SUBROUTINE smd_vel_read0(un,natms,vvv)

      implicit none

      integer un,natms
      double precision vvv(3,natms)
c
      integer i

      do i=1,natms

       read(un,*) vvv(1,i),vvv(2,i),vvv(3,i)

      enddo

      return

      END

      SUBROUTINE smd_vel_write0(un,natms,vvv)

      implicit none

      integer un,natms
      double precision vvv(3,natms)
c
      integer i

      do i=1,natms

       write(un,*) vvv(1,i),vvv(2,i),vvv(3,i)

      enddo

      return

      END

      SUBROUTINE smd_vel_ke0(natms,ke,vvv,typmass,atmtype)

      implicit none

      integer natms
      double precision targetke
      double precision vvv(3,natms)
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
     $          +typmass(iatm)*((vvv(1,i)**2+vvv(2,i)**2+vvv(3,i)**2))

      enddo

      ke=0.5*ke

      write(*,*) "kinetic energy",ke
      return

      END

      subroutine smd_velocfile_input(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_velocfile_input"

      result = .true.
      call smd_rtdb_get_string("smd:veloc:input",1,
     >                           filename,result)

      end
c $Id$
