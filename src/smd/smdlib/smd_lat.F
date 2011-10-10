      subroutine smd_lat_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_lat
      character*32 tag,pname
      logical result

      pname = "smd_lat_init_system"
c
      tag = "lattice"
      call smd_system_get_component(sp_lat,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_lat_init(sp_lat,result)
      if(.not.result) then
       tag = "lattice"
       call smd_system_unset_component(tag)
      end if
c

      return
      end

      subroutine smd_lat_init(namespace,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      logical result
c
      integer rtdb
      character*32 pname
      character*80 tag
      integer i_lc,i_lrc,i_lfc
      double precision vol
c
      pname = "smd_lat_init"
c
c      write(*,*) "in "//pname
c
      call smd_rtdb_get_handle(rtdb)
c
c     check if there is any lattice in rtdb
      call smd_lat_rtdb_check(rtdb,result)
      if(.not.result) then
        call util_warning(
     >       pname//'no lattice found in rtdb',0,0)
        return
      end if
c
      call smd_namespace_create(namespace)
c
c     create lattice data structures
c     ------------------------------
      tag = "lat:fconst"
      call smd_data_create_get(namespace,tag,2,MT_DBL,i_lfc)
      tag = "lat:cell"
      call smd_data_create_get(namespace,tag,9,MT_DBL,i_lc)
      tag = "lat:rcell"
      call smd_data_create_get(namespace,tag,9,MT_DBL,i_lrc)
      
      call smd_lat_rtdb_read(rtdb,dbl_mb(i_lc))
      call smd_lat_invrt(dbl_mb(i_lc),dbl_mb(i_lrc))
      call smd_latt_vol(dbl_mb(i_lc),vol)
      dbl_mb(i_lfc) = vol
      return
      end

      subroutine smd_lat_rtdb_check(rtdb,olatt)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
      logical olatt
c
      double precision latt(3,3)
      character*32 pname
      character*80 tag
      double precision a(3)
      integer i
c
      pname = "smd_lat_rtdb_read"
c
c      write(*,*) "in "//pname
c
      olatt = .true.
      tag="smd:lat_a"
      call smd_rtdb_get_dbl(tag,3,a(1),olatt)
c      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
c     >      olatt=.false.

      return
      end

      subroutine smd_lat_rtdb_read(rtdb,latt)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      double precision latt(3,3)
      integer rtdb
c
      character*32 pname
      character*80 tag
      double precision a(3)
      integer i
      logical result
c
      pname = "smd_lat_rtdb_read"
c
c      write(*,*) "in "//pname
c
      tag="smd:lat_a"
      call smd_rtdb_get_dbl(tag,3,a(1),result)
      if (.not.result) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

c      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
c     >      call errquit(pname//'failed to get'//tag,0,
c     >       RTDB_ERR)
      do i=1,3
       latt(i,1)=a(i)
      end do
      tag="smd:lat_b"
      call smd_rtdb_get_dbl(tag,3,a(1),result)
      if (.not.result) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

c      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
c     >      call errquit(pname//'failed to get'//tag,0,
c     >       RTDB_ERR)
      do i=1,3
       latt(i,2)=a(i)
      end do
      tag="smd:lat_c"
      call smd_rtdb_get_dbl(tag,3,a(1),result)
      if (.not.result) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

c      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
c     >      call errquit(pname//'failed to get'//tag,0,
c     >       RTDB_ERR)
      do i=1,3
       latt(i,3)=a(i)
      end do
      return
      end

      subroutine smd_lat_invrt(latt,rlatt)
      implicit none
      double precision  latt(3,3),rlatt(3,3)
c
      double precision  det

      rlatt(1,1)=latt(2,2)*latt(3,3)-latt(3,2)*latt(2,3)
      rlatt(2,1)=latt(3,1)*latt(2,3)-latt(2,1)*latt(3,3)
      rlatt(3,1)=latt(2,1)*latt(3,2)-latt(3,1)*latt(2,2)
      rlatt(1,2)=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      rlatt(2,2)=latt(1,1)*latt(3,3)-latt(3,1)*latt(1,3)
      rlatt(3,2)=latt(3,1)*latt(1,2)-latt(1,1)*latt(3,2)
      rlatt(1,3)=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)
      rlatt(2,3)=latt(2,1)*latt(1,3)-latt(1,1)*latt(2,3)
      rlatt(3,3)=latt(1,1)*latt(2,2)-latt(2,1)*latt(1,2)
      
      det=latt(1,1)*rlatt(1,1)+latt(1,2)*rlatt(2,1)+latt(1,3)*rlatt(3,1)
      if(abs(det).gt.0.d0)det=1.d0/det
      
      rlatt(1,1)=det*rlatt(1,1)
      rlatt(2,1)=det*rlatt(2,1)
      rlatt(3,1)=det*rlatt(3,1)
      rlatt(1,2)=det*rlatt(1,2)
      rlatt(2,2)=det*rlatt(2,2)
      rlatt(3,2)=det*rlatt(3,2)
      rlatt(1,3)=det*rlatt(1,3)
      rlatt(2,3)=det*rlatt(2,3)
      rlatt(3,3)=det*rlatt(3,3)

      return

      end

      subroutine smd_latt_vol(latt,vol)
      implicit none
      real*8 x,y,z,latt,vol

      dimension latt(3,3)

      x=latt(2,2)*latt(3,3)-latt(2,3)*latt(2,3)
      y=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      z=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)

      vol=abs(latt(1,1)*x+latt(2,1)*y+latt(3,1)*z)

      return

      END

      subroutine smd_latt_get_vol(vol)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      double precision vol
c
      character*72 tag
      character*30 pname
      integer i_fconst
      logical result

      pname = "smd_latt_vol"
      tag = "lat:fconst"
      call smd_get_ind(tag,i_fconst,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting ntot '//tag,0, RTDB_ERR)
      vol = dbl_mb(i_fconst)

      return
      end

      subroutine smd_lat_rebox_alt(n,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer n
      double precision c(n,3)
c
      character*32 sp_lattice
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_lat_rebox"
c
c     get lattice params if any
c     -------------------------
      call smd_system_get_component(sp_lattice,"lattice",result)
      if(.not.result) then
        call util_warning(
     >       pname//'skipping reboxing as there is no lattice ',0,0)
        return
      end if

      tag = "lat:cell"
      call smd_data_get_index(sp_lattice,tag,i_lc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "lat:rcell"
      call smd_data_get_index(sp_lattice,tag,i_lrc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_util_rebox_alt(n,
     >                    dbl_mb(i_lc),
     >                    dbl_mb(i_lrc),
     >                    c)

      return
      end

      subroutine smd_lat_rebox(n,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer n
      double precision c(3,n)
c
      character*32 sp_lattice
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_lat_rebox"
c
c     get lattice params if any
c     -------------------------
      call smd_system_get_component(sp_lattice,"lattice",result)
      if(.not.result) then
        call util_warning(
     >       pname//'skipping reboxing as there is no lattice ',0,0)
        return
      end if

      tag = "lat:cell"
      call smd_data_get_index(sp_lattice,tag,i_lc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "lat:rcell"
      call smd_data_get_index(sp_lattice,tag,i_lrc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_util_rebox(n,
     >                    dbl_mb(i_lc),
     >                    dbl_mb(i_lrc),
     >                    c)

10    continue
      return
      end
C>    Folds coordinates back into the box
C>    based on center of coordinates
C>    so that "molecule" remains intact
C>    @param[in] n number of atoms
C>    @param[in,out] c[3,n] coordinates
      subroutine smd_lat_rebox_by_cog(n,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer n
      double precision c(3,n)
c
      character*32 sp_lattice
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_lat_rebox1"
c
c     get lattice params if any
c     -------------------------
      call smd_system_get_component(sp_lattice,"lattice",result)
      if(.not.result) then
        call util_warning(
     >       pname//'skipping reboxing as there is no lattice ',0,0)
        return
      end if

      tag = "lat:cell"
      call smd_data_get_index(sp_lattice,tag,i_lc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "lat:rcell"
      call smd_data_get_index(sp_lattice,tag,i_lrc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_util_rebox_by_cog(n,
     >                    dbl_mb(i_lc),
     >                    dbl_mb(i_lrc),
     >                    c)

10    continue
      return
      end

C>    Folds coordinates back into the box
C>    while preserving molecules/segments
C>    @param[in] nt number of atoms
C>    @param[in] n number of segments
C>    @param[in] in[n] starting index of given segment in coordinate array
C>    @param[in] nn[n] number of atoms in given segment
C>    @param[in,out] c[3,nt] coordinates
      subroutine smd_lat_rebox_segments(nt,n,in,nn,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer nt
      integer n
      integer in(n)
      integer nn(n)
      double precision c(3,nt)
c
      character*32 sp_lattice
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_lat_rebox1"
c
c     get lattice params if any
c     -------------------------
      call smd_system_get_component(sp_lattice,"lattice",result)
      if(.not.result) then
        call util_warning(
     >       pname//'skipping reboxing as there is no lattice ',0,0)
        return
      end if

      tag = "lat:cell"
      call smd_data_get_index(sp_lattice,tag,i_lc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "lat:rcell"
      call smd_data_get_index(sp_lattice,tag,i_lrc,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_util_rebox_segments(nt,n,in,nn,
     >                    dbl_mb(i_lc),
     >                    dbl_mb(i_lrc),
     >                    c)

10    continue
      return
      end
c $Id$
