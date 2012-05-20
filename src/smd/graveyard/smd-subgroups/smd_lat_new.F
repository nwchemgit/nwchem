      subroutine smd_lat_init(namespace,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
              call FCAT_smd_lat(1)
c     
              call FCAT_smd_lat(2)
      character*(*) namespace
              call FCAT_smd_lat(3)
      logical result
              call FCAT_smd_lat(4)
c
              call FCAT_smd_lat(5)
      integer rtdb
              call FCAT_smd_lat(6)
      character*32 pname
              call FCAT_smd_lat(7)
      character*80 tag
              call FCAT_smd_lat(8)
      integer i_lc,i_lrc,i_lfc
              call FCAT_smd_lat(9)
      double precision vol
              call FCAT_smd_lat(10)
c
              call FCAT_smd_lat(11)
      pname = "smd_lat_init"
              call FCAT_smd_lat(12)
c
              call FCAT_smd_lat(13)
      write(*,*) "in "//pname
              call FCAT_smd_lat(14)
c
              call FCAT_smd_lat(15)
      call smd_rtdb_get_handle(rtdb)
              call FCAT_smd_lat(16)
c
              call FCAT_smd_lat(17)
c     check if there is any lattice in rtdb
              call FCAT_smd_lat(18)
      call smd_lat_rtdb_check(rtdb,result)
              call FCAT_smd_lat(19)
      if(.not.result) then
              call FCAT_smd_lat(20)
        call util_warning(
              call FCAT_smd_lat(21)
     >       pname//'no lattice found in rtdb',0,0)
              call FCAT_smd_lat(22)
        return
              call FCAT_smd_lat(23)
      end if
              call FCAT_smd_lat(24)
c
              call FCAT_smd_lat(25)
      call smd_namespace_create(namespace)
              call FCAT_smd_lat(26)
c
              call FCAT_smd_lat(27)
c     create lattice data structures
              call FCAT_smd_lat(28)
c     ------------------------------
              call FCAT_smd_lat(29)
      call smd_data_create(namespace,"lat:fconst",2,MT_DBL)
              call FCAT_smd_lat(30)
      call smd_data_create(namespace,"lat:cell",9,MT_DBL)
              call FCAT_smd_lat(31)
      call smd_data_create(namespace,"lat:rcell",9,MT_DBL)
              call FCAT_smd_lat(32)
c
              call FCAT_smd_lat(33)
c     get memory pointers
              call FCAT_smd_lat(34)
c     -------------------
              call FCAT_smd_lat(35)
      tag = "lat:cell"
              call FCAT_smd_lat(36)
      call smd_data_get_index(namespace,tag,i_lc,result)
              call FCAT_smd_lat(37)
      if(.not. result) 
              call FCAT_smd_lat(38)
     >  call errquit(
              call FCAT_smd_lat(39)
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
              call FCAT_smd_lat(40)
      tag = "lat:rcell"
              call FCAT_smd_lat(41)
      call smd_data_get_index(namespace,tag,i_lrc,result)
              call FCAT_smd_lat(42)
      if(.not. result) 
              call FCAT_smd_lat(43)
     >  call errquit(
              call FCAT_smd_lat(44)
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
              call FCAT_smd_lat(45)
      tag = "lat:fconst"
              call FCAT_smd_lat(46)
      call smd_data_get_index(namespace,tag,i_lfc,result)
              call FCAT_smd_lat(47)
      if(.not. result) 
              call FCAT_smd_lat(48)
     >  call errquit(
              call FCAT_smd_lat(49)
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      
              call FCAT_smd_lat(50)
      call smd_lat_rtdb_read(rtdb,dbl_mb(i_lc))
              call FCAT_smd_lat(51)
      call smd_lat_invrt(dbl_mb(i_lc),dbl_mb(i_lrc))
              call FCAT_smd_lat(52)
      call smd_latt_vol(dbl_mb(i_lc),vol)
              call FCAT_smd_lat(53)
      dbl_mb(i_lfc) = vol
              call FCAT_smd_lat(54)
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
              call FCAT_smd_lat(55)
c     
              call FCAT_smd_lat(56)
      integer rtdb
              call FCAT_smd_lat(57)
      logical olatt
              call FCAT_smd_lat(58)
c
              call FCAT_smd_lat(59)
      double precision latt(3,3)
              call FCAT_smd_lat(60)
      character*32 pname
              call FCAT_smd_lat(61)
      character*80 tag
              call FCAT_smd_lat(62)
      double precision a(3)
              call FCAT_smd_lat(63)
      integer i
              call FCAT_smd_lat(64)
c
              call FCAT_smd_lat(65)
      pname = "smd_lat_rtdb_read"
              call FCAT_smd_lat(66)
c
              call FCAT_smd_lat(67)
      write(*,*) "in "//pname
              call FCAT_smd_lat(68)
c
              call FCAT_smd_lat(69)
      olatt = .true.
              call FCAT_smd_lat(70)
      tag="smd:lat_a"
              call FCAT_smd_lat(71)
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
              call FCAT_smd_lat(72)
     >      olatt=.false.

              call FCAT_smd_lat(73)
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
              call FCAT_smd_lat(74)
c     
              call FCAT_smd_lat(75)
      double precision latt(3,3)
              call FCAT_smd_lat(76)
      integer rtdb
              call FCAT_smd_lat(77)
c
              call FCAT_smd_lat(78)
      character*32 pname
              call FCAT_smd_lat(79)
      character*80 tag
              call FCAT_smd_lat(80)
      double precision a(3)
              call FCAT_smd_lat(81)
      integer i
              call FCAT_smd_lat(82)
c
              call FCAT_smd_lat(83)
      pname = "smd_lat_rtdb_read"
              call FCAT_smd_lat(84)
c
              call FCAT_smd_lat(85)
      write(*,*) "in "//pname
              call FCAT_smd_lat(86)
c
              call FCAT_smd_lat(87)
      tag="smd:lat_a"
              call FCAT_smd_lat(88)
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
              call FCAT_smd_lat(89)
     >      call errquit(pname//'failed to get'//tag,0,
              call FCAT_smd_lat(90)
     >       RTDB_ERR)
              call FCAT_smd_lat(91)
      do i=1,3
              call FCAT_smd_lat(92)
       latt(i,1)=a(i)
              call FCAT_smd_lat(93)
      end do
              call FCAT_smd_lat(94)
      tag="smd:lat_b"
              call FCAT_smd_lat(95)
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
              call FCAT_smd_lat(96)
     >      call errquit(pname//'failed to get'//tag,0,
              call FCAT_smd_lat(97)
     >       RTDB_ERR)
              call FCAT_smd_lat(98)
      do i=1,3
              call FCAT_smd_lat(99)
       latt(i,2)=a(i)
              call FCAT_smd_lat(100)
      end do
              call FCAT_smd_lat(101)
      tag="smd:lat_c"
              call FCAT_smd_lat(102)
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,a(1))) 
              call FCAT_smd_lat(103)
     >      call errquit(pname//'failed to get'//tag,0,
              call FCAT_smd_lat(104)
     >       RTDB_ERR)
              call FCAT_smd_lat(105)
      do i=1,3
              call FCAT_smd_lat(106)
       latt(i,3)=a(i)
              call FCAT_smd_lat(107)
      end do
              call FCAT_smd_lat(108)
      return
      end

      subroutine smd_lat_invrt(latt,rlatt)
      implicit none
      double precision  latt(3,3),rlatt(3,3)
              call FCAT_smd_lat(109)
c
              call FCAT_smd_lat(110)
      double precision  det

              call FCAT_smd_lat(111)
      rlatt(1,1)=latt(2,2)*latt(3,3)-latt(3,2)*latt(2,3)
              call FCAT_smd_lat(112)
      rlatt(2,1)=latt(3,1)*latt(2,3)-latt(2,1)*latt(3,3)
              call FCAT_smd_lat(113)
      rlatt(3,1)=latt(2,1)*latt(3,2)-latt(3,1)*latt(2,2)
              call FCAT_smd_lat(114)
      rlatt(1,2)=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
              call FCAT_smd_lat(115)
      rlatt(2,2)=latt(1,1)*latt(3,3)-latt(3,1)*latt(1,3)
              call FCAT_smd_lat(116)
      rlatt(3,2)=latt(3,1)*latt(1,2)-latt(1,1)*latt(3,2)
              call FCAT_smd_lat(117)
      rlatt(1,3)=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)
              call FCAT_smd_lat(118)
      rlatt(2,3)=latt(2,1)*latt(1,3)-latt(1,1)*latt(2,3)
              call FCAT_smd_lat(119)
      rlatt(3,3)=latt(1,1)*latt(2,2)-latt(2,1)*latt(1,2)
      
              call FCAT_smd_lat(120)
      det=latt(1,1)*rlatt(1,1)+latt(1,2)*rlatt(2,1)+latt(1,3)*rlatt(3,1)
              call FCAT_smd_lat(121)
      if(abs(det).gt.0.d0)det=1.d0/det
      
              call FCAT_smd_lat(122)
      rlatt(1,1)=det*rlatt(1,1)
              call FCAT_smd_lat(123)
      rlatt(2,1)=det*rlatt(2,1)
              call FCAT_smd_lat(124)
      rlatt(3,1)=det*rlatt(3,1)
              call FCAT_smd_lat(125)
      rlatt(1,2)=det*rlatt(1,2)
              call FCAT_smd_lat(126)
      rlatt(2,2)=det*rlatt(2,2)
              call FCAT_smd_lat(127)
      rlatt(3,2)=det*rlatt(3,2)
              call FCAT_smd_lat(128)
      rlatt(1,3)=det*rlatt(1,3)
              call FCAT_smd_lat(129)
      rlatt(2,3)=det*rlatt(2,3)
              call FCAT_smd_lat(130)
      rlatt(3,3)=det*rlatt(3,3)

              call FCAT_smd_lat(131)
      return

      end

      subroutine smd_latt_vol(latt,vol)
      implicit none
      real*8 x,y,z,latt,vol

              call FCAT_smd_lat(132)
      dimension latt(3,3)

              call FCAT_smd_lat(133)
      x=latt(2,2)*latt(3,3)-latt(2,3)*latt(2,3)
              call FCAT_smd_lat(134)
      y=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
              call FCAT_smd_lat(135)
      z=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)

              call FCAT_smd_lat(136)
      vol=abs(latt(1,1)*x+latt(2,1)*y+latt(3,1)*z)

              call FCAT_smd_lat(137)
      return

      END

      subroutine smd_latt_get_vol(namespace,vol)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
              call FCAT_smd_lat(138)
c     
              call FCAT_smd_lat(139)
      character*(*) namespace
              call FCAT_smd_lat(140)
      double precision vol
              call FCAT_smd_lat(141)
c
              call FCAT_smd_lat(142)
      character*72 tag
              call FCAT_smd_lat(143)
      character*30 pname
              call FCAT_smd_lat(144)
      integer i_fconst
              call FCAT_smd_lat(145)
      logical result

              call FCAT_smd_lat(146)
      pname = "smd_latt_vol"
              call FCAT_smd_lat(147)
      tag = "lat:fconst"
              call FCAT_smd_lat(148)
      call smd_data_get_index(namespace,tag,i_fconst,result)
              call FCAT_smd_lat(149)
      if(.not. result) 
              call FCAT_smd_lat(150)
     >  call errquit(
              call FCAT_smd_lat(151)
     >       pname//'error getting ntot '//tag,0, RTDB_ERR)
              call FCAT_smd_lat(152)
      vol = dbl_mb(i_fconst)

              call FCAT_smd_lat(153)
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
              call FCAT_smd_lat(154)
c     
              call FCAT_smd_lat(155)
      integer n
              call FCAT_smd_lat(156)
      double precision c(n,3)
              call FCAT_smd_lat(157)
c
              call FCAT_smd_lat(158)
      character*32 sp_lattice
              call FCAT_smd_lat(159)
c
              call FCAT_smd_lat(160)
      character*72 tag
              call FCAT_smd_lat(161)
      character*30 pname
              call FCAT_smd_lat(162)
      integer na
              call FCAT_smd_lat(163)
      integer i_c,i_lrc,i_lc
              call FCAT_smd_lat(164)
      logical result

              call FCAT_smd_lat(165)
      pname = "smd_lat_rebox"
              call FCAT_smd_lat(166)
c
              call FCAT_smd_lat(167)
c     get lattice params if any
              call FCAT_smd_lat(168)
c     -------------------------
              call FCAT_smd_lat(169)
      call smd_system_get_component(sp_lattice,"lattice",result)
              call FCAT_smd_lat(170)
      if(.not.result) then
              call FCAT_smd_lat(171)
        call util_warning(
              call FCAT_smd_lat(172)
     >       pname//'skipping reboxing as there is no lattice ',0,0)
              call FCAT_smd_lat(173)
        return
              call FCAT_smd_lat(174)
      end if

              call FCAT_smd_lat(175)
      tag = "lat:cell"
              call FCAT_smd_lat(176)
      call smd_data_get_index(sp_lattice,tag,i_lc,result)
              call FCAT_smd_lat(177)
      if(.not. result)
              call FCAT_smd_lat(178)
     >  call errquit(
              call FCAT_smd_lat(179)
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

              call FCAT_smd_lat(180)
      tag = "lat:rcell"
              call FCAT_smd_lat(181)
      call smd_data_get_index(sp_lattice,tag,i_lrc,result)
              call FCAT_smd_lat(182)
      if(.not. result)
              call FCAT_smd_lat(183)
     >  call errquit(
              call FCAT_smd_lat(184)
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

              call FCAT_smd_lat(185)
      call smd_util_rebox(n,
              call FCAT_smd_lat(186)
     >                    dbl_mb(i_lc),
              call FCAT_smd_lat(187)
     >                    dbl_mb(i_lrc),
              call FCAT_smd_lat(188)
     >                    c)

              call FCAT_smd_lat(189)
      return
      end
       module FCAT_smd_lat_mod
        double precision,dimension (190):: &
     &   FCAT_smd_lat_counter = 0
        integer :: FCAT_smd_lat_nline = 189
       end module FCAT_smd_lat_mod
       subroutine FCAT_smd_lat(n)
         use FCAT_smd_lat_mod
         integer :: n
         FCAT_smd_lat_counter(n) = &
     &    FCAT_smd_lat_counter(n)  + 1
         if (FCAT_smd_lat_counter(n) == 1) then
           write(*,"(a,i10)") "FCAT_smd_lat_",n
         endif
       end subroutine FCAT_smd_lat
       subroutine FCAT_smd_lat_rep()
          use FCAT_smd_lat_mod
          integer :: i
          do i = 1, FCAT_smd_lat_nline
             write(*,"(a,i10,i10)") &
     &      "FCAT_smd_lat_count",i, &
     &      int(FCAT_smd_lat_counter(i)+0.1)
          end do
       end subroutine FCAT_smd_lat_rep
c $Id$
