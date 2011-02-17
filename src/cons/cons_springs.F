c $I#d: cons.F,v 1.1 2004/01/28 01:30:59 marat Exp $
      subroutine cons_unload_hbonds()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"
      integer nhb,h_rhb,h_khb,h_ijhb

      h_ijhb = cons_get_h_hbond_id()
      h_khb  = cons_get_h_hbond_k()
      h_rhb  = cons_get_h_hbond_r()
c
c     unload harmonic constraints if any
c
      call cons_get_hbond_nhb(nhb)
      if(nhb.gt.0) then
      if (.not. ma_free_heap(h_rhb) ) call errquit(
     &    'cons_bond_input: unable to free h_rhb',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_khb) ) call errquit(
     &    'cons_bond_input: unable to free h_khb',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_ijhb) ) call errquit(
     &    'cons_bond_input: unable to free h_ijhb',
     &    0, MA_ERR)
      call cons_set_hbond_nhb(0)
      end if
      end
c
      subroutine cons_load_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      integer nhb
      integer i_rhb,i_khb,i_ijhb
      integer h_rhb,h_khb,h_ijhb
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k
      integer ma_type,ma_n
c
      call cons_unload_hbonds()
c
      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,tag_n,mt_int,1,nhb)
c     
c     exit if no harm bonds
      if(.not.status .or. nhb.eq.0) then
c            call errquit(
c     >     'cons_load_hbonds: unable to get number of harm bonds:',
c     >      nhb, MA_ERR)
      return
      end if

      if ( .not. rtdb_ma_get( rtdb, tag_id,ma_type, ma_n, 
     &      h_ijhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_get_index(h_ijhb, 
     &      i_ijhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_k,ma_type, ma_n, 
     &      h_khb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_get_index(h_khb, 
     &      i_khb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_r,ma_type, ma_n, 
     &      h_rhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_get_index(h_rhb, 
     &      i_rhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
     
      
       call cons_set_hbond_nhb(nhb)
       call cons_set_h_hbond_id(h_ijhb)
       call cons_set_h_hbond_k(h_khb)
       call cons_set_h_hbond_r(h_rhb)

      end
c
      subroutine cons_delete_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      status =              rtdb_delete(rtdb,tag_n)
      status = status .and. rtdb_delete(rtdb,tag_id)
      status = status .and. rtdb_delete(rtdb,tag_k)
      status = status .and. rtdb_delete(rtdb,tag_r)

      end
c
      subroutine cons_add_hbond_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      logical status
      logical oprint
      integer inb
      integer iat,jat
      double precision r
      double precision r0,k
      double precision e,enb
      integer i_c,i_rhb,i_khb,i_ijhb
      integer nhb
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbond_nhb(nhb)
      i_c    = cons_get_i_c()
      i_ijhb = cons_get_i_hbond_id()
      i_khb  = cons_get_i_hbond_k()
      i_rhb  = cons_get_i_hbond_r()

      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)

       call cons_spring_energy(k,r0,
     >         dbl_mb(i_c+(iat-1)*3),
     >         dbl_mb(i_c+(jat-1)*3),
     >         r,
     >         enb)

       if(oprint) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,enb
       end if
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(*,*) "cons energy simulation"
        energy=e
      else 
        energy = energy + e
      end if
c
      call util_print_pop()
      return
      end
c      
      subroutine cons_add_hdihed_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      double precision gx(*)
c      
      logical status
      logical oprint
      integer inb
      integer iat,jat,kat,lat
      double precision r
      double precision r0,k
      double precision e,enb
      double precision f1(3),f2(3),f3(3),f4(3)
      integer i_c,i_r,i_k,i_id
      integer nc
      integer idm
      integer i
c
      idm = 4
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hdihed_n(nc)
      i_c    = cons_get_i_c()
      i_id = cons_get_i_hdihed_id()
      i_k  = cons_get_i_hdihed_k()
      i_r  = cons_get_i_hdihed_r()

      e=0.0d0
      do inb=1,nc
       iat=int_mb(i_id+idm*(inb-1))
       jat=int_mb(i_id+idm*(inb-1)+1)
       kat=int_mb(i_id+idm*(inb-1)+2)
       lat=int_mb(i_id+idm*(inb-1)+3)
       r0 =dbl_mb(i_r+inb-1)
       k  =dbl_mb(i_k+inb-1)

       call cons_dihed_force(k,r0,
     >         dbl_mb(i_c+(iat-1)*3),
     >         dbl_mb(i_c+(jat-1)*3),
     >         dbl_mb(i_c+(kat-1)*3),
     >         dbl_mb(i_c+(lat-1)*3),
     >         r,
     >         enb,
     >         f1,
     >         f2,
     >         f3,
     >         f4)

c       if(oprint) then
        write(6,*)"adding dihedral# ",inb
        write(6,*)" dihedral parameters (i,j,k,l,K,r0):",
     >         iat,jat,kat,lat,k,r0
        write(6,*)"  dihedral anlge and energy    :",r,enb
        write(6,*)"dihedral forces:"
        write(6,*) "f1:",(f1(i),i=1,3)
        write(6,*) "f2:",(f2(i),i=1,3)
        write(6,*) "f3:",(f3(i),i=1,3)
        write(6,*) "f4:",(f4(i),i=1,3)
c       end if
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(*,*) "cons energy simulation"
        do i=1,3
          gx((iat-1)*3+i)=f1(i)
          gx((jat-1)*3+i)=f2(i)
          gx((kat-1)*3+i)=f3(i)
          gx((lat-1)*3+i)=f4(i)
        end do
        energy=e
      else 
        do i=1,3
          gx((iat-1)*3+i)=gx((iat-1)*3+i)+f1(i)
          gx((jat-1)*3+i)=gx((jat-1)*3+i)+f2(i)
          gx((kat-1)*3+i)=gx((kat-1)*3+i)+f3(i)
          gx((lat-1)*3+i)=gx((lat-1)*3+i)+f4(i)
        end do
         energy = energy + e
      end if
c
      call util_print_pop()
      return
      end
c      
      subroutine cons_add_hdihed_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      logical status
      logical oprint
      integer inb
      integer iat,jat,kat,lat
      double precision r
      double precision r0,k
      double precision e,enb
      integer i_c,i_r,i_k,i_id
      integer nc
      integer idm
c
      idm = 4
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hdihed_n(nc)
      i_c    = cons_get_i_c()
      i_id = cons_get_i_hdihed_id()
      i_k  = cons_get_i_hdihed_k()
      i_r  = cons_get_i_hdihed_r()

      e=0.0d0
      do inb=1,nc
       iat=int_mb(i_id+idm*(inb-1))
       jat=int_mb(i_id+idm*(inb-1)+1)
       kat=int_mb(i_id+idm*(inb-1)+2)
       lat=int_mb(i_id+idm*(inb-1)+3)
       r0 =dbl_mb(i_r+inb-1)
       k  =dbl_mb(i_k+inb-1)

       call cons_dihed_energy(k,r0,
     >         dbl_mb(i_c+(iat-1)*3),
     >         dbl_mb(i_c+(jat-1)*3),
     >         dbl_mb(i_c+(kat-1)*3),
     >         dbl_mb(i_c+(lat-1)*3),
     >         r,
     >         enb)

       if(oprint) then
        write(6,*)"adding dihedral# ",inb
        write(6,*)" dihedral parameters (i,j,k,l,K,r0):",
     >         iat,jat,kat,lat,k,r0
        write(6,*)"  dihedral anlge and energy    :",r,enb
       end if
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(*,*) "cons energy simulation"
        energy=e
      else 
        energy = energy + e
      end if
c
      call util_print_pop()
      return
      end
c
      subroutine cons_add_hbond_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical status
      logical oprint
      integer i
      integer inb,nhb
      integer iat,jat
      double precision r
      double precision r0,k
      double precision e,f(3)
      integer i_c,i_rhb,i_khb,i_ijhb
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbond_nhb(nhb)
      i_c    = cons_get_i_c()
      i_ijhb = cons_get_i_hbond_id()
      i_khb  = cons_get_i_hbond_k()
      i_rhb  = cons_get_i_hbond_r()

c
      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)
c
       call cons_spring_force(k,r0,
     >           dbl_mb(i_c+(iat-1)*3),
     >           dbl_mb(i_c+(jat-1)*3),
     >           r,e,f)

       if(oprint) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,e
        write(6,*)"  spring forces               :",(f(i),i=1,3)
        write(6,*)"  spring deriv                :",
     >            -2.0*k*(r-r0)
       end if

      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
       do i=1,3
         gx((iat-1)*3+i)=f(i)
         gx((jat-1)*3+i)=-f(i)
       end do
       energy = e
      else 
       do i=1,3
         gx((iat-1)*3+i)=gx((iat-1)*3+i)+f(i)
         gx((jat-1)*3+i)=gx((jat-1)*3+i)-f(i)
       end do
       energy = energy+e
      end if

      end do
      call util_print_pop()

      return
      end  
c
      subroutine cons_spring_force(k,r0,r1,r2,r,energy,f)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,r0
      double precision r1(3)
      double precision r2(3)
      double precision r
      double precision energy
      double precision f(3)
c
      integer i
c
       r=(r1(1)-r2(1))**2+
     >   (r1(2)-r2(2))**2+
     >   (r1(3)-r2(3))**2
       r=sqrt(r)
       energy=k*(r-r0)**2

       do i=1,3
         f(i)=2*k*(r-r0)*
     &     (r1(i)-r2(i))/r
       end do

      return
      end

      subroutine cons_spring_energy(k,r0,r1,r2,r,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,r0
      double precision r1(3)
      double precision r2(3)
      double precision r
      double precision energy
c
      integer i
c
       r=(r1(1)-r2(1))**2+
     >   (r1(2)-r2(2))**2+
     >   (r1(3)-r2(3))**2
       r=sqrt(r)
       energy=k*(r-r0)**2
      return
      end

      subroutine cons_dihed_energy(k,phi0,r1,r2,r3,r4,phi,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,phi0
      double precision r1(3)
      double precision r2(3)
      double precision r3(3)
      double precision r4(3)
      double precision phi
      double precision energy
c
      integer i
c
c     calculate dihedral angle
c     ------------------------
      call cons_dihed(r1,r2,r3,r4,phi)
      write(*,*) "dihedral angle is : ",57.2957795*phi
      energy = 0.5d0*k*(phi-phi0)**2
      write(*,*) "dihedral energy is : ",energy
      return
      end

      subroutine cons_dihed_force(k,phi0,r1,r2,r3,r4,
     >                            phi,energy,f1,f2,f3,f4)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,phi0
      double precision r1(3)
      double precision r2(3)
      double precision r3(3)
      double precision r4(3)
      double precision phi
      double precision energy
      double precision f1(3)
      double precision f2(3)
      double precision f3(3)
      double precision f4(3)
c
      integer i
      double precision a
c
c     calculate dihedral angle
c     ------------------------
      call cons_dihed(r1,r2,r3,r4,phi)
      write(*,*) "dihedral angle is : ",57.2957795*phi
      energy = 0.5d0*k*(phi-phi0)**2
      write(*,*) "dihedral energy is : ",energy
      call cons_dihed_deriv(r1,r2,r3,r4,f1,f2,f3,f4)
      a = k*(phi-phi0)
      do i=1,3
        f1(i) = f1(i)*a
        f2(i) = f2(i)*a
        f3(i) = f3(i)*a
        f4(i) = f4(i)*a
      end do
      return
      end

      subroutine cons_load_hdihed(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      integer nhd
      integer i_rhd,i_khd,i_ijhd
      integer h_rhd,h_khd,h_ijhd
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k
      integer ma_type,ma_n
c
      call cons_unload_hdihed()
c
      call cons_hdihed_id_tag(namespace,tag_id)
      call cons_hdihed_n_tag(namespace,tag_n)
      call cons_hdihed_k_tag(namespace,tag_k)
      call cons_hdihed_r_tag(namespace,tag_r)
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,tag_n,mt_int,1,nhd)
c     
c     exit if no harm dihedrals
      if(.not.status .or. nhd.eq.0) then
c            call errquit(
c     >     'cons_load_hdiheds: unable to get number of harm bonds:',
c     >      nhd, MA_ERR)
      return
      end if

      if ( .not. rtdb_ma_get( rtdb, tag_id,ma_type, ma_n, 
     &      h_ijhd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)
      if ( .not. ma_get_index(h_ijhd, 
     &      i_ijhd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_k,ma_type, ma_n, 
     &      h_khd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)
      if ( .not. ma_get_index(h_khd, 
     &      i_khd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_r,ma_type, ma_n, 
     &      h_rhd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)
      if ( .not. ma_get_index(h_rhd, 
     &      i_rhd) ) call errquit(
     &      'cons_load_hdiheds: unable to allocate cons scratch space',
     &      2*nhd, MA_ERR)
     
      
       call cons_set_hdihed_n(nhd)
       call cons_set_h_hdihed_id(h_ijhd)
       call cons_set_h_hdihed_k(h_khd)
       call cons_set_h_hdihed_r(h_rhd)

      end
c
      subroutine cons_unload_hdihed()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"
      integer nhd,h_rhd,h_khd,h_idhd

      h_idhd = cons_get_h_hdihed_id()
      h_khd  = cons_get_h_hdihed_k()
      h_rhd  = cons_get_h_hdihed_r()
c
c     unload harmonic constraints if any
c
      call cons_get_hdihed_n(nhd)
      if(nhd.gt.0) then
      if (.not. ma_free_heap(h_rhd) ) call errquit(
     &    'cons_bond_input: unable to free h_rhd',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_khd) ) call errquit(
     &    'cons_bond_input: unable to free h_khd',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_idhd) ) call errquit(
     &    'cons_bond_input: unable to free h_ijhd',
     &    0, MA_ERR)
      call cons_set_hdihed_n(0)
      end if
      end
