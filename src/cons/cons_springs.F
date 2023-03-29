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
#include "stdio.fh"
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
        write(luout,*)"adding spring # ",inb
        write(luout,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(luout,*)"  spring length and energy    :",r,enb
       end if
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(luout,*) "cons energy simulation"
        energy=e
      else 
        energy = energy + e
      end if
c
      call util_print_pop()
      return
      end
c 

      subroutine cons_unload_hbondings()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"
      integer nhc
      integer h_hbondings_n0
      integer h_hbondings_indx
      integer h_hbondings_coef
      integer h_hbondings_k0
      integer h_hbondings_gamma0

      h_hbondings_n0     = cons_get_h_hbondings_n0()
      h_hbondings_indx   = cons_get_h_hbondings_indx()
      h_hbondings_coef   = cons_get_h_hbondings_coef()
      h_hbondings_k0     = cons_get_h_hbondings_k0()
      h_hbondings_gamma0 = cons_get_h_hbondings_gamma0()
c
c     unload harmonic constraints if any
c
      call cons_get_hbondings_nhc(nhc)
      if(nhc.gt.0) then
      if (.not. ma_free_heap(h_hbondings_n0) ) call errquit(
     &    'cons_bonding_input: unable to free h_hbondings_n0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_hbondings_indx) ) call errquit(
     &    'cons_bonding_input: unable to free h_hbondings_indx',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_hbondings_coef) ) call errquit(
     &    'cons_bonding_input: unable to free h_hbondings_coef',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_hbondings_k0) ) call errquit(
     &    'cons_bonding_input: unable to free h_hbondings_k0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_hbondings_gamma0) ) call errquit(
     &    'cons_bonding_input: unable to free h_hbondings_gamma0',
     &    0, MA_ERR)
      call cons_set_hbondings_nhc(0)
      end if
      end

c     
      subroutine cons_load_hbondings(namespace,rtdb)
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
      integer nhc
      integer h_hbondings_n0,i_hbondings_n0
      integer h_hbondings_indx,i_hbondings_indx
      integer h_hbondings_coef,i_hbondings_coef
      integer h_hbondings_k0,i_hbondings_k0
      integer h_hbondings_gamma0,i_hbondings_gamma0
      character*255 tag_n
      character*255 tag_n0
      character*255 tag_indx
      character*255 tag_coef
      character*255 tag_k0
      character*255 tag_gamma0
      integer ma_type,ma_n
c
      call cons_unload_hbondings()
c
      call cons_hbondings_n_tag(namespace,tag_n)

      call cons_hbondings_n0_tag(namespace,tag_n0)
      call cons_hbondings_indx_tag(namespace,tag_indx)
      call cons_hbondings_coef_tag(namespace,tag_coef)
      call cons_hbondings_k0_tag(namespace,tag_k0)
      call cons_hbondings_gamma0_tag(namespace,tag_gamma0)
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,tag_n,mt_int,1,nhc)
c     
c     exit if no harm bonds
      if(.not.status .or. nhc.eq.0) then
c            call errquit(
c     >     'cons_load_hbonds: unable to get number of harm bonds:',
c     >      nhb, MA_ERR)
      return
      end if

      if ( .not. rtdb_ma_get( rtdb, tag_n0,ma_type, ma_n,
     &      h_hbondings_n0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)
      if ( .not. ma_get_index(h_hbondings_n0,
     &      i_hbondings_n0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_indx,ma_type, ma_n,
     &      h_hbondings_indx) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)
      if ( .not. ma_get_index(h_hbondings_indx,
     &      i_hbondings_indx) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_coef,ma_type, ma_n,
     &      h_hbondings_coef) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)
      if ( .not. ma_get_index(h_hbondings_coef,
     &      i_hbondings_coef) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_k0,ma_type, ma_n,
     &      h_hbondings_k0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)
      if ( .not. ma_get_index(h_hbondings_k0,
     &      i_hbondings_k0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_gamma0,ma_type, ma_n,
     &      h_hbondings_gamma0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)
      if ( .not. ma_get_index(h_hbondings_gamma0,
     &      i_hbondings_gamma0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhc, MA_ERR)


       call cons_set_hbondings_nhc(nhc)
       call cons_set_hbondings_n0(h_hbondings_n0)
       call cons_set_hbondings_indx(h_hbondings_indx)
       call cons_set_hbondings_coef(h_hbondings_coef)
       call cons_set_hbondings_k0(h_hbondings_k0)
       call cons_set_hbondings_gamma0(h_hbondings_gamma0)

      return
      end

c
      subroutine cons_delete_hbondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      character*255 tag_n
      character*255 tag_n0
      character*255 tag_indx
      character*255 tag_coef
      character*255 tag_k0
      character*255 tag_gamma0

      call cons_hbondings_n_tag(namespace,tag_n)
      call cons_hbondings_n0_tag(namespace,tag_n0)
      call cons_hbondings_indx_tag(namespace,tag_indx)
      call cons_hbondings_coef_tag(namespace,tag_coef)
      call cons_hbondings_k0_tag(namespace,tag_k0)
      call cons_hbondings_k0_tag(namespace,tag_gamma0)

      status =              rtdb_delete(rtdb,tag_n)
      status = status .and. rtdb_delete(rtdb,tag_n0)
      status = status .and. rtdb_delete(rtdb,tag_indx)
      status = status .and. rtdb_delete(rtdb,tag_coef)
      status = status .and. rtdb_delete(rtdb,tag_k0)
      status = status .and. rtdb_delete(rtdb,tag_gamma0)

      end


c
      subroutine cons_add_hbondings_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      logical oprint
      integer i,ii,n,j
      double precision K0,gamma0,gamm
      double precision e,enb
      integer i_c,i_n,i_indx,i_coef,i_k0,i_gamma0
      integer nhc

c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbondings_nhc(nhc)
      i_c      = cons_get_i_c()
      i_n      = cons_get_i_hbondings_n0()
      i_indx   = cons_get_i_hbondings_indx()
      i_coef   = cons_get_i_hbondings_coef()
      i_k0     = cons_get_i_hbondings_k0()
      i_gamma0 = cons_get_i_hbondings_gamma0()

      e=0.0d0
      ii = 0
      do i=1,nhc
         n      = int_mb(i_n+i-1)
         K0     = dbl_mb(i_k0+i-1)
         gamma0 = dbl_mb(i_gamma0+i-1)

         call cons_bondings_energy(n,int_mb(i_indx+2*ii),
     >                               dbl_mb(i_coef+ii),
     >                               K0,gamma0,dbl_mb(i_c),gamm,enb)
       if(oprint) then
        write(luout,*)
        write(luout,'(A,I4)') "Adding bondings spring # ",i
        write(luout,'(A,F18.9,E14.6)') 
     >             "   - spring parameters (K0,gamma0):",
     >             K0,gamma0
        write(luout,'(A,E18.6)') 
     >  "   - gamma                        :",gamm
        write(luout,'(A,F18.9)') 
     >  "   - spring energy                :",enb
        write(luout,'(A,F18.9)') 
     >  "   - previous/base energy         :",energy
        write(luout,'(A,F18.9)') 
     >  "   - constrained energy           :",energy+enb
        write(luout,*) "  - coefficient index1 index2"
        do j=1,n
           write(luout,'(F16.6,I7,I7)') dbl_mb(i_coef+ii+j-1),
     >                           int_mb(i_indx+2*ii+2*(j-1)),
     >                           int_mb(i_indx+2*ii+2*(j-1)+1)
        end do
        write(luout,*)
       end if
       ii = ii + n
       e=e+enb
      end do
      energy = energy + e
c
      call util_print_pop()
      return
      end


c     
      subroutine cons_add_hbondings_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical oprint
      integer i,j,n,ii
      integer nhc
      double precision gamma0,K0,gamm
      double precision e,enb
      integer i_c,i_n,i_indx,i_coef,i_k0,i_gamma0
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbondings_nhc(nhc)
      i_c      = cons_get_i_c()
      i_n      = cons_get_i_hbondings_n0()
      i_indx   = cons_get_i_hbondings_indx()
      i_coef   = cons_get_i_hbondings_coef()
      i_k0     = cons_get_i_hbondings_k0()
      i_gamma0 = cons_get_i_hbondings_gamma0()


      e=0.0d0
      ii = 0
      do i=1,nhc
         n      = int_mb(i_n+i-1)
         K0     = dbl_mb(i_k0+i-1)
         gamma0 = dbl_mb(i_gamma0+i-1)

         call cons_bondings_force(n,int_mb(i_indx+2*ii),
     >                               dbl_mb(i_coef+ii),
     >                               K0,gamma0,dbl_mb(i_c),gamm,enb,gx)
       if(oprint) then
        write(luout,*)
        write(luout,'(A,I4)') "Adding bondings spring # ",i
        write(luout,'(A,F18.9,E14.6)')
     >             "   - spring parameters (K0,gamma0):",
     >             K0,gamma0
        write(luout,'(A,E18.6)')
     >  "   - gamma                        :",gamm
        write(luout,'(A,F18.9)')
     >  "   - spring energy                :",enb
        write(luout,'(A,F18.9)')
     >  "   - previous/base energy         :",energy
        write(luout,'(A,F18.9)')
     >  "   - constrained energy           :",energy+enb
        write(luout,*) "  - coefficient index1 index2"
        do j=1,n
           write(luout,'(F16.6,I7,I7)') dbl_mb(i_coef+ii+j-1),
     >                           int_mb(i_indx+2*ii+2*(j-1)),
     >                           int_mb(i_indx+2*ii+2*(j-1)+1)
        end do
        write(luout,*)
       end if

       ii = ii + n
       e=e+enb
      end do
      energy = energy + e
c
      call util_print_pop()
      return
      end



      subroutine cons_unload_pbondings()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"
      integer nhp
      integer h_pbondings_n0
      integer h_pbondings_indx
      integer h_pbondings_coef
      integer h_pbondings_k0
      integer h_pbondings_gcut0
      integer h_pbondings_gamma0
      integer h_pbondings_c0

      h_pbondings_n0     = cons_get_h_pbondings_n0()
      h_pbondings_indx   = cons_get_h_pbondings_indx()
      h_pbondings_coef   = cons_get_h_pbondings_coef()
      h_pbondings_k0     = cons_get_h_pbondings_k0()
      h_pbondings_gcut0  = cons_get_h_pbondings_gcut0()
      h_pbondings_gamma0 = cons_get_h_pbondings_gamma0()
      h_pbondings_c0     = cons_get_h_pbondings_c0()
c
c     unload harmonic constraints if any
c
      call cons_get_pbondings_nhp(nhp)
      if(nhp.gt.0) then
      if (.not. ma_free_heap(h_pbondings_n0) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_n0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_indx) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_indx',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_coef) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_coef',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_k0) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_k0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_gcut0) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_gcut0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_gamma0) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_gamma0',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_pbondings_c0) ) call errquit(
     &    'cons_bonding_input: unable to free h_pbondings_c0',
     &    0, MA_ERR)
      call cons_set_pbondings_nhp(0)
      end if
      end

c     
      subroutine cons_load_pbondings(namespace,rtdb)
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
      integer nhp
      integer h_pbondings_n0,i_pbondings_n0
      integer h_pbondings_indx,i_pbondings_indx
      integer h_pbondings_coef,i_pbondings_coef
      integer h_pbondings_k0,i_pbondings_k0
      integer h_pbondings_gcut0,i_pbondings_gcut0
      integer h_pbondings_gamma0,i_pbondings_gamma0
      integer h_pbondings_c0,i_pbondings_c0
      character*255 tag_n
      character*255 tag_n0
      character*255 tag_indx
      character*255 tag_coef
      character*255 tag_k0
      character*255 tag_gcut0
      character*255 tag_gamma0
      character*255 tag_c0
      integer ma_type,ma_n
c
      call cons_unload_pbondings()
c
      call cons_pbondings_n_tag(namespace,tag_n)

      call cons_pbondings_n0_tag(namespace,tag_n0)
      call cons_pbondings_indx_tag(namespace,tag_indx)
      call cons_pbondings_coef_tag(namespace,tag_coef)
      call cons_pbondings_k0_tag(namespace,tag_k0)
      call cons_pbondings_gcut0_tag(namespace,tag_gcut0)
      call cons_pbondings_gamma0_tag(namespace,tag_gamma0)
      call cons_pbondings_c0_tag(namespace,tag_c0)
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,tag_n,mt_int,1,nhp)
c     
c     exit if no harm bonds
      if(.not.status .or. nhp.eq.0) then
c            call errquit(
c     >     'cons_load_hbonds: unable to get number of harm bonds:',
c     >      nhb, MA_ERR)
      return
      end if

      if ( .not. rtdb_ma_get( rtdb, tag_n0,ma_type, ma_n,
     &      h_pbondings_n0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_n0,
     &      i_pbondings_n0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_indx,ma_type, ma_n,
     &      h_pbondings_indx) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_indx,
     &      i_pbondings_indx) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_coef,ma_type, ma_n,
     &      h_pbondings_coef) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_coef,
     &      i_pbondings_coef) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_k0,ma_type, ma_n,
     &      h_pbondings_k0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_k0,
     &      i_pbondings_k0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_gcut0,ma_type, ma_n,
     &      h_pbondings_gcut0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_gcut0,
     &      i_pbondings_gcut0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_gamma0,ma_type, ma_n,
     &      h_pbondings_gamma0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_gamma0,
     &      i_pbondings_gamma0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

      if ( .not. rtdb_ma_get( rtdb, tag_c0,ma_type, ma_n,
     &      h_pbondings_c0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)
      if ( .not. ma_get_index(h_pbondings_c0,
     &      i_pbondings_c0) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhp, MA_ERR)

       call cons_set_pbondings_nhp(nhp)
       call cons_set_pbondings_n0(h_pbondings_n0)
       call cons_set_pbondings_indx(h_pbondings_indx)
       call cons_set_pbondings_coef(h_pbondings_coef)
       call cons_set_pbondings_k0(h_pbondings_k0)
       call cons_set_pbondings_gcut0(h_pbondings_gcut0)
       call cons_set_pbondings_gamma0(h_pbondings_gamma0)
       call cons_set_pbondings_c0(h_pbondings_c0)

      return
      end

c
      subroutine cons_delete_pbondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      character*255 tag_n
      character*255 tag_n0
      character*255 tag_indx
      character*255 tag_coef
      character*255 tag_k0
      character*255 tag_gcut0
      character*255 tag_gamma0
      character*255 tag_c0

      call cons_pbondings_n_tag(namespace,tag_n)
      call cons_pbondings_n0_tag(namespace,tag_n0)
      call cons_pbondings_indx_tag(namespace,tag_indx)
      call cons_pbondings_coef_tag(namespace,tag_coef)
      call cons_pbondings_k0_tag(namespace,tag_k0)
      call cons_pbondings_gcut0_tag(namespace,tag_gcut0)
      call cons_pbondings_gamma0_tag(namespace,tag_gamma0)
      call cons_pbondings_c0_tag(namespace,tag_c0)

      status =              rtdb_delete(rtdb,tag_n)
      status = status .and. rtdb_delete(rtdb,tag_n0)
      status = status .and. rtdb_delete(rtdb,tag_indx)
      status = status .and. rtdb_delete(rtdb,tag_coef)
      status = status .and. rtdb_delete(rtdb,tag_k0)
      status = status .and. rtdb_delete(rtdb,tag_gcut0)
      status = status .and. rtdb_delete(rtdb,tag_gamma0)
      status = status .and. rtdb_delete(rtdb,tag_c0)

      end


c
      subroutine cons_add_pbondings_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      logical oprint
      integer i,ii,n,j
      double precision K0,gcut0,gamma0,c0,gamm
      double precision e,enb
      integer i_c,i_n,i_indx,i_coef,i_k0,i_gcut0,i_gamma0,i_c0
      integer nhp

c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_pbondings_nhp(nhp)
      i_c      = cons_get_i_c()
      i_n      = cons_get_i_pbondings_n0()
      i_indx   = cons_get_i_pbondings_indx()
      i_coef   = cons_get_i_pbondings_coef()
      i_k0     = cons_get_i_pbondings_k0()
      i_gcut0  = cons_get_i_pbondings_gcut0()
      i_gamma0 = cons_get_i_pbondings_gamma0()
      i_c0     = cons_get_i_pbondings_c0()

      e=0.0d0
      ii = 0
      do i=1,nhp
         n      = int_mb(i_n+i-1)
         K0     = dbl_mb(i_k0+i-1)
         gcut0  = dbl_mb(i_gcut0+i-1)
         gamma0 = dbl_mb(i_gamma0+i-1)
         c0     = dbl_mb(i_c0+i-1)

         call cons_pbondings_energy(n,int_mb(i_indx+2*ii),
     >                          dbl_mb(i_coef+ii),
     >                          K0,gcut0,gamma0,c0,dbl_mb(i_c),gamm,enb)
       if(oprint) then
        write(luout,*)
        write(luout,'(A,I4)') "Adding bondings spring # ",i
        write(luout,'(A,F18.9,F18.9,E14.6,F18.9)') 
     >             "   - spring parameters (K0,gcut0,gamma0,c0):",
     >             K0,gcut0,gamma0,c0
        write(luout,'(A,E18.6)') 
     >  "   - gamma                        :",gamm
        write(luout,'(A,F18.9)') 
     >  "   - spring energy                :",enb
        write(luout,'(A,F18.9)') 
     >  "   - previous/base energy         :",energy
        write(luout,'(A,F18.9)') 
     >  "   - constrained energy           :",energy+enb
        write(luout,*) "  - coefficient index1 index2"
        do j=1,n
           write(luout,'(F16.6,I7,I7)') dbl_mb(i_coef+ii+j-1),
     >                           int_mb(i_indx+2*ii+2*(j-1)),
     >                           int_mb(i_indx+2*ii+2*(j-1)+1)
        end do
        write(luout,*)
       end if
       ii = ii + n
       e=e+enb
      end do
      energy = energy + e
c
      call util_print_pop()
      return
      end


c     
      subroutine cons_add_pbondings_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical oprint
      integer i,j,n,ii
      integer nhp
      double precision gamma0,K0,gcut0,c0,gamm
      double precision e,enb
      integer i_c,i_n,i_indx,i_coef,i_k0,i_gcut0,i_gamma0,i_c0
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_pbondings_nhp(nhp)
      i_c      = cons_get_i_c()
      i_n      = cons_get_i_pbondings_n0()
      i_indx   = cons_get_i_pbondings_indx()
      i_coef   = cons_get_i_pbondings_coef()
      i_k0     = cons_get_i_pbondings_k0()
      i_gcut0  = cons_get_i_pbondings_gcut0()
      i_gamma0 = cons_get_i_pbondings_gamma0()
      i_c0     = cons_get_i_pbondings_c0()


      e=0.0d0
      ii = 0
      do i=1,nhp
         n      = int_mb(i_n+i-1)
         K0     = dbl_mb(i_k0+i-1)
         gcut0  = dbl_mb(i_gcut0+i-1)
         gamma0 = dbl_mb(i_gamma0+i-1)
         c0     = dbl_mb(i_c0+i-1)

         call cons_pbondings_force(n,int_mb(i_indx+2*ii),
     >                       dbl_mb(i_coef+ii),
     >                       K0,gcut0,gamma0,c0,dbl_mb(i_c),gamm,enb,gx)
       if(oprint) then
        write(luout,*)
        write(luout,'(A,I4)') "Adding bondings penalty # ",i
        write(luout,'(A,F18.9,F18.9,E14.6,F18.9)')
     >             "   - spring parameters (K0,gcut0,gamma0,c0):",
     >             K0,gcut0,gamma0,c0
        write(luout,'(A,E18.6)')
     >  "   - gamma                        :",gamm
        write(luout,'(A,F18.9)')
     >  "   - penalty energy               :",enb
        write(luout,'(A,F18.9)')
     >  "   - previous/base energy         :",energy
        write(luout,'(A,F18.9)')
     >  "   - penalized energy             :",energy+enb
        write(luout,*) "  - coefficient index1 index2"
        do j=1,n
           write(luout,'(F16.6,I7,I7)') dbl_mb(i_coef+ii+j-1),
     >                           int_mb(i_indx+2*ii+2*(j-1)),
     >                           int_mb(i_indx+2*ii+2*(j-1)+1)
        end do
        write(luout,*)
       end if

       ii = ii + n
       e=e+enb
      end do
      energy = energy + e
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
#include "stdio.fh"
#include "cons.fh"
#include "cons_data.fh"
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

      if(oprint) then  
        write(luout,*)"Processing dihedral constraints"
        write(luout,1004)
      end if
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

       if(oprint) then
        write(luout,1005) iat,jat,kat,lat,k,
     >                 r0*rad_to_deg,r*rad_to_deg,
     >                 enb,f1(1),f2(1),f3(1),f4(1)
        write(luout,1006) f1(2),f2(2),f3(2),f4(2)
        write(luout,1006) f1(3),f2(3),f3(3),f4(3)
c        write(6,1002) "dihedral forces (i,j,k,l):"
c        write(6,1001) "f_i", iat,(f1(i),i=1,3)
c        write(6,1001) "f_j", jat,(f2(i),i=1,3)
c        write(6,1001) "f_k", kat,(f3(i),i=1,3)
c        write(6,1001) "f_l", lat,(f4(i),i=1,3)
       end if
        do i=1,3
          gx((iat-1)*3+i)=gx((iat-1)*3+i)+f1(i)
          gx((jat-1)*3+i)=gx((jat-1)*3+i)+f2(i)
          gx((kat-1)*3+i)=gx((kat-1)*3+i)+f3(i)
          gx((lat-1)*3+i)=gx((lat-1)*3+i)+f4(i)
        end do
        energy = energy + enb
      end do
c
      call util_print_pop()
1002    FORMAT(1X,A,/,T6,"atom",T13,"fx",T24,"fy",T35,"fz")
1001    FORMAT(1X,A,1X,I4,3(1X,F10.6))
1004   FORMAT(T5,"i",4X,"j",4X,"k",4X,"l",T23,"Kphi",T31,"phi0",
     > T39,"phi",T47,"Energy",T55,"f1",T62,"f2",T69"f3",T76,"f4",/,
     > T5,76("_"))
1005   FORMAT( 4(1X,I4),3(1X,F8.3),1X,5(F7.3))
1006   FORMAT( T53,4F7.3)
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
#include "stdio.fh"
#include "cons.fh"
#include "cons_data.fh"
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

      if(oprint) then  
        write(luout,*)"Processing dihedral constraints"
        write(luout,1004)
      end if
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
         write(luout,1005) iat,jat,kat,lat,k,
     >                 r0*rad_to_deg,r*rad_to_deg,
     >                 enb
       end if

       energy=energy+enb
      end do
c
      call util_print_pop()
      return
1004  FORMAT(T5,"i",4X,"j",4X,"k",4X,"l",T23,"Kphi",T31,"phi0",
     >        T39,"phi",T47,"Energy",/,
     >        T5,49("_"))
1005  FORMAT( 4(1X,I4),3(1X,F7.3),1X,G12.6)
      end
c
      subroutine cons_add_hbond_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
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
        write(luout,*)"adding spring # ",inb
        write(luout,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(luout,*)"  spring length and energy    :",r,e
        write(luout,*)"  spring forces               :",(f(i),i=1,3)
        write(luout,*)"  spring deriv                :",
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
      double precision r,x,y,z
      double precision energy
      double precision f(3)
c
      integer i
c
      x = r1(1) - r2(1)
      y = r1(2) - r2(2)
      z = r1(3) - r2(3)
      call cons_lattice_min_diff(x,y,z)
      r=sqrt(x*x + y*y + z*z)
      energy=k*(r-r0)**2

      f(1)=2*k*(r-r0)*(x/r)
      f(2)=2*k*(r-r0)*(y/r)
      f(3)=2*k*(r-r0)*(z/r)

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
      double precision r,x,y,z
      double precision energy
c
      integer i
c
      x =  r1(1) - r2(1) 
      y =  r1(2) - r2(2) 
      z =  r1(3) - r2(3) 
      call cons_lattice_min_diff(x,y,z)
      r=sqrt(x*x + y*y + z*z)
      energy=k*(r-r0)**2
      return
      end


      subroutine cons_bondings_energy(n,indx,coef,K0,gamma0,rion,
     >                                gamm,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer n,indx(*)
      double precision coef(*),K0,gamma0
      double precision rion(3,*)
      double precision gamm,energy
c
      integer i,i0,i1
      double precision x,y,z,dist
c
      gamm = 0.0d0
      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i)
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         gamm = gamm +  coef(i)*dist
      end do
      energy = K0*(gamm-gamma0)**2

      return
      end

      subroutine cons_bondings_force(n,indx,coef,K0,gamma0,rion,
     >                               gamm,energy,f)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer n,indx(*)
      double precision coef(*),K0,gamma0
      double precision rion(3,*)
      double precision gamm,energy,f(3,*)
c
      integer i,ii,i0,i1
      double precision x,y,z,dist,denergy
c
      gamm = 0.0d0
      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i)
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         gamm = gamm +  coef(i)*dist
      end do
      energy  = K0*(gamm-gamma0)**2
      denergy = 2.0d0*K0*(gamm-gamma0)

      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i)
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         f(1,i0) = f(1,i0) + denergy*coef(i)*(x/dist)
         f(2,i0) = f(2,i0) + denergy*coef(i)*(y/dist)
         f(3,i0) = f(3,i0) + denergy*coef(i)*(z/dist)
         f(1,i1) = f(1,i1) - denergy*coef(i)*(x/dist)
         f(2,i1) = f(2,i1) - denergy*coef(i)*(y/dist)
         f(3,i1) = f(3,i1) - denergy*coef(i)*(z/dist)
      end do

      return
      end



      subroutine cons_pbondings_energy(n,indx,coef,K0,gcut0,gamma0,c0,
     >                                 rion,gamm,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer n,indx(*)
      double precision coef(*),K0,gcut0,gamma0,c0
      double precision rion(3,*)
      double precision gamm,energy
c
      integer i,i0,i1
      double precision x,y,z,dist
c
      gamm = 0.0d0
      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i) 
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         gamm = gamm +  coef(i)*dist
      end do
      energy = K0*exp(-((gamm-gamma0)/gcut0)**c0)

      return
      end



      subroutine cons_pbondings_force(n,indx,coef,K0,gcut0,gamma0,c0,
     >                                rion,gamm,energy,f)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer n,indx(*)
      double precision coef(*),K0,gcut0,gamma0,c0
      double precision rion(3,*)
      double precision gamm,energy,f(3,*)
c
      integer i,ii,i0,i1
      double precision x,y,z,dist,denergy
c
      gamm = 0.0d0
      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i)
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         gamm = gamm +  coef(i)*dist
      end do
      energy  = K0*exp(-((gamm-gamma0)/gcut0)**c0)
      denergy = -(c0*energy*(gamm-gamma0)**(c0-1.0d0)) / (gcut0**c0)

      do i=1,n
         i0 = indx(2*i-1)
         i1 = indx(2*i)
         x = rion(1,i0) - rion(1,i1)
         y = rion(2,i0) - rion(2,i1)
         z = rion(3,i0) - rion(3,i1)
         call cons_lattice_min_diff(x,y,z)
         dist = dsqrt(x*x + y*y + z*z)
         f(1,i0) = f(1,i0) + denergy*coef(i)*(x/dist)
         f(2,i0) = f(2,i0) + denergy*coef(i)*(y/dist)
         f(3,i0) = f(3,i0) + denergy*coef(i)*(z/dist)
         f(1,i1) = f(1,i1) - denergy*coef(i)*(x/dist)
         f(2,i1) = f(2,i1) - denergy*coef(i)*(y/dist)
         f(3,i1) = f(3,i1) - denergy*coef(i)*(z/dist)
      end do

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
      double precision diff
      double precision pi
c
      integer i
c
c     calculate dihedral angle
c     ------------------------
      call cons_dihed(r1,r2,r3,r4,phi,"rads")

      pi = acos(-1d0)
      diff = phi - phi0
      if(diff.gt.pi) then
        diff = diff - 2d0*pi
      elseif(diff.lt.-pi) then
        diff = diff + 2d0*pi
      endif
      energy = 0.5d0*k*diff**2
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
      double precision diff
      double precision pi
c
c     calculate dihedral angle
c     ------------------------
      call cons_dihed(r1,r2,r3,r4,phi,"rads")
      pi = acos(-1d0)
      diff = phi-phi0
      if(diff.gt.pi) then
        diff = diff - 2d0*pi
      elseif(diff.lt.-pi) then
        diff = diff + 2d0*pi
      endif

      energy = 0.5d0*k*diff**2
      call cons_dihed_deriv(r1,r2,r3,r4,f1,f2,f3,f4,"rads")
      a = k*diff
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
c $Id$
