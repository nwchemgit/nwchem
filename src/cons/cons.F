c $I#d: cons.F,v 1.1 2004/01/28 01:30:59 marat Exp $
      subroutine cons_load_geom(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_data.fh"
      integer rtdb
c
      integer geom
      integer nat
      logical status
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_init: no geometry ',72, RTDB_ERR)
c
c     get cart coordinates
c     --------------------
      status=geom_ncent(geom,nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if ( .not. ma_alloc_get( MT_DBL, nat*3, 'cons coord',
     &      h_c, i_c) ) call errquit(
     &      'input_set: unable to allocate cons coord',
     &      nat*3, MA_ERR)
      status=geom_cart_coords_get(geom, dbl_mb(i_c))
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
 
      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)

      return
      end
c
      subroutine cons_unload_geom()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons_data.fh"

      if ( .not. ma_free_heap(h_c)) call errquit(
     &      'input_set: unable to allocate cons coord',
     &      0, MA_ERR)

      return
      end

      subroutine cons_unload_hbonds()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_data.fh"
c     
c
c     unload harmonic constraints
c
      if (.not. ma_pop_stack(h_rhb) ) call errquit(
     &    'cons_bond_input: unable to free h_rhb',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_khb) ) call errquit(
     &    'cons_bond_input: unable to free h_khb',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_ijhb) ) call errquit(
     &    'cons_bond_input: unable to free h_ijhb',
     &    0, MA_ERR)

      end

      subroutine cons_load_hbonds(nhb,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_data.fh"
c     
      integer nhb
      integer rtdb
c
      logical status
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb)
      if(.not.status .or. nhb.eq.0) call errquit(
     >     'input_set: unable to get prior harm bonds',
     >      nhb, MA_ERR)
 
      if ( .not. ma_push_get( MT_INT, 2*nhb, 'h_ijhb',
     &      h_ijhb, i_ijhb) ) call errquit(
     &      'input_set: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_khb',
     &      h_khb, i_khb) ) call errquit(
     &      'input_set: unable to allocate cons scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_rhb',
     &      h_rhb, i_rhb) ) call errquit(
     &      'input_set: unable to allocate cons scratch space',
     &      nhb, MA_ERR)

        status = rtdb_get(rtdb,"cons:hbond_ij",
     >     mt_int,2*nhb,int_mb(i_ijhb)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_kij",
     >     mt_dbl,nhb,dbl_mb(i_khb)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_rij",
     >     mt_dbl,nhb,dbl_mb(i_rhb)) 
        if(.not.status) call errquit(
     >     'input_set: unable to get prior harm bonds',
     >      nhb, MA_ERR)

c      if(rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb)) then
c        status = rtdb_ma_get(rtdb,"cons:hbond_ij",
c     >     ma_type,ma_size,h_ijhb) 
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get cons:hbond_n',
c     >      nhb, RTDB_ERR)
c        status = ma_get_index(h_ijhb,i_ijhb)
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get ijhb',
c     >      nhb, RTDB_ERR)
c        status = rtdb_ma_get(rtdb,"cons:hbond_kij",
c     >     ma_type,ma_size,h_khb) 
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get h_khb',
c     >      nhb, RTDB_ERR)
c        status = ma_get_index(h_khb,i_khb)
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get i_khb',
c     >      nhb, RTDB_ERR)
c        status = rtdb_ma_get(rtdb,"cons:hbond_rij",
c     >     ma_type,ma_size,h_rhb) 
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get h_rhb',
c     >      nhb, RTDB_ERR)
c        status = ma_get_index(h_rhb,i_rhb)
c        if(.not.status) call errquit(
c     >     'cons_init: unable to get i_irhb',
c     >      nhb, RTDB_ERR)
c      end if
      end
c
      subroutine cons_load_fixed_atoms(nfp,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cons_data.fh"
c     
      integer nfp
      integer rtdb
c
      logical status
c
c     load fixed atoms
c
      status = rtdb_get(rtdb,"cons:fpos_i",
     >     mt_int,nfp,int_mb(i_fpos)) 
        if(.not.status) call errquit(
     >     'input_set: unable to get prior fixed atoms',
     >      nfp, MA_ERR)
      end
c
      subroutine cons_delete_hbonds(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      integer rtdb
c
      logical status
      status =              rtdb_delete(rtdb,"cons:hbond_n")
      status = status .and. rtdb_delete(rtdb,"cons:hbond_ij")
      status = status .and. rtdb_delete(rtdb,"cons:hbond_kij")
      status = status .and. rtdb_delete(rtdb,"cons:hbond_rij")
      if(.not.status) call errquit(
     >     'input_set: unable to delete prior harm bonds',
     >      0, MA_ERR)
      end
c
      subroutine cons_delete_fixed_atoms(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      integer rtdb
c
      logical status
      status = rtdb_delete(rtdb,"cons:fpos_i")
      if(.not.status) call errquit(
     >     'input_set: unable to delete prior fixed atoms',
     >      0, MA_ERR)
      end
c
      subroutine cons_add_hbond_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      double precision energy
      logical pdefault
      logical status
      integer i
      integer inb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e,enb
      integer nhb
c
      pdefault=ga_nodeid().eq.0
c
      call cons_load_hbonds(nhb,rtdb)
c
      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)
c
       r=0.0d0
       do i=1,3
         ciat(i)=dbl_mb(i_c+(iat-1)*3+i-1)
         cjat(i)=dbl_mb(i_c+(jat-1)*3+i-1)
         r=r+(ciat(i)-cjat(i))**2
       end do
       r=sqrt(r)
       enb=k*(r-r0)**2
       if(pdefault) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,enb
       end if
c
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(*,*) "cons energy simulation"
        energy=e
      else 
        energy = energy + e
      end if
c
      call cons_unload_hbonds(rtdb)
      return
      end
c
      subroutine cons_add_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      double precision energy
      double precision gx(*)
      logical status
      integer nhb
c
      if(.not. rtdb_get(rtdb, 'cons:active', mt_log, 1, status)) then
        return
      else if (.not.status) then
        return
      end if
c
      call cons_load_geom(rtdb)
c
      call util_print_centered(6,
     >     "NWChem Constraints Module",
     >     36,.true.)
c
      if(rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb)) then
        if(nhb .gt. 0 ) then
           call cons_add_hbond_egrad(rtdb,energy,gx)
        end if
      end if
c
      call cons_unload_geom()
c
      return
      end
c
      subroutine cons_add_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      double precision energy
      logical status
      integer nhb
c
      if(.not. rtdb_get(rtdb, 'cons:active', mt_log, 1, status)) then
        return
      else if (.not.status) then
        return
      end if
      call util_print_centered(6,
     >     "NWChem Constraints Module",
     >     36,.true.)
c
      call cons_load_geom(rtdb)
c
      if(rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb)) then
        if(nhb .gt. 0 ) then
           call cons_add_hbond_energy(rtdb,energy)
        end if
      end if
c
      call cons_unload_geom()
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
#include "cons_data.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical status
      logical pdefault
      integer i
      integer inb,nhb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e,f(3)
c
      pdefault=ga_nodeid().eq.0
c
      call cons_load_hbonds(nhb,rtdb)
c
      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)
c
       r=0.0d0
       do i=1,3
         ciat(i)=dbl_mb(i_c+(iat-1)*3+i-1)
         cjat(i)=dbl_mb(i_c+(jat-1)*3+i-1)
         r=r+(ciat(i)-cjat(i))**2
       end do
       r=sqrt(r)
c
       e=k*(r-r0)**2
       if(pdefault) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,e
       end if
c
       do i=1,3
         f(i)=2*k*(r-r0)*
     &     (ciat(i)-cjat(i))/r
       end do

       if(pdefault) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,e
        write(6,*)"  spring forces               :",(f(i),i=1,3)
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
         if(i.eq.1) write(6,*) "cons gradient",f
       end do
       energy = energy+e
      end if

      end do
      call cons_unload_hbonds(rtdb)

      return
      end  
c


