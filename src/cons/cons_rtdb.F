c $Id$
      subroutine cons_create(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      call cons_set_namespace(namespace)
      call cons_set_irtdb(rtdb)

c     return if space already exists
      if(cons_inquire(namespace,rtdb)) return

      call cons_tag(namespace,tag)

c     contraints will be active by default
      if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit('cons_create', 0,
     >        RTDB_ERR)
      
c     no fixed atoms by default
      call cons_fxd_atoms_n_tag(namespace,tag)
      if(.not.rtdb_put(rtdb,tag,mt_int,1,0)) 
     >      call errquit('cons_create', 0,
     >        rtdb_err)

c     no springs by default
      call cons_hbond_n_tag(namespace,tag)
      if(.not.rtdb_put(rtdb,tag,mt_int,1,0)) 
     >      call errquit('cons_create', 0,
     >        rtdb_err)


      end

      subroutine cons_activate(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      if(.not.cons_inquire(namespace,rtdb)) then
           call errquit('nonexistent constraints', 0,
     >        0)

      end if

      if(.not. rtdb_cput(rtdb,"constraints" , 1, namespace)) then
         call errquit('qmmm: Failed constraints',0,0)
      end if
      

      end

      subroutine cons_get_namespace(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
#include "cons_data.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      if(.not. rtdb_cget(rtdb,"constraints" , 1, namespace)) then
         namespace=cons_default
      end if
      

      end

      function cons_inquire(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_inquire
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      integer n

      call cons_fxd_atoms_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return

      call cons_hbond_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return

      call cons_hbondings_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return


      call cons_pbondings_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return

      call cons_hdihed_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return

      
      return
      end

      subroutine cons_add_fixed_atoms(namespace,nat,id,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
      integer nat
      integer id(nat)
      integer rtdb
c
      logical status
      integer i,j
      integer n0
      integer n
      integer h_tmp_id,i_tmp_id 
      
      character*255 tag_id
      character*255 tag_n


      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      call cons_sort_i(nat,id)
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n0)) then
         n0 = 0
      end if
      n = n0 + nat

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       if(n0 .gt. 0) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,n0,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
        status =  rtdb_delete(rtdb,tag_id)
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to delete prior fixed atoms',
     >      0, MA_ERR)
       end if
       do i=1,n0
        do j=1,nat
           if(int_mb(i_tmp_id + i-1).eq.id(j)) then
             id(j)=-1
           end if
        end do
       end do
       
       j=0
       do i=1,nat
         if(id(i) .ne. -1) then
          int_mb(i_tmp_id + n0+j)=id(i)
          j=j+1
         end if
       end do
       n=n0+j

       status = rtdb_put(rtdb,tag_n,mt_int,1,n)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,n,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior harm bonds',
     >      n, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_add_spring(namespace,id1,id2,k,r,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c     
      character*(*) namespace
      integer id1,id2
      double precision k
      double precision r
      integer rtdb
c
      logical status
      integer i
      integer nhb0
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      integer id01,id02
      double precision k0
      double precision r0
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb0)) then
         nhb0 = 0
      end if
      nhb = nhb0 + 1

      if ( .not. ma_push_get( MT_INT, 2*nhb, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)

       if(nhb0 .ge. 1) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,2*nhb0,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nhb0,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nhb0,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       status =              rtdb_delete(rtdb,tag_n)
       status = status .and. rtdb_delete(rtdb,tag_id)
       status = status .and. rtdb_delete(rtdb,tag_k)
       status = status .and. rtdb_delete(rtdb,tag_r)
       if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bonds',
     >      0, MA_ERR)
       end if
       
      status = .true.
      do i=1,nhb0
       id01 = int_mb(i_tmp_id+2*(i-1))
       id02 = int_mb(i_tmp_id+2*(i-1)+1)
       k0 = dbl_mb(i_tmp_k+i-1)
       r0 = dbl_mb(i_tmp_r+i-1)
       if(id01.eq.id1 .and. id02.eq.id2) then
         write(luout,*) "cons_add_spring: replacing existing spring"
         dbl_mb(i_tmp_k+i-1)=k
         dbl_mb(i_tmp_r+i-1)=r
         status = .false.
       end if
      end do
        
      if(status) then
       int_mb(i_tmp_id + 2*nhb0)     = id1
       int_mb(i_tmp_id + 2*nhb0 + 1) = id2
       dbl_mb(i_tmp_k + nhb0)=k
       dbl_mb(i_tmp_r + nhb0)=r
      else
       nhb=nhb0
      end if



       status = rtdb_put(rtdb,tag_n,mt_int,1,nhb)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_put(rtdb,tag_k,
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_put(rtdb,tag_r,
     >     mt_dbl,nhb,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_add_dihed(namespace,id,kd,r,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c     
      character*(*) namespace
      integer id(4)
      double precision kd
      double precision r
      integer rtdb
c
      logical status
      integer i
      integer nc0
      integer nc
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      integer id01,id02
      double precision k0
      double precision r0
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k
      integer idm
      logical oduplicate
      integer k
c
c     4 indices for dihedral
      idm = 4
c
c     rtdb tags
c
      call cons_hdihed_id_tag(namespace,tag_id)
      call cons_hdihed_n_tag(namespace,tag_n)
      call cons_hdihed_k_tag(namespace,tag_k)
      call cons_hdihed_r_tag(namespace,tag_r)
      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nc0)) then
         nc0 = 0
      end if
      nc = nc0 + 1

      if ( .not. ma_push_get( MT_INT, idm*nc, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nc, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nc, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nc, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nc, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nc, MA_ERR)

       if(nc0 .ge. 1) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,idm*nc0,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nc0,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nc0,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nc, MA_ERR)
       status =              rtdb_delete(rtdb,tag_n)
       status = status .and. rtdb_delete(rtdb,tag_id)
       status = status .and. rtdb_delete(rtdb,tag_k)
       status = status .and. rtdb_delete(rtdb,tag_r)
       if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bonds',
     >      0, MA_ERR)
       end if
       
c     check for duplicate
      oduplicate = .false.
      do i=1,nc0
       do k=1,idm
         if(id(k).ne.int_mb(i_tmp_id+idm*(i-1)+k-1)) goto 10
       end do
       write(luout,*) "cons_add_dihedral: replacing existing spring"
       dbl_mb(i_tmp_k+i-1)=kd
       dbl_mb(i_tmp_r+i-1)=r
       oduplicate = .true.
       goto 12
10     continue    
      end do
12    continue
c      
      if(oduplicate) then
       nc=nc0
      else
       do k=1,idm
         int_mb(i_tmp_id+idm*(i-1)+k-1)=id(k)
       end do
       dbl_mb(i_tmp_k + nc0)=kd
       dbl_mb(i_tmp_r + nc0)=r
      end if

       status = rtdb_put(rtdb,tag_n,mt_int,1,nc)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,idm*nc,int_mb(i_tmp_id)) 
        status = status . and. rtdb_put(rtdb,tag_k,
     >     mt_dbl,nc,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_put(rtdb,tag_r,
     >     mt_dbl,nc,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nc, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

*     ************************************
*     *                                  *
*     *      cons_add_bondings           *
*     *                                  *
*     ************************************

      subroutine cons_add_bondings(namespace,n,indx,coef,K0,gamma0,rtdb)
      implicit none
      character*(*) namespace
      integer n,indx(*)
      double precision coef(*),K0,gamma0
      integer rtdb

#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical status
      integer i
      integer nhc0,mm
      integer nhc
      integer h_tmp_n0,i_tmp_n0
      integer h_tmp_indx,i_tmp_indx
      integer h_tmp_coef,i_tmp_coef
      integer h_tmp_k0,i_tmp_k0
      integer h_tmp_gamma0,i_tmp_gamma0
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
      call cons_hbondings_gamma0_tag(namespace,tag_gamma0)


      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhc0)) then
         nhc0 = 0
      end if
      nhc = nhc0 + 1
      mm = 0

      if(nhc0.eq.0) then
         if (.not. ma_push_get( MT_INT, nhc, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)
         if (.not. ma_push_get( MT_INT, nhc*2*n, 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc*2*n, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhc*n, 'h_tmp_coef',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc*n, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhc, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhc, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)

      else
         if ( .not. ma_push_get( MT_INT, nhc, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)
         status = rtdb_get(rtdb,tag_n0,mt_int,nhc0,int_mb(i_tmp_n0))
         mm = 0
         do i=1,nhc0
            mm = mm + int_mb(i_tmp_n0+i-1)
         end do
         if ( .not. ma_push_get( MT_INT, 2*(mm+n), 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*(mm+n), MA_ERR)
         if ( .not. ma_push_get( MT_DBL, (mm+n), 'h_tmp_indx',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      (mm+n), MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhc, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhc, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhc, MA_ERR)

         status = status . and. rtdb_get(rtdb,tag_indx,
     >     mt_int,2*mm,int_mb(i_tmp_indx))
         status = status . and. rtdb_get(rtdb,tag_coef,
     >     mt_dbl,mm,dbl_mb(i_tmp_coef))
         status = status . and. rtdb_get(rtdb,tag_k0,
     >     mt_dbl,nhc0,dbl_mb(i_tmp_k0))
         status = status . and. rtdb_get(rtdb,tag_gamma0,
     >     mt_dbl,nhc0,dbl_mb(i_tmp_gamma0))
         if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bondings',
     >      nhc, MA_ERR)

         status =              rtdb_delete(rtdb,tag_n)
         status = status .and. rtdb_delete(rtdb,tag_n0)
         status = status .and. rtdb_delete(rtdb,tag_indx)
         status = status .and. rtdb_delete(rtdb,tag_coef)
         status = status .and. rtdb_delete(rtdb,tag_k0)
         status = status .and. rtdb_delete(rtdb,tag_gamma0)
         if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bondingss',
     >      0, MA_ERR)
      end if

      int_mb(i_tmp_n0+nhc-1) = n
      do i=1,2*n
         int_mb(i_tmp_indx+2*mm+i-1) = indx(i)
      end do
      do i=1,n
         dbl_mb(i_tmp_coef+mm+i-1) = coef(i)
      end do
      dbl_mb(i_tmp_k0+nhc-1) = K0
      dbl_mb(i_tmp_gamma0+nhc-1) = gamma0

       status = rtdb_put(rtdb,tag_n,mt_int,1,nhc)
       status = status.and. rtdb_put(rtdb,tag_n0,
     >     mt_int,nhc,int_mb(i_tmp_n0))
       status = status.and. rtdb_put(rtdb,tag_indx,
     >     mt_int,2*(mm+n),int_mb(i_tmp_indx))
       status = status.and. rtdb_put(rtdb,tag_coef,
     >     mt_dbl,(mm+n),dbl_mb(i_tmp_coef))
       status = status.and. rtdb_put(rtdb,tag_k0,
     >     mt_dbl,nhc,dbl_mb(i_tmp_k0))
       status = status.and. rtdb_put(rtdb,tag_gamma0,
     >     mt_dbl,nhc,dbl_mb(i_tmp_gamma0))

       if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bondings',
     >      nhc, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_gamma0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_gamma0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_k0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_coef) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_coef',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_indx) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_indx',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_n0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_n0',
     &    0, MA_ERR)

      return
      end

      subroutine cons_unfix_atoms(namespace,nat,id,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
      integer nat
      integer id(nat)
      integer rtdb
c
      logical status
      integer i,j
      integer n0
      integer h_tmp_id,i_tmp_id 

      character*255 tag_id
      character*255 tag_n

      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n0)) then
         n0=0
      end if
      if(n0.eq.0) then
         return
      end if

      if ( .not. ma_push_get( MT_INT, n0, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      n0, MA_ERR)

        status = rtdb_get(rtdb,tag_id,
     >     mt_int,n0,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
        status =  rtdb_delete(rtdb,tag_id)
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to delete prior fixed atoms',
     >      0, MA_ERR)

       do i=1,n0
        do j=1,nat
           if(int_mb(i_tmp_id + i-1).eq.id(j)) then
             int_mb(i_tmp_id + i-1)=-1
           end if
        end do
       end do

       call cons_sort_i(n0,int_mb(i_tmp_id))
       
       j=0
       do i=1,n0
         if(int_mb(i_tmp_id + i-1).eq.-1) then
           j=j+1
         end if
       end do


       status = rtdb_put(rtdb,tag_n,mt_int,1,n0-j)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,n0-j,int_mb(i_tmp_id+j)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior harm bonds',
     >      n0, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      integer i
      integer n
      integer h_tmp_id,i_tmp_id 

      character*255 tag_id
      character*255 tag_n

      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n)) then
         n = 0
      end if
      if (n .eq. 0 ) then
         if(ga_nodeid().eq.0)
     >     write(luout,*) "no fixed atoms found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_info_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       status = rtdb_get(rtdb,tag_id,
     >     mt_int,n,int_mb(i_tmp_id)) 
       if(.not.status) call errquit(
     >     'cons_info_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
       
       if(ga_nodeid().eq.0) then
          write(luout,'(" fixed atoms",T31,":",(T32,4(I4,2X),:))')
     >         (int_mb(i_tmp_id +i-1),i=1,n)
          write(luout,'(/)')
       end if

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_springs(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) namespace
      integer rtdb
c     
      logical status
      integer i
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      character*255 tag_id
      character*255 tag_n
      character*255 tag_k
      character*255 tag_r

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) then
         nhb = 0
      end if
       
      if (nhb .eq. 0 ) then
         if(ga_nodeid().eq.0)    
     >     write(luout,*) "no springs found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, 2*nhb, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)

       if(nhb .ge. 1) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nhb,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       end if

     
      if(ga_nodeid().eq.0) then   
      do i=1,nhb
      write(luout,FMT=200) int_mb(i_tmp_id+2*(i-1)),
     >           int_mb(i_tmp_id+2*(i-1)+1),
     >           dbl_mb(i_tmp_k+i-1),
     >           dbl_mb(i_tmp_r+i-1)
      end do
      end if
200   FORMAT(" spring bond (i,j,k,r0)",T29,I6,2X,I6,2X,F6.3,2X,F6.3)
      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_id',
     &    0, MA_ERR)
      end


      subroutine cons_info_bondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) namespace
      integer rtdb
c     
      logical status
      integer i,ii,mm,j,n
      integer nhc
      integer h_tmp_n0,i_tmp_n0
      integer h_tmp_indx,i_tmp_indx
      integer h_tmp_coef,i_tmp_coef
      integer h_tmp_k0,i_tmp_k0
      integer h_tmp_gamma0,i_tmp_gamma0
      double precision K0, gamma0
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
      call cons_hbondings_gamma0_tag(namespace,tag_gamma0)


      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhc)) then
         nhc = 0
      end if

      if (nhc .eq. 0 ) then
         if(ga_nodeid().eq.0)
     >     write(luout,*) "no bondings springs found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, nhc, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_info_bondings: unable to allocate scratch space',
     &      nhc, MA_ERR)
      status = rtdb_get(rtdb,tag_n0,mt_int,nhc,int_mb(i_tmp_n0))
      if(.not.status) call errquit(
     >     'cons_info_bondings: unable to get n0',
     >      nhc, MA_ERR)
      mm = 0
      do i=1,nhc
         mm = mm + int_mb(i_tmp_n0+i-1)
      end do
      if ( .not. ma_push_get( MT_INT, 2*mm, 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_info_bondings: unable to allocate scratch space',
     &      mm, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, mm, 'h_tmp_coef',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_info_bondings: unable to allocate scratch space',
     &      mm, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhc, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_info_bondings: unable to allocate scratch space',
     &      nhc, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhc, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_info_bondings: unable to allocate scratch space',
     &      nhc, MA_ERR)

      status = rtdb_get(rtdb,tag_indx,mt_int,2*mm,int_mb(i_tmp_indx))
      status = status.and.
     >         rtdb_get(rtdb,tag_coef,mt_dbl,mm,dbl_mb(i_tmp_coef))
      status = status.and.
     >         rtdb_get(rtdb,tag_k0,mt_dbl,nhc,dbl_mb(i_tmp_k0))
      status = status.and.
     >         rtdb_get(rtdb,tag_gamma0,mt_dbl,nhc,dbl_mb(i_tmp_gamma0))
      if(.not.status) call errquit(
     >     'cons_info_bondings: unable to get prior harm bondings',
     >      nhc, MA_ERR)
     
      if(ga_nodeid().eq.0) then
         ii = 0
         do i=1,nhc
            n      = int_mb(i_tmp_n0+i-1)
            K0     = dbl_mb(i_tmp_k0+i-1)
            gamma0 = dbl_mb(i_tmp_gamma0+i-1)

            write(luout,*)
            write(luout,'(A,I5)') " Bondings Spring #",i
            write(luout,'(A,F18.9,E14.6)')
     >             "   - parameters (K0,gamma0):",
     >             K0,gamma0
            write(luout,*) "  - coefficient index1 index2"
            do j=1,n
              write(luout,'(F16.6,I7,I7)') dbl_mb(i_tmp_coef+ii+j-1),
     >                           int_mb(i_tmp_indx+2*ii+2*(j-1)),
     >                           int_mb(i_tmp_indx+2*ii+2*(j-1)+1)
            end do
            write(luout,*)
            ii = ii + n
         end do
      end if

      if (.not. ma_pop_stack(h_tmp_gamma0) ) call errquit(
     &    'cons_info_bondings: unable to free h_tmp_gamma0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k0) ) call errquit(
     &    'cons_info_bondings: unable to free h_tmp_k0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_coef) ) call errquit(
     &    'cons_info_bondings: unable to free h_tmp_coef',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_indx) ) call errquit(
     &    'cons_info_bondings: unable to free h_tmp_indx',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_n0) ) call errquit(
     &    'cons_info_bondings: unable to free h_tmp_n0',
     &    0, MA_ERR)
      end








      subroutine cons_info_dihed(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) namespace
      integer rtdb
c     
      logical status
      integer i
      integer nc
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      character*255 tag_id
      character*255 tag_n
      character*255 tag_k
      character*255 tag_r
      integer idm

      idm = 4

      call cons_hdihed_id_tag(namespace,tag_id)
      call cons_hdihed_n_tag(namespace,tag_n)
      call cons_hdihed_k_tag(namespace,tag_k)
      call cons_hdihed_r_tag(namespace,tag_r)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nc)) then
         nc = 0
      end if
       
      if (nc .eq. 0 ) then
         if(ga_nodeid().eq.0)    
     >     write(luout,*) "no dihedral springs found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, idm*nc, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nc, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nc, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nc, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nc, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nc, MA_ERR)

       if(nc .ge. 1) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,idm*nc,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nc,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nc,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nc, MA_ERR)
       end if

     
      if(ga_nodeid().eq.0) then   
      do i=1,nc
      write(luout,FMT=200) int_mb(i_tmp_id+2*(i-1)),
     >           int_mb(i_tmp_id+2*(i-1)+1),
     >           int_mb(i_tmp_id+2*(i-1)+2),
     >           int_mb(i_tmp_id+2*(i-1)+3),
     >           dbl_mb(i_tmp_k+i-1),
     >           dbl_mb(i_tmp_r+i-1)
      end do
      end if
200   FORMAT(" spring dihedral (i,j,k,r0)",T29,4(I6,2X),F6.3,2X,F6.3)
      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_id',
     &    0, MA_ERR)
      end

      subroutine cons_print_gen_info(namespace,rtdb)
      implicit none
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      character*(*) namespace
      integer rtdb

      if(ga_nodeid().eq.0) then
        call cons_print_header()
        call util_print_centered(luout,
     >     "General Information",
     >     1,.true.)
      end if

      call cons_info_fixed_atoms(namespace,rtdb)
      call cons_info_springs(namespace,rtdb)
      call cons_info_bondings(namespace,rtdb)
      call cons_info_pbondings(namespace,rtdb)
      call cons_info_dihed(namespace,rtdb)

      end

      function cons_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_fixed_atoms
      character*(*) namespace
      integer rtdb
c
      integer nfa
      character*255 tag_n

      call cons_fxd_atoms_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nfa)) nfa=0
     
      if(nfa.eq.0) then
        cons_fixed_atoms=.false.
      else
        cons_fixed_atoms=.true.
      end if

      return
      end

      function cons_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_hbonds
      character*(*) namespace
      integer rtdb
c
      integer nhb
      character*255 tag_n

      call cons_hbond_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0
     
      if(nhb.eq.0) then
        cons_hbonds=.false.
      else
        cons_hbonds=.true.
      end if

      return
      end


      logical function cons_hbondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      character*(*) namespace
      integer rtdb
c
      integer nhb
      character*255 tag_n

      call cons_hbondings_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0
      
      if(nhb.eq.0) then
        cons_hbondings=.false.
      else
        cons_hbondings=.true.
      end if

      return
      end


      logical function cons_pbondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      character*(*) namespace
      integer rtdb
c
      integer nhb
      character*255 tag_n

      call cons_pbondings_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0

      if(nhb.eq.0) then
        cons_pbondings=.false.
      else
        cons_pbondings=.true.
      end if

      return
      end





      function cons_hdihed(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_hdihed
      character*(*) namespace
      integer rtdb
c
      integer nhb
      character*255 tag_n

      call cons_hdihed_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0
     
      if(nhb.eq.0) then
        cons_hdihed=.false.
      else
        cons_hdihed=.true.
      end if

      return
      end

      function cons_status(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_status
      character*(*) namespace
      integer rtdb
c
      character*255 tag

      call cons_status_tag(namespace,tag)

      if(.not.rtdb_get(rtdb,tag,mt_log,1,cons_status)) then
        cons_status=.true.
      end if
     
      return
      end

      subroutine cons_set_status(namespace,cons_status,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_status
      character*(*) namespace
      integer rtdb
c
      character*255 tag

      call cons_status_tag(namespace,tag)

      if(.not.rtdb_put(rtdb,tag,mt_log,1,cons_status)) then
        call errquit("cons_set_status",0,RTDB_ERR)
      end if
     
      return
      end

      subroutine cons_delete_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag_id
      character*255 tag_n
c
      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)
c
      ignore = rtdb_delete(rtdb,tag_n)
      ignore = rtdb_delete(rtdb,tag_id)
c
      end



      subroutine cons_info_pbondings(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) namespace
      integer rtdb
c     
      logical status
      integer i,ii,mm,j,n
      integer nhp
      integer h_tmp_n0,i_tmp_n0
      integer h_tmp_indx,i_tmp_indx
      integer h_tmp_coef,i_tmp_coef
      integer h_tmp_k0,i_tmp_k0
      integer h_tmp_gcut0,i_tmp_gcut0
      integer h_tmp_gamma0,i_tmp_gamma0
      integer h_tmp_c0,i_tmp_c0
      double precision K0, gcut0,gamma0,c0
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


      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhp)) then
         nhp = 0
      end if

      if (nhp .eq. 0 ) then
         if(ga_nodeid().eq.0)
     >     write(luout,*) "no bondings springs found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, nhp, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      nhp, MA_ERR)
      status = rtdb_get(rtdb,tag_n0,mt_int,nhp,int_mb(i_tmp_n0))
      if(.not.status) call errquit(
     >     'cons_info_pbondings: unable to get n0',
     >      nhp, MA_ERR)
      mm = 0
      do i=1,nhp
         mm = mm + int_mb(i_tmp_n0+i-1)
      end do
      if ( .not. ma_push_get( MT_INT, 2*mm, 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      mm, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, mm, 'h_tmp_coef',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      mm, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      nhp, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_gcut0',
     &      h_tmp_gcut0, i_tmp_gcut0) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      nhp, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      nhp, MA_ERR)

      if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_c0',
     &      h_tmp_c0, i_tmp_c0) ) call errquit(
     &      'cons_info_pbondings: unable to allocate scratch space',
     &      nhp, MA_ERR)

      status = rtdb_get(rtdb,tag_indx,mt_int,2*mm,int_mb(i_tmp_indx))
      status = status.and.
     >         rtdb_get(rtdb,tag_coef,mt_dbl,mm,dbl_mb(i_tmp_coef))
      status = status.and.
     >         rtdb_get(rtdb,tag_k0,mt_dbl,nhp,dbl_mb(i_tmp_k0))
      status = status.and.
     >         rtdb_get(rtdb,tag_gcut0,mt_dbl,nhp,dbl_mb(i_tmp_gcut0))
      status = status.and.
     >         rtdb_get(rtdb,tag_gamma0,mt_dbl,nhp,dbl_mb(i_tmp_gamma0))
      status = status.and.
     >         rtdb_get(rtdb,tag_c0,mt_dbl,nhp,dbl_mb(i_tmp_c0))

      if(.not.status) call errquit(
     >     'cons_info_pbondings: unable to get prior harm bondings',
     >      nhp, MA_ERR)
     
      if(ga_nodeid().eq.0) then
         ii = 0
         do i=1,nhp
            n      = int_mb(i_tmp_n0+i-1)
            K0     = dbl_mb(i_tmp_k0+i-1)
            gcut0  = dbl_mb(i_tmp_gcut0+i-1)
            gamma0 = dbl_mb(i_tmp_gamma0+i-1)
            c0     = dbl_mb(i_tmp_c0+i-1)

            write(luout,*)
            write(luout,'(A,I5)') " Bondings Penalty #",i
            write(luout,'(A,F18.9,F18.9,E14.6,F18.9)')
     >             "   - parameters (K0,gcut0,gamma0,c0):",
     >             K0,gcut0,gamma0,c0
            write(luout,*) "  - coefficient index1 index2"
            do j=1,n
              write(luout,'(F16.6,I7,I7)') dbl_mb(i_tmp_coef+ii+j-1),
     >                           int_mb(i_tmp_indx+2*ii+2*(j-1)),
     >                           int_mb(i_tmp_indx+2*ii+2*(j-1)+1)
            end do
            write(luout,*)
            ii = ii + n
         end do
      end if

      if (.not. ma_pop_stack(h_tmp_c0) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_c0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_gamma0) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_gamma0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_gcut0) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_gcut0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k0) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_k0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_coef) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_coef',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_indx) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_indx',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_n0) ) call errquit(
     &    'cons_info_pbondings: unable to free h_tmp_n0',
     &    0, MA_ERR)
      end





*     ************************************
*     *                                  *
*     *      cons_add_pbondings          *
*     *                                  *
*     ************************************

      subroutine cons_add_pbondings(namespace,n,indx,coef,
     >                              K0,gcut0,gamma0,c0,rtdb)
      implicit none
      character*(*) namespace
      integer n,indx(*)
      double precision coef(*),K0,gcut0,gamma0,c0
      integer rtdb

#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical status
      integer i
      integer nhp0,mm
      integer nhp
      integer h_tmp_n0,i_tmp_n0
      integer h_tmp_indx,i_tmp_indx
      integer h_tmp_coef,i_tmp_coef
      integer h_tmp_k0,i_tmp_k0
      integer h_tmp_gcut0,i_tmp_gcut0
      integer h_tmp_gamma0,i_tmp_gamma0
      integer h_tmp_c0,i_tmp_c0
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


      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhp0)) then
         nhp0 = 0
      end if
      nhp = nhp0 + 1
      mm = 0

      if(nhp0.eq.0) then
         if (.not. ma_push_get( MT_INT, nhp, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if (.not. ma_push_get( MT_INT, nhp*2*n, 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp*2*n, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhp*n, 'h_tmp_coef',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp*n, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhp, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhp, 'h_tmp_gcut0',
     &      h_tmp_gcut0, i_tmp_gcut0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhp, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if (.not. ma_push_get( MT_DBL, nhp, 'h_tmp_c0',
     &      h_tmp_c0, i_tmp_c0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)

      else
         if ( .not. ma_push_get( MT_INT, nhp, 'h_tmp_n0',
     &      h_tmp_n0, i_tmp_n0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         status = rtdb_get(rtdb,tag_n0,mt_int,nhp0,int_mb(i_tmp_n0))
         mm = 0
         do i=1,nhp0
            mm = mm + int_mb(i_tmp_n0+i-1)
         end do
         if ( .not. ma_push_get( MT_INT, 2*(mm+n), 'h_tmp_indx',
     &      h_tmp_indx, i_tmp_indx) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*(mm+n), MA_ERR)
         if ( .not. ma_push_get( MT_DBL, (mm+n), 'h_tmp_indx',
     &      h_tmp_coef, i_tmp_coef) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      (mm+n), MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_k0',
     &      h_tmp_k0, i_tmp_k0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_gcut0',
     &      h_tmp_gcut0, i_tmp_gcut0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_gamma0',
     &      h_tmp_gamma0, i_tmp_gamma0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)
         if ( .not. ma_push_get( MT_DBL, nhp, 'h_tmp_c0',
     &      h_tmp_c0, i_tmp_c0) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhp, MA_ERR)

         status = status . and. rtdb_get(rtdb,tag_indx,
     >     mt_int,2*mm,int_mb(i_tmp_indx))
         status = status . and. rtdb_get(rtdb,tag_coef,
     >     mt_dbl,mm,dbl_mb(i_tmp_coef))
         status = status . and. rtdb_get(rtdb,tag_k0,
     >     mt_dbl,nhp0,dbl_mb(i_tmp_k0))
         status = status . and. rtdb_get(rtdb,tag_gcut0,
     >     mt_dbl,nhp0,dbl_mb(i_tmp_gcut0))
         status = status . and. rtdb_get(rtdb,tag_gamma0,
     >     mt_dbl,nhp0,dbl_mb(i_tmp_gamma0))
         status = status . and. rtdb_get(rtdb,tag_c0,
     >     mt_dbl,nhp0,dbl_mb(i_tmp_c0))

         if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bondings',
     >      nhp, MA_ERR)

         status =              rtdb_delete(rtdb,tag_n)
         status = status .and. rtdb_delete(rtdb,tag_n0)
         status = status .and. rtdb_delete(rtdb,tag_indx)
         status = status .and. rtdb_delete(rtdb,tag_coef)
         status = status .and. rtdb_delete(rtdb,tag_k0)
         status = status .and. rtdb_delete(rtdb,tag_gcut0)
         status = status .and. rtdb_delete(rtdb,tag_gamma0)
         status = status .and. rtdb_delete(rtdb,tag_c0)
         if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bondingss',
     >      0, MA_ERR)
      end if

      int_mb(i_tmp_n0+nhp-1) = n
      do i=1,2*n
         int_mb(i_tmp_indx+2*mm+i-1) = indx(i)
      end do
      do i=1,n
         dbl_mb(i_tmp_coef+mm+i-1) = coef(i)
      end do
      dbl_mb(i_tmp_k0+nhp-1) = K0
      dbl_mb(i_tmp_gcut0+nhp-1) = gcut0
      dbl_mb(i_tmp_gamma0+nhp-1) = gamma0

       status = rtdb_put(rtdb,tag_n,mt_int,1,nhp)
       status = status.and. rtdb_put(rtdb,tag_n0,
     >     mt_int,nhp,int_mb(i_tmp_n0))
       status = status.and. rtdb_put(rtdb,tag_indx,
     >     mt_int,2*(mm+n),int_mb(i_tmp_indx))
       status = status.and. rtdb_put(rtdb,tag_coef,
     >     mt_dbl,(mm+n),dbl_mb(i_tmp_coef))
       status = status.and. rtdb_put(rtdb,tag_k0,
     >     mt_dbl,nhp,dbl_mb(i_tmp_k0))
       status = status.and. rtdb_put(rtdb,tag_gcut0,
     >     mt_dbl,nhp,dbl_mb(i_tmp_gcut0))
       status = status.and. rtdb_put(rtdb,tag_gamma0,
     >     mt_dbl,nhp,dbl_mb(i_tmp_gamma0))
       status = status.and. rtdb_put(rtdb,tag_c0,
     >     mt_dbl,nhp,dbl_mb(i_tmp_c0))

       if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bondings',
     >      nhp, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_c0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_c0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_gamma0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_gamma0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_gcut0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_gcut0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_k0',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_coef) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_coef',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_indx) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_indx',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_n0) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_n0',
     &    0, MA_ERR)

      return
      end
