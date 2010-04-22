      subroutine smd_geom_monomer_create(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "bq.fh"
#include "smd_const_data.fh"
c
      integer im
c     
      character*32 tag,pname
      logical result
      integer i_ir,i_ta,i_q,i_c,i_m,i_qf
      integer i_ct,h_ct
      integer i_cd,h_cd
      integer i,j
      integer nt,nm,i0
      integer rtdb
      integer igeom
      integer bq_handle
      integer nbq,i_qbq,i_cbq
      double precision scale
      double precision cm(3)
      logical nobq
c
      logical geom_cart_set_gen1
      external geom_cart_set_gen1
      pname = "smd_geom_monomer_create"
c
      call smd_rtdb_get_handle(rtdb)
c     
      tag = "coords"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "charge"
      call smd_get_ind(tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:name"
      call smd_get_ind(tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:formal_charge"
      call smd_get_ind(tag,i_qf,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
      tag = "mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c    
      call smd_geom_monomer_index_size(im,nt,int_mb(i_ir),nm,i0)
c
      call smd_geom_find_neighbors(i0)
c
      write(*,*) "monomer size",i0,nm
c
      if(.not.ma_push_get(mt_dbl,3*nm,'tmp',h_cd,i_cd))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)
 
c
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)
 
c
c     shift box so that monomer is in the center
c     -----------------------------------------
      do i=1,3*nt
        dbl_mb(i_ct+i-1) = dbl_mb(i_c+i-1)
      end do

c      call smd_coords_print_pdb("old.pdb")

      call smd_geom_recenter(nt,nm,i0,dbl_mb(i_cd),
     >                       dbl_mb(i_ct))

c      call smd_coords_print_pdb1("new.pdb",dbl_mb(i_ct))

      if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set_gen1(igeom,i0,nm,nt,smd_string_size,
     + byte_mb(i_ta),
     + dbl_mb(i_ct),dbl_mb(i_qf)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nm,dbl_mb(i_m+i0-1)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
c
      if (.not. geom_print_xyz(igeom,97))
     + call errquit('qmmmsi geom_print_xyz',0, GEOM_ERR)
c
      if(.not.geom_rtdb_store(rtdb,igeom,"geometry"))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c      call smd_coords_print(55)
c
c     now do the charges
c
      call smd_rtdb_get_log("smd:nobq",1,nobq,result)
      if(.not.result) nobq=.false.
      write(*,*) "nobq is",nobq
      if(.not.nobq) then
        nbq = nt - nm
        if(.not.bq_create("monomer charges",bq_handle))
     +   call errquit(pname//'Failed bq_create',0,CALC_ERR)

        if(.not.bq_alloc(bq_handle,nbq))
     +   call errquit(pname//'Failed bq_create',0,CALC_ERR)


        if(.not.bq_index_coord(bq_handle,i_cbq))
     +   call errquit(pname//'Failed bq_index',0,CALC_ERR)

        if(.not.bq_index_charge(bq_handle,i_qbq))
     +   call errquit(pname//'Failed bq_index',0,CALC_ERR)

        scale = 1.88972598858d0

        do i=1,i0-1
           do j=1,3
              dbl_mb(i_cbq+3*(i-1)+j-1)=dbl_mb(i_ct+3*(i-1)+j-1)*scale
           end do
           dbl_mb(i_qbq+i-1)=dbl_mb(i_q+i-1)
        end do

        do i=i0,nbq
           do j=1,3
              dbl_mb(i_cbq+3*(i-1)+j-1)=dbl_mb(i_ct+3*(i-1)+j-1+3*nm)
     >                                  *scale
           end do
           dbl_mb(i_qbq+i-1)=dbl_mb(i_q+nm+i-1)
        end do


        call bq_print_xyz(bq_handle,98)

        if(.not.bq_activate(bq_handle))
     +   call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c

c        call bq_print_info(bq_handle)
      end if

      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_cd))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
c


      return
      end
c
      subroutine smd_geom_monomer_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "bq.fh"
c
      integer bq_handle

      character*32 pname

      pname = "smd_geom_end"
      call smd_message_post("in "//pname)
      if(bq_get_active(bq_handle)) then
      call util_message("found bq, destroying it")

      if(.not.bq_destroy(bq_handle))
     + call errquit(pname//'Failed bq_destroy',0,CALC_ERR)

      end if

      call smd_rtdb_pop()
      call smd_message_post("out "//pname)
      return
      end
c
      subroutine smd_geom_recenter(nt,nm,i0,cm,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer nt
      integer nm
      integer i0
      double precision cm(3,nm)
      double precision c(3,nt)
c
      integer i,j
      double precision cshift(3)
      double precision com(3)
c
c     assemble a complete residue
c     ---------------------------
      do j=1,3
         cshift(j)=c(j,i0)
      end do

      do i=1,nm
      do j=1,3
         cm(j,i)=c(j,i0+i-1)-cshift(j)
      end do
      end do
 
      call smd_lat_rebox(nm,cm) 

      do i=1,nm
      do j=1,3
         cm(j,i)=cm(j,i)+cshift(j)
      end do
      end do
c
c     shift the entire box so that
c     monomer geom. center is at the 
c     center
c     ------------------------------ 
      com(1) = 0.0d0
      com(2) = 0.0d0
      com(3) = 0.0d0
      
      do i=1,nm
      do j=1,3
         com(j)=com(j)+cm(j,i)
      end do
      end do

      com(1) = com(1)/nm
      com(2) = com(2)/nm
      com(3) = com(3)/nm

      do i=1,nt
      do j=1,3
         c(j,i)=c(j,i)-com(j)
      end do
      end do
 
      call smd_lat_rebox(nt,c) 

c
      return
c
  
      end
c
      subroutine smd_geom_push_charges(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c
      integer im
c     
      character*32 sp_coords,sp_atom,sp_charge,sp_mass
      character*32 tag,pname
      logical result
      integer i_ir,i_ta,i_q,i_c,i_m,i_qf
      integer i_ct,h_ct
      integer i,j
      integer nt,nm,i0,nqm
      integer rtdb
      integer igeom
c
      logical geom_cart_set1
      external geom_cart_set1
      logical geom_rtdb_ncent
      external geom_rtdb_ncent
c
      pname = "smd_geom_push_charges"
c
      call smd_rtdb_get_handle(rtdb)
c     
      tag = "coords"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "charge"
      call smd_get_ind(tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:name"
      call smd_get_ind(tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:formal_charge"
      call smd_get_ind(tag,i_qf,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
      tag = "mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c    
      call smd_geom_monomer_index_size(im,nt,int_mb(i_ir),nm,i0)
c
      if(.not.ma_push_get(mt_dbl,3*nm,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)

      if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

c
c     pick the best available charge geometry
c     ---------------------------------------
      tag = "cresp:geometry"
      result = geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      if(.not.result) then
        tag = "resp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.result) then
        tag = "cesp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.result) then
        tag = "esp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if

      if(.not.result) then
         call errquit(pname//'cannot find charges')
      end if

      if(nqm.ne.nm)
     + call errquit(pname//'mismatch in sizes'//tag,0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,igeom,tag(1:inp_strlen(tag))))
     + call errquit(pname//'Failed to load geometry'//tag,0, GEOM_ERR)

      if(.not.geom_cart_get_charges(igeom,nm,dbl_mb(i_q+i0-1)))
     + call errquit(pname//'Failed to get charges',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

c
      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
c
      return
      end
cc      
      subroutine smd_geom_push_gradients(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer im
c     
      character*32 tag,pname
      logical result
      integer i_ir,i_f
      integer i_g,h_g
      integer i,j
      integer nt,nm,i0,nqm
      integer rtdb
      integer ma_type,n
c
      pname = "smd_geom_push_gradients"
c
      call smd_rtdb_get_handle(rtdb)
c     
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "force"
      call smd_get_ind(tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c    
      call smd_geom_monomer_index_size(im,nt,int_mb(i_ir),nm,i0)
c
      if (rtdb_ma_get(rtdb, 'task:gradient', ma_type,
     $        n, h_g)) then

        if (.not. ma_get_index(h_g, i_g))
     $           call errquit(pname,h_g,
     &       MA_ERR)
       end if

      do i=1,n
        write(*,*) "gradient ref",dbl_mb(i_g+i-1)
      end do
      write(*,*) "i0,nm=",i0,nm
      do i=1,3*nm
        dbl_mb(i_f+3*(i0-1)+i-1) = dbl_mb(i_g+i-1)
        write(*,*) "gradient",dbl_mb(i_f+i0+i-1)
      end do
      call smd_force_print(89)
      
      return
      end
cc
      subroutine smd_geom_monomer_index_size(im,nt,ir,nm,i0)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer im
      integer rtdb
      integer nt
      integer ir(nt)
      integer nm
      integer i0
c
      integer i,j
c
c     count number of atoms in residue im
c     contigious storage is presumed
c     -----------------------------------
      nm = 0
      do i=1,nt
       j = i
       if(im.eq.ir(i)) then
100     continue
        nm = nm + 1
        j = j + 1
        if(im.eq.ir(j)) then
         goto 100
        else
         goto 101 
        end if
       end if
      end do
101   continue
      i0=i
      return
c
  
      end

      subroutine smd_geom_find_neighbors(icr)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      integer icr
c
      character*30 pname
      character*32 tag
      integer i_f,i_c,i_rlat
      integer i_lat
      integer i_l,i_p
      integer i_ir
      integer i_inb,h_nb
      integer nl,na,nb
      double precision rcutsq
      logical result
      integer i
      integer i_tmp
      integer nr

      pname = "smd_geom_find_neighbors"
c
c     get verlet list data
c     --------------------
      tag = "vlist:pointer"
      call smd_get_ind(tag,i_p,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "vlist:list"
      call smd_get_ind_dim(tag,i_l,nl,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "vlist:distances"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c      
      tag = "nres"
      call smd_get_ind(tag,i_tmp,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      nr = int_mb(i_tmp)
      write(*,*) "number of residues is", nr
c
c     get cutoff data
c     ---------------
      call smd_cutoff_get_rcut(rcutsq)
      rcutsq = rcutsq*rcutsq
c      
      if(.not.ma_push_get(mt_int,nr,'tmp inb',h_nb,i_inb))
     + call errquit(pname//'Failed to allocate memory for tmp inb',
     + na, MA_ERR)


      call    smd_geom_find_neighbors_gen(icr,
     >                                  na,
     >                                  nr,
     >                                  nl,
     >                                  rcutsq,
     >                                  int_mb(i_ir),
     >                                  dbl_mb(i_c),
     >                                  int_mb(i_p),
     >                                  int_mb(i_l),
     >                                  nb,
     >                                  int_mb(i_inb))

      if(.not.ma_pop_stack(h_nb))
     & call errquit(pname//'Failed to deallocate stack h_nb',na,
     &       MA_ERR)

      return

      END
      subroutine smd_geom_find_neighbors_gen(icr,
     >                           na,
     >                           nr,
     >                           nl,
     >                           rcutsq,
     >                           ir,
     >                           ccc,
     >                           point,
     >                           list,
     >                           nb,
     >                           inb)

      implicit none

#include "smd_const_data.fh"


c       icr        [in] residue number to find neighbors
c       na         [in] number of atoms
c       nr         [in] number of residues
c       nl         [in] size of verlet list
c       rcutsq     [in] cutoff radius 
c       ir(na)     [in] residue index
c       point(na)  [in] verlet pointer
c       list(nl)   [in] verlet list
c       ccc(3,nl)  [in] verlet pair distance vector
c       nb         [out] number of neighbors
c       inb(na)    [out] indices of neighbors

      integer icr
      integer na
      integer nr
      integer nl
      double precision rcutsq
      integer ir(na)
      integer point(na)
      integer list(nl)
      integer nb
      double precision ccc(3,nl)
      integer inb(na)
c
      integer i,j,k,jnab
      integer jbeg,jend
      integer ii
      integer ir1,ir2,nt

      double precision  rsq

      double precision x,y,z
      character*30 pname

      pname = "smd_geom_find_neighbors_gen"

      call util_message("in "//pname)

      do i=1,nr
        inb(i) = 0
      end do

      nb = 0
      do i=1,na-1

       ir1 = ir(i)
       jbeg=point(i)
       jend=point(i+1)-1

      do jnab=jbeg,jend

       j=list(jnab)
       ir2 = ir(j)
       if(ir1.eq.icr) then
         inb(ir2) = 1
         nb = nb + 1
       else if(ir2.eq.icr) then
         inb(ir1) = 1
         nb = nb +1
       end if
      end do
      end do

      write(79,*) "neighbors of residue",icr
      do i=1,nr
      if(inb(i).eq.1) then
        write(79,*) i
      end if
      end do
      return

      end
c      
      subroutine smd_geom_find_neighbors_gen1(icr,
     >                           na,
     >                           nl,
     >                           rcutsq,
     >                           ir,
     >                           ccc,
     >                           point,
     >                           list,
     >                           nb,
     >                           inb)

      implicit none

#include "smd_const_data.fh"


c       icr        [in] residue number to find neighbors
c       na         [in] number of atoms
c       nl         [in] size of verlet list
c       rcutsq     [in] cutoff radius 
c       ir(na)     [in] residue index
c       point(na)  [in] verlet pointer
c       list(nl)   [in] verlet list
c       ccc(3,nl)  [in] verlet pair distance vector
c       nb         [out] number of neighbors
c       inb(na)    [out] indices of neighbors

      integer icr
      integer na
      integer nl
      double precision rcutsq
      integer ir(na)
      integer point(na)
      integer list(nl)
      integer nb
      double precision ccc(3,nl)
      integer inb(na)
c
      integer i,j,k,jnab
      integer jbeg,jend
      integer ii
      integer ir1,ir2,nt

      double precision  rsq

      double precision x,y,z
      character*30 pname

      pname = "smd_geom_find_neighbors_gen"

      call util_message("in "//pname)

      do i=1,na
        write(*,*) "residue index",i,ir(i)
      end do
      ii = 0
      do i=1,na-1

       ir1 = ir(i)
       jbeg=point(i)
       jend=point(i+1)-1
       write(*,*) "jbeg,jend,i",jbeg,jend,i

      do jnab=jbeg,jend

       j=list(jnab)
       ir2 = ir(j)
       if(ir1.eq.icr.or.ir2.eq.icr) then

         x=ccc(1,jnab)
         y=ccc(2,jnab)
         z=ccc(3,jnab)

         rsq=x*x+y*y+z*z

         if(rsq.lt.rcutsq) then

          ii = ii + 1
          inb(ii) = ir(j)
          write(*,*) "setting neighbor index",ii,j,inb(ii),ir(j)
         endif
       end if


      end do
      end do

      nb = ii
      do i=1,nb
        write(*,*) "neighbor is",i,inb(i)
      end do
      return

      end
