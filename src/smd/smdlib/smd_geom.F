      subroutine smd_geom_init_system0()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(smd_string_size) sp_geom,sp_atom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer h_c,i_c,h_ct,i_ct
      integer i_nb
      integer nt
      integer ng
      integer i,j
      integer i_pr,i_mr
      integer i_im,i_nm
      integer i_c1,i_ms1,i_t1,i_qf1
      integer i_sht
      integer i_qf,i_m,i_ta
      integer rtdb

      pname = "smd_geom_init_system0"

c      write(*,*) "in "//pname
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
      call util_message(pname//":0")
      call smd_namespace_create(sp_geom)
c      write(*,*) "out of "//pname

      return
      end

      subroutine smd_geom_init_system(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer im
c     
      character*(smd_string_size) sp_geom,sp_atom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer h_c,i_c,h_ct,i_ct
      integer i_nb
      integer nt
      integer ng
      integer i,j
      integer i_pr,i_mr
      integer i_im,i_nm
      integer i_c1,i_ms1,i_t1,i_qf1
      integer i_sht
      integer i_qf,i_m,i_ta
      integer rtdb

      pname = "smd_geom_init_system"

      call util_message("in "//pname)
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      call smd_namespace_create(sp_geom)
c      call util_message(pname//":1")
c
      tag = "atom:res_pointer"
      call smd_get_ind_dim(tag,i_pr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:res_size"
      call smd_get_ind_dim(tag,i_mr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c
c     determine neighbors around im
c     -----------------------------
      if(.not.ma_push_get(mt_int,nr,'tmp i',h_itmp,i_itmp))
     + call errquit(pname//'Failed to allocate memory for h_itmp',0,0)
c
      nb = nr
      call smd_group_message("1: "//pname)
      call smd_geom_find_neighbors(im,nb,int_mb(i_itmp))
      tag = "neighbor_index"
      call smd_data_put_int(sp_geom,tag,nb,int_mb(i_itmp))
      call smd_data_get_index(sp_geom,tag,i_nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
      call smd_group_message("2: "//pname)
c
c     determine position and size (atom based) of the cluster residues
c
      tag = "position"
      call smd_data_create_get(sp_geom,tag,nb,MT_INT,i_im)
      tag = "size"
      call smd_data_create_get(sp_geom,tag,nb,MT_INT,i_nm)
      do i=1,nb
        j=int_mb(i_nb+i-1)
        int_mb(i_im+i-1)=int_mb(i_pr+j-1)
        int_mb(i_nm+i-1)=int_mb(i_mr+j-1)
      end do
      call smd_group_message("3: "//pname)
c
c     total number of atoms in the geometry
c     ------------------------------------
      ng = 0
      do i=1,nb
        j=int_mb(i_nb+i-1)
        ng=ng+int_mb(i_mr+j-1)
      end do
      tag = "ntot"
      call smd_data_put_int(sp_geom,tag,1,ng)
c      
      if(.not.ma_pop_stack(h_itmp))
     & call errquit(pname//'Failed to deallocate stack h_itmp',0,
     &       MA_ERR)
      call smd_group_message("4: "//pname)
c     
c     get copy of coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
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
c    
      call smd_geom_create()
      call smd_group_message("6: "//pname)

      call util_message("out "//pname)

      return
      end
c
      subroutine smd_geom_create()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "smd_const_data.fh"
c
      logical geom_cart_set_gen1
      external geom_cart_set_gen1
c      
      integer rtdb,m,n
c
      integer i,j
      integer i0,i1,n0,k
      character*30 pname
      character*(smd_string_size) tag
      character*(smd_string_size) sp_geom
      integer igeom
      integer nact,i_act,h_act
      double precision shift(3)
      integer i_t1,h_t1
      integer i_q1,h_q1
      integer i_c1,h_c1
      integer i_ms1,h_ms1
      integer n1,ns
      integer i_c,i_qf,i_m,i_ta
      integer i_im,i_nm
      integer ioff
      logical result
c
      pname = "smd_geom_create_0"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,m,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,m,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c    
c
      ns = smd_string_size
c     
c     get copy of coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,n,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
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
c c
c     size
c     ----
      n1=0
      do i=1,m
        n1 = n1 + int_mb(i_nm+i-1)
      end do
c      
      if(.not.ma_push_get(mt_dbl,3*n1,'tmp c',h_c1,i_c1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_dbl,n1,'tmp q',h_q1,i_q1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_dbl,n1,'tmp m',h_ms1,i_ms1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_byte,n1*smd_string_size,'tmp t',h_t1,i_t1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c
      call smd_geom_cog(n,int_mb(i_nm),int_mb(i_im),
     +                  dbl_mb(i_c),shift)
      tag = "shift vector"
      call smd_data_put_dbl(sp_geom,tag,3,shift)
      j=0
      do i=1,m
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       do i1=1,n0
         j = j+1
         ioff = i0+i1-1
         do k=1,3
           dbl_mb(i_c1+(j-1)*3+k-1)=
     +         dbl_mb(i_c+ioff*3+k-1)-shift(k)
         end do
         dbl_mb(i_ms1+j-1)=dbl_mb(i_m+ioff)
         dbl_mb(i_q1+j-1)=dbl_mb(i_qf+ioff)
         do k=1,ns
           byte_mb(i_t1+(j-1)*ns+k-1)=byte_mb(i_ta+ioff*ns+k-1)
         end do
       end do
      end do
      call smd_lat_rebox(n1,dbl_mb(i_c1)) 

       i0=1
       if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set_gen1(igeom,i0,n1,n1,smd_string_size,
     + byte_mb(i_t1),dbl_mb(i_c1),dbl_mb(i_q1)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,n1,dbl_mb(i_ms1)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
c
      if(ga_nodeid().eq.0) then
      if (.not. geom_print_xyz(igeom,105))
     + call errquit('qmmmsi geom_print_xyz',0, GEOM_ERR)
      end if
c
      call ga_sync()
c
      call smd_rtdb_get_handle(rtdb)
      if(.not.geom_rtdb_store(rtdb,igeom,"geometry"))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
      nact = int_mb(i_nm)
      if(.not.ma_push_get(mt_int,nact,'tmp',h_act,i_act))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nact, MA_ERR)
      do i=1,nact
       int_mb(i_act+i-1) = i
      end do
      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >       mt_int,nact,int_mb(i_act)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if(.not.ma_pop_stack(h_act))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_t1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_ms1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)

      return
c
  
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
      subroutine smd_geom_cog(nt,nm,i0,c,cog)
c     finds center of geometry for a give segment
c     nt      [in]  total number of atoms
c     nm      [in]  number of atoms in the segment
c     i0      [in]  starting index for the segment
c     c[3,nt] [in]  coordinates of the overall system
c     cog[3]  [out] center of mass coordinates
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
      double precision c(3,nt)
      double precision cog(3)
c
      integer i,j
      double precision cshift(3)
      integer i_ftmp,h_ftmp
      character*30 pname
c
      pname = "smd_geom_cog"
c      
      if(.not.ma_push_get(mt_dbl,3*nm,'tmp i',h_ftmp,i_ftmp))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c
c     put 1st atom in the middle
c     ---------------------------
      do j=1,3
         cshift(j)=c(j,i0)
      end do

      do i=1,nm
      do j=1,3
         dbl_mb(i_ftmp+(i-1)*3+j-1)=
     +         c(j,i0+i-1)-cshift(j)
      end do
      end do
 
c     fold in other coordinates potentially
c     cutoff by the box boundary
      call smd_lat_rebox(nm,dbl_mb(i_ftmp)) 
c     calculate cog for shifted segment
c     --------------------------------- 
      cog(1) = 0.0d0
      cog(2) = 0.0d0
      cog(3) = 0.0d0
      
      do i=1,nm
      do j=1,3
         cog(j)=cog(j)+dbl_mb(i_ftmp+(i-1)*3+j-1)
      end do
      end do

      cog(1) = cog(1)/nm
      cog(2) = cog(2)/nm
      cog(3) = cog(3)/nm
c
c     correct for the shift
c     ---------------------
      do j=1,3
        cog(j) = cog(j)+cshift(j)
      end do
c
      if(.not.ma_pop_stack(h_ftmp))
     & call errquit(pname//'Failed to deallocate stack h_ftmp',0,
     &       MA_ERR)
c
      return
c
  
      end
c
      subroutine smd_geom_get_charges()
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
      character*32 sp_geom
      character*32 tag,pname
      logical result
      integer i_ir,i_ta,i_q,i_c,i_m,i_qf
      integer i_ct,h_ct
      integer i,j,i0,i1,n0
      integer nt,nm,nqm
      integer rtdb
      integer igeom
      integer i_q1,h_q1
      integer ng
      integer nb
      integer i_nm,i_im
c
      logical geom_cart_set1
      external geom_cart_set1
      logical geom_rtdb_ncent
      external geom_rtdb_ncent
c
      pname = "smd_geom_get_charges"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      call smd_rtdb_get_handle(rtdb)
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c      
c      tag = "charge"
c      call smd_data_get_ind_dim(sp_geom,tag,i_q1,ng,result)
c      if(.not. result)  call errquit(
c     >       pname//'error getting index for'//tag,0,0 )
c
      tag="ntot"
      call smd_data_get_int(sp_geom,tag,1,ng,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c    
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

      if(nqm.ne.ng)
     + call errquit(pname//'mismatch in sizes'//tag,0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,igeom,tag(1:inp_strlen(tag))))
     + call errquit(pname//'Failed to load geometry'//tag,0, GEOM_ERR)

      if(.not.ma_push_get(mt_dbl,ng,'tmp q',h_q1,i_q1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)

      if(.not.geom_cart_get_charges(igeom,ng,dbl_mb(i_q1)))
     + call errquit(pname//'Failed to get charges',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)
c
      j=0
      do i=1,nb
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       do i1=1,n0
         j = j+1
         dbl_mb(i_q+i0+i1-1) = dbl_mb(i_q1+j-1)
       end do
      end do

      if(.not.ma_pop_stack(h_q1))
     & call errquit(pname//'Failed to deallocate stack ',0,
     &       MA_ERR)
c      
      return
      end
c      
      subroutine smd_geom_get_gradients()
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
      character*32 tag,pname,sp_geom
      logical result
      integer i_f
      integer i_g1,h_g1
      integer i,j,i1,i0
      integer ng
      integer rtdb
      integer i_im,i_nm
      integer nb,k,n0
c
      pname = "smd_geom_get_gradients"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag="ntot"
      call smd_data_get_int(sp_geom,tag,1,ng,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
c      tag = "gradients"
c      call smd_data_get_ind_dim(sp_geom,tag,i_g1,nt,result)
c      if(.not. result)  call errquit(
c     >       pname//'error getting index for'//tag,0,0 )
c      nt = nt/3
c
      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c    
c
      call smd_rtdb_get_handle(rtdb)
c
      tag = "force"
      call smd_get_ind(tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c    
      if(.not.ma_push_get(mt_dbl,3*ng,'tmp g',h_g1,i_g1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c      
      if (.not.rtdb_get(rtdb, 'task:gradient', mt_dbl,
     $        3*ng, dbl_mb(i_g1))) 
     $           call errquit(pname,0,RTDB_ERR)

      j=0
      do i=1,nb
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       do i1=1,n0
         j = j+1
         do k=1,3
           dbl_mb(i_f+3*(i0+i1-1)+k-1) = dbl_mb(i_g1+3*(j-1)+k-1)
         end do
       end do
      end do


      if(.not.ma_pop_stack(h_g1))
     & call errquit(pname//'Failed to deallocate stack ',0,
     &       MA_ERR)

      return
      end
c
      subroutine smd_geom_get_coords()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "smd_const_data.fh"
c
      integer im
c     
      character*32 tag,pname,sp_geom
      logical result
      integer i_c
      integer i_c1,h_c1
      integer i,j,i1,i0
      integer nt
      integer rtdb
      integer i_im,i_nm
      integer nb,k,n0
      integer igeom
      double precision sht(3)
c
      pname = "smd_geom_get_gradients"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag="ntot"
      call smd_data_get_int(sp_geom,tag,1,nt,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c    
c
      call smd_rtdb_get_handle(rtdb)
c
      tag = "coords"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c     
      if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
c      
      if(.not.geom_rtdb_load(rtdb,igeom,"geometry"))
     + call errquit(pname//'Failed to load geometry'//tag,0, GEOM_ERR)
c    
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp g',h_c1,i_c1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c    
      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c1)))
     $           call errquit(pname,0,RTDB_ERR)
      call dscal(3*nt,cau2ang,dbl_mb(i_c1),1)
c
      tag = "shift vector"
      call smd_data_get_dbl(sp_geom,tag,3,sht,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
       
      do i=1,nt
        do k=1,3
          dbl_mb(i_c1+3*(i-1)+k-1)=dbl_mb(i_c1+3*(i-1)+k-1)+
     >                             sht(k)
        end do
      end do
      call smd_lat_rebox(nt,dbl_mb(i_c1)) 
c
      j=0
      do i=1,nb
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       do i1=1,n0
         j = j+1
         do k=1,3
           dbl_mb(i_c+3*(i0+i1-1)+k-1) = dbl_mb(i_c1+3*(j-1)+k-1)
         end do
       end do
      end do

      if(.not.ma_pop_stack(h_c1))
     & call errquit(pname//'Failed to deallocate stack ',0,
     &       MA_ERR)

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

      subroutine smd_geom_find_neighbors(icr,nb,inb)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      integer icr
      integer nb
      integer inb(nb)
c
      character*30 pname
      character*32 tag
      integer i_f,i_c,i_rlat
      integer i_lat
      integer i_l,i_p
      integer i_ir
      integer i_inb,h_nb
      integer nl,na
      double precision rcutsq
      logical result
      integer i,j
      integer i_tmp
      integer nr,nrmax

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
c
c     get cutoff data
c     ---------------
      call smd_cutoff_get_rcut(rcutsq)
      rcutsq = rcutsq*rcutsq
c      
      if(.not.ma_push_get(mt_int,nr,'tmp inb',h_nb,i_inb))
     + call errquit(pname//'Failed to allocate memory for tmp inb',
     + na, MA_ERR)
c      
      call    smd_geom_find_neighbors_gen(icr,
     >                                  na,
     >                                  nr,
     >                                  nl,
     >                                  rcutsq,
     >                                  int_mb(i_ir),
     >                                  dbl_mb(i_c),
     >                                  int_mb(i_p),
     >                                  int_mb(i_l),
     >                                  int_mb(i_inb))

      j = 1
      inb(j) = icr
      do i=1,nr
       if(int_mb(i_inb+i-1).ne.0) then
         j = j+1
         inb(j) = i
       end if
      end do
      nb = j

c      j = 0
c      do i=1,nr
c       if(int_mb(i_inb+i-1).ne.0) then
c         j = j+1
c         int_mb(i_inb+j-1) = i
c       end if
c      end do
c      nb=j+1
c      int_mb(i_inb+nb-1) = int_mb(i_inb)
c      int_mb(i_inb) = icr 
c
c      write(*,*) "test :"//pname
c      do i=1,nb
c         write(*,*) int_mb(i_inb+i-1),inb(i)
c      end do



      if(.not.ma_pop_stack(h_nb))
     & call errquit(pname//'Failed to deallocate stack h_nb',na,
     &       MA_ERR)

      return

      END
c
      subroutine smd_geom_find_neighbors_gen(icr,
     >                           na,
     >                           nr,
     >                           nl,
     >                           rcutsq,
     >                           ir,
     >                           ccc,
     >                           point,
     >                           list,
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
c       inb(nr)    [out] indices of neighbors

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
      integer inb(nr)
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

      call util_message("1: "//pname)
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

c      write(79,*) "neighbors of residue",icr
c      do i=1,nr
c      if(inb(i).eq.1) then
c        write(79,*) i
c      end if
c      end do
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

      ii = 0
      do i=1,na-1

       ir1 = ir(i)
       jbeg=point(i)
       jend=point(i+1)-1

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
         endif
       end if


      end do
      end do

      nb = ii
c      do i=1,nb
c        write(*,*) "neighbor is",i,inb(i)
c      end do
      return

      end

