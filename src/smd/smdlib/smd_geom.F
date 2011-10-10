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

      subroutine smd_geom_repair()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "stdio.fh"
      integer im
      logical oinclude
      character*(smd_string_size) sp_geom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer i_nb
      integer nt
      integer ng
      integer i,j,ir
      integer i_pr,i_mr
      integer i_im,i_nm
      integer i_c
      character*64 msg

      pname = "smd_geom_init_system"
c
c     get locations of residues in the atomic arrays
c     ----------------------------------------------
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
c     get copy of coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      do ir=1,nr
        call smd_geom_rebuild(nt,int_mb(i_mr+ir-1),
     >                   int_mb(i_pr+ir-1),dbl_mb(i_c)) 
      end do
      call smd_coords_print_pdb("repaired.pdb")
      return
      end

      subroutine smd_geom_prepare()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
#include "smd_system.fh"
#include "stdio.fh"
      integer im
      logical oinclude
      character*(smd_string_size) sp_geom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer i_nb
      integer nt
      integer ng
      integer i,j
      integer i_pr,i_mr,i_xr
      integer i_im,i_nm
      double precision shift(3)
      integer i0,n0
      character*64 msg
      logical operiodic
      integer i_c

      pname = "smd_geom_init_system"
c
      operiodic = smd_system_periodic()
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c     get locations of residues in the atomic arrays
c     ----------------------------------------------
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
      tag = "atom:res_type"
      call smd_get_ind_dim(tag,i_xr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     determine cluster residues
c     --------------------------
      tag = "icenter"
      call smd_data_get_int(sp_geom,tag,1,im,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting central residue'//tag,0, RTDB_ERR)
c
c     check if attempting to compute nonexisting residue
      if(im.gt.nr)   
     >  call errquit(
     >       pname//'residue index out of bounds',0, RTDB_ERR)
c
c
      tag = "oinclude"
      call smd_data_get_log(sp_geom,tag,1,oinclude,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting include flag'//tag,0, RTDB_ERR)

      if(.not.ma_push_get(mt_int,nr,'tmp i',h_itmp,i_itmp))
     + call errquit(pname//'Failed to allocate memory for h_itmp',0,0)
c
      nb = nr
      call smd_geom_find_neighbors(oinclude,im,nb,int_mb(i_itmp))
c      if(ga_nodeid().eq.0) then
c        write(luout,*) "cluster information",
c     +    (int_mb(i_itmp+i-1),i=1,nb)
c        call util_flush(luout)
c      end if
      tag = "neighbor_index"
      call smd_data_put_int(sp_geom,tag,nb,int_mb(i_itmp))
      call smd_data_get_index(sp_geom,tag,i_nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
c     set residue type to QM for cluter segments
c     ------------------------------------------
      do i=1,nb
        j=int_mb(i_nb+i-1)
        int_mb(i_xr+j-1) = i_tqm
      end do
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
c     determine translation vector to put 1st residue in the center
c     -------------------------------------------------------------
c     
c     get copy of coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      if(operiodic) then
        i0=int_mb(i_pr+im-1)
        n0=int_mb(i_mr+im-1)
        call smd_geom_shift(n0,
     +                  dbl_mb(i_c+3*(i0-1)),shift)
      else
        do i=1,3
          shift(i) = 0.0d0
        end do
      end if
      tag = "shift vector"
      call smd_data_put_dbl(sp_geom,tag,3,shift)
c
c      
      if(.not.ma_pop_stack(h_itmp))
     & call errquit(pname//'Failed to deallocate stack h_itmp',0,
     &       MA_ERR)
c    
c     create geometry
c     --------------
      call smd_geom_create()

      call util_message("out "//pname)

      return
      end

      subroutine smd_geom_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
#include "stdio.fh"
      integer im
      logical oinclude
      character*(smd_string_size) sp_geom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer i_nb
      integer nt
      integer ng
      integer i,j
      integer i_pr,i_mr,i_xr
      integer i_im,i_nm
      character*64 msg

      pname = "smd_geom_init_system"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "atom:res_type"
      call smd_get_ind_dim(tag,i_xr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     reset residue type to MM
c     ------------------------
      do i=1,nr
        int_mb(i_xr+i-1) = i_tmm
      end do
      call util_message("out "//pname)

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
#include "stdio.fh"
      integer im
c     
      character*(smd_string_size) sp_geom
      character*(smd_string_size) tag,pname
      logical result
      integer h_itmp,i_itmp
      integer nb,nr
      integer i_nb
      integer nt
      integer ng
      integer i,j
      integer i_pr,i_mr
      integer i_im,i_nm
      character*64 msg

      pname = "smd_geom_init_system"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get locations of residues in the atomic arrays
c     ----------------------------------------------
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
c     check if attempting to compute nonexisting residue
      if(im.gt.nr)   
     >  call errquit(
     >       pname//'residue index out of bounds',0, RTDB_ERR)
c
c     determine neighbors around im
c     -----------------------------
      if(.not.ma_push_get(mt_int,nr,'tmp i',h_itmp,i_itmp))
     + call errquit(pname//'Failed to allocate memory for h_itmp',0,0)
c
      nb = nr
      call smd_geom_find_neighbors(.true.,im,nb,int_mb(i_itmp))
c      if(ga_nodeid().eq.0) then
c        write(luout,*) "cluster information",
c     +    (int_mb(i_itmp+i-1),i=1,nb)
c        call util_flush(luout)
c      end if
      tag = "neighbor_index"
      call smd_data_put_int(sp_geom,tag,nb,int_mb(i_itmp))
      call smd_data_get_index(sp_geom,tag,i_nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
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
c    
c     create geometry
c     --------------
      call smd_geom_create()

      call util_message("out "//pname)

      return
      end
C>    Creates geometry object and
C>    stores it in rtdb
C>    Input parameters are passed silently
C>    through the smd objects (geom,coords,atom)
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
#include "smd_system.fh"
#include "stdio.fh"
c
      logical geom_cart_set_gen1
      external geom_cart_set_gen1
c      
      integer rtdb,m,n,ncent
c
      integer i,j
      integer i0,i1,n0,k
      character*20 pname
      character*64 msg
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
      integer igroup
      integer un
      logical operiodic
      integer im
      integer i_gim,i_gnm
c
      pname = "smd_geom_create"
c      call smd_geom_repair()
c      stop

c
      call smd_task_get_counter(im)
      write(msg,'(A20," 1 ",I3)') pname,im
c      call smd_coords_check(msg)
      operiodic = smd_system_periodic()
c
      call smd_group_id_get(igroup)
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
      if(.not.ma_push_get(mt_dbl,3*n1,pname//'tmp c',h_c1,i_c1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_dbl,n1,'tmp q',h_q1,i_q1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_dbl,n1,'tmp m',h_ms1,i_ms1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
      if(.not.ma_push_get(mt_byte,n1*smd_string_size,'tmp t',h_t1,i_t1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c
c     determine translation vector to put 1st residue in the center
c     -------------------------------------------------------------
c      if(operiodic) then
c        i0=int_mb(i_im+i-1)-1
c        call smd_geom_shift(int_mb(i_nm),
c     +                  dbl_mb(i_c+3*()),shift)
c      else
c        do i=1,3
c          shift(i) = 0.0d0
c        end do
c      end if
c      tag = "shift vector"
c      call smd_data_put_dbl(sp_geom,tag,3,shift)
      tag = "shift vector"
      call smd_data_get_dbl(sp_geom,tag,3,shift,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting shift vector '//tag,0, RTDB_ERR)

      tag = "res-index"
      call smd_data_create_get(sp_geom,tag,m,MT_INT,i_gim)
      tag = "res-size"
      call smd_data_create_get(sp_geom,tag,m,MT_INT,i_gnm)
c
c     copy geometry of the region into temp buffers
c     ---------------------------------------------
      j=0
      do i=1,m
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       int_mb(i_gim+i-1)=j+1
       int_mb(i_gnm+i-1)=n0
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
c      if(operiodic)
c     + call smd_lat_rebox(n1,dbl_mb(i_c1)) 
      if(operiodic)
     + call smd_lat_rebox_segments(n1,m,int_mb(i_gim),int_mb(i_gnm),
     +                     dbl_mb(i_c1)) 

C      write(*,*) "before shift"
C      do i=1,n1
C        write(*,*) (dbl_mb(i_c1+(i-1)*3+k-1),k=1,3)
C      end do
      do i=1,n1
         do k=1,3
           dbl_mb(i_c1+(i-1)*3+k-1)=
     +         dbl_mb(i_c1+(i-1)*3+k-1)+shift(k)
         end do
       end do
C
C      write(*,*) "after shift"
C      do i=1,n1
C        write(*,*) (dbl_mb(i_c1+(i-1)*3+k-1),k=1,3)
C      end do
c
c
c     create geometry
c     ---------------
       if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      i0=1
      if(.not.geom_cart_set_gen1(igeom,i0,n1,n1,smd_string_size,
     + byte_mb(i_t1),dbl_mb(i_c1),dbl_mb(i_q1)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,n1,dbl_mb(i_ms1)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
      if(.not.geom_ncent(igeom, ncent))
     $           call errquit(pname,0,RTDB_ERR)

c
      if(ga_nodeid().eq.0) then
      call smd_file_open(un,"tmp.xyz")
      if (.not. geom_print_xyz(igeom,un))
     + call errquit('qmmmsi geom_print_xyz',0, GEOM_ERR)
      call util_flush(un)
      close(un)
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
c     set active atoms for optimizations (this really should not be here)
c     -----------------------------------------------------------------
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

      if(.not.ma_chop_stack(h_c1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)


c      write(msg,'(A20," 2 ",I3)') pname,im
c      call smd_coords_check(msg)
      return
c
  
      end
c
      subroutine smd_geom_rebuild(nt,nm,i0,c)
c     rebuilds potentially cut segment
c     nt      [in]  total number of atoms
c     nm      [in]  number of atoms in the segment
c     i0      [in]  starting index for the segment
c     c[3,nt] [in]  coordinates of the overall system
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
      do i=1,nm
      do j=1,3
         c(j,i0+i-1)=dbl_mb(i_ftmp+(i-1)*3+j-1)+cshift(j)
      end do
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
      subroutine smd_geom_cog1(nt,c,cog)
c     finds center of geometry for a give segment
c     nt      [in]  number of atoms
c     c[3,nt] [in]  coordinates 
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
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp i',h_ftmp,i_ftmp))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c
c     put 1st atom in the middle
c     ---------------------------
      do j=1,3
         cshift(j)=c(j,1)
      end do

      do i=1,nt
      do j=1,3
         dbl_mb(i_ftmp+(i-1)*3+j-1)=
     +         c(j,i)-cshift(j)
      end do
      end do
 
c     fold in other coordinates potentially
c     cutoff by the box boundary
      call smd_lat_rebox(nt,dbl_mb(i_ftmp)) 
c     calculate cog for shifted segment
c     --------------------------------- 
      cog(1) = 0.0d0
      cog(2) = 0.0d0
      cog(3) = 0.0d0
      
      do i=1,nt
      do j=1,3
         cog(j)=cog(j)+dbl_mb(i_ftmp+(i-1)*3+j-1)
      end do
      end do

      cog(1) = cog(1)/nt
      cog(2) = cog(2)/nt
      cog(3) = cog(3)/nt
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

      subroutine smd_geom_shift(nt,c,shift)
c     finds center of geometry for a give segment
c     nt      [in]  total number of atoms
c     c[3,nt] [in]  coordinates of the overall system
c     shift[3]  [out] center of mass coordinates
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
      double precision shift(3)
      character*72 tag,sp_lat
c
      integer i,j
      double precision cshift(3)
      integer i_ftmp,h_ftmp
      character*30 pname
      double precision a(3,3),dx(3)
      logical result
c
      pname = "smd_geom_shift"
c
      tag = "lattice"
      call smd_system_get_component(sp_lat,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
      tag = "lat:cell"
      call smd_data_get_dbl(sp_lat,tag,9,a,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no lat:cell'//tag,0,0)

      shift(1) = 0.0d0
      shift(2) = 0.0d0
      shift(3) = 0.0d0
      
      do i=1,nt
      do j=1,3
         shift(j)=shift(j)+c(j,i)
      end do
      end do

      shift(1) = shift(1)/nt
      shift(2) = shift(2)/nt
      shift(3) = shift(3)/nt
c
      return
c
  
      end

      subroutine smd_geom_coc(nt,nm,i0,q,c,cog)
c     finds center of charge   for a give segment
c     nt      [in]  total number of atoms
c     nm      [in]  number of atoms in the segment
c     i0      [in]  starting index for the segment
c     q[nt]   [in]  charges  of the overall system
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
      double precision q(nt)
      double precision c(3,nt)
      double precision cog(3)
c
      integer i,j
      double precision cshift(3)
      integer i_ftmp,h_ftmp
      character*30 pname
      double precision qtot
c
      pname = "smd_geom_coc"
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
      
      qtot = 0.0d0
      do i=1,nm
      do j=1,3
         cog(j)=cog(j)+q(i0+i-1)*dbl_mb(i_ftmp+(i-1)*3+j-1)
      end do
      qtot = qtot + q(i0+i-1)
      end do

      cog(1) = cog(1)/qtot
      cog(2) = cog(2)/qtot
      cog(3) = cog(3)/qtot
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
C>    Updates local charge arrays with computed
C>    ESP charges
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
      character*32 sp_global
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
      integer ga_handle, ilo,ihi
c
      logical geom_cart_set1
      external geom_cart_set1
      logical geom_rtdb_ncent
      external geom_rtdb_ncent
c
      pname = "smd_geom_get_charges"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
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
c      This code updates global charge array.
c      Orginally
c      it was meant to serve as replacement for local charge update
c      as done below, keeping local charge array pristine
cc
cc     get global array
cc     ---------------
      call smd_data_get_int(sp_global,"charges",1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
      do i=1,nb
       i0=int_mb(i_im+i-1)
       n0=int_mb(i_nm+i-1)
       ilo = i0
       ihi = (i0-1)+n0
      if(ga_nodeid().eq.0) then
        call nga_put(ga_handle,ilo,ihi, dbl_mb(i_q1), 1)
      end if
      end do
cc      call smd_global_print_file("charges","charge.tmp")
      
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
      if(.not.ma_push_get(mt_dbl,3*ng,'tmp h_g1',h_g1,i_g1))
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
      subroutine smd_geom_get_dipole()
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
      character*32 tag,pname,sp_geom,theory,sp_atom
      character*255 key
      logical result
      integer i,j,i1,i0
      integer ng
      integer rtdb
      integer i_nb
      integer nb,k,n0
      double precision dipole(3)
      integer ihi,ilo,ic,nr
c
      pname = "smd_geom_get_dipole"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c    
c
      call smd_rtdb_get_handle(rtdb)
c
c     get the theory 
c     --------------
      if (.not.(rtdb_cget(rtdb,'task:theory',1,theory)))
     +     call errquit('rtdb_cget task:theory failed',0,0)
c
       key = theory(1:inp_strlen(theory))//":dipole"
         if (.not.rtdb_get(rtdb,key,mt_dbl,3,
     +            dipole))
     +            call errquit(key,
     +            0,RTDB_ERR)
         if(ga_nodeid().eq.0) 
     +     write(*,*) "dipole is",dipole(1),dipole(2),
     +                 dipole(3)
c      
      tag = "neighbor_index"
      call smd_data_get_ind_dim(sp_geom,tag,i_nb,nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      if(nb.gt.1) 
     >  call errquit(
     >       pname//'not a monomer',0, RTDB_ERR)
c
      ic = int_mb(i_nb)
      ilo=3*(ic-1)+1
      ihi=3*(ic-1)+3
      tag = "res_dipoles"
      write(*,*) "ic,ilo,ihi",ic,ilo,ihi
      call smd_data_put_dbl1(sp_atom,tag,ilo,ihi,dipole)
      call smd_data_namespace_print(sp_atom,6)
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
#include "stdio.fh"
#include "smd_system.fh"
c
      integer im
c     
      character*32 tag,sp_geom
      character*20 pname
      logical result
      integer i_c
      integer i_c1,h_c1
      integer i,j,i1,i0
      integer nt
      integer rtdb
      integer i_im,i_nm,i_nb
      integer nb,k,n0
      integer igeom
      double precision sht(3)
      integer ncent
      character*64 msg
      integer i_gim,i_gnm
      logical operiodic
c
      pname = "smd_geom_get_coords"
c
      operiodic = smd_system_periodic()
c
      call smd_task_get_counter(im)
c
c      write(msg,'(A20," 1 ",I3)') pname,im
c      call smd_coords_check(msg)
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
      tag = "neighbor_index"
      call smd_data_get_index(sp_geom,tag,i_nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "res-index"
      call smd_data_get_index(sp_geom,tag,i_gim,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "res-size"
      call smd_data_get_index(sp_geom,tag,i_gnm,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
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
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp h_c1',h_c1,i_c1))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c    
      if(.not.geom_ncent(igeom, ncent))
     $           call errquit(pname,0,RTDB_ERR)
      if(ncent.ne.nt) then
        call smd_task_get_counter(im)
        write(luout,*) "node",ga_nodeid()
        write(luout,*) "ncent,nt,im",ncent,nt,im
        write(luout,*) "current neighbor index",
     +       (int_mb(i_nb+i-1),i=1,nb)
c
        call util_flush(luout)
        call errquit(
     >       pname//' ncent ne nt',im,0 )
      end if

      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c1)))
     $           call errquit(pname,0,RTDB_ERR)
      call dscal(3*nt,cau2ang,dbl_mb(i_c1),1)
c
      tag = "shift vector"
      call smd_data_get_dbl(sp_geom,tag,3,sht,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
       
c      do i=1,nt
c        do k=1,3
c          dbl_mb(i_c1+3*(i-1)+k-1)=dbl_mb(i_c1+3*(i-1)+k-1)+
c     >                             sht(k)
c        end do
c      end do
      if(operiodic)
     + call smd_lat_rebox_segments(nt,nb,int_mb(i_gim),int_mb(i_gnm),
     +                     dbl_mb(i_c1)) 
c      call smd_lat_rebox(nt,dbl_mb(i_c1)) 
      
c
      j=0
      do i=1,nb
       i0=int_mb(i_im+i-1)-1
       n0=int_mb(i_nm+i-1)
       do i1=1,n0
         j = j+1
         if(j.gt.nt) 
     +     call errquit(
     >       pname//'out of bounds j nt',j,nt )
         do k=1,3
           dbl_mb(i_c+3*(i0+i1-1)+k-1) = dbl_mb(i_c1+3*(j-1)+k-1)
         end do
       end do
      end do

c      write(msg,'(A20," 2 ",I3)') pname,im
c      call smd_coords_check(msg)
       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

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

      subroutine smd_geom_find_neighbors(oinclude,icr,nb,inb)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
#include "stdio.fh"
c     if oinclude central residue is included
      logical oinclude
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

      if(oinclude) then
        j = 1
        inb(j) = icr
      else
        j = 0
      end if
      do i=1,nr
       if(int_mb(i_inb+i-1).ne.0) then
         j = j+1
         inb(j) = i
       end if
      end do
      nb = j


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

      return

      end

c $Id$
