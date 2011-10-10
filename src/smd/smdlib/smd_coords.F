      subroutine smd_coords_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_coords,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_coords_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_coords_init(sp_coords)
      
      call smd_coords_read(sp_coords)

      call smd_coords_print_pdb("old.pdb")

      return
      end

      subroutine smd_coords_init(sp_coords)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_coords
c
      character*32 pname
      integer na
c
      pname = "smd_coords_init"
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
c     create coords data structures
c     ---------------------------
      call smd_namespace_create(sp_coords)
      call smd_data_create(sp_coords,"coords",3*na,MT_DBL)

      return
      end

      subroutine smd_coords_read(sp_coords)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_coords
      integer rtdb
c
      character*32 pname
      character*72 tag
      integer i_c,i_c1
      integer na
      logical result
c
      pname = "smd_coords_read"
c
c      write(*,*) "in "//pname
c
c     fill in coordinates from pdb file if any
c     ----------------------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c1,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
cc
      na = na/3

      call smd_coordfile_read_coords1(na,
     +                          dbl_mb(i_c1))

      return
      end

C>    Folds coordinates back into the box
C>    while preserving molecules/segments
      subroutine smd_coords_rebox_segments()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      integer i_pr,i_mr,nr,ir
      integer i0,n0
      logical result

      pname = "smd_coords_rebox"
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
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3

      call smd_lat_rebox_segments(na,nr,int_mb(i_pr),int_mb(i_mr),
     >                    dbl_mb(i_c))


      return
      end

      subroutine smd_coords_rebox()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_coords_rebox"

      if(ga_nodeid().eq.0) write(*,*) "NO REBOXING"
      return
c
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3

      call smd_lat_rebox(na,
     >                    dbl_mb(i_c))

      return
      end

      subroutine smd_coords_check(msg)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c 
      character*(*) msg
c     
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c
      integer i_c1,h_c1
      logical result
      double precision diff,diff0
      integer i
      character*64 buffer

      pname = "smd_coords_check"

c
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3

      if(.not.ma_push_get(mt_dbl,3*na,pname//'tmp c',h_c1,i_c1))
     + call errquit(pname//'Failed to allocate memory for h_c1',0,0)

      if(ga_nodeid().eq.0) then
      call ga_brdcst(msg_qmmm_misc,dbl_mb(i_c),
     >     3*na*ma_sizeof(mt_int,1,mt_byte),0)
      else
      call ga_brdcst(msg_qmmm_misc,dbl_mb(i_c1),
     >     3*na*ma_sizeof(mt_int,1,mt_byte),0)
      end if

      diff0=0.0d0
      buffer = msg
      if(ga_nodeid().ne.0) then
      do i=1,3*na
         diff=dbl_mb(i_c+i-1)-dbl_mb(i_c1+i-1)
         diff = abs(diff)
         if(diff.gt.diff0) diff0=diff
         if(diff.gt.1e-4) then
           write(*,*) "big diff in ",
     +            buffer(1:inp_strlen(buffer)),diff0,
     +            ga_nodeid()
           call errquit(buffer,i,0)
         end if
      end do
      end if
      
      write(*,*) "coord check is ok in ",
     +            buffer(1:inp_strlen(buffer)),diff0,
     +            ga_nodeid()
      if(.not.ma_chop_stack(h_c1))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
      return
      end

      subroutine smd_coords_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer un
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_coords_print"

c
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3

      write(un,*) "printing coord" 
      call smd_util_print_force_array1(un,na,
     >                           dbl_mb(i_c))  

      return
      end

      subroutine smd_coords_update()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "smd_system.fh"
#include "smd_const_data.fh"
c     
      character*72 sp_vel
      character*72 sp_coords
      character*72 sp_mass
      character*72 sp_force
      character*72 sp_shakelist
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_f
      integer i_m,i_c,i_v,i_c1
      logical result
      integer rtdb
      double precision ekin, tstep
      logical oshake
      integer i_is1,i_is2,i_ds,ns
      integer h_ncc,i_ncc
      integer h_nvv,i_nvv
      integer h_dcc,i_dcc
      integer h_nrij,i_nrij
      integer h_orij,i_orij
c
      pname = "smd_coords_update"
c
c      write(*,*) "in "//pname
c
      call  smd_rtdb_get_handle(rtdb)
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
      tag = "mass"
      call smd_get_ind(tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get force array
c     ---------------
      tag = "force"
      call smd_get_ind(tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

c
c     get coord array
c     ---------------
      tag = "coords"
      call smd_get_ind(tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

c
c     get time step
c     -------------
      tag="smd:step"
        if (.not.rtdb_get(rtdb,tag,mt_dbl,1,tstep)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
 

      oshake = smd_system_shake()

      if(oshake) then

        tag = "shake:i1"
        call smd_get_ind(tag,i_is1,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        tag = "shake:i2"
        call smd_get_ind(tag,i_is2,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        tag = "shake:distance"
        call smd_get_ind_dim(tag,i_ds,ns,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        if(.not.ma_push_get(mt_dbl,na*3,'i_ncc',h_ncc,i_ncc))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,na*3,'i_nvv',h_nvv,i_nvv))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,na*3,'i_dcc',h_dcc,i_dcc))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,ns*3,'i_nrij',h_nrij,i_nrij))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,ns*3,'i_orij',h_orij,i_orij))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

      end if
      if(oshake) then
       call smd_leapf_shake(na,
     >                           ns,
     >                           tstep,
     >                           ekin,
     >                           dbl_mb(i_m),
     >                           int_mb(i_is1),
     >                           int_mb(i_is2),
     >                           dbl_mb(i_ds),
     >                           dbl_mb(i_ncc),
     >                           dbl_mb(i_nvv),
     >                           dbl_mb(i_dcc),
     >                           dbl_mb(i_nrij),
     >                           dbl_mb(i_orij),
     >                           dbl_mb(i_f),
     >                           dbl_mb(i_v),
     >                           dbl_mb(i_c))


      else
      call smd_leapf(na,
     >                 tstep,
     >                 ekin,
     >                 dbl_mb(i_m),
     >                 dbl_mb(i_f),
     >                 dbl_mb(i_v),
     >                 dbl_mb(i_c))

      end if

      call smd_energy_set_component("kinetic",ekin/convfct2)
      call smd_coords_rebox_segments()

      
      if(oshake) then

        if(.not.ma_pop_stack(h_orij))
     &   call errquit(pname//'Failed to deallocate stack h_orij',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_nrij))
     &   call errquit(pname//'Failed to deallocate stack h_nrij',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_dcc))
     &   call errquit(pname//'Failed to deallocate stack h_dcc',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_nvv))
     &   call errquit(pname//'Failed to deallocate stack h_nvv',0,
     &         MA_ERR)


        if(.not.ma_pop_stack(h_ncc))
     &   call errquit(pname//'Failed to deallocate stack h_ncc',0,
     &         MA_ERR)

      end if

      return
      end

      subroutine smd_coords_print_pdb(fname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      character*(*) fname
c     
      integer un
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer nta
      integer i,i0,j
      integer i_c,i_ta,i_tr,i_ir
      logical result
c
      pname = "smd_coords_print"
c
      if(.not.util_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
      open(unit=un,status="unknown",form="formatted",file=fname)
c
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3
c
      tag = "atom:name"
      call smd_get_ind(tag,i_ta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind(tag,i_ir,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resname"
      call smd_get_ind(tag,i_tr,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
 
      tag = "atom:name"
      call smd_get_size(tag,nta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting size for '//tag,0, RTDB_ERR)
c
      nta=nta/na
      do i=1,na
       i0=nta*(i-1)
       write(un,FMT=9000)
     >           i,
     >           (byte_mb(i_ta+i0+j-1),j=1,4),
     >           (byte_mb(i_tr+i0+j-1),j=1,3),
     >            int_mb(i_ir+i-1),
     >            (dbl_mb(i_c+3*(i-1)+j-1),j=1,3)

      end do
      close(un)
9000  FORMAT("ATOM",T7,I5,T13,4A,T18,3A,T23,
     >       I4,T31,F8.3,T39,F8.3,T47,F8.3)

      return
      end

C>    Generates lattice replicated pdb (-1,0,1)
C     in all three directions.
C>    for current coordinates
C>    Used to check reboxing routines
      subroutine smd_coords_print_pdb_repl(fname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      character*(*) fname
c     
      integer un
      character*32 sp_coords
c
      character*72 tag,sp_lat
      character*30 pname
      integer na
      integer nta
      integer i,i0,j
      integer i_c,i_ta,i_tr,i_ir
      logical result
      double precision a(3,3),dx(3)
      integer k1,k2,k3,ir
      integer ia
      integer nr
c
      pname = "smd_coords_print"
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
     >       pname//'no component '//tag,0,0)
c
      if(.not.util_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
      open(unit=un,status="unknown",form="formatted",file=fname)
c
c     get atomic coordinates
c     ----------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      na = na/3
c
      call smd_system_nres(nr)
c
      tag = "atom:name"
      call smd_get_ind(tag,i_ta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind(tag,i_ir,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resname"
      call smd_get_ind(tag,i_tr,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
 
      tag = "atom:name"
      call smd_get_size(tag,nta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting size for '//tag,0, RTDB_ERR)
c
      nta=nta/na
      ia = 0
      do k1=-1,1
      do k2=-1,1
      do k3=-1,1
       dx(1) = k1*a(1,1)+k2*a(1,2)+k3*a(1,3)
       dx(2) = k1*a(2,1)+k2*a(2,2)+k3*a(2,3)
       dx(3) = k1*a(3,1)+k2*a(3,2)+k3*a(3,3)
      do i=1,na
       i0=nta*(i-1)
       ia = ia + 1
       write(un,FMT=9000)
     >           ia,
     >           (byte_mb(i_ta+i0+j-1),j=1,4),
     >           (byte_mb(i_tr+i0+j-1),j=1,3),
     >            int_mb(i_ir+i-1)+(k1+k2+k3+3)*nr,
     >            dbl_mb(i_c+3*(i-1))+dx(1),
     >            dbl_mb(i_c+3*(i-1)+1)+dx(2),
     >            dbl_mb(i_c+3*(i-1)+2)+dx(3),dx(1),dx(2),dx(3)

      end do
      end do
      end do
      end do
      close(un)
9000  FORMAT("ATOM",T7,I5,T13,4A,T18,3A,T23,
     >       I4,T31,F8.3,T39,F8.3,T47,F8.3,3X,3(F8.3))

c
      call smd_coords_rebox_segments()
      call smd_coords_print_pdb("reboxed-new.pdb")
c

      return
      end

      subroutine smd_coords_save(oscratch)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      logical oscratch
c     
      character*256 fname
      character*30 pname
c
      pname = "smd_coords_save"
c
      if(ga_nodeid().eq.0) then
      call util_file_name("pdb", oscratch, .false., fname) 
      call smd_coords_print_pdb(fname)
      write(*,*) "saved coord file into",
     +    ga_pgroup_nodeid(ga_pgroup_get_world()),
     +    fname
      end if
      return
      end

      subroutine smd_coords_save_numbered(ic,oscratch)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      integer ic
      logical oscratch
c     
      character*256 fname
      character*30 pname,stub
c
      pname = "smd_coords_save"
c
      if(ga_nodeid().eq.0) then
      write(stub,'(I4.4,A)') ic,".pdb"
      call util_file_name(stub, oscratch, .false., fname) 
      call smd_coords_print_pdb(fname)
      write(*,*) "saved coord file into",
     +    ga_pgroup_nodeid(ga_pgroup_get_world()),
     +    fname
      end if
      return
      end

C>    Prints PDB file with coordinates different
C>    from those in coord object. 
C>    No check on array bounds is performed which
C>    may be unsafe.
C>    @param[in] c(*) coordinates to be printed 
      subroutine smd_coords_print_pdb_c(fname,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      character*(*) fname
      double precision c(*)
c     
      integer un
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer nta
      integer i,i0,j
      integer i_ta,i_tr,i_ir
      logical result
c
      pname = "smd_coords_print"
c
      if(.not.util_get_io_unit(un)) 
     >   call errquit("cannot get file number",0,0)
      open(unit=un,status="unknown",form="formatted",file=fname)
c
      tag = "atom:name"
      call smd_get_ind(tag,i_ta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind(tag,i_ir,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:resname"
      call smd_get_ind(tag,i_tr,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "atom:name"
      call smd_get_size(tag,nta,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting size for '//tag,0, RTDB_ERR)
c 
      nta=nta/na

      do i=1,na
       i0=nta*(i-1)
       write(un,FMT=9000)
     >           i,
     >           (byte_mb(i_ta+i0+j-1),j=1,4),
     >           (byte_mb(i_tr+i0+j-1),j=1,3),
     >            int_mb(i_ir+i-1),
     >            c(i),
     >            c(i+na),
     >            c(i+2*na)

      end do
      close(un)
9000  FORMAT("ATOM",T7,I5,T13,4A,T18,3A,T23,
     >       I4,T31,F8.3,T39,F8.3,T47,F8.3)

      return
      end

      subroutine smd_coord_write_default()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(255) filename
      logical result
      call smd_coordfile_output(filename,result)
c
      if(result) then
        call smd_coords_print_pdb(filename) 
      end if
      return
      end
c $Id$
