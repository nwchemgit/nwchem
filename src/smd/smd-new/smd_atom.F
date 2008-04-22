      subroutine smd_atom_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 pname
      character*80 tag
      integer nt
      integer i_t,i_tr,i_ir,i_qf
      integer i_iconst
      logical result
      integer srtdb
c
      pname = "smd_atom_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms from pdb file
c     ---------------------------------------
      call smd_coordfile_natoms(nt)
#ifdef SMD_DEBUG
      write(*,*) "number of atoms",nt
#endif

c      write(*,*) "in 1"//pname
      if(nt.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
       call smd_srtdb_get_handle(srtdb)
       call smd_data_create(srtdb,"atom:resid",nt,MT_INT)
cc
cc     create atom namespace
cc     ---------------------
c      call smd_namespace_create(namespace)
cc      write(*,*) "in 2"//pname
cc
cc     create atom data structures
cc     ---------------------------
c      call smd_data_create(namespace,"atom:iconst",1,MT_INT)
c      call smd_data_create(namespace,"atom:name",16*nt,MT_BYTE)
cc      call smd_data_create(namespace,"atom:formal_charge",nt,MT_DBL)
c      call smd_data_create(namespace,"atom:resname",16*nt,MT_BYTE)
c      call smd_data_create(namespace,"atom:resid",nt,MT_INT)
c
cc
cc    don't ask
cc    --------
c      tag = "atom:iconst"
c      call smd_data_get_index(namespace,tag,i_iconst,result)
c      if(.not. result) 
c     >  call errquit(
c     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c      int_mb(i_iconst) = nt
cc
cc     fill in names from pdb file
cc     ---------------------------
c      tag = "atom:name"
c      call smd_data_get_index(namespace,tag,i_t,result)
c      if(.not. result) 
c     >  call errquit(
c     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c      tag = "atom:resname"
c      call smd_data_get_index(namespace,tag,i_tr,result)
c      if(.not. result) 
c     >  call errquit(
c     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c      tag = "atom:resid"
c      call smd_data_get_index(namespace,tag,i_ir,result)
c      if(.not. result) 
c     >  call errquit(
c     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
cc      tag = "atom:formal_charge"
cc      call smd_data_get_index(namespace,tag,i_qf,result)
cc      if(.not. result) 
cc     >  call errquit(
cc     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c
cc      write(*,*) "in 3"//pname,ga_nodeid()
c      call smd_coordfile_read_atomres(nt,
c     +                         byte_mb(i_t),
c     +                         byte_mb(i_tr),
c     +                         int_mb(i_ir))
c
cc      write(*,*) "in 4"//pname,ga_nodeid()
cc      call util_flush(6)
cc      if(.not.geom_tag_to_charge(nt,byte_mb(i_t),dbl_mb(i_qf)))
cc     >  call errquit(
cc     >       pname//'error setting formal charge',0, RTDB_ERR)
cc      write(*,*) "out "//pname,ga_nodeid()
cc      call util_flush(6)
c      return
      end

      subroutine smd_atom_end(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
c
      return
      end

