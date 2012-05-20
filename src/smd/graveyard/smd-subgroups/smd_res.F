      subroutine smd_res_init(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer rtdb
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer nr,nt
      integer i_s,i_t,i_ind
      logical result
c
      pname = "smd_res_init"
c
      write(*,*) "in "//pname
c
      if(.not.rtdb_cget(rtdb,'smd:coordfile',1,filename))
     >  call errquit(
     >       pname//'error reading coordinate file',0, RTDB_ERR)
      write(*,*) "found coordinates",filename
c
c     get number of residues
c     ----------------------
      call smd_pdb_nres(filename,nr)
      write(*,*) "number of res",nr
c
c     get number of atoms
c     -------------------
      call smd_atom_ntot(namespace,nt)
      write(*,*) "number of atoms",nt
c
c     create residue data structures
c     ------------------------------
      call smd_data_create(namespace,"res:name",16*nr,MT_BYTE)
      call smd_data_create(namespace,"res:size",nr,MT_INT)
      call smd_data_create(namespace,"res:",nr,MT_INT)

      tag = "atom:res"
      call smd_data_get_index(namespace,tag,i_ind,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "res:size"
      call smd_data_get_index(namespace,tag,i_s,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "res:name"
      call smd_data_get_index(namespace,tag,i_t,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)


      call smd_pdb_read_res(filename,nt,nr,
     +                         byte_mb(i_t),
     +                         int_mb(i_ind),
     +                         int_mb(i_s))

      return
      end
c $Id$
