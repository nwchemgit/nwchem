      subroutine smd_atom_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_atom_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_atom_init(sp_atom)

      return
      end

      subroutine smd_atom_init(namespace)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
#include "geom.fh"
c     
      character*(*) namespace
c
      character*32 pname
      character*80 tag
      integer i
      integer nt,ns
      integer i_t,i_tr,i_ir,i_qf
      integer i_iconst
      logical result
      integer nr,nrmax
      integer i_sr
      integer i_mr,i_pr
      integer i_xr
c
      logical geom_tag_to_charge_gen
      external geom_tag_to_charge_gen
c
      pname = "smd_atom_init"
c
      call util_message("in "//pname)
c
c     get total number of atoms from pdb file
c     ---------------------------------------
      call smd_coordfile_natoms(nt)
      if(nt.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c
c     create atom namespace
c     ---------------------
      call smd_namespace_create(namespace)
c
c     create atom data structures
c     ---------------------------
      tag = "atom:iconst"
      call smd_data_create_get(namespace,tag,1,MT_INT,i_iconst)
      tag = "atom:name"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_t)
      tag  = "atom:resname"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_tr)
      tag = "atom:resid"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_sr)
      tag = "atom:resid0"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ir)
      tag = "atom:formal_charge"
      call smd_data_create_get(namespace,tag,nt,MT_DBL,i_qf)

c
c    don't ask
c    --------
      int_mb(i_iconst) = nt
c
c     fill in names from pdb file
c     ---------------------------
      call smd_coordfile_read_atomres(nt,
     +                         byte_mb(i_t),
     +                         byte_mb(i_tr),
     +                         int_mb(i_ir))

      tag = "atom:name"
      call smd_get_size(tag,ns,result)
      if(.not.result) call errquit(
     >       pname//'error getting size',0, RTDB_ERR)

      ns=ns/nt
c
      if(.not.geom_tag_to_charge_gen(nt,ns,
     >                        byte_mb(i_t),dbl_mb(i_qf)))
     >  call errquit(
     >       pname//'error setting formal charge',0, RTDB_ERR)
c
c     find/store total number of residues
c     -----------------------------
      call smd_atom_compute_nres(nr,nt,int_mb(i_ir))
      tag  = "nres"
      call smd_data_put_int(namespace,tag,1,nr)
c
c     determine res size and boundaries
c     ---------------------------------
      tag = "atom:res_size"
      call smd_data_create_get(namespace,tag,nr,MT_INT,i_mr)
      tag = "atom:res_pointer"
      call smd_data_create_get(namespace,tag,nr,MT_INT,i_pr)
      call smd_atom_compute_res(nt,nr,int_mb(i_ir),int_mb(i_sr),
     +                          int_mb(i_mr),int_mb(i_pr))
c
      tag = "res_dipoles"
      call smd_data_create(namespace,tag,3*nr,MT_DBL)
c
c     define residue type (QM, Bq, ...)
      tag = "atom:res_type"
      call smd_data_create_get(namespace,tag,nr,MT_INT,i_xr)
      call ifill(nr,i_tmm,int_mb(i_xr),1)
c
      call util_message("out "//pname)

      return
      end
c      
      subroutine smd_atom_dipole_read(oparallel)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "geom.fh"
#include "msgids.fh"
c     
      logical oparallel
c
      character*32 pname
      character*80 tag,filename
      integer i
      integer i_fd,nr,un
      logical result
      integer k
c
      pname = "smd_atom_dipole_read"
c
      call util_message("in "//pname)
      tag = "res_dipoles"
      call smd_get_ind_dim(tag,i_fd,nr,result)
      nr=nr/3
      if(.not.result) 
     + call errquit(pname//'Failed to get'//tag,
     + 0, 0)
      call smd_rtdb_get_string("smd:dipole:input",1,
     >                           filename,result)
      if(.not.result) return
c
      if(ga_nodeid().eq.0) then
      call smd_file_open(un,filename)
      do i=1,nr
        read(un,*) (dbl_mb(i_fd+(i-1)*3+k-1),k=1,3)
      end do
      close(un)
      end if
      if(oparallel) then
          call ga_brdcst(msg_smd,dbl_mb(i_fd),
     >     3*nr*ma_sizeof(mt_dbl,1,mt_byte),0)
      end if
      call util_message("out "//pname)

      return
      end

      subroutine smd_atom_dipole_write()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "geom.fh"
c     
c
      character*32 pname
      character*80 tag,filename
      integer i
      integer i_fd,nr,un
      logical result
      integer k
c
      pname = "smd_atom_dipole_read"
c
      call util_message("in "//pname)
      tag = "res_dipoles"
      call smd_get_ind_dim(tag,i_fd,nr,result)
      nr=nr/3
      if(.not.result) 
     + call errquit(pname//'Failed to get'//tag,
     + 0, 0)
      call smd_rtdb_get_string("smd:dipole:output",1,
     >                           filename,result)
      if(.not.result) return
c
      call smd_file_open(un,filename)
      do i=1,nr
        write(un,*) (dbl_mb(i_fd+(i-1)*3+k-1),k=1,3)
      end do
      call util_message("out "//pname)

      return
      end
c      
      subroutine smd_atom_compute_res(nt,nr,ir,is,mr,pr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer nt
      integer nr
      integer ir(nt)
      integer is(nt)
      integer mr(nr)
      integer pr(nr)
      integer i0
c
      integer i,j
      integer h_it,i_it
      integer irc
      character*30 pname
c
      pname = "smd_atom_compute_res_seq"
c
      if(.not.ma_push_get(mt_log,nr,'tmp',h_it,i_it))
     + call errquit(pname//'Failed to allocate memory for tmp',
     + nr, MA_ERR)
c
      do i=1,nr
        log_mb(i_it+i-1)=.false.
        pr(i) = 0
        mr(i) = 0
      end do

      irc=1
      is(1)=1
      pr(1)=1
      mr(1)=1
      do i=2,nt
        if(ir(i).ne.ir(i-1)) then
          if(log_mb(i_it+ir(i)-1))
     +      call errquit(pname//"discont. residue",
     +      0,0)

          log_mb(i_it+ir(i)-1)=.true.
          irc=irc+1
          pr(irc)=i
        end if        
        is(i) = irc
        mr(irc) = mr(irc)+1
      end do
c      
      if(.not.ma_pop_stack(h_it))
     & call errquit(pname//'Failed to deallocate stack h_it',0,
     &       MA_ERR)
c      
      end subroutine
c
      subroutine smd_atom_compute_res_seq(nt,ir,is)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer nt
      integer ir(nt)
      integer is(nt)
      integer i0
c
      integer i,j,nrmax
      integer h_it,i_it
      integer irc
      character*30 pname
c
      pname = "smd_atom_compute_res_seq"
c
c     find maximum residue number
c     ---------------------------
      nrmax = 0
      do i=1,nt
        if(ir(i).ge.nrmax) 
     +     nrmax = ir(i)
      end do
c
c
      if(.not.ma_push_get(mt_log,nrmax,'tmp',h_it,i_it))
     + call errquit(pname//'Failed to allocate memory for tmp',
     + nrmax, MA_ERR)
c
      do i=1,nrmax
        log_mb(i_it+i-1)=.false.
      end do

      irc=1
      is(1)=1
      do i=2,nt
        if(ir(i).ne.ir(i-1)) then
          if(log_mb(i_it+ir(i)-1))
     +      call errquit(pname//"discont. residue",
     +      0,0)

          log_mb(i_it+ir(i)-1)=.true.
          irc=irc+1
        end if        
        is(i) = irc
      end do
c      
      if(.not.ma_pop_stack(h_it))
     & call errquit(pname//'Failed to deallocate stack h_it',0,
     &       MA_ERR)
c      
      end subroutine
c
      subroutine smd_atom_compute_nres(nr,nt,ir)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer nr
      integer nt
      integer ir(nt)
      integer i0
c
      integer i,j,nrmax
      integer h_it,i_it
      character*30 pname
c
      pname = "smd_atom_compute_nres"
c
c     find maximum residue number
c     ---------------------------
      nrmax = 0
      do i=1,nt
        if(ir(i).ge.nrmax) 
     +     nrmax = ir(i)
      end do
c
      if(.not.ma_push_get(mt_int,nrmax,'tmp',h_it,i_it))
     + call errquit(pname//'Failed to allocate memory for tmp',
     + nrmax, MA_ERR)
c
      do i=1,nrmax
        int_mb(i_it+i-1)=0
      end do
c
      int_mb(i_it+ir(1)-1)=1
      do i=2,nt
        if(ir(i).ne.ir(i-1)) int_mb(i_it+ir(i)-1)=1
      end do
c
      nr = 0
      do i=1,nrmax
        nr = nr+int_mb(i_it+i-1)
      end do
c      
      if(.not.ma_pop_stack(h_it))
     & call errquit(pname//'Failed to deallocate stack h_it',0,
     &       MA_ERR)
c      
      end subroutine
c      
      subroutine smd_atom_get_dipoles()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c
      character*32 tag,pname
      logical result
      integer i_fd,nr
      integer rtdb
c
      pname = "smd_atom_get_dipoles"
c
      call smd_rtdb_get_handle(rtdb)
c
      tag = "res_dipoles"
      call smd_get_ind_dim(tag,i_fd,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      
      return
      end
ccc      
ccc      
c      subroutine smd_atom_compute_res_size(nt,ir,nr,sr)
c      implicit none
c#include "errquit.fh"
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "util.fh"
c#include "global.fh"
cc
c      integer nt
c      integer ir(nt)
c      integer nr
c      integer sr(nr)
c      integer i0
cc
c      integer i,j
cc
cc     count number of atoms in residue im
cc     contigious storage is presumed
cc     -----------------------------------
c      nm = 0
c      do i=1,nt
c       j = i
c       if(im.eq.ir(i)) then
c100     continue
c        nm = nm + 1
c        j = j + 1
c        if(im.eq.ir(j)) then
c         goto 100
c        else
c         goto 101 
c        end if
c       end if
c      end do
c101   continue
c      i0=i
c      return
cc
c  
c      end
c
      subroutine smd_atom_ntot(ntot)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer ntot
c
      character*72 tag
      character*30 pname
      integer i_iconst
      logical result

      pname = "smd_atom_ntot"
      tag = "atom:iconst"
      call smd_get_ind(tag,i_iconst,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting ntot '//tag,0, RTDB_ERR)
      ntot = int_mb(i_iconst)

      return
      end

      subroutine smd_atom_nres(namespace,nr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer nr
c
      character*72 tag
      character*30 pname
      integer i_ir
      logical result
      integer nt
      integer ir,ir0,i

      pname = "smd_atom_nres"
      tag = "atom:resid"
      call smd_data_get_index(namespace,tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting '//tag,0, RTDB_ERR)
      call smd_data_get_dim(namespace,tag,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting size '//tag,0, RTDB_ERR)

      nr = 0
      ir0 = 0
      do i=1,nt
       ir = int_mb(i_ir+i-1)
       if(ir0.ne.ir) then
          ir0=ir
          nr = nr + 1
       end if
      end do
       
      return
      end

      subroutine smd_atom_ntot1(namespace,ntot)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer ntot
c
      character*72 tag
      character*30 pname
      logical result

      pname = "smd_atom_ntot"
      tag = "atom:xyz"
      call smd_data_get_dim(namespace,tag,ntot,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting size for '//tag,0, RTDB_ERR)
      ntot = ntot/3

      return
      end

c $Id$
