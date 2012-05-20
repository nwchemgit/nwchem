
      subroutine smd_exlist_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_exlist,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_exlist_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "excl_list"
      call smd_system_get_component(sp_exlist,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_exlist_init(sp_exlist,sp_atom,result)
c
      if(.not.result) then
       tag = "excl_list"
       call smd_system_unset_component(tag)
      end if

      return
      end

      subroutine smd_exlist_init(sp_exlist,sp_atom,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_exlist
      character*(*) sp_atom
      logical result
c
      character*32 pname
      character*80 tag
      integer np,nl
      integer i_p,i_l,i_ir
      integer h_l
      integer i_list
      integer i
c
      pname = "smd_exlist_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms
c     -------------------------
      call smd_atom_ntot(sp_atom,np)
c
c     gestimate the size of pair list
c     ------------------------------
      nl =  min( np*200, ma_inquire_avail(MT_INT))
c
c     create pointer memory
c     ---------------------
      call smd_namespace_create(sp_exlist)
      tag = "exlist:pointer"
      call smd_data_create(sp_exlist,"exlist:pointer",np,MT_INT)
      call smd_data_get_index(sp_exlist,tag,i_p,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
c    create temporary scratch array for list since
c    we do not know the size yet
c    ---------------------------------------------
      if(.not.ma_push_get(mt_int,nl,'tmp l',h_l,i_l))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nl, MA_ERR)

c
c     get atom residue index for exclusion criteria
c     ---------------------------------------------
      tag = "atom:resid"
      call smd_data_get_index(sp_atom,tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)


      call smd_exlist_set(nl,
     +                    int_mb(i_l),
     +                    np,
     +                    int_mb(i_p),
     +                    int_mb(i_ir))

c
c     create list memory
c     nl now contains the actual size
c     -------------------------------

      if(nl.eq.0) then
        result = .false.
        call smd_namespace_destroy(sp_exlist)
        goto 200
      end if
c
c      write(*,*) "excluded size list",nl
      tag = "exlist:list"
      call smd_data_create(sp_exlist,tag,nl,MT_INT)
      call smd_data_get_index(sp_exlist,tag,i_list,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)


      do i=1,nl
       int_mb(i_list+i-1) = int_mb(i_l+i-1)
      end do 

      if(.not.ma_pop_stack(h_l))
     & call errquit(pname//'Failed to deallocate stack h_l',nl,
     &       MA_ERR)

200   continue
      return
      end
c
c
      subroutine smd_exlist_set(nl,list,natms,point,resid)
c
c     constructs excluded pairt list
c     point(i) is an index to list() array
c     that contains all the pairs of atom i
c     In other words all the atoms paired with atom i
c     are contained in list(point(i)),...,list(point(i+1)-1)
c     
      implicit none
#include "errquit.fh"
      integer nl,natms
      integer list(nl)
      integer point(natms)
      integer resid(natms)
c
      integer i,j
      integer nlist

      character*30 pname
 
      pname = "smd_exlist_set"

      nlist = 0
      do i=1,natms-1

       point(i)=nlist+1

       do j=i+1,natms

         if(resid(i).eq.resid(j))then

          nlist=nlist+1

          if(nlist.gt.nl)  
     +     call errquit(pname//'exceeded list dimensions',
     +     nl, 0)

          list(nlist)=j

         endif

       enddo
      enddo
c
c     we need to set this to mark
c     the end of pair list index for natms-1
c     it would be used as point(natms)-1
c     -------------------------------------
      point(natms)=nlist+1

      nl = nlist

      return

      END
c $Id$
