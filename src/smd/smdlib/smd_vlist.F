      subroutine smd_vlist_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_vlist,sp_exlist,sp_coords
      character*32 tag,pname
      logical result
      logical oexlist

      pname = "smd_vlist_init_system"
c
      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      oexlist = .true.
      tag = "excl_list"
      call smd_system_get_component(sp_exlist,tag,result)
      if(.not.result) oexlist = .false.

      tag = "verlet_list"
      call smd_system_get_component(sp_vlist,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_vlist_init(oexlist,sp_vlist,sp_exlist,sp_coords)
c
      return
      end
c
      subroutine smd_vlist_init(oexlist,sp_vlist,sp_exlist,sp_coord)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_vlist
      character*(*) sp_exlist
      character*(*) sp_coord
      logical oexlist
c
      character*32 pname
      character*80 tag
      integer np,nl
      integer i_p,i_l,i_c
      integer h_l
      integer i_list,i_clist
      integer i
      integer i_xp,i_xl
      integer i_cl,h_cl
      integer i_ct,h_ct
      double precision rc2
      integer nlb
      logical result
c
      pname = "smd_vlist_init"
c
c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      np = np/3
c
c     get excluded list information
c     -----------------------------
      if(oexlist) then
      tag = "exlist:pointer"
      call smd_get_ind(tag,i_xp,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      tag = "exlist:list"
      call smd_get_ind(tag,i_xl,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      end if
c
c     gestimate the size of pair list
c     ------------------------------
      nl =  min( 7*np*200, ma_inquire_avail(MT_INT))
      nl = nl/7
c
c     create pointer memory
c     ---------------------
      call smd_namespace_create(sp_vlist)
      tag = "vlist:pointer"
      call smd_data_create_get(sp_vlist,"vlist:pointer",np,MT_INT,i_p)
c
c    create temporary scratch array for list since
c    we do not know the size yet
c    ---------------------------------------------
      if(.not.ma_push_get(mt_int,nl,'tmp l',h_l,i_l))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nl, MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*nl,'tmp cl',h_cl,i_cl))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nl, MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*np,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nl, MA_ERR)


      call smd_cutoff_get_rcut_verlet(rc2)
      rc2=rc2*rc2

      if(oexlist) then
      call smd_vlist_set(np,
     +                   nl,
     +                   rc2,
     +                   dbl_mb(i_c),
     +                   int_mb(i_xp),
     +                   int_mb(i_xl),
     +                   int_mb(i_p),
     +                   int_mb(i_l),
     +                   dbl_mb(i_cl),
     +                   dbl_mb(i_ct),
     +                   result)

      else

      call smd_vlist_set1(np,
     +                   nl,
     +                   rc2,
     +                   dbl_mb(i_c),
     +                   int_mb(i_p),
     +                   int_mb(i_l),
     +                   dbl_mb(i_cl),
     +                   dbl_mb(i_ct),
     +                   result)


      end if
c
c     create list memory
c     nl now contains the actual size
c     we will buffer it though to allow for possible expansion
c     --------------------------------------------------------
      nlb = 500
      tag = "vlist:list"
      call smd_data_create_get(sp_vlist,tag,nl+nlb,MT_INT,i_list)

      tag = "vlist:distances"
      call smd_data_create_get(sp_vlist,tag,3*(nl+nlb),MT_DBL,i_clist)

      tag = "vlist:displacement"
      call smd_data_create(sp_vlist,tag,3*np,MT_DBL)
 
      do i=1,3*nl
       dbl_mb(i_clist+i-1) = dbl_mb(i_cl+i-1)
      end do 

      do i=1,nl
        int_mb(i_list+i-1)=int_mb(i_l+i-1)
      end do

      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nl,
     &       MA_ERR)

     
      if(.not.ma_pop_stack(h_cl))
     & call errquit(pname//'Failed to deallocate stack h_l',nl,
     &       MA_ERR)


      if(.not.ma_pop_stack(h_l))
     & call errquit(pname//'Failed to deallocate stack h_l',nl,
     &       MA_ERR)


      return
      end
c
      subroutine smd_vlist_update(olist,ocoord)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      logical olist
      logical ocoord
c
      character*32 sp_vlist
      character*32 sp_exlist
      character*32 sp_coord
c
      character*32 pname
      character*80 tag
      integer np,nl
      integer i_p,i_l,i_c
      integer h_l
      integer i_list,i_clist
      integer i
      integer i_xp,i_xl
      double precision rc2
      integer nlb
      logical result
      integer i_ct,h_ct
      logical oexlist 
c
      pname = "smd_vlist_update"
c
c
      if((.not.olist).and.(.not.ocoord)) return
c
c     get components
c     --------------
      tag = "coordinates"
      call smd_system_get_component(sp_coord,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      oexlist = .true.
      tag = "excl_list"
      call smd_system_get_component(sp_exlist,tag,result)
      if(.not.result) oexlist =.false.

      tag = "verlet_list"
      call smd_system_get_component(sp_vlist,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c
      if(.not.olist) then
        call smd_vlist_update_coord(sp_vlist,sp_coord)
        goto 200
      end if
c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      np = np/3 
c
c     get excluded list information
c     -----------------------------
      if(oexlist) then
      tag = "exlist:pointer"
      call smd_get_ind(tag,i_xp,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      tag = "exlist:list"
      call smd_get_ind(tag,i_xl,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      end if
c
c     get verlet list information
c     ---------------------------
      tag = "vlist:pointer"
      call smd_get_ind(tag,i_p,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "vlist:list"
      call smd_get_ind_dim(tag,i_list,nl,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
 
      tag = "vlist:distances"
      call smd_get_ind(tag,i_clist,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c    create temporary scratch array for list
c    ---------------------------------------------

      if(.not.ma_push_get(mt_dbl,3*np,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nl, MA_ERR)


      call smd_cutoff_get_rcut_verlet(rc2)
      rc2=rc2*rc2

      if(oexlist) then
      call smd_vlist_set(np,
     +                   nl,
     +                   rc2,
     +                   dbl_mb(i_c),
     +                   int_mb(i_xp),
     +                   int_mb(i_xl),
     +                   int_mb(i_p),
     +                   int_mb(i_list),
     +                   dbl_mb(i_clist),
     +                   dbl_mb(i_ct),
     +                   result)

      else

      call smd_vlist_set1(np,
     +                   nl,
     +                   rc2,
     +                   dbl_mb(i_c),
     +                   int_mb(i_p),
     +                   int_mb(i_list),
     +                   dbl_mb(i_clist),
     +                   dbl_mb(i_ct),
     +                   result)


      end if
       
      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nl,
     &       MA_ERR)


200   continue
      return
      end
c
      subroutine smd_vlist_update_coord(sp_vlist,sp_coord)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_vlist
      character*(*) sp_coord
c
      character*32 pname
      character*80 tag
      integer np,nl
      integer i_p,i_c
      integer i_list,i_clist
      logical result
      integer i
c
      pname = "smd_vlist_init"
c
c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      np = np/3
c
c     get list information
c     --------------------
      tag = "vlist:pointer"
      call smd_get_ind(tag,i_p,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "vlist:list"
      call smd_get_ind(tag,i_list,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)


      tag = "vlist:distances"
      call smd_get_ind_dim(tag,i_clist,nl,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      nl = nl/3

      call smd_vlist_update_coord0(np,nl,
     +                   dbl_mb(i_c),
     +                   int_mb(i_p),
     +                   int_mb(i_list),
     +                   dbl_mb(i_clist))

      call smd_lat_rebox(nl,dbl_mb(i_clist))

      return
      end
c
      subroutine smd_vlist_update_coord0(np,nl,c,point,
     >                         list,cl)
c
c     update coordinates in verlet pair list
c     point(i) is an index to list() array
c     that contains all the pairs of atom i
c     In other words all the atoms paired with atom i
c     are contained in list(point(i)),...,list(point(i+1)-1)
c     np     [in]  number of atoms (which is also size of pointer array)
c     nl     [in]  list size
c     c      [in]  coordinates
c     point  [in] verlet pointer
c     list   [in] verlet list
c     cl     [out] list of vectors (ri-rj)
      implicit none
#include "errquit.fh"
      integer np,nl
      double precision c(3,np)
      integer point(np)
      integer list(nl)
      double precision cl(3,nl)
c
      integer i,j
      integer nlist
      integer jnab,jbeg,jend

      character*30 pname
 
      pname = "smd_vlist_update_coord0"
      nlist=0
      do i=1,np-1
       jbeg=point(i)
       jend=point(i+1)-1

       if(jbeg.le.jend)then

        do jnab=jbeg,jend

         j=list(jnab)

         nlist = nlist + 1
         if(nlist.gt.nl)
     >       call errquit(
     >       pname//'out of bounds',0, RTDB_ERR)

         cl(1,nlist)=c(1,i)-c(1,j)
         cl(2,nlist)=c(2,i)-c(2,j)
         cl(3,nlist)=c(3,i)-c(3,j)

        enddo

       end if
      enddo


200   continue
      return

      END
c
      subroutine smd_vlist_set(np,nl,vcutsq,c,xp,xl,point,
     >                         list,cl,ct,result)
c
c     constructs verlet pairt list
c     point(i) is an index to list() array
c     that contains all the pairs of atom i
c     In other words all the atoms paired with atom i
c     are contained in list(point(i)),...,list(point(i+1)-1)
c     np     [in]  number of atoms (which is also size of pointer array)
c     nl     [inout] size of verlet list 
c     c      [in]  coordinates
c     xp     [in]  excluded pointer 
c     xl     [in]  excluded list
c     point  [out] verlet pointer
c     list   [out] verlet list
c     cl     [out] list of vectors (ri-rj)
c     result [out] status of subroutine
      implicit none
#include "errquit.fh"
      integer np
      integer nl
      double precision vcutsq
      double precision c(3,np)
      integer xp(np)
      integer xl(*)
      integer point(np)
      integer list(nl)
      double precision cl(3,nl)
      double precision ct(3,np)
      logical result
c
      integer i,j,k
      integer nlist
      integer eatm
      double precision rij(3),rijsq,cc(1,3)

      character*30 pname
 
      pname = "smd_vlist_set"

      result = .false.
      nlist=0

      do i=1,np-1

       point(i)=nlist+1
       if(xp(i).ne.xp(i+1))eatm=xp(i)

       k = 0
       do j=i+1,np

         k = k +1
         ct(1,k)=c(1,i)-c(1,j)
         ct(2,k)=c(2,i)-c(2,j)
         ct(3,k)=c(3,i)-c(3,j)

       end do


       call smd_lat_rebox(np,ct)


       k = 0
       do j=i+1,np

        k = k + 1
        if((xp(i).ne.xp(i+1)).and.(xl(eatm).eq.j))then

         eatm=min(eatm+1,(xp(i+1)-1))

        else

         rij(1)=ct(1,k)
         rij(2)=ct(2,k)
         rij(3)=ct(3,k)

         rijsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)


         if(rijsq.lt.vcutsq)then

          nlist=nlist+1

          if(nlist.gt.nl)then
           result = .false.
           goto 200
          endif

          list(nlist)=j
          cl(1,nlist)=rij(1)
          cl(2,nlist)=rij(2)
          cl(3,nlist)=rij(3)

         endif

        endif

       enddo

      enddo

      point(np)=nlist+1

      nl = nlist

      result = .true.

200   continue
      return

      END

      subroutine smd_vlist_set1(np,nl,vcutsq,c,point,
     >                         list,cl,ct,result)
c
c     constructs verlet pairt list
c     point(i) is an index to list() array
c     that contains all the pairs of atom i
c     In other words all the atoms paired with atom i
c     are contained in list(point(i)),...,list(point(i+1)-1)
c     np     [in]  number of atoms (which is also size of pointer array)
c     nl     [inout] size of verlet list 
c     c      [in]  coordinates
c     xp     [in]  excluded pointer 
c     xl     [in]  excluded list
c     point  [out] verlet pointer
c     list   [out] verlet list
c     cl     [out] list of vectors (ri-rj)
c     result [out] status of subroutine
      implicit none
#include "errquit.fh"
      integer np
      integer nl
      double precision vcutsq
      double precision c(3,np)
      integer point(np)
      integer list(nl)
      double precision cl(3,nl)
      double precision ct(3,np)
      logical result
c
      integer i,j,k
      integer nlist
      integer eatm
      double precision rij(3),rijsq,cc(1,3)

      character*30 pname
 
      pname = "smd_vlist_set"

      result = .false.
      nlist=0

      do i=1,np-1

       point(i)=nlist+1

       k = 0
       do j=i+1,np

         k = k +1
         ct(1,k)=c(1,i)-c(1,j)
         ct(2,k)=c(2,i)-c(2,j)
         ct(3,k)=c(3,i)-c(3,j)

       end do


       call smd_lat_rebox(np,ct)


       k = 0
       do j=i+1,np

        k = k + 1

         rij(1)=ct(1,k)
         rij(2)=ct(2,k)
         rij(3)=ct(3,k)

         rijsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)


         if(rijsq.lt.vcutsq)then

          nlist=nlist+1

          if(nlist.gt.nl)then
           result = .false.
           goto 200
          endif

          list(nlist)=j
          cl(1,nlist)=rij(1)
          cl(2,nlist)=rij(2)
          cl(3,nlist)=rij(3)

         endif

       enddo

      enddo

      point(np)=nlist+1

      nl = nlist

      result = .true.

200   continue
      return

      END

      SUBROUTINE smd_vlist_test(lupdate)
      implicit none
#include "smd_system.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c 
      logical lupdate
c
      integer na
      integer i_v,i_vd
      double precision t
      character*72 sp_vel
      character*72 sp_vlist
      logical result
      character*72 tag
      character*30 pname
      double precision vcut,rcut
c
      pname = "smd_vlist_test"
c
c     get velocity aray
c     -----------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "vel"
      call smd_get_ind_dim(tag,i_v,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3
c
c     get time step
c     ------------
      if(.not.smd_system_tstep(t)) 
     >  call errquit(
     >       pname//'no time step ',0,0)

c
c     get verlet displacement array
c     -----------------------------
      tag = "vlist:displacement"
      call smd_get_ind(tag,i_vd,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
 
       call smd_cutoff_get_rcut(rcut)
       call smd_cutoff_get_rcut_verlet(vcut)

       call smd_vlist_test0(na,t,rcut,vcut,
     >     dbl_mb(i_vd),dbl_mb(i_v),lupdate)

      return

      END

      SUBROUTINE smd_vlist_test0
     +    (natms,tstep,rcut,vcut,ivv,vvv,lupdate)

      implicit none
c
      integer natms
      double precision rcut,vcut
      double precision tstep
      logical lupdate
      double precision ivv(3,natms)
      double precision vvv(3,natms)
c
      integer i,exceed

      double precision  tstepsq
      double precision  dispmax,dispxsq,dispysq,dispzsq,disprsq 

      logical lnew

      data lnew/.true./

      save lnew

      tstepsq=tstep**2

      if(lnew)then

       lupdate=.true.
       lnew=.false.

       do i=1,natms

        ivv(1,i)=0.0
        ivv(2,i)=0.0
        ivv(3,i)=0.0

       enddo

      else

       lupdate=.false.

       dispmax=((vcut-rcut)/2.0)**2

       do i=1,natms

        ivv(1,i)=ivv(1,i)+vvv(1,i)
        ivv(2,i)=ivv(2,i)+vvv(2,i)
        ivv(3,i)=ivv(3,i)+vvv(3,i)

       enddo

       exceed=0

       do i=1,natms

        dispxsq=ivv(1,i)**2
        dispysq=ivv(2,i)**2
        dispzsq=ivv(3,i)**2
        disprsq=tstepsq*(dispxsq+dispysq+dispzsq)
        if(disprsq.gt.dispmax) then
          exceed=exceed+1
        end if
        if(exceed.ge.2)lupdate=.true.
       enddo

       if(lupdate)then

        do i=1,natms

         ivv(1,i)=0
         ivv(2,i)=0
         ivv(3,i)=0

        enddo

       endif

      endif

      return

      END
c $Id$
