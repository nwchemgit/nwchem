
      block data initial_mm_bond_coords_data
      implicit none
#include "mm_bond_coords_data.fh"

      data nbtot /0/
      data nbqm  /0/
      data nbmm  /0/

      end

      subroutine mm_bond_coords_init()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
  
      integer n      
      integer fn
      character*180 filename
      character*180 message
      character*30 pname
      character*2  card
      character*2  tag
      double precision  coords(3)
      integer indx
      integer nbqm0, nbmm0
 
      pname = "mm_bond_coords_init"   
            
c     deallocate all previous allocated arrays just in case
      call mm_bond_coords_end()

      if(.not.ma_alloc_get(MT_LOG,nqm,"add qm to bond list",
     &         h_addqm,i_addqm))  goto 911

      if(.not.ma_alloc_get(MT_LOG,nmm,"add mm to bond list",
     &         h_addmm,i_addmm))  goto 911

      call mm_bond_coords_map(nbqm0,nqm,log_mb(i_addqm),
     &     int_mb(i_iqm))

      call mm_bond_coords_map(nbmm0,nmm,log_mb(i_addmm),
     &     int_mb(i_imm))

      call mm_bond_coords_allocate(nbqm0,nbmm0)
      call mm_bond_coords_load()

      if(.not.ma_free_heap(h_addqm))
     & call errquit('mm: 
     &              Failed to deallocate stack addqm',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_addmm))
     & call errquit('mm: 
     &              Failed to deallocate stack addmm',nmm,
     &       MA_ERR)

      return
911   call errquit("error "//trim(message),0,
     >        -1)

      end 
 
      subroutine mm_bond_coords_map(nb,n,addtyp,indx_typ)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
      
      integer nb
      integer n
      logical addtyp(n)
      integer indx_typ(n)
c     
      integer i, j 
      character*32 pname
      integer indx_i,indx_j,indx_k,indx_l
     
      nb = 0
      addtyp(1:n) = .false.

c     loop over nbon
      do i=1,n
        do j=1,nbon
          indx_i = int_mb(i_ibon+j-1) 
          indx_j = int_mb(i_jbon+j-1)
          if(indx_i.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if

          if(indx_j.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if
        end do
      end do

c     loop over nangl 
      do i=1,n
        do j=1,nang
          indx_i = int_mb(i_iang+j-1)
          indx_j = int_mb(i_jang+j-1)
          indx_k = int_mb(i_kang+j-1)
          if(indx_i.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if
            
          if(indx_j.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if

          if(indx_k.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if
        end do
      end do

c     loop over ndih      
      do i=1,n
        do j=1,ndih
          indx_i = int_mb(i_idih+j-1)
          indx_j = int_mb(i_jdih+j-1)
          indx_k = int_mb(i_kdih+j-1)
          indx_l = int_mb(i_ldih+j-1)
          if(indx_i.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if

          if(indx_j.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if

          if(indx_k.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if

          if(indx_l.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nb = nb + 1
            end if
          end if
        end do
      end do
       
      end
 
      subroutine mm_bond_coords_allocate(nbqm0,nbmm0)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_coords_data.fh"
      integer nbqm0, nbmm0

      character*180 message
      character*30 pname

      integer nbtot0   
   
      pname = "mm_bond_coords_allocate"
      nbtot0 = nbqm0 + nbmm0
      if(nbtot0.ne.nbtot) then
        call mm_bond_coords_end()
        if(.not.ma_alloc_get(mt_dbl,3*nbtot0, 
     +                      "mm bonded coords",
     +                       h_br,i_br)) goto 911
        if(.not.ma_alloc_get(mt_byte,16*nbtot0,
     +                      'mm bonded tags',
     +                      h_bt,i_bt))  goto 911
        if(.not.ma_alloc_get(mt_int,nbtot0,
     +                      'mm bonded indices',
     +                      h_indx,i_indx)) goto 911
        if(.not.ma_alloc_get(mt_log,nbtot0,
     +                      'mm bonded indices',
     +                      h_lqm,i_lqm)) goto 911
      end if
      nbtot = nbtot0 
      nbqm  = nbqm0
      nbmm  = nbmm0

      call dfill(3*nbtot,0.0d0,dbl_mb(i_br),1)
      call ifill(nbtot,0.0,int_mb(i_indx),1)
      call mm_coords_tags_reset(nbtot,byte_mb(i_bt))
     
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_bond_coords_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_coords_data.fh"
c
      character*30 pname
      pname = "mm_bond_coords_end"
c
      if(nbtot.gt.0) then
          if (.not.ma_free_heap(h_bt))   goto 911
          if (.not.ma_free_heap(h_br))   goto 911
          if (.not.ma_free_heap(h_indx)) goto 911
          if (.not.ma_free_heap(h_lqm)) goto 911
          nbtot = 0
          nbqm  = 0
          nbmm = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return

      end

      subroutine mm_bond_coords_load()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      
      character*30  pname
      character*180 filename
      character*2   card
            
      integer n,fn 
      integer i
      integer i_btmm
      integer indx_qm, indx_mm
      logical addatom
      double precision coord(3)

      pname = "mm_bond_coords_load"
      
      i_btmm = i_bt+16*nbqm

      log_mb(i_lqm:i_lqm+nbtot-1) = .false.
      n = 0

      do i=1,nqm
        addatom = log_mb(i_addqm+i-1)
        indx_qm = int_mb(i_iqm+i-1)
        coord(1) = dbl_mb(i_rqm+3*(i-1))
        coord(2) = dbl_mb(i_rqm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rqm+3*(i-1)+2)
        if(addatom) then
          n = n + 1
          log_mb(i_lqm+n-1) = .true.
          int_mb(i_indx+n-1) = indx_qm
          dbl_mb(i_br+3*(n-1))   = coord(1)
          dbl_mb(i_br+3*(n-1)+1) = coord(2)
          dbl_mb(i_br+3*(n-1)+2) = coord(3)
        end if
        if(n == nbqm) exit
      end do
      
      n = 0

      do i=1,nmm
        addatom = log_mb(i_addmm+i-1)
        indx_mm = int_mb(i_imm+i-1)
        coord(1) = dbl_mb(i_rmm+3*(i-1))
        coord(2) = dbl_mb(i_rmm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rmm+3*(i-1)+2)
        if(addatom) then
          n = n + 1
          int_mb(i_indx+nbqm+n-1) = indx_mm
          dbl_mb(i_br+3*nbqm+3*(n-1))   = coord(1)
          dbl_mb(i_br+3*nbqm+3*(n-1)+1) = coord(2)
          dbl_mb(i_br+3*nbqm+3*(n-1)+2) = coord(3)
        end if
        if(n == nbmm) exit
      end do

      call mm_coords_tags_set(nbqm,"QM",byte_mb(i_bt))
      call mm_coords_tags_set(nbmm,"MM",byte_mb(i_btmm))
c     call mm_bond_coords_test()

      return
911   call errquit("error "//trim(pname),0,
     >        -1)
      return
      end

      subroutine mm_bond_qmcoords_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      integer n,fn
      integer geom
      integer nat
      logical status
      character*30 message
      character*30 pname
      double precision scale
      integer i, j
      integer iqm, ibqm
      integer i_ctmp, h_ctmp

      pname = "mm_bond_qmcoords_load"
c     write(*,*) pname

c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_load_geom',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_load_geom',0, RTDB_ERR)

c     get cart coordinates
c     --------------------
      status=geom_ncent(geom,nat)

      if(.not.ma_push_get(mt_dbl,3*nat,'ctmp',h_ctmp,i_ctmp))
     & call errquit( pname//'Failed to allocate memory for ctmp',
     & 3*nat, MA_ERR)

      if(.not. geom_cart_coords_get(geom,dbl_mb(i_ctmp))) 
     &    call errquit("mm:geom_cart_coords_get",0,0)

      call util_convert_units("au","angstrom",scale)
      call dscal(3*nat, scale,dbl_mb(i_ctmp),1)
  
      
      do i=1,nact
        iqm = int_mb(i_iact+i-1)
        do j=1,nbtot
          ibqm = int_mb(i_indx+j-1)
          if(ibqm.eq.iqm) then
            dbl_mb(i_br+3*(j-1)) = dbl_mb(i_ctmp+3*(i-1))
            dbl_mb(i_br+3*(j-1)+1) = dbl_mb(i_ctmp+3*(i-1)+1)
            dbl_mb(i_br+3*(j-1)+2) = dbl_mb(i_ctmp+3*(i-1)+2)
          end if
        end do
      end do

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('mm: 
     &              Failed to deallocate stack c_tmp',nat,
     &       MA_ERR)

      if(.not.geom_destroy(geom)) 
     &    goto 911

      return    
 
911   call errquit("error "//trim(pname),0,-1)

      end

      subroutine mm_bond_coords_test()
      implicit none
#include "mm_bond_coords_data.fh"
#include "mafdecls.fh"

      integer i,j

      do i=1,nbtot
         write(6,'(16A1,1X,3(1X,F10.6),I5)') 
     $        (byte_mb(i_bt+16*(i-1)+j),j=0,15),
     $        (dbl_mb(i_br+3*(i-1)+j),j=0,2),
     $        (int_mb(i_indx+i-1))
      end do

      end
