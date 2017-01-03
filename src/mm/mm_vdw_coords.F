      subroutine mm_vdw_coords_init()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_data.fh"
#include "mm_vdw_coords_data.fh"
      integer i,j
      integer nvtot0,nvqm0,nvmm0
      integer nv14tot0,nv14qm0,nv14mm0
      character*30 pname     
 
      pname="mm_vdw_coords_init"
c     write(*,*) pname

      if(.not.ma_alloc_get(mt_log,nqm,"add qm to vdw coords list",
     &   h_advqm, i_advqm)) goto 911      

      if(.not.ma_alloc_get(mt_log,nmm,"add mm to vdw coords list",
     &   h_advmm, i_advmm)) goto 911      

      if(.not.ma_alloc_get(mt_log,nqm,"add qm to vdw14 coords list",
     &   h_advqm14, i_advqm14)) goto 911      

      if(.not.ma_alloc_get(mt_log,nmm,"add mm to vdw14 coords list",
     &   h_advmm14, i_advmm14)) goto 911      

      call mm_vdw_coords_map(nvqm0,nqm,log_mb(i_advqm),
     &   int_mb(i_iqm),nvdw,int_mb(i_ivdw),int_mb(i_jvdw))

      call mm_vdw_coords_map(nvmm0,nmm,log_mb(i_advmm),
     &   int_mb(i_imm),nvdw,int_mb(i_ivdw),int_mb(i_jvdw))

      call mm_vdw_coords_map(nv14qm0,nqm,log_mb(i_advqm14),
     &   int_mb(i_iqm),nvdw14,int_mb(i_ivdw14),int_mb(i_jvdw14))

      call mm_vdw_coords_map(nv14mm0,nmm,log_mb(i_advmm14),
     &   int_mb(i_imm),nvdw14,int_mb(i_ivdw14),int_mb(i_jvdw14))

      call mm_vdw_coords_allocate(nvqm0,nvmm0)
      call mm_vdw14_coords_allocate(nv14qm0,nv14mm0)

      call mm_vdw_coords_load(nvtot,nvqm,nvmm,dbl_mb(i_rvdw),
     &     int_mb(i_icvdw),log_mb(i_lvqm),
     &     log_mb(i_advqm),log_mb(i_advmm),.false.)

      call mm_vdw_coords_load(nv14tot,nv14qm,nv14mm,dbl_mb(i_rvdw14),
     &     int_mb(i_icvdw14),log_mb(i_lv14qm),
     &     log_mb(i_advqm14),log_mb(i_advmm14),.true.)
       
      if(.not.ma_free_heap(h_advqm))
     & call errquit('mm: 
     &              Failed to deallocate stack advqm',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_advmm))
     & call errquit('mm: 
     &              Failed to deallocate stack advmm',nmm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_advqm14))
     & call errquit('mm: 
     &              Failed to deallocate stack advqm14',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_advmm14))
     & call errquit('mm: 
     &              Failed to deallocate stack advqm14',nmm,
     &       MA_ERR)

      return
911   call errquit("error "//trim(pname),0,
     >        -1)
    
      end
 
      subroutine mm_vdw_coords_map(nvt,nt,addtyp,
     &           indx_typ,nij,i_i,i_j)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer nvt
      integer nt
      logical addtyp(nt)
      integer indx_typ(nt)
      integer nij
      integer i_i(nij),i_j(nij)

      integer i,j
      integer indx_i, indx_j
      
      nvt = 0
      addtyp(1:nt) = .false.      
 
      do i=1,nt
        do j=1,nij
          indx_i = i_i(j)
          indx_j = i_j(j)
          if(indx_i.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nvt = nvt + 1
            end if
          end if

          if(indx_j.eq.indx_typ(i)) then
            if(.not.addtyp(i)) then
              addtyp(i) = .true.
              nvt = nvt + 1
            end if
          end if
        end do
      end do

      end

      subroutine mm_vdw_coords_allocate(nvqm0,nvmm0)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_vdw_coords_data.fh"
      integer nvqm0, nvmm0

      character*180 message
      character*30 pname
          
      integer nvtot0

      pname = "mm_vdw_coords_allocate"
      nvtot0 = nvqm0 + nvmm0
      if(nvtot0.ne.nvtot) then
        call mm_vdw_coords_end()
        if(.not.ma_alloc_get(mt_dbl,3*nvtot0,
     &                        "mm vdw coords",
     &                        h_rvdw,i_rvdw)) goto 911   
        if(.not.ma_alloc_get(mt_int,nvtot0,
     &                        "mm vdw indices",
     &                        h_icvdw,i_icvdw)) goto 911   
        if(.not.ma_alloc_get(mt_log,nvtot0,
     &                        "mm vdw quantum flags",
     &                        h_lvqm,i_lvqm)) goto 911   
      end if
      nvtot = nvtot0
      nvqm = nvqm0
      nvmm = nvmm0
      call dfill(3*nvtot,0.0d0,dbl_mb(i_rvdw),1)
      call ifill(nvtot,0.0,int_mb(i_icvdw),1)

      return
911   call errquit("error "//trim(pname),0,-1)
      end

      subroutine mm_vdw14_coords_allocate(nv14qm0,nv14mm0)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_vdw_coords_data.fh"
      integer nv14qm0, nv14mm0

      character*180 message
      character*30 pname
          
      integer nv14tot0

      pname = "mm_vdw14_coords_allocate"
      nv14tot0 = nv14qm0 + nv14mm0
      if(nv14tot0.ne.nv14tot) then
        call mm_vdw14_coords_end()
        if(.not.ma_alloc_get(mt_dbl,3*nv14tot0,
     &                        "mm vdw14 coords",
     &                        h_rvdw14,i_rvdw14)) goto 911   
        if(.not.ma_alloc_get(mt_int,nv14tot0,
     &                        "mm vdw14 indices",
     &                        h_icvdw14,i_icvdw14)) goto 911   
        if(.not.ma_alloc_get(mt_dbl,nv14tot0,
     &                        "mm vdw14 charges",
     &                        h_chgmm14,i_chgmm14)) goto 911   
        if(.not.ma_alloc_get(mt_log,nv14tot0,
     &                        "mm vdw14 quantum flags",
     &                        h_lv14qm,i_lv14qm)) goto 911   
      end if
      nv14tot = nv14tot0
      nv14qm = nv14qm0
      nv14mm = nv14mm0
      call dfill(3*nv14tot,0.0d0,dbl_mb(i_rvdw14),1)
      call dfill(nv14tot,0.0d0,dbl_mb(i_chgmm14),1)
      call ifill(nv14tot,0.0,int_mb(i_icvdw14),1)

      return
911   call errquit("error "//trim(pname),0,-1)
      end

      subroutine mm_vdw_coords_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_coords_data.fh"

      character*30 pname
      pname = "mm_vdw_coords_end"

      if(nvtot.gt.0) then
          if (.not.ma_free_heap(h_rvdw))   goto 911
          if (.not.ma_free_heap(h_icvdw))   goto 911
          if (.not.ma_free_heap(h_lvqm))   goto 911
          nvtot = 0
          nvqm = 0
          nvmm = 0
      end if

      return
911   call errquit("error "//trim(pname),0,-1)

      end

      subroutine mm_vdw14_coords_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_coords_data.fh"

      character*30 pname
      pname = "mm_vdw14_coords_end"

      if(nv14tot.gt.0) then
          if (.not.ma_free_heap(h_rvdw14))   goto 911
          if (.not.ma_free_heap(h_icvdw14))   goto 911
          if (.not.ma_free_heap(h_chgmm14))   goto 911
          if (.not.ma_free_heap(h_lv14qm))   goto 911
          nv14tot = 0
          nv14qm = 0
          nv14mm = 0
      end if

      return
911   call errquit("error "//trim(pname),0,-1)

      end

      subroutine mm_vdw_coords_deallocate()
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
        call mm_vdw_coords_end()
        call mm_vdw14_coords_end()
 
      end 

      subroutine mm_vdw_coords_load(nvt,nvq,nvm,c,ic,lvqm,
     +                              advqm,advmm,load_mmchg)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_coords_data.fh"
      integer nvt,nvq,nvm
      double precision c(3*nvt)
      integer ic(nvt)
      logical lvqm(nvt)
      logical advqm(nqm)
      logical advmm(nmm)
      logical load_mmchg

      character*30 pname
      integer n
      integer i
      integer indx_qm, indx_mm
      logical addatom
      double precision coord(3)

      pname = "mm_vdw_coords_load"

      lvqm(1:nvt) = .false.
      n = 0 

      do i=1,nqm
        addatom = advqm(i)
        indx_qm = int_mb(i_iqm+i-1)
        coord(1) = dbl_mb(i_rqm+3*(i-1))
        coord(2) = dbl_mb(i_rqm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rqm+3*(i-1)+2)
        if(addatom) then
          n = n + 1
          lvqm(n) = .true.
          ic(n) = indx_qm
          c(3*(n-1)+1)   = coord(1)
          c(3*(n-1)+2)   = coord(2)
          c(3*(n-1)+3)   = coord(3)
        end if
        if(n == nvq) exit
      end do

      n = 0     

      do i=1,nmm
        addatom = advmm(i)
        indx_mm = int_mb(i_imm+i-1)
        coord(1) = dbl_mb(i_rmm+3*(i-1))
        coord(2) = dbl_mb(i_rmm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rmm+3*(i-1)+2)
        if(addatom) then
          n = n + 1
          ic(nvq+n) = indx_mm
          c(3*(nvq+n-1)+1)   = coord(1)
          c(3*(nvq+n-1)+2)   = coord(2)
          c(3*(nvq+n-1)+3)   = coord(3)
          if(load_mmchg) then
            dbl_mb(i_chgmm14+nvq+n-1) = dbl_mb(i_chgmm+i-1)
          end if
        end if
        if(n == nvm) exit
      end do

c     call mm_vdw_coords_test(nvt,c,ic)
      return
911   call errquit("error "//trim(pname),0,
     >        -1)
      return
      end

      subroutine mm_vdw_qmcoords_load(rtdb)
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
#include "mm_vdw_coords_data.fh"
      integer rtdb      

      integer n,fn
      integer geom
      integer nat
      logical status
      character*30 message
      character*30 pname
      double precision scale

c     local variables
      character*16 t(nqm)
      integer i, j
      integer iqm, ivqm
      integer i_ctmp, h_ctmp

      pname = "mm_vdw_qmcoords_load"
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
      call dscal(3*nact, scale,dbl_mb(i_ctmp),1)

      do i=1,nact
        iqm = int_mb(i_iact+i-1)
        do j=1,nvtot
          ivqm = int_mb(i_icvdw+j-1)
          if(ivqm.eq.iqm) then
            dbl_mb(i_rvdw+3*(j-1)) = dbl_mb(i_ctmp+3*(i-1))
            dbl_mb(i_rvdw+3*(j-1)+1) = dbl_mb(i_ctmp+3*(i-1)+1)
            dbl_mb(i_rvdw+3*(j-1)+2) = dbl_mb(i_ctmp+3*(i-1)+2)
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

      subroutine mm_vdw14_qmcoords_load(rtdb)
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
#include "mm_vdw_coords_data.fh"
      integer rtdb      

      integer n,fn
      integer geom
      integer nat
      logical status
      character*30 message
      character*30 pname
      double precision scale

c     local variables
      character*16 t(nqm)
      integer i, j
      integer iqm, ivqm
      integer i_ctmp, h_ctmp

      pname = "mm_vdw14_qmcoords_load"
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
      call dscal(3*nact, scale,dbl_mb(i_ctmp),1)

      do i=1,nact
        iqm = int_mb(i_iact+i-1)
        do j=1,nv14tot
          ivqm = int_mb(i_icvdw14+j-1)
          if(ivqm.eq.iqm) then
            dbl_mb(i_rvdw14+3*(j-1)) = dbl_mb(i_ctmp+3*(i-1))
            dbl_mb(i_rvdw14+3*(j-1)+1) = dbl_mb(i_ctmp+3*(i-1)+1)
            dbl_mb(i_rvdw14+3*(j-1)+2) = dbl_mb(i_ctmp+3*(i-1)+2)
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


      subroutine mm_vdw_coords_test(nvt,c,ic)
      implicit none
#include "mm_vdw_coords_data.fh"
#include "mafdecls.fh"
      integer nvt
      double precision c(nvt)
      integer ic(nvt)

      integer i,j
      do i=1,nvt
         write(6,'(1X,3(1X,F10.6),I5)') 
     $        (c(3*(i-1)+j),j=1,3),
     $        ic(i)
      end do

      end
