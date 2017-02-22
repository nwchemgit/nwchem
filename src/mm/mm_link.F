c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      subroutine mm_links_init(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "geom.fh"
#include "mm_link_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb

      character*32 pname
      integer i, j, k     
      integer lnk_indx, mm_indx      
      logical geom_tag_to_atn
      external geom_tag_to_atn

      pname = "mm_links_init"
c     write(*,*) pname

      if(.not.rtdb_get(rtdb,'mm:qmlink',mt_log,1,qmlink))
     >    qmlink = .true.

c     deallocate all previous allocated arrays just in case
      call mm_links_end()

      call mm_find_links()

      nqml = nqm + nlink
      
      if(.not.ma_alloc_get(mt_dbl,3*nqml,'qmlink coords',h_rqml,i_rqml))
     & call errquit(pname//'Failed to allocate heap',3*nqml,
     &       MA_ERR)

      if(.not.ma_alloc_get(mt_int,nqml,'qmlink indices',h_iqml,i_iqml))
     & call errquit(pname//'Failed to allocate heap',nqml,
     &       MA_ERR)

      if(.not.ma_alloc_get(mt_byte,16*nqml,'qmlink tag',h_tqml,i_tqml))
     & call errquit(pname//'Failed to allocate heap',16*nqml,
     &       MA_ERR)

      if(nlink.gt.0) then
        if(.not.ma_alloc_get(mt_dbl,nlink,'qmlink mm charge',
     &                       h_lnkchg,i_lnkchg))
     &   call errquit(pname//'Failed to allocate heap',nlink,
     &         MA_ERR)
      
        if(.not.ma_alloc_get(mt_int,nlink,'qmlink mm atn',
     &                       h_lnkatn,i_lnkatn))
     &   call errquit(pname//'Failed to allocate heap',nlink,
     &         MA_ERR)
      end if

      do i=1,nqm
        dbl_mb(i_rqml+3*(i-1))   = dbl_mb(i_rqm+3*(i-1))
        dbl_mb(i_rqml+3*(i-1)+1) = dbl_mb(i_rqm+3*(i-1)+1)
        dbl_mb(i_rqml+3*(i-1)+2) = dbl_mb(i_rqm+3*(i-1)+2)
        int_mb(i_iqml+i-1) = int_mb(i_iqm+i-1)
        call mm_coords_tag_set(byte_mb(i_tqm+16*(i-1)),
     &                         byte_mb(i_tqml+16*(i-1)))
      end do

      do i=1,nlink
        lnk_indx = int_mb(i_lb+nlink+i-1) 
        int_mb(i_iqml+nqm+i-1) = lnk_indx
        do j=1,nmm
          mm_indx = int_mb(i_imm+j-1)
          if(mm_indx.eq.lnk_indx) then
            dbl_mb(i_rqml+3*(nqm+i-1))     = dbl_mb(i_rmm+3*(j-1))
            dbl_mb(i_rqml+3*(nqm+i-1)+1)   = dbl_mb(i_rmm+3*(j-1)+1)
            dbl_mb(i_rqml+3*(nqm+i-1)+2)   = dbl_mb(i_rmm+3*(j-1)+2)
            dbl_mb(i_lnkchg+i-1)           = dbl_mb(i_chgmm+j-1)
            call mm_coords_tag_set(byte_mb(i_tmm+16*(j-1)),
     &                             byte_mb(i_tqml+16*(nqm+i-1)))
            exit
          end if
        end do
      end do

      if(nlink.gt.0) then
        if(.not.geom_tag_to_atn(nlink,byte_mb(i_tqml+16*nqm),
     &                          int_mb(i_lnkatn)))
     &  call errquit('mm: Failed to get lnkatn from tag',nlink,MA_ERR)
      end if

      call mm_link_excluded_list_init()
 
      end

      subroutine mm_link_excluded_list_init()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mm_link_data.fh"
#include "mm_vdw_data.fh"
#include "mm_bond_data.fh"

      integer n0
      integer i, j
      integer lnk_indx
      character*30 pname
      integer nnb_i
      integer h_excld, i_excld

      pname = "mm_link_excluded_list_init"

      if(nlink.eq.0) return

      xc_maxatm = nbon + nang + ndih + nvdw14
 
      if(.not.ma_alloc_get(mt_int,nlink,'no. excl atoms',h_nxc,i_nxc))
     >  call errquit('error allocating no. excl atoms',nlink, MA_ERR)

      if(.not.ma_alloc_get(mt_int,nlink*xc_maxatm,'link excluded list',
     >                     h_lxc,i_lxc))
     >  call errquit('error allocating link excluded list',
     >                xc_maxatm*nlink, MA_ERR)

      if(.not.ma_push_get(mt_int,xc_maxatm,'excl_atoms',
     >                    h_excld,i_excld))
     >   call errquit('mm: Failed to push stack excl_atoms',
     >                 xc_maxatm,MA_ERR)

      call ifill(xc_maxatm*nlink,0.0,int_mb(i_lxc),1)
      call ifill(nlink,0.0,int_mb(i_nxc),1)

      do i = 1, nlink
        lnk_indx = int_mb(i_lb+nlink+i-1)
        nnb_i = 0
        call ifill(xc_maxatm,0.0,int_mb(i_excld),1)
        call mm_link_find_excluded_atoms(lnk_indx,
     >                  int_mb(i_excld),
     >                  nnb_i)
        do j=1,nnb_i
          int_mb(i_lxc+(i-1)*xc_maxatm+j-1) = int_mb(i_excld+j-1)
        end do
        int_mb(i_nxc+i-1) = nnb_i
      end do

      if(.not.ma_pop_stack(h_excld))
     >   call errquit('mm: Failed to pop stack excl_atoms',
     >                 xc_maxatm, MA_ERR)

      end

      subroutine mm_link_find_excluded_atoms(indx_i,i_excl_atom,nnb_i)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_bond_data.fh"
#include "mm_vdw_data.fh"

      integer indx_i
      integer i_excl_atom(*)
      integer nnb_i       
 
      integer i, j, k
      integer ii, jj, kk, ll
      character*30 pname

      pname = "mm_link_find_excluded_atoms"
c     write(*,*) pname
      
c     Find excluded atoms for bond
      do i = 1, nbon
        ii = int_mb(i_ibon+i-1)
        jj = int_mb(i_jbon+i-1)
        if(indx_i.eq.ii) then
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
        else if(indx_i.eq.jj) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
        end if
      end do
      
c     Find excluded atoms for angle
      do i = 1, nang
        ii = int_mb(i_iang+i-1)
        jj = int_mb(i_jang+i-1)
        kk = int_mb(i_kang+i-1)
        if(indx_i.eq.ii) then
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(kk,i_excl_atom,nnb_i)
        else if(indx_i.eq.jj) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(kk,i_excl_atom,nnb_i)
        else if(indx_i.eq.kk) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
         end if
       end do

c     Find excluded atoms for dihedral
      do i = 1, ndih
        ii = int_mb(i_idih+i-1)
        jj = int_mb(i_jdih+i-1)
        kk = int_mb(i_kdih+i-1)
        ll = int_mb(i_ldih+i-1)
        if(indx_i.eq.ii) then
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(kk,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(ll,i_excl_atom,nnb_i)
        else if(indx_i.eq.jj) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(kk,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(ll,i_excl_atom,nnb_i)
        else if(indx_i.eq.kk) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(ll,i_excl_atom,nnb_i)
        else if(indx_i.eq.ll) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
          call mm_link_add_excluded_list(kk,i_excl_atom,nnb_i)
        end if
      end do

c     Find excluded atoms for 1-4 VDW
      do i = 1, nvdw14
        ii = int_mb(i_ivdw14+i-1)
        jj = int_mb(i_jvdw14+i-1)
        if(indx_i.eq.ii) then
          call mm_link_add_excluded_list(jj,i_excl_atom,nnb_i)
        else if(indx_i.eq.jj) then
          call mm_link_add_excluded_list(ii,i_excl_atom,nnb_i)
        end if
      end do

      end 
 
      subroutine mm_link_add_excluded_list(indx,excl_list,n)
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
      integer indx
      integer excl_list(*)
      integer n

      integer i, ii
      integer cnt
      character*30 pname
      logical add_indx
  
      pname = "mm_link_add_excluded_list"
c     write(*,*) pname

      add_indx = .true.
       
      if(n.eq.0) then
        n = n + 1 
        excl_list(n) = indx
      else
        cnt = n
        do i=1,cnt
          ii = excl_list(i)
          if(ii.eq.indx) then
            add_indx = .false.
            exit
          end if
        end do
c       if not already in the list, add it
        if(add_indx) then
          n = n + 1
          excl_list(n) = indx
        end if
      end if

      end

      subroutine mm_links_bq_scaled_add_egrad(rtdb,e,n,g)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_vdw_data.fh"
#include "mm_vdw_coords_data.fh"
      integer rtdb
      integer n
      double precision e
      double precision g(3,n)    

      integer i, j, k, m
      integer i1, i2
      integer i_i, i_j
      integer i_lnk, j_lnk, lnk_indx, indx
      integer iqm, imm, cnt
      double precision scee
      double precision r1(3), r2(3)
      double precision dx, dy, dz, dr, dr2, dr3
      double precision gx, gy, gz
      double precision e_el14, eel14
      double precision scale
      double precision z1, z2
      logical i_is_lnk, j_is_lnk
      logical i_is_qm, j_is_qm
    
      character*30 pname

      parameter(scee = 0.833333d0)

      pname = "mm_links_bq_scaled_add_egrad"
c     write(*,*) pname

      if(nlink.eq.0) return

      call util_convert_units("ang","au",scale)

      e_el14 = 0.0d0
      eel14  = 0.0d0
      do i = 1,nvdw14
        i1 = int_mb(i_ivdw14+i-1)
        i2 = int_mb(i_jvdw14+i-1)
c       scee = dbl_mb(i_scee+i-1)

        i_is_lnk = .false.
        j_is_lnk = .false.
        do j = 1, nlink
          lnk_indx = int_mb(i_lb+nlink+j-1)
          if(i1.eq.lnk_indx) then
            i_is_lnk = .true. 
            i_lnk = j
          end if
          if(i2.eq.lnk_indx) then
            j_is_lnk = .true. 
            j_lnk = j
          end if
        end do

        i_is_qm = .false.
        j_is_qm = .false.
        do j = 1, nqm
          iqm = int_mb(i_iqm+j-1)
          if(i1.eq.iqm) then
            i_is_qm = .true. 
          end if
          if(i2.eq.iqm) then
            j_is_qm= .true. 
          end if
        end do

        if((i_is_lnk.and..not.j_is_qm).or.
     >     (j_is_lnk.and..not.i_is_qm)) then
          do k = 1,nv14tot
            indx = int_mb(i_icvdw14+k-1)
            if(i1.eq.indx) i_i = k
            if(i2.eq.indx) i_j = k
          end do

          z1 = dbl_mb(i_chgmm14+i_i-1)
          z2 = dbl_mb(i_chgmm14+i_j-1)

          r1(1) = dbl_mb(i_rvdw14+(i_i-1)*3)*scale
          r1(2) = dbl_mb(i_rvdw14+(i_i-1)*3+1)*scale
          r1(3) = dbl_mb(i_rvdw14+(i_i-1)*3+2)*scale
          r2(1) = dbl_mb(i_rvdw14+(i_j-1)*3)*scale
          r2(2) = dbl_mb(i_rvdw14+(i_j-1)*3+1)*scale
          r2(3) = dbl_mb(i_rvdw14+(i_j-1)*3+2)*scale

          dx = r1(1) - r2(1)
          dy = r1(2) - r2(2)
          dz = r1(3) - r2(3)

          dr2 = dx*dx + dy*dy + dz*dz
          dr  = sqrt(dr2)
          dr3 = dr*dr*dr

          gx   =  z1*z2*dx/dr3*scee
          gy   =  z1*z2*dy/dr3*scee
          gz   =  z1*z2*dz/dr3*scee
          eel14 = z1*z2/dr*scee

          e_el14 = e_el14 + eel14

          if(i_is_lnk) then
            g(1,nqm+i_lnk) = g(1,nqm+i_lnk) - gx
            g(2,nqm+i_lnk) = g(2,nqm+i_lnk) - gy
            g(3,nqm+i_lnk) = g(3,nqm+i_lnk) - gz
          end if
          if(j_is_lnk) then
            g(1,nqm+j_lnk) = g(1,nqm+j_lnk) + gx
            g(2,nqm+j_lnk) = g(2,nqm+j_lnk) + gy
            g(3,nqm+j_lnk) = g(3,nqm+j_lnk) + gz
          end if
        end if
      end do 

      e = e + e_el14

      end

      subroutine mm_links_bq_scaled_add_energy(rtdb,e)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_vdw_data.fh"
#include "mm_vdw_coords_data.fh"
      integer rtdb
      double precision e

      integer i, j, k, m
      integer i1, i2
      integer i_i, i_j
      integer lnk_indx, i_lnk, j_lnk, indx
      integer iqm, imm, cnt
      double precision scee
      double precision r1(3), r2(3)
      double precision dx, dy, dz, dr, dr2
      double precision e_el14, eel14
      double precision scale
      double precision z1, z2
      logical i_is_lnk, j_is_lnk
      logical i_is_qm, j_is_qm
    
      character*30 pname

      parameter(scee = 0.833333d0)

      pname = "mm_links_bq_scaled_add_energy"
c     write(*,*) pname

      if(nlink.eq.0) return

      call util_convert_units("ang","au",scale)

      e_el14 = 0.0d0
      eel14  = 0.0d0
      do i = 1,nvdw14
        i1 = int_mb(i_ivdw14+i-1)
        i2 = int_mb(i_jvdw14+i-1)
c       scee = dbl_mb(i_scee+i-1)

        i_is_lnk = .false.
        j_is_lnk = .false.
        do j = 1, nlink
          lnk_indx = int_mb(i_lb+nlink+j-1)
          if(i1.eq.lnk_indx) then
            i_is_lnk = .true. 
            i_lnk = j
          end if
          if(i2.eq.lnk_indx) then
            j_is_lnk = .true. 
            j_lnk = j
          end if
        end do

        i_is_qm = .false.
        j_is_qm = .false.
        do j = 1, nqm
          iqm = int_mb(i_iqm+j-1)
          if(i1.eq.iqm) then
            i_is_qm = .true. 
          end if
          if(i2.eq.iqm) then
            j_is_qm= .true. 
          end if
        end do

        if((i_is_lnk.and..not.j_is_qm).or.
     >     (j_is_lnk.and..not.i_is_qm)) then
          do k = 1,nv14tot
            indx = int_mb(i_icvdw14+k-1)
            if(i1.eq.indx) i_i = k
            if(i2.eq.indx) i_j = k
          end do

          z1 = dbl_mb(i_chgmm14+i_i-1)
          z2 = dbl_mb(i_chgmm14+i_j-1)

          r1(1) = dbl_mb(i_rvdw14+(i_i-1)*3)*scale
          r1(2) = dbl_mb(i_rvdw14+(i_i-1)*3+1)*scale
          r1(3) = dbl_mb(i_rvdw14+(i_i-1)*3+2)*scale
          r2(1) = dbl_mb(i_rvdw14+(i_j-1)*3)*scale
          r2(2) = dbl_mb(i_rvdw14+(i_j-1)*3+1)*scale
          r2(3) = dbl_mb(i_rvdw14+(i_j-1)*3+2)*scale

          dx = r1(1) - r2(1)
          dy = r1(2) - r2(2)
          dz = r1(3) - r2(3)

          dr2 = dx*dx + dy*dy + dz*dz
          dr  = sqrt(dr2)

          eel14 = z1*z2/dr*scee

          e_el14 = e_el14 + eel14

        end if
      end do 

      e = e + e_el14

      end

      subroutine mm_links_bq_add_egrad(rtdb,e,n,g)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_vdw_data.fh"
#include "mm_bond_data.fh"
      integer rtdb
      integer n
      double precision e
      double precision g(3,n)    

      character*30 pname
      integer i, j, k
      integer nn
      integer indx_i, indx_j, i_excl
      integer nnb_i
      integer h_xcl, i_xcl
      double precision cut, cut2
      double precision dx, dy, dz, dr, dr2, dr3
      double precision gx, gy, gz
      double precision e_lbq, elbq
      double precision scale
      double precision lnk_coord(3)  
      double precision bq_coord(3)  
      double precision zlnk, zbq  
      logical j_xcl

      pname = "mm_links_bq_add_egrad"
c     write(*,*) pname
 
      if(nlink.eq.0) return

c     cutoff in angstrom
      if(.not.rtdb_get(rtdb,'mm:cutoff',mt_dbl,1,cut))
     >   cut = 999.0d0
c     write(*,*) 'cutoff = ', cut, 'A'
c     
      call util_convert_units("ang","au",scale)

      cut = cut*scale
      cut2 = cut*cut

      e_lbq = 0.0d0
      elbq = 0.0d0
      do i=1,nlink
        indx_i = int_mb(i_lb+nlink+i-1)  
        nnb_i = int_mb(i_nxc+i-1)
        do j=1,nmm
          indx_j = int_mb(i_imm+j-1)
          j_xcl = .false.
          do k=1,nnb_i
            i_excl = int_mb(i_lxc+(i-1)*xc_maxatm+k-1)
            if(indx_j.eq.i_excl) then
              j_xcl = .true.
              exit
            end if
          end do
          if(indx_i.ne.indx_j.and..not.j_xcl) then
            lnk_coord(1) = dbl_mb(i_rqml+(nqm+i-1)*3)*scale
            lnk_coord(2) = dbl_mb(i_rqml+(nqm+i-1)*3+1)*scale 
            lnk_coord(3) = dbl_mb(i_rqml+(nqm+i-1)*3+2)*scale
            bq_coord(1)  = dbl_mb(i_rmm+(j-1)*3)*scale 
            bq_coord(2)  = dbl_mb(i_rmm+(j-1)*3+1)*scale 
            bq_coord(3)  = dbl_mb(i_rmm+(j-1)*3+2)*scale 

            dx = lnk_coord(1)-bq_coord(1)
            dy = lnk_coord(2)-bq_coord(2)
            dz = lnk_coord(3)-bq_coord(3)

            dr2 = dx*dx + dy*dy + dz*dz
c           
            if(dr2 < cut2) then
              zlnk = dbl_mb(i_lnkchg+i-1)
              zbq  = dbl_mb(i_chgmm+j-1)
              dr   = sqrt(dr2)
              dr3  = dr*dr*dr

              gx   =  zlnk*zbq*dx/dr3
              gy   =  zlnk*zbq*dy/dr3
              gz   =  zlnk*zbq*dz/dr3
              elbq = zlnk*zbq/dr   
          
              g(1,nqm+i) = g(1,nqm+i) - gx
              g(2,nqm+i) = g(2,nqm+i) - gy
              g(3,nqm+i) = g(3,nqm+i) - gz
              e_lbq = e_lbq + elbq
            end if
          end if           
        end do
      end do

      e = e + e_lbq

      end

      subroutine mm_links_bq_add_energy(rtdb,e)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_bond_data.fh"
#include "mm_vdw_data.fh"
      integer rtdb
      double precision e

      character*30 pname
      integer i, j, k
      integer nn
      integer indx_i, indx_j, i_excl
      integer h_xcl, i_xcl
      integer nnb_i
      double precision cut, cut2
      double precision dx, dy, dz, dr, dr2
      double precision e_lbq, elbq
      double precision scale
      double precision lnk_coord(3)  
      double precision bq_coord(3)  
      double precision zlnk, zbq  
      logical j_xcl

      pname = "mm_links_bq_add_energy"
c     write(*,*) pname
 
      if(nlink.eq.0) return

c     cutoff in angstrom
      if(.not.rtdb_get(rtdb,'mm:cutoff',mt_dbl,1,cut))
     >   cut = 999.0d0
c     write(*,*) 'cutoff = ', cut, 'A'
      
      call util_convert_units("angstrom","au",scale)

      cut = cut*scale
      cut2 = cut*cut

      e_lbq = 0.0d0
      elbq = 0.0d0
      do i=1,nlink
        indx_i = int_mb(i_lb+nlink+i-1)  
        nnb_i = int_mb(i_nxc+i-1)
        do j=1,nmm
          indx_j = int_mb(i_imm+j-1)
          j_xcl = .false.
          do k=1,nnb_i
            i_excl = int_mb(i_lxc+(i-1)*xc_maxatm+k-1)
            if(indx_j.eq.i_excl) then
              j_xcl = .true.
              exit
            end if
          end do
          if(indx_i.ne.indx_j.and..not.j_xcl) then
            lnk_coord(1) = dbl_mb(i_rqml+(nqm+i-1)*3)*scale
            lnk_coord(2) = dbl_mb(i_rqml+(nqm+i-1)*3+1)*scale 
            lnk_coord(3) = dbl_mb(i_rqml+(nqm+i-1)*3+2)*scale
            bq_coord(1)  = dbl_mb(i_rmm+(j-1)*3)*scale 
            bq_coord(2)  = dbl_mb(i_rmm+(j-1)*3+1)*scale 
            bq_coord(3)  = dbl_mb(i_rmm+(j-1)*3+2)*scale 

            dx = lnk_coord(1)-bq_coord(1)
            dy = lnk_coord(2)-bq_coord(2)
            dz = lnk_coord(3)-bq_coord(3)

            dr2 = dx*dx + dy*dy + dz*dz
            
            if(dr2 < cut2) then
              zlnk = dbl_mb(i_lnkchg+i-1)
              zbq  = dbl_mb(i_chgmm+j-1)
              dr   = sqrt(dr2)
              elbq = zlnk*zbq/dr   
              e_lbq = e_lbq + elbq
            end if           
          end if
        end do
      end do

      e = e + e_lbq

      end

      subroutine mm_links_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_link_data.fh"

      character*30 pname
      pname = "mm_links_end"

      if(nqml.gt.0) then
          if (.not.ma_free_heap(h_rqml))     goto 911
          if (.not.ma_free_heap(h_iqml))     goto 911
          if (.not.ma_free_heap(h_tqml))     goto 911
          nqml = 0
      end if
c 
      if(nlink.gt.0) then
          if (.not.ma_free_heap(h_lb))       goto 911
          if (.not.ma_free_heap(h_lnkchg))   goto 911
          if (.not.ma_free_heap(h_lnkatn))   goto 911
          if (.not.ma_free_heap(h_nxc))      goto 911
          if (.not.ma_free_heap(h_lxc))      goto 911
          nlink = 0
      end if

      return
911   call errquit("error "//trim(pname),0,-1)

      end

      subroutine mm_find_links()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_link_data.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"

      integer i,j
      integer indx_i, indx_j
      integer ind_qm, ind_l
      integer indx_crd
      integer cnt
      character*16 tag
      character*32 pname
  
      pname = "qmmm_init_hlinks"
c     write(*,*) pname

      nlink = 0
      do i = 1, nbon
        indx_i = int_mb(i_ibon+i-1)
        indx_j = int_mb(i_jbon+i-1)
        do j=1,nqm
          indx_crd = int_mb(i_iqm+j-1)
          if(indx_i.eq.indx_crd.or.indx_j.eq.indx_crd)
     &        nlink = nlink + 1
        end do
      end do

      if(nlink.gt.0) then
        if(.not.ma_alloc_get(mt_int,2*nlink,'link bond index',
     &                       h_lb,i_lb))
     &   call errquit(pname//'Failed to allocate heap',nlink,
     &        MA_ERR)
      end if

      cnt = 0
      do i=1,nbon
        indx_i = int_mb(i_ibon+i-1)
        indx_j = int_mb(i_jbon+i-1)
        do j=1,nqm
          indx_crd = int_mb(i_iqm+j-1)
          if(indx_i.eq.indx_crd) then 
            cnt = cnt + 1
            int_mb(i_lb+cnt-1) = indx_i
            int_mb(i_lb+nlink+cnt-1) = indx_j
          else if(indx_j.eq.indx_crd) then
            cnt = cnt + 1
            int_mb(i_lb+cnt-1) = indx_j
            int_mb(i_lb+nlink+cnt-1) = indx_i
          end if
        end do
      end do
  
      end 

      subroutine mm_links_adjust(n, ifg, atn, t, c, q)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "mm_link_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"

      integer n
      integer ifg(n)
      integer atn(n)
      character*16 t(n)
      double precision c(3, n) 
      double precision q(n) 
    
c     local variables     
      integer i,j,k,m
      character*32 pname
      double precision a
      double precision rq,rl,rh
      integer icl,icq,icrd

      pname = "mm_links_adjust"     
c     write(*,*) pname 
 
      a = 0.709d0
      do i=1,nlink
        icq = int_mb(i_lb+i-1)
        icl = int_mb(i_lb+nlink+i-1)
        do k=1,3
          do m=1,nqml
            icrd = int_mb(i_iqml+m-1)
            if(icq.eq.icrd) rq = dbl_mb(i_rqml+3*(m-1)+k-1)
            if(icl.eq.icrd) rl = dbl_mb(i_rqml+3*(m-1)+k-1)
          end do
          rh = rq - (rq-rl)*a
          c(k,nqm+i) = rh
        end do
        q(nqm+i) = 1.0d0
        atn(nqm+i) = int_mb(i_lnkatn+i-1)
        t(nqm+i) = "H_L"
      end do

      end

      subroutine mm_link_ebq_add_grad(rtdb,n,g)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
#include "mm_bq_data.fh"
      integer rtdb
      integer n
      double precision g(3,n)
      
      integer bq_ncent
      integer i, j
      integer h_c, i_c
      integer h_g, i_g
      character*30 pname
      double precision scale
 
      pname = "mm_link_ebq_add_grad"
c     write(*,*) pname

      if(.not.qmlink.or.nlink.eq.0) return 

      bq_ncent = nlink 

      if(.not.ma_push_get(mt_dbl,3*bq_ncent,'mmlink bq c',h_c,i_c))
     & call errquit(pname//': Failed to allocate memory for c',
     & 3*bq_ncent, MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*bq_ncent,'mmlink bq g',h_g,i_g))
     & call errquit(pname//': Failed to allocate memory for g',
     & 3*bq_ncent, MA_ERR)

      call dfill(3*bq_ncent,0.0d0,dbl_mb(i_g),1)
#ifdef USE_BLASMOD
      call dcopy( 3*nlink,dbl_mb(i_c),1,dbl_mb(i_rqml+nqm*3),1)
#else
      do i=1,nlink
        dbl_mb(i_c+(i-1)*3)   = dbl_mb(i_rqml+(nqm+i-1)*3)
        dbl_mb(i_c+(i-1)*3+1) = dbl_mb(i_rqml+(nqm+i-1)*3+1)
        dbl_mb(i_c+(i-1)*3+2) = dbl_mb(i_rqml+(nqm+i-1)*3+2)
      end do
#endif

      call util_convert_units("angstrom","au",scale)
      call dscal(3*bq_ncent, scale,dbl_mb(i_c),1)

      call qmmm_bq_forces(rtdb,bq_ncent,
     &                     dbl_mb(i_c),
     &                     dbl_mb(i_lnkchg),
     &                     dbl_mb(i_g))
#ifdef USE_BLASMOD
      call daxpy(3*nlink,1d0,dbl_mb(i_g),1,g(1,nqm+1))
#else
CDIR$ NOVECTOR
      do i=1,nlink
        g(1,nqm+i) = g(1,nqm+i) + dbl_mb(i_g+(i-1)*3)
        g(2,nqm+i) = g(2,nqm+i) + dbl_mb(i_g+(i-1)*3+1)
        g(3,nqm+i) = g(3,nqm+i) + dbl_mb(i_g+(i-1)*3+2)
      end do
#endif

      if(.not.ma_pop_stack(h_g))
     & call errquit('mm: Failed to deallocate stack g_all',3*bq_ncent,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('mm: Failed to deallocate stack c_all',3*bq_ncent,
     &       MA_ERR)

      end

      subroutine mm_link_update_bq_coords(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "mm_link_data.fh"
#include "mm_coords_data.fh"
#include "mm_bq_data.fh"
      integer rtdb     
 
      logical  bq_ncenter
      external bq_ncenter
c
      logical  bq_index_coord
      external bq_index_coord
c
      logical  bq_index_charge
      external bq_index_charge
c
      logical  bq_destroy
      external bq_destroy
c
      logical  bq_get_active
      external bq_get_active

      character*32 pname
      integer bq_handle
      integer bq_ncent
      integer ncent
      integer i,j
      integer indx_lnk, indx_bq
      integer i_cbq,i_qbq
      double precision scale
      double precision coord(3)

      pname = "mm_link_update_bq_coords"
c     write(*,*) pname

      if(.not.qmlink.or.nlink.eq.0) return
c   
c     get bq information
c     ------------------
      if(.not.bq_get_active(bq_handle))
     >   call errquit(pname//':no active bq handle',0,0)
      if(.not.bq_ncenter(bq_handle,bq_ncent))
     >   call errquit(pname//':no bq centers',0,0)
      if(.not.bq_index_coord(bq_handle,i_cbq))
     >   call errquit(pname//':no bq coords',0,0)

      call util_convert_units("angstrom","au",scale)

      do i=1,nlink
        indx_lnk = int_mb(i_lb+nlink+i-1)
        do j=1, bq_ncent
          indx_bq = int_mb(i_ibq+j-1)
          if(indx_bq.eq.indx_lnk) then
            coord(1) = dbl_mb(i_rqml+(nqm+i-1)*3)*scale
            coord(2) = dbl_mb(i_rqml+(nqm+i-1)*3+1)*scale
            coord(3) = dbl_mb(i_rqml+(nqm+i-1)*3+2)*scale
            dbl_mb(i_cbq+(j-1)*3)   = coord(1)
            dbl_mb(i_cbq+(j-1)*3+1) = coord(2)
            dbl_mb(i_cbq+(j-1)*3+2) = coord(3)
            exit
          end if
        end do
      end do
     
      end 

      subroutine mm_links_adjust_forces(n, ai, g)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "mm_link_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer n
      integer ai(n)
      double precision g(3,n)

c     local variables
      integer i,j,k,m
      character*32 pname
      double precision a
      double precision gq,gl,gh
      integer il, iq      

      pname = "mm_links_adjust_forces"
c     write(*,*) pname

      a = 0.709d0
 
      do i=1,nlink
        il = 0
        iq = 0
        do j=1,nqml
          if(ai(j).eq.int_mb(i_lb+i-1)) iq = j
          if(ai(j).eq.int_mb(i_lb+nlink+i-1)) il = j
        end do
        if(iq.ne.0 .and. il.ne.0) then
          do k=1,3
             gq  = g(k,iq)
             gl  = g(k,il)
             g(k,il) = a*gl
             g(k,iq) = gq+(1.0d0-a)*gl
          end do
        end if
      end do

      end 
