c $Id$

*     *********************************
*     *                               *
*     *         nwpw_hartree_init     *
*     *                               *
*     *********************************
      subroutine nwpw_hartree_init(nkatm,
     >                        nprj,nbasis,psp_type,lmax,
     >                        nprj_max0,l_prj,m_prj,b_prj,
     <                        hartree_tag)
      implicit none
      integer nkatm
      integer nprj(*),nbasis(*),psp_type(*),lmax(*)
      integer nprj_max0
      integer l_prj(nprj_max0,*),m_prj(nprj_max0,*),b_prj(nprj_max0,*)
      integer hartree_tag


#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_hartree.fh"

c     **** local variables ****
      logical ok
      integer ii,ia,nsize
      integer ic
      integer l,m,li,mi,lj,mj,bi,bj,iprj,jprj,lm
      integer li1,mi1,lj1,mj1,bi1,bj1,iprj1,jprj1
      integer i_p,i_t
      integer i_tlm,i_plm,i_rlm
      integer matr_ptr
      integer nb
      real*8  tmp_theta,cs_theta,tmp_phi,tmp_gaunt

c     **** external functions ****
      integer  psi_data_get_chnk
      external psi_data_get_chnk
      real*8   rtheta_lm,drtheta_lm,nwpw_gaunt
      external rtheta_lm,drtheta_lm,nwpw_gaunt

      call nwpw_timing_start(4)
      ok = .true.

*      ***** set up index arrays for nwpw_hartree_solve2 *****
      ok = BA_alloc_get(mt_int,nkatm,"nindx_hartree",
     >                  nindx_hartree(2),nindx_hartree(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nkatm,"shift_hartree",
     >                  shift_hartree(2),shift_hartree(1))
      if (.not.ok)
     >call errquit("nwpw_hartree_init:error allocating work arrays",0,0)

      nsize = 0
      do ia=1,nkatm
         nb = nbasis(ia)
         matr_ptr = psi_data_get_chnk(hartree_tag(ia))
         int_mb(shift_hartree(1)+ia-1) = nsize
         do jprj = 1,nprj(ia)
            lj = l_prj(jprj,ia)
            mj = m_prj(jprj,ia)
            bj = b_prj(jprj,ia)
            do iprj = 1,nprj(ia)
               li = l_prj(iprj,ia)
               mi = m_prj(iprj,ia)
               bi = b_prj(iprj,ia)
               do jprj1 = 1,nprj(ia)
                  lj1 = l_prj(jprj1,ia)
                  mj1 = m_prj(jprj1,ia)
                  bj1 = b_prj(jprj1,ia)
                  do iprj1 = 1,nprj(ia)
                     li1 = l_prj(iprj1,ia)
                     mi1 = m_prj(iprj1,ia)
                     bi1 = b_prj(iprj1,ia)
                     do l=0,2*lmax(ia)
                        do m=-l,l
                           tmp_gaunt = 
     >                        nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *nwpw_gaunt(.false.,l,m,li1,mi1,lj1,mj1)
     >                       *dbl_mb(matr_ptr
     >                              + (bi-1) 
     >                              + (bj-1)*nb
     >                              + (bi1-1)*nb*nb
     >                              + (bj1-1)*nb*nb*nb
     >                              + l*nb*nb*nb*nb)
                           if (dabs(tmp_gaunt).gt.1.0d-11) then
                              nsize = nsize + 1
                           end if
                        end do
                     end do

                  end do
               end do

            end do
         end do
         int_mb(nindx_hartree(1)+ia-1)= nsize
     >                                - int_mb(shift_hartree(1)+ia-1)
      end do

      ok = BA_alloc_get(mt_int,nsize,"iprj_hartree",
     >                 iprj_hartree(2),iprj_hartree(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"jprj_hartree",
     >                 jprj_hartree(2),jprj_hartree(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"iprj1_hartree",
     >                 iprj1_hartree(2),iprj1_hartree(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"jprj1_hartree",
     >                 jprj1_hartree(2),jprj1_hartree(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"coeff_hartree",
     >                 coeff_hartree(2),coeff_hartree(1))
      if (.not.ok)
     >call errquit("nwpw_hartree_init:error allocating work arrays",0,0)


      nsize = 0
      do ia=1,nkatm
         nb = nbasis(ia)
         matr_ptr = psi_data_get_chnk(hartree_tag(ia))
         do jprj = 1,nprj(ia)
            lj = l_prj(jprj,ia)
            mj = m_prj(jprj,ia)
            bj = b_prj(jprj,ia)
            do iprj = 1,nprj(ia)
               li = l_prj(iprj,ia)
               mi = m_prj(iprj,ia)
               bi = b_prj(iprj,ia)
               do jprj1 = 1,nprj(ia)
                  lj1 = l_prj(jprj1,ia)
                  mj1 = m_prj(jprj1,ia)
                  bj1 = b_prj(jprj1,ia)
                  do iprj1 = 1,nprj(ia)
                     li1 = l_prj(iprj1,ia)
                     mi1 = m_prj(iprj1,ia)
                     bi1 = b_prj(iprj1,ia)
                     do l=0,2*lmax(ia)
                        do m=-l,l
                           tmp_gaunt =
     >                        nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *nwpw_gaunt(.false.,l,m,li1,mi1,lj1,mj1)
     >                       *dbl_mb(matr_ptr
     >                              + (bi-1)
     >                              + (bj-1)*nb
     >                              + (bi1-1)*nb*nb
     >                              + (bj1-1)*nb*nb*nb
     >                              + l*nb*nb*nb*nb)
                           if (dabs(tmp_gaunt).gt.1.0d-11) then
                              int_mb(iprj_hartree(1)+nsize)  = iprj
                              int_mb(jprj_hartree(1)+nsize)  = jprj
                              int_mb(iprj1_hartree(1)+nsize) = iprj1
                              int_mb(jprj1_hartree(1)+nsize) = jprj1
                              nsize = nsize + 1
                           end if
                        end do
                     end do

                  end do
               end do

            end do
         end do
      end do

      call nwpw_timing_end(4)
      return
      end 

*     ********************************************
*     *                                          *
*     *             nwpw_hartree_solve           *
*     *                                          *
*     ********************************************
      subroutine nwpw_hartree_solve(ia,
     >                              ispin,ne,nprj,sw1,sw2)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_hartree.fh"

      integer shift

      call nwpw_timing_start(4)
      call nwpw_timing_start(22)

*     **** hartree potential  non-local matrix elements ****
      shift = int_mb(shift_hartree(1)+ia-1)
      call nwpw_hartree_sw1tosw2(ispin,ne,nprj,sw1,sw2,
     >                           int_mb(nindx_hartree(1)+ia-1),
     >                           int_mb(iprj_hartree(1)+shift),
     >                           int_mb(jprj_hartree(1)+shift),
     >                           int_mb(iprj1_hartree(1)+shift),
     >                           int_mb(jprj1_hartree(1)+shift),
     >                           dbl_mb(coeff_hartree(1)+shift))

      call nwpw_timing_end(4)
      call nwpw_timing_end(22)

      return
      end 


*     *****************************************
*     *                                       *
*     *             nwpw_hartree_end          *
*     *                                       *
*     *****************************************
      subroutine nwpw_xc_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_hartree.fh"

      logical ok
      
      call nwpw_timing_start(4)
      ok = .true.
      ok = ok.and.BA_free_heap(coeff_hartree(2))
      ok = ok.and.BA_free_heap(jprj_hartree(2))
      ok = ok.and.BA_free_heap(iprj_hartree(2))
      ok = ok.and.BA_free_heap(jprj1_hartree(2))
      ok = ok.and.BA_free_heap(iprj1_hartree(2))
      ok = ok.and.BA_free_heap(shift_hartree(2))
      ok = ok.and.BA_free_heap(nindx_hartree(2))

      if (.not.ok)
     > call errquit("nwpw_hartree_end: error freeing heap",0,MA_ERR)

      call nwpw_timing_end(4)
      return
      end 

c     *********************************************
c     *                                           *
c     *           nwpw_hartree_sw1tosw2           *
c     *                                           *
c     *********************************************
      subroutine nwpw_hartree_sw1tosw2(ispin,ne,nprj,sw1,sw2,
     >                                 nindx,
     >                                 iprj_hartree, jprj_hartree,
     >                                 iprj1_hartree,jprj1_hartree,
     >                                 coeff_hartree)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer n1dgrid,nbasis,lmax2
      integer nindx,iprj_hartree(*),jprj_hartree(*)
      integer iprj1_hartree(*),jprj1_hartree(*)
      real*8  coeff_hartree(*)

      integer n,i,iprj,jprj,iprj1,jprj1
      real*8  coeff,w,scal

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      scal = 1.0d0/dsqrt(lattice_omega())
      scal = scal/lattice_omega()

*     ***init to zero***
      do i=1,nindx
         iprj  = iprj_hartree(i)
         jprj  = jprj_hartree(i)
         iprj1 = iprj1_hartree(i)
         jprj1 = jprj1_hartree(i)
         coeff = coeff_hartree(i)*scal
         w = 0.0d0
         do n=1,ne(1)+ne(2)
            w = w + sw1(n,iprj1)*sw1(n,jprj1)
         end do
         do n=1,ne(1)+ne(2)
            sw2(n,iprj) = sw2(n,iprj) + coeff*w*sw1(n,jprj)
         end do 
      end do
      call nwpw_timing_end(21)
      return
      end 
