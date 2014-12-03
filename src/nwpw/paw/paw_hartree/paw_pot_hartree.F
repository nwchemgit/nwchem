c
c $Id$
c

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_pot_hartree()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_ma.fh"
#include "paw_geom.fh"

      !*** local variables ***
      logical ok
      integer in,ia

      ! integer paw_pot_hartree_size

  
c     !*** calculate the size and indexing ***
      ok = my_alloc(mt_int,ion_nion(),"i_paw_pot_hartree",
     >                                     i_paw_pot_hartree)
      if (.not.ok) 
     > call errquit('init_paw_pot_hartree:out of heap memory',0,0)
      paw_pot_hartree_size = 0
      do in=1,ion_nion()
         int_mb(i_paw_pot_hartree(1) + in - 1) = paw_pot_hartree_size
         ia = ion_katm(in)
         paw_pot_hartree_size = paw_pot_hartree_size
     >                        + paw_proj_nbasis(ia)**2
      end do


      ok = my_alloc(mt_dcpl,paw_pot_hartree_size,
     >               "paw_pot_hartree",paw_pot_hartree)
      if (.not.ok)
     > call errquit("init_paw_pot_hartree:out of heap memory",0,1)
  
      
      end !subroutine init_paw_pot_hartree

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  end_paw_pot_hartree()
      implicit none

#include "paw_basis.fh"
#include "bafdecls.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_pot_hartree)
      ok = ok.and.my_dealloc(paw_pot_hartree)
      if (.not.ok)
     > call errquit("end_paw_pot_hartree: error freeing heap",0,0)

      return
      end !subroutine end_paw_pot_hartree

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!************************************************
      subroutine paw_pot_hartree_solve()
      
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_kin_data.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_matrix_hartree_data.fh"

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      double precision tmp_gaunt
      integer in,ia
      integer i_w,i_w0
      integer i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer ip,jp,lip,ljp,mip,mjp
      integer nb,nb2,nb3,nb4
      integer l,m
      integer mult_l
      integer i_hp,i_hp0
      integer nilm,njlm,niplm,njplm
      complex*16 tmp_hartree_pot

      
      !*** zero out paw_pot_hartree ***
      call dcopy(2*paw_pot_hartree_size,
     >           0.0d0,0,
     >           dcpl_mb(paw_pot_hartree(1)),1)

      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        nb2 = nb*nb
        nb3 = nb2*nb
        nb4 = nb3*nb
        mult_l = paw_basis_mult_l(ia)
        i_mtr0 = int_mb(i_paw_hartree_matrix(1) + ia - 1)
        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
        i_hp0 = int_mb(i_paw_pot_hartree(1) + in - 1)
        do l=0,mult_l

        nilm = 0
        do i=1,nb
        li = paw_basis_orb_l(i,ia)

        njlm = 0
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)

        niplm = 0
        do ip=1,nb
        lip = paw_basis_orb_l(ip,ia)

        njplm = 0
        do jp=1,nb
        ljp = paw_basis_orb_l(jp,ia)

        if( (l.le.(li+lj))    .and. 
     >      (l.ge.ABS(li-lj)) .and. 
     >      (l.le.(lip+ljp))  .and. 
     >      (l.GE.ABS(lip-ljp)) ) then

        do m=-l,l
          do mi=-li,li
          do mj=-lj,lj
          if ( m.eq.(mi-mj)) then

            do mip=-lip,lip
            do mjp=-ljp,ljp
              if (m.eq.(mjp-mip)) then

c                !**** indicies in wrong order ****
                i_w   = i_w0-1+(njplm+ljp+mjp+1)+
     >                   (niplm+lip+mip)*paw_proj_nbasis(ia)
                i_hp  = i_hp0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
c                !**** indicies in wrong order ****

c                 i_w   = i_w0-1+(niplm+lip+mip+1)+
c     >                   (njplm+ljp+mjp)*paw_proj_nbasis(ia)  ! overlap weights index
c                 i_hp  = i_hp0-1+(nilm+li+mi+1)+
c     >                   (njlm+lj+mj)*paw_proj_nbasis(ia)   ! G(alpha,beta) index


                 i_mtr = i_mtr0-1+i+(j-1)*nb+(ip-1)*nb2+(jp-1)*nb3+l*nb4

                 tmp_gaunt = gen_gaunt_coeff(l,m,li,mi,lj,mj)*
     >                  gen_gaunt_coeff(l,m,ljp,mjp,lip,mip)

                  tmp_hartree_pot = dbl_mb(paw_hartree_matrix(1)+i_mtr)*
     >                  dcpl_mb(paw_ovlp_w(1)+i_w)*tmp_gaunt

                  dcpl_mb(paw_pot_hartree(1)+i_hp) = 
     >                  dcpl_mb(paw_pot_hartree(1)+i_hp)+
     >                  tmp_hartree_pot

!           write(21,*) "i,mi,j,mj,in,ip,mip,jp,mjp",
!     >                i,mi,j,mj,in,ip,mip,jp,mjp
!           write(21,'(A10,E16.8)') "gaunt=",tmp_gaunt
!                              
!      write(21,'(A10,E16.8)')"matr=",dbl_mb(paw_hartree_matrix(1)+i_mtr)
!      write(21,'(A10,2E16.8)')"w=",dcpl_mb(paw_ovlp_w(1)+i_w)
!      write(21,'(A10,2E16.8)')"tmp_hp=",tmp_hartree_pot
!      write(21,'(A10,2E16.8)')"hp=",dcpl_mb(paw_pot_hartree(1)+i_hp)
           
 
              end if
            end do !mjp
            end do !mip
          end if
          end do !mi
          end do !mj
        end do !m
        end if

        njplm = njplm + 2*ljp+1
        end do !jp

        niplm = niplm + 2*lip+1
        end do !ip

        njlm = njlm + 2*lj+1
        end do !j

        nilm = nilm + 2*li+1
        end do !i  
        end do !l
      end do !in
 
!*     *** print hartree potential ***
!      write(26,*)paw_pot_hartree(3)
!      do in=1,ion_nion()
!        ia = ion_katm(in)
!        i_hp0 = int_mb(i_paw_pot_hartree(1) + in - 1)
!        nb = paw_basis_nbasis(ia)
!        nilm = 0
!        do i=1,nb
!          li = paw_basis_orb_l(i,ia)
!          do mi=-li,li
!            njlm = 0
!            do j=1,nb
!              lj = paw_basis_orb_l(j,ia)
!                do mj=-lj,lj
!                 i_hp  = i_hp0-1+(njlm+lj+mj+1)+
!     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
!
!                  write(26,*) i,mi,j,mj,in,
!     >                     dcpl_mb(paw_pot_hartree(1)+i_hp)
!                end do
!                njlm=njlm + 2*lj+1
!             end do
!           end do   
!           nilm = nilm + 2*li+1
!         end do
!      end do
       
      return
      end



      subroutine paw_hartree_pot_ptr(ptr)
      implicit none
      integer ptr

#include "paw_pot_hartree_data.fh"

      ptr = paw_pot_hartree(1)

      return
      end



