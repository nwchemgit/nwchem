c
c $Id: paw_pot_hartree.F,v 1.1 2003-02-13 20:23:04 marat Exp $
c

!**************************************************
!
!	name:
!
!	purpose 
!
!	created:	7/9/2002
!**************************************************
      subroutine  init_paw_pot_hartree()
      implicit none

#include "mafdecls.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_ma.fh"
#include "paw_geom.fh"

      !*** local variables ***
      logical ok
      integer in,ia
      integer paw_pot_hartree_size

  
c     !*** calculate the size and indexing ***
      ok = my_alloc(mt_int,ion_nion(),"i_paw_pot_hartree",
     >                                     i_paw_pot_hartree)
      if (.not.ok) call errquit('init_paw_pot_hartree: alloc heap',1)
      paw_pot_hartree_size = 0
      do in=1,ion_nion()
         int_mb(i_paw_pot_hartree(1) + in - 1) = paw_pot_hartree_size
         ia = ion_katm(in)
         paw_pot_hartree_size = paw_pot_hartree_size
     >                + paw_proj_nbasis(ia)**2*ion_natm(ia)
      end do


      ok = my_alloc(mt_dcpl,paw_pot_hartree_size,
     >               "paw_pot_hartree",paw_pot_hartree)
      if (.not.ok)
     > call errquit("init_paw_pot_hartree: error allocating heap",1)
  
      write(*,*) "initialized paw_hartree_pot"
      end !subroutine init_paw_pot_hartree

!**************************************************
!
!	name:
!
!	purpose 
!
!	created:	7/9/2002
!**************************************************
      subroutine  end_paw_pot_hartree()
      implicit none

#include "paw_basis.fh"
#include "mafdecls.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_pot_hartree)
      ok = ok.and.my_dealloc(paw_pot_hartree)
      if (.not.ok)
     > call errquit("end_paw_pot_hartree: error freeing heap",0)

      return
      end !subroutine end_paw_pot_hartree

!**************************************************
!
!	name:
!
!	purpose 
!
!	created:	7/9/2002
!************************************************
      subroutine paw_pot_hartree_solve()
      
      implicit none

#include "mafdecls.fh"
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
      integer in,ia
      integer i_w,i_w1,i_w0
      integer i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer ip,jp,lip,ljp,mip,mjp
      integer nb,nb2,nb3,nb4
      integer l,m
      integer mult_l
      integer i_hp,i_hp0
      integer nilm,njlm,niplm,njplm
      complex tmp_hartree_pot

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
        IF( l.LE.(li+lj) .AND. l.GE.ABS(li-lj) .and. l.LE.(lip+ljp)
     >                            .AND. l.GE.ABS(lip-ljp) ) THEN

        do m=-l,l
          do mi=-li,li
          do mj=-lj,lj
          if ( m.eq.(mi-mj)) then
            i_w1 = 0
            do mip=-lip,lip
            do mjp=-ljp,ljp
              if (m.eq.(mjp-mip)) then
                 i_w   = i_w0-1+(njplm+ljp+mjp+1)+
     >                   (niplm+lip+mip)*paw_proj_nbasis(in)

                 i_hp  = i_hp0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(in)

                 i_mtr = i_mtr0-1+i+(j-1)*nb+(ip-1)*nb2+(jp-1)*nb3+l*nb4

                  tmp_hartree_pot = dbl_mb(paw_hartree_matrix(1)+i_mtr)*
     >                  gen_gaunt_coeff(l,m,li,mi,lj,mj)*
     >                  gen_gaunt_coeff(l,m,ljp,mjp,lip,mip)*
     >                  dcpl_mb(paw_ovlp_w(1)+i_w)

                  dcpl_mb(paw_pot_hartree(1)+i_hp) = 
     >                  dcpl_mb(paw_pot_hartree(1)+i_hp)+
     >                  tmp_hartree_pot

c      write(*,*) "hpot:i,mi,j,mj",i,mi,j,mj
c      write(*,*) "hpot:ip,mip,jp,mjp",ip,mip,jp,mjp
c      write(*,*) "hpot:matrix",dbl_mb(paw_hartree_matrix(1)+i_mtr)
c      write(*,*) "hpot:g",gen_gaunt_coeff(l,m,li,mi,lj,mj)*
c     >                     gen_gaunt_coeff(l,m,ljp,mjp,lip,mip)
c      write(*,*) "hpot:ovlp",dcpl_mb(paw_ovlp_w(1)+i_w)
c      write(*,*) "hpot:tmp_pot", tmp_hartree_pot
c      write(*,*) "hpot:pot",dcpl_mb(paw_pot_hartree(1)+i_hp),i_hp

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
c        write(*,*) "hpot:i,mi,j,mj",i,mi,j,mj,
c     >              dcpl_mb(paw_pot_hartree(1)+i_hp)
c        write(*,*) "hpot:pot",dcpl_mb(paw_pot_hartree(1)+i_hp)
        njlm = njlm + 2*lj+1
        end do !j
        nilm = nilm + 2*li+1
        end do !i  
        end do !l
      end do !in
 
*     *** print hartree potential ***
      write(*,*) "hartree pot"
      do in=1,ion_nion()
        ia = ion_katm(in)
        i_hp0 = int_mb(i_paw_pot_hartree(1) + in - 1)
        nb = paw_basis_nbasis(ia)
        nilm = 0
        do i=1,nb
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
            njlm = 0
            do j=1,nb
              lj = paw_basis_orb_l(j,ia)
                do mj=-lj,lj
                 i_hp  = i_hp0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(in)

                  write(*,*) i,mi,j,mj,dcpl_mb(paw_pot_hartree(1)+i_hp),
     >                     i_hp
                end do
                njlm=njlm + 2*lj+1
             end do
           end do   
           nilm = nilm + 2*li+1
         end do
      end do
       
      return
      end subroutine
