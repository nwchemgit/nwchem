*
* $Id$
*

*     ***********************************************************
*     *															*
*     *   		   		lcao_diis_dn library					*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*
*

*     ***********************************
*     *									*
*     *	       lcao_diis_dn_init		*	
*     *									*
*     ***********************************
      subroutine lcao_diis_dn_init(maxsize0)
      implicit none
      integer maxsize0

#include "bafdecls.fh"
#include "lcao_diis_dn.fh"


      maxsize = maxsize0

*     **** initialize dnlist ****
      call lcao_dnlist_init(maxsize)

      return
      end

*     ***********************************
*     *                                 *
*     *        lcao_diis_dn_end        *
*     *                                 *
*     ***********************************
      subroutine lcao_diis_dn_end()
      implicit none

#include "bafdecls.fh"
#include "lcao_diis_dn.fh"


      call lcao_dnlist_end()

      return
      end



*     ***********************************
*     *									*
*     *	       lcao_diis_dn_generate	*	
*     *									*
*     ***********************************
      subroutine lcao_diis_dn_generate(alpha,ispin,n2ft3d,
     >                                 rho_new,
     >                                 rho_out,
     >                                 failed)
      implicit none
#include "errquit.fh"
      real*8  alpha
      integer ispin,n2ft3d
      real*8  rho_new(n2ft3d,ispin)
      real*8  rho_out(n2ft3d,ispin)
      logical failed

#include "bafdecls.fh"
#include "global.fh"
#include "lcao_diis_dn.fh"

*     ***** local variables ****
      logical value
      integer i,j,size,count,indx,indxt
      integer rho_indx
      integer err_indx1(2),err_indx2(2)
      integer d(2),A(2),ierr
       integer nx,ny,nz
      real*8  scal1,scal2,dv
      real*8  sum,sumt
      integer  lcao_dnlist_size
      real*8   dsum,lattice_omega
      external lcao_dnlist_size
      external dsum,lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()


      size = lcao_dnlist_size()

*     **** create new item in dnlist ***
      if (size.lt.maxsize) then
        call lcao_dnlist_new(ispin*n2ft3d,rho_indx)

*     **** shift dnlist ****
      else
         call lcao_dnlist_shift()
         count = size
         call lcao_dnlist_get(count,rho_indx)
         size = size - 1

      end if

*     **** add rho_new to dnlist ***
      call dcopy((ispin*n2ft3d),rho_new,1,dbl_mb(rho_indx),1)



*     **** exit early *** 
      if (size.lt.2) then
         call dscal(ispin*n2ft3d,(1.0d0-alpha),rho_out,1)
         call daxpy(ispin*n2ft3d,(alpha),rho_new,1,rho_out,1)
         failed = .false.
         return
      end if

*     **** calculate A and B matrices ***
      value = BA_push_get(mt_dbl,size*size,'A',A(2),A(1))
      value = BA_push_get(mt_dbl,ispin*n2ft3d,'err_indx1',
     >                    err_indx1(2),err_indx1(1))
      value = BA_push_get(mt_dbl,ispin*n2ft3d,'err_indx2',
     >                    err_indx2(2),err_indx2(1))

      do i=1,size 
        call lcao_dnlist_get(i+1,rho_indx)
        call dcopy(ispin*n2ft3d,
     >             dbl_mb(rho_indx),1,
     >             dbl_mb(err_indx1(1)),1) 
        call lcao_dnlist_get(i,rho_indx)
        call daxpy(ispin*n2ft3d,(-1.0d0),
     >             dbl_mb(rho_indx),1,
     >             dbl_mb(err_indx1(1)),1) 
        do j=1,i
           call lcao_dnlist_get(j+1,rho_indx)
           call dcopy(ispin*n2ft3d,
     >                dbl_mb(rho_indx),1,
     >                dbl_mb(err_indx2(1)),1) 
           call lcao_dnlist_get(j,  rho_indx)
           call daxpy(ispin*n2ft3d,(-1.0d0),
     >                dbl_mb(rho_indx),1,
     >                dbl_mb(err_indx2(1)),1) 
           call D3dB_rr_dot(1,dbl_mb(err_indx1(1)),
     >                        dbl_mb(err_indx2(1)),sum)
           if (ispin.gt.1) then
             call D3dB_rr_dot(1,dbl_mb(err_indx1(1)+n2ft3d),
     >                          dbl_mb(err_indx2(1)+n2ft3d),sumt)
             sum = sum + sumt
           end if

            sum = sum*dv
            indx  = (i-1) + (j-1)*size
            indxt = (j-1) + (i-1)*size
            dbl_mb(A(1)+indx)  = sum
            dbl_mb(A(1)+indxt) = sum
        end do
      end do
      value = BA_pop_stack(err_indx2(2))
      value = BA_pop_stack(err_indx1(2))

c     write(*,*) "Amatrix:"
c     call lcao_write_matrix_seq(6,size,dbl_mb(A(1)))


*     **** calculate inverse(A) matrix ***

*     **** cholesky factorize ****
      call dpotrf( 'L', size, dbl_mb(A(1)), size, ierr)

*     **** construct inverse of original matrix    ****
*     **** from its lower Cholesky decomposition L ****
      call dpotri( 'L', size, dbl_mb(A(1)), size, ierr)

c     write(*,*) "inv A:"
c     call lcao_write_matrix_seq(6,size,dbl_mb(A(1)))


*     **** calculate d(i) coefficients ****
      failed = .false.
      value = BA_push_get(mt_dbl,(size),'d',d(2),d(1))
      if (.not. value)
     > call errquit('lcao_diis_dn_generate: out of stack memory',0,
     &       MA_ERR)
      sumt = 0.0d0
      do i=1,size
         indx = (i-1)*size     
         sum = dsum(size,dbl_mb(A(1)+indx),1) 
         dbl_mb(d(1)+i-1) = sum
         if (sum.lt.0.0d0) failed = .true.
         sumt = sumt+sum
      end do

      do i=1,size
         dbl_mb(d(1)+i-1) = dbl_mb(d(1)+i-1)/sumt
      end do

      if (.not.failed) then
*     **** calculate diis density ****
      call dcopy(ispin*n2ft3d,0.0d0,0,rho_out,1)
      sumt = 0.0d0
      do i=1,size
         sum = dbl_mb(d(1)+i-1)
         sumt = sumt + sum
         write(*,*) "d:",i,sum
         call lcao_dnlist_get(i,rho_indx)
         call daxpy(ispin*n2ft3d,sum,
     >              dbl_mb(rho_indx),1,
     >              rho_out,1)
      end do
      write(*,*) "sum(d(*)):",sumt
      end if

      value = BA_pop_stack(d(2))
      value = value.and.BA_pop_stack(A(2))
      if (.not. value)
     > call errquit('lcao_diis_dn_generate: popping stack memory',0,
     &       MA_ERR)

      return
      end

