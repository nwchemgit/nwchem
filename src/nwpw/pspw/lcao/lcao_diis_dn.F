*
* $Id: lcao_diis_dn.F,v 1.1 2001-12-04 02:39:30 bylaska Exp $
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

#include "mafdecls.fh"
#include "lcao_diis_dn.fh"

*     **** local variables ****
      logical value


      maxsize = maxsize0

*     ***** allocate memory from heap memory ****
      value = MA_alloc_get(mt_dbl,(maxsize),'dnerror',
     >                     dnerror(2),dnerror(1))
      if (.not. value)
     > call errquit('lcao_diis_dn_init: out of heap memory',0)

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

#include "mafdecls.fh"
#include "lcao_diis_dn.fh"

*     **** local variables ****
      logical value

      call lcao_dnlist_end()
      value = MA_free_heap(dnerror(2))
      if (.not. value)
     > call errquit('lcao_diis_dn_end: freeing heap memory',0)

      return
      end



*     ***********************************
*     *									*
*     *	       lcao_diis_dn_generate	*	
*     *									*
*     ***********************************
      subroutine lcao_diis_dn_generate(ispin,n2ft3d,rho_new,
     >                                 rho_error,
     >                                 rho_out)
      implicit none
      integer ispin,n2ft3d
      real*8  rho_new(n2ft3d,ispin)
      real*8  rho_error
      real*8  rho_out(n2ft3d,ispin)

#include "mafdecls.fh"
#include "global.fh"
#include "lcao_diis_dn.fh"

*     ***** local variables ****
      logical value
      integer i,j,size,count,indx,indxt
      integer rho_indx,d(2),A(2),B(2),ipiv(2),ierr
      real*8  sum,sumt
      integer  lcao_dnlist_size
      real*8   dsum
      external lcao_dnlist_size
      external dsum


      size = lcao_dnlist_size()


*     **** add rho_new to dnlist ***
      if (size.lt.maxsize) then
        call lcao_dnlist_new(ispin*n2ft3d,rho_indx)

*     **** shift dnlist and dnerror and add rho_new to dnlist ***
      else
         call lcao_dnlist_shift()
         count = size
         call lcao_dnlist_get(count,rho_indx)

*        **** shift dnerror****
         sum = dbl_mb(dnerror(1)) 
         do i=1,size-1
             dbl_mb(dnerror(1)+i-1)=dbl_mb(dnerror(1)+i)
         end do
         dbl_mb(dnerror(1)+size-1)=sum

      end if
      call dcopy((ispin*n2ft3d),rho_new,1,dbl_mb(rho_indx),1)


*     **** put rho_error in dnerror list ***
      size = lcao_dnlist_size()
      dbl_mb(dnerror(1)+size-1) = rho_error


*     **** calculate A and B matrices ***
      value = MA_push_get(mt_dbl,size*size,'B',B(2),B(1))
      value = value.and.
     >        MA_push_get(mt_dbl,size*size,'A',A(2),A(1))

      do i=1,size 
        do j=1,i
            sum = dbl_mb(dnerror(1)+i-1)*dbl_mb(dnerror(1)+j-1)
            indx  = (j-1)*size+(i-1)
            indxt = (i-1)*size+(j-1)
            dbl_mb(A(1)+indx)  = sum
            dbl_mb(A(1)+indxt) = sum

            sum = 0.0d0
            if (i.eq.j) sum = 1.0d0
            dbl_mb(B(1)+indx)  = sum
            dbl_mb(B(1)+indxt) = sum
        end do
      end do
      write(*,*) "Amatrix:"
      call lcao_write_matrix_seq(6,size,dbl_mb(A(1)))
      write(*,*) "Bmatrix:"
      call lcao_write_matrix_seq(6,size,dbl_mb(B(1)))


*     **** calculate inverse(A) matrix ***
      value = MA_push_get(mt_int,size,'ipiv',ipiv(2),ipiv(1))
      if (.not. value) call errquit('out of stack memory',0)
      do i=1,size
         int_mb(ipiv(1)+i-1) = 0
      end do

      call DGESV(size,size,dbl_mb(A(1)),size,
     >                      int_mb(ipiv(1)),
     >                      dbl_mb(B(1)),size,
     >                      ierr)
      write(*,*) "DGESV INFO:",ierr
      value = MA_pop_stack(ipiv(2))
      value = value.and.MA_pop_stack(A(2))
      if (.not. value) call errquit('popping stack memory',0)
      write(*,*) "inv A:"
      call lcao_write_matrix_seq(6,size,dbl_mb(B(1)))



*     **** calculate d(i) coefficients ****
      value = MA_push_get(mt_dbl,(size),'d',d(2),d(1))
      if (.not. value)
     > call errquit('lcao_diis_dn_generate: out of stack memory',0)
      sumt = 0.0d0
      do i=1,size
         indx = (i-1)*size     
         sum = dsum(size,dbl_mb(B(1)+indx),1) 
         dbl_mb(d(1)+i-1) = sum
         sumt = sumt+sum
      end do

      do i=1,size
         dbl_mb(d(1)+i-1) = dbl_mb(d(1)+i-1)/sumt
      end do

*     **** calculate diis density ****
      call dcopy(ispin*n2ft3d,0.0d0,0,rho_out,1)
      sumt = 0.0d0
      do i=1,size
         sum = dbl_mb(d(1)+i-1)
         sumt = sumt + sum
         write(*,*) "d:",i,sum
         call daxpy(ispin*n2ft3d,sum,
     >              dbl_mb(rho_indx),1,
     >              rho_out,1)
      end do
      write(*,*) "sum(d(*)):",sumt

      value = MA_pop_stack(d(2))
      value = value.and.MA_pop_stack(B(2))
      if (.not. value)
     > call errquit('lcao_diis_dn_generate: popping stack memory',0)

      return
      end

