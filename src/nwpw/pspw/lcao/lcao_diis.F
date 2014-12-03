*
* $Id$
*

*     ***********************************************************
*     *															*
*     *   		   		lcao_diis library					*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*
*

*     ***********************************
*     *									*
*     *	       lcao_diis_init		*	
*     *									*
*     ***********************************
      subroutine lcao_diis_init(maxsize0,ispin0,nbasis0)
      implicit none
      integer maxsize0
      integer ispin0
      integer nbasis0

#include "bafdecls.fh"
#include "lcao_diis.fh"


      maxsize = maxsize0
      ispin   = ispin0
      nbasis  = nbasis0

*     **** initialize mlist ****
      call lcao_mlist_init(maxsize,ispin,nbasis)

      return
      end

*     ***********************************
*     *                                 *
*     *        lcao_diis_end 	        *
*     *                                 *
*     ***********************************
      subroutine lcao_diis_end()
      implicit none

#include "bafdecls.fh"
#include "lcao_diis.fh"


      call lcao_mlist_end()

      return
      end



*     ***********************************
*     *									*
*     *	       lcao_diis_generate		*	
*     *									*
*     ***********************************
      subroutine lcao_diis_generate(alpha,
     >                              vhxc_new,
     >                              vhxc_out,
     >                              failed)
      implicit none
#include "errquit.fh"
      real*8  alpha
      integer vhxc_new(*)
      integer vhxc_out(*)
      logical failed

#include "bafdecls.fh"
#include "global.fh"
#include "lcao_diis.fh"

*     ***** local variables ****
      logical value
      integer i,j,size,count,indx,indxt
      integer vhxc_indx(2),vhxc_indx2(2)
      integer err_indx1(2),err_indx2(2)
      integer d(2),A(2),ierr
      real*8  scal
      real*8  sum,sumt
      integer  lcao_mlist_size
      real*8   dsum,lattice_omega
      external lcao_mlist_size
      external dsum,lattice_omega

      scal = 1.0d0/dble(nbasis*nbasis)
      

      size = lcao_mlist_size()

*     **** create new item in mlist ***
      if (size.lt.maxsize) then
        call lcao_mlist_new(vhxc_indx)

*     **** shift mlist ****
      else
         call lcao_mlist_shift()
         count = size
         call lcao_mlist_get(count,vhxc_indx)
         size = size - 1

      end if
      call ga_copy(vhxc_new(1),vhxc_indx(1))
      if (ispin.gt.1)
     > call ga_copy(vhxc_new(2),vhxc_indx(2))

*     **** exit early *** 
      if (size.lt.3) then
         call ga_add(alpha,vhxc_new(1),
     >               (1.0d0-alpha),vhxc_out(1),
     >                vhxc_out(1))
         if (ispin.gt.1)
     >    call ga_add(alpha,vhxc_new(2),
     >                (1.0d0-alpha),vhxc_out(2),
     >                 vhxc_out(2))
         failed = .false.
         return
      end if

*     **** calculate A  matrices ***
      value = BA_push_get(mt_dbl,size*size,'A',A(2),A(1))
      if (.not.value)
     > call errquit('lcao_diis: out of memory',0, MA_ERR)

      value = ga_create(mt_dbl,nbasis,nbasis,
     >                  'err_indx1',nbasis,0,err_indx1(1))
      if (ispin.gt.1)
     > value = value.and.
     >         ga_create(mt_dbl,nbasis,nbasis,
     >                   'errt',nbasis,0,err_indx1(2))

      value = value.and.ga_create(mt_dbl,nbasis,nbasis,
     >                  'err_indx1',nbasis,0,err_indx2(1))
      if (ispin.gt.1)
     > value = value.and.
     >         ga_create(mt_dbl,nbasis,nbasis,
     >                   'errt',nbasis,0,err_indx2(2))

      if (.not.value)
     > call errquit('lcao_diis: out of ga memory',0, MA_ERR)


      do i=1,size 
        call lcao_mlist_get(i+1,vhxc_indx)
        call lcao_mlist_get(i,  vhxc_indx2)
        call ga_add((1.0d0), vhxc_indx(1),
     >              (-1.0d0),vhxc_indx2(1),
     >              err_indx1(1))
        if (ispin.gt.1)
     >  call ga_add((1.0d0), vhxc_indx(2),
     >              (-1.0d0),vhxc_indx2(2),
     >              err_indx1(2))

        do j=1,i
          call lcao_mlist_get(j+1,vhxc_indx)
          call lcao_mlist_get(j,  vhxc_indx2)
          call ga_add((1.0d0), vhxc_indx(1),
     >                (-1.0d0),vhxc_indx2(1),
     >                err_indx2(1))
          if (ispin.gt.1)
     >    call ga_add((1.0d0), vhxc_indx(2),
     >                (-1.0d0),vhxc_indx2(2),
     >                err_indx2(2))

            sum = ga_ddot(err_indx1(1),err_indx2(1))
            if (ispin.gt.1)
     >        sum = sum + ga_ddot(err_indx1(2),err_indx2(2))
            sum = sum*scal
            indx  = (i-1) + (j-1)*size
            indxt = (j-1) + (i-1)*size
            dbl_mb(A(1)+indx)  = sum
            dbl_mb(A(1)+indxt) = sum
        end do
      end do
      value = ga_destroy(err_indx1(1))
      value = ga_destroy(err_indx2(1))
      if (ispin.gt.1)
     > value = ga_destroy(err_indx1(2))
      if (ispin.gt.1)
     > value = ga_destroy(err_indx2(2))

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
      call ga_zero(vhxc_out(1))
      if (ispin.gt.0)
     >  call ga_zero(vhxc_out(2))
      sumt = 0.0d0
      do i=1,size
         sum = dbl_mb(d(1)+i-1)
         sumt = sumt + sum
         write(*,*) "d:",i,sum
         call lcao_mlist_get(i,vhxc_indx)
         call ga_add(sum,vhxc_indx(1),
     >               (1.0d0),vhxc_out(1),
     >                vhxc_out(1))
         if (ispin.gt.1)
     >   call ga_add(sum,vhxc_indx(2),
     >               (1.0d0),vhxc_out(2),
     >                vhxc_out(2))
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

