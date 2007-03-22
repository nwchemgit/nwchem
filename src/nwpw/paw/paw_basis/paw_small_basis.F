*
* $Id: paw_small_basis.F,v 1.1 2007-03-22 20:46:20 bylaska Exp $
*

!**************************************************
!
!	Name: init_paw_small_basis
!
!	Purpose 
!
!	Created:	7/11/2002
!**************************************************
      subroutine init_paw_small_basis()
      implicit none
      
#include "paw_basis.fh"
#include "mafdecls.fh"
#include "errquit.fh"
  
*     !*** local variables ***
      integer it, ir, nt,ngrid
      integer l_nb,l_ng
      integer i_nb,i_ng
      integer i_lo,i_hi,ii,mult_l
      integer  max_mult_l
      integer max_ic
      logical ok
      double precision log_amesh,amesh
      character*255 fname

c      !*** external functions ***
      integer ion_nkatm
      external ion_nkatm
      
*     *** check for existence of _basis files ***
      nt = ion_nkatm()

*     **** allocating memory from stack ****
      ok = MA_push_get(MT_INT,nt, "nb_tmp_array", l_nb,i_nb)
      ok = ok.and.
     >     MA_push_get(MT_INT,nt, "ng_tmp_array", l_ng,i_ng)
      if (.not.ok) 
     > call errquit('init_paw_small_basis:out of stack memory',0,MA_ERR)
  
      do it = 1, nt
         ngrid = paw_basis_ngrid(it)
         write(*,*) "small basis, ngrid=",it,ngrid
      end do


*     **** deallocating memory from stack ****
      ok =        MA_pop_stack(l_ng)
      ok = ok.and.MA_pop_stack(l_nb)
      if (.not.ok) 
     > call errquit('init_paw_small_basis:error popping stack',1,MA_ERR)

      return
      end 

