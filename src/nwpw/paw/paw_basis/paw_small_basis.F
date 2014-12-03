*
* $Id$
*

!**************************************************
!
!       Name: init_paw_small_basis
!
!       Purpose 
!
!       Created:        7/11/2002
!**************************************************
      subroutine init_paw_small_basis()
      implicit none
      
#include "paw_basis.fh"
#include "bafdecls.fh"
#include "errquit.fh"
  
*     !*** local variables ***
      integer it, ir, nt,nn
      integer l_nb,l_ng
      integer i_nb,i_ng
      logical ok

c      !*** external functions ***
      integer ion_nkatm
      external ion_nkatm
      
*     *** check for existence of _basis files ***
      nt = ion_nkatm()

*     **** allocating memory from stack ****
      ok = BA_push_get(mt_int,nt,"nb_tmp_array",l_nb,i_nb)
      ok = ok.and.
     >     BA_push_get(mt_int,nt,"ng_tmp_array",l_ng,i_ng)
      if (.not.ok) 
     > call errquit('init_paw_small_basis:out of stack memory',0,MA_ERR)

  
      do it = 1, nt
         nn = paw_basis_ngrid(it)
      end do


*     **** deallocating memory from stack ****
      ok =        BA_pop_stack(l_ng)
      ok = ok.and.BA_pop_stack(l_nb)
      if (.not.ok) 
     > call errquit('init_paw_small_basis:error popping stack',1,MA_ERR)

      return
      end 

