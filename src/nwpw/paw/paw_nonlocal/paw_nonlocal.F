*
* $Id: paw_nonlocal.F,v 1.1 2003-02-26 02:13:16 bylaska Exp $
*

!**************************************************
!
!	Name: paw_nonlocal_init
!
!	Purpose: initializes the paw nonlocal operator
!
!	Created:	2/25/2003
!**************************************************
      subroutine paw_nonlocal_init(ispin0,ne)
      implicit none
      integer ispin0
      integer ne(2)

#include "mafdecls.fh"
#include "paw_proj.fh"
#include "paw_nonlocal_data.fh" 
#include "paw_geom.fh" 
#include "paw_ma.fh" 
#include "paw_basis.fh" 
  
*     **** local variables ****
      logical value
      integer ii,ia,nion

      nion = ion_nion()

      !*** calculate total number of (n,l,m) projectors  - Yet Again ***
      total_nbasis = 0
      do ia=1,ion_nkatm()
         total_nbasis = total_nbasis 
     >                + paw_proj_nbasis(ia)*ion_natm(ia)
      end do

      !*** calculate the size and indexing for Gop  - Yet Again ***
      value = my_alloc(mt_int,nion,"i_paw_Gop",i_paw_Gop)
      if (.not.value) call errquit('paw_nonlocal_init: alloc heap',1)

      paw_Gop_size = 0
      do ii=1,nion
         int_mb(i_paw_Gop(1)+ii-1) = paw_Gop_size
         ia = ion_katm(ii)
         paw_Gop_size = paw_Gop_size
     >                + paw_proj_nbasis(ia)**2
      end do


      !*** allocate space for the non-local operator ***
      !-- total number of electron states
      ispin = ispin0
      neall = ne(1) + ne(2)
      ne_u=ne(1)
      ne_d=ne(2)

      !*** space for the non-local operator ***
      value = value.and.MA_alloc_get(mt_dcpl,
     >            (total_nbasis*neall),
     >            'paw_nl_coeff',paw_nl_coeff(2),paw_nl_coeff(1))
      value = value.and.my_alloc(mt_dcpl,paw_Gop_size,
     >                           "paw_Gop_u",paw_Gop_u)
      value = value.and.my_alloc(mt_dcpl,paw_Gop_size,
     >                           "paw_Gop_d",paw_Gop_d)
      if (.not.value) call errquit('paw_nonlocal_init: alloc heap',1)


      return
      end

!**************************************************
!
!	Name: paw_nonlocal_end
!
!	Purpose: removes space used by the paw nonlocal
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_nonlocal_end()
      implicit none      

#include "mafdecls.fh"
#include "paw_ma.fh" 
#include "paw_nonlocal_data.fh"

      !*** local variables ***
      logical value
      
      value = value.and.my_dealloc(i_paw_Gop)
      value = value.and.my_dealloc(paw_Gop_u)
      value = value.and.my_dealloc(paw_Gop_d)
      value = value.and.my_dealloc(paw_nl_coeff)
      if (.not.value) call errquit('paw_nonlocal_end: dealloc heap',0)
      return
      end


!**************************************************
!
!	Name: paw_nonlocal_gen_nl_coeff
!
!	Purpose: 
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_nonlocal_gen_nl_coeff()
      implicit none      

#include "mafdecls.fh"
#include "paw_ma.fh" 
#include "paw_nonlocal_data.fh"

      return
      end
