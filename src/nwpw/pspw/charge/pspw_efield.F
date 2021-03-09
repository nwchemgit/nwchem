c
c $Id$
c
  
***********************************************************
*                                                         *
*                   PSPW-Efield module                    *
*                                                         *
***********************************************************

*     **********************************
*     *	                               *
*     *          pspw_Efield_init      *
*     *                                *
*     **********************************
 
      subroutine pspw_Efield_init()
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

      return
      end



*     **********************************
*     *	                               *
*     *        pspw_Efield_end         *
*     *                                *
*     **********************************
      subroutine pspw_Efield_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

      return
      end


*     **********************************
*     *	                               *
*     *         pspw_Efield_found      *
*     *                                *
*     **********************************
      logical function pspw_Efield_found()
      implicit none

#include "pspw_efield.fh"
 

      pspw_efield_found = efield_on
      return
      end

