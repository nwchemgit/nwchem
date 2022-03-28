c
C> \ingroup cosmo
C> @{
C>
C> \file cosmo_tidy.F
C> Tidy up after COSMO
C>
C> \brief Tidy up COSMO data structures
C>
C> This operation currently involves saving the current COSMO charges
C> on the RTDB, and destroying the COSMO geometry instance.
c
      subroutine cosmo_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "prop.fh"
#include "geom.fh"
#include "bq.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb !< [Input] the RTDB handle
c
      call util_print_pop
      if (.not. bq_destroy(cosmo_bq_efc)) call errquit
     $     ('cosmo_tidy: bq destroy failed', 0, GEOM_ERR)
      if (.not. bq_destroy(cosmo_bq_invnuc)) call errquit
     $     ('cosmo_tidy: bq destroy failed', 0, GEOM_ERR)
      return
      end
C>
C> @}
c $Id$
