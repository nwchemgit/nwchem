C
C     Saves some DFT matrices to file.  This is mainly used for the
C     real-time TDDFT code but might be useful elsewhere.
C
C     Uses dmat_io() in rt_tddft/rtutils
C
C     It currently prints the following matrices:
C     - converged density matrix in AO basis (alpha, beta for UKS)
C     
      subroutine dft_save_mats (g_dens)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
      
C     == Inputs ==
      integer g_dens(2)
      
C     == Variables ==
      integer me
      logical, external :: dmat_io_dump
      character(255) fname

      me = ga_nodeid()

C
C     Density matrix
C
      if (ipol .eq. 1) then
         call util_file_name ("dens_rks",
     $        .false., .false., fname)
         if (.not. dmat_io_dump (g_dens(1), fname))
     $        call errquit ("Failed to dump"//trim(fname),0,0)
         if (me.eq.0) then
            write (luout, *) "Saved RKS density matrix to file"
         endif
         
      elseif (ipol .eq. 2) then
         call util_file_name ("dens_uks_alpha",
     $        .false., .false., fname)
         if (.not. dmat_io_dump (g_dens(1), fname))
     $        call errquit ("Failed to dump"//trim(fname),0,0)
         call util_file_name ("dens_uks_beta",
     $        .false., .false., fname)
         if (.not. dmat_io_dump (g_dens(2), fname))
     $        call errquit ("Failed to dump"//trim(fname),0,0)
         if (me.eq.0) then
            write (luout, *) "Saved UKS density matrices to file"
         endif
      else
         call errquit ("Bad ipol value",0,0)
      endif
      
      end subroutine
