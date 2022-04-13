      subroutine dft_quickguess(ipol,nbf_ao,g_dens,oskel,
     N     ntotel,nclosed,
     D     densmat)
      implicit none
#include "global.fh"      
#include "errquit.fh"
      integer ipol
      integer nbf_ao
      logical oskel
      integer ntotel,nclosed
      integer g_dens(*)
      logical densmat ![in]
c
      character*256 file_mat1
c      
      if(densmat) then
         call util_file_name('densmat', .false.,.false.,file_mat1)
         if(ga_nodeid().eq.0)
     W        write(6,*) ' reading DM from file ',file_mat1
         call ao_1prdm_read(nbf_ao,g_dens(1),file_mat1)
      endif

      if(ipol.eq.2) then
         if(densmat) then
            call errquit(' code not ready for open-shell ',0,CAPMIS_ERR)
            call util_file_name('densmat', .false.,.false.,file_mat1)
            if(ga_nodeid().eq.0)
     W           write(6,*) ' reading DM from file ',file_mat1
            call ao_1prdm_read(nbf_ao,g_dens(2),file_mat1)
         else
            call ga_copy(g_dens(1),g_dens(2))
            call ga_dscal(g_dens(1),dble(ntotel-nclosed)/(ntotel))
            call ga_dscal(g_dens(2),dble(nclosed)/(ntotel))
         endif
      endif
      if(oskel) then
        call ga_symmetrize(g_dens(1))
        if(ipol.eq.2) call ga_symmetrize(g_dens(2))
      endif
      return
      end
