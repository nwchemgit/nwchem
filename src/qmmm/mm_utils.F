       subroutine mm_set_irtdb(irtdb0)
         implicit none
#include "md_common.fh"
         integer irtdb0
         
         irtdb=irtdb0

       end

       subroutine mm_get_nsaloc(nsaloc0)
         implicit none
#include "md_common.fh"
         integer nsaloc0
         
         nsaloc0=nsaloc

       end

       subroutine mm_get_msa(msa0)
         implicit none
#include "md_common.fh"
         integer msa0
         
         msa0=msa

       end

       subroutine mm_get_psdt(psdt)
         implicit none
#include "md_common.fh"
         integer psdt
         
         psdt = i_is+(lsdyn-1)*msa

       end

       subroutine mm_get_psl(psl)
         implicit none
#include "md_common.fh"
         integer psl
         
         psl = i_is+(lsgan-1)*msa

       end

       subroutine mm_get_psat(psat)
         implicit none
#include "md_common.fh"
         integer psat
         
         psat = i_is+(lsatt-1)*msa

       end

       subroutine mm_get_pxs(pxs)
         implicit none
#include "md_common.fh"
         integer pxs
         
         pxs = i_xs

       end


