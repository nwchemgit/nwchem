       subroutine mm_set_irtdb(irtdb0)
         implicit none
#include "md_common.fh"
         integer irtdb0
         
         irtdb=irtdb0

       end

       subroutine mm_get_lquant(lquant0)
         implicit none
#include "md_common.fh"
         integer lquant0
         
         lquant0 = lquant

       end

       subroutine mm_get_lqlink(lqlink0)
         implicit none
#include "md_common.fh"
         integer lqlink0
         
         lqlink0 = lqlink

       end

       subroutine mm_get_nsaloc(nsaloc0)
         implicit none
#include "md_common.fh"
         integer nsaloc0
         
         nsaloc0=nsaloc

       end

       subroutine mm_get_nwmloc(nwmloc0)
         implicit none
#include "md_common.fh"
         integer nwmloc0
         
         nwmloc0=nwmloc

       end

       subroutine mm_get_msa(msa0)
         implicit none
#include "md_common.fh"
         integer msa0
         
         msa0=msa

       end

       subroutine mm_get_mwm(mwm0)
         implicit none
#include "md_common.fh"
         integer mwm0
         
         mwm0=mwm

       end

       subroutine mm_get_mwa(mwa0)
         implicit none
#include "md_common.fh"
         integer mwa0
         
         mwa0=mwa

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

       subroutine mm_get_psq(psq)
         implicit none
#include "md_common.fh"
         integer psq
         
         psq = i_is+(lsct1-1)*msa

       end

       subroutine mm_get_psz(psz)
         implicit none
#include "md_common.fh"
         integer psz
         
         psz = i_isz

       end

       subroutine mm_get_pwz(pwz)
         implicit none
#include "md_common.fh"
         integer pwz
         
         pwz = i_iwz

       end

       subroutine mm_get_pxs(pxs)
         implicit none
#include "md_common.fh"
         integer pxs
         
         pxs = i_xs

       end

       subroutine mm_get_pxw(pxw)
         implicit none
#include "md_common.fh"
         integer pxw
         
         pxw = i_xw

       end

       subroutine mm_get_pfs(pfs)
         implicit none
#include "md_common.fh"
         integer pfs
         
         pfs = i_fs

       end


