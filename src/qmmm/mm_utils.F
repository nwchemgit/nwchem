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

       subroutine mm_get_pwl(pwl)
         implicit none
#include "md_common.fh"
         integer pwl
         
         pwl = i_iw+(lwgmn-1)*mwm

       end

       subroutine mm_get_psl(psl)
         implicit none
#include "md_common.fh"
         integer psl
         
         psl = i_is+(lsgan-1)*msa

       end

       subroutine mm_get_psg(psg)
         implicit none
#include "md_common.fh"
         integer psg
         
         psg = i_is+(lsgrp-1)*msa

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

       subroutine mm_get_solute_coord(i,r)
         implicit none
#include "md_common.fh"
#include "mafdecls.fh"
         integer i
         double precision r(3) 
         integer j
         
         do j=1,3
           r(j) = dbl_mb(i_xs+(i-1)+(j-1)*msa)
         end do

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

       subroutine mm_get_pfw(pfw)
         implicit none
#include "md_common.fh"
         integer pfw
         
         pfw = i_fw

       end

      subroutine mm_clear_bqzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
c        local variables
         integer psz

         call mm_get_psz(psz)

         int_mb(psz+i-1)=0
         
       end

      subroutine mm_clear_bqwzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
c        local variables
         integer pwz

         call mm_get_pwz(pwz)

         int_mb(pwz+i-1)=0
         
       end

       function mm_quant(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_quant
c        local variables
         integer psdt

         call mm_get_psdt(psdt)
         mm_quant = iand(int_mb(psdt+i-1),lquant).ne.0
         
       end

       function mm_quant_link(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_quant_link
c        local variables
         integer psdt

         call mm_get_psdt(psdt)
         mm_quant_link = iand(int_mb(psdt+i-1),lquant).ne.0
         mm_quant_link = mm_quant_link .or. 
     >                   iand(int_mb(psdt+i-1),lqlink).ne.0
          
       end

       function mm_link(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_link
c        local variables
         integer psdt

         call mm_get_psdt(psdt)
         mm_link = iand(int_mb(psdt+i-1),lqlink).ne.0
         
       end

       function mm_lqsoft(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_lqsoft
c        local variables
         integer psdt,lqsoft

         lqsoft = 64
         call mm_get_psdt(psdt)
         mm_lqsoft = iand(int_mb(psdt+i-1),lqsoft).ne.0
         
       end       

      function mm_bqwzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_bqwzone
c        local variables
         integer pwz

         call mm_get_pwz(pwz)

         mm_bqwzone = int_mb(pwz+i-1).gt.0
         
       end

      function mm_bqzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_bqzone
c        local variables
         integer psdt,psz

         call mm_get_psdt(psdt)
         call mm_get_psz(psz)

         mm_bqzone = iand(int_mb(psdt+i-1),lquant).eq.0
         mm_bqzone = mm_bqzone .and. iand(int_mb(psdt+i-1),lqlink).eq.0
         mm_bqzone = mm_bqzone .and. int_mb(psz+i-1).gt.0
         
       end

      function mm_fixed(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_fixed
c        local variables
         integer psdt

         call mm_get_psdt(psdt)

         mm_fixed = iand(int_mb(psdt+i-1),mfixed).eq.lfixed
         
       end

      function mm_identity(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_identity

         mm_identity = .true.
         
       end

      function mm_link_bqzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_link_bqzone
c        local variables
         integer psdt,psz

         call mm_get_psdt(psdt)
         call mm_get_psz(psz)

         mm_link_bqzone = iand(int_mb(psdt+i-1),lquant).eq.0
         mm_link_bqzone = mm_link_bqzone .and. int_mb(psz+i-1).gt.0
         
       end

      function mm_element(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         integer mm_element
c        local variables
         integer psat
         integer cf_element
         external cf_element

         call mm_get_psat(psat)

         mm_element = cf_element(int_mb(psat+i-1))
         
       end

      function mm_solute_charge(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         double precision  mm_solute_charge
c        local variables
         double precision cf_charge
         external cf_charge

         mm_solute_charge= cf_charge(int_mb(i_is+(lsct1-1)*msa+i-1))
         
       end

       subroutine mm_print_psl_psz()
         implicit none
#include "md_common.fh"
#include "mafdecls.fh"
         integer ns,psz,psl
         integer i
         
         call  mm_get_psz(psz)
         call  mm_get_psl(psl)
         call mm_get_nsaloc(ns)
         
         do i=1,ns
c           write(*,*) int_mb(psl+i-1),
c     >                int_mb(psz+i-1)
         end do

       end


