       subroutine mm_set_irtdb(irtdb0)
         implicit none
#include "md_common.fh"
         integer irtdb0
         
         irtdb=irtdb0

       end

       subroutine mm_set_sd_maxiter(n)
         implicit none
#include "md_common.fh"
         integer n
         
         msdit = n 

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

       subroutine mm_get_mdacq(mdacq0)
         implicit none
#include "md_common.fh"
         integer mdacq0
         
         mdacq0 = mdacq

       end

       subroutine mm_set_mdacq(mdacq0)
         implicit none
#include "md_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
         integer mdacq0
         
         mdacq = mdacq0
         ldacq = mdacq0
     
        if(.not.rtdb_put(irtdb,'md:mdacq',mt_int,1,mdacq0))
     + call md_abort('md_input: rtdb_put failed',0)
        if(.not.rtdb_put(irtdb,'md:ldacq',mt_int,1,mdacq0))
     + call md_abort('md_input: rtdb_put failed',0)
   
       end

       subroutine mm_set_high_precision()
         implicit none
#include "md_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
         nprec = 2
         if(.not.rtdb_put(irtdb,'md:nprec',mt_int,1,nprec))
     + call md_abort('md_input: rtdb_put failed',0)

       end

       subroutine mm_get_mequi(mequi0)
         implicit none
#include "md_common.fh"
         integer mequi0
         
         mequi0 = mequi

       end

       subroutine mm_set_mequi(mequi0)
         implicit none
#include "md_common.fh"
         integer mequi0
         
         mequi = mequi0

       end

       subroutine mm_set_nfc(nf,nfs)
         implicit none
#include "rtdb.fh"
#include "md_common.fh"
#include "mafdecls.fh"
         integer nf,nfs

         if(.not.rtdb_put(irtdb,'md:nfcoor',mt_int,1,nf))
     +    call md_abort('md_rtdbin: rtdb_get failed',50)
         if(.not.rtdb_put(irtdb,'md:nfscoo',mt_int,1,nfs))
     +    call md_abort('md_rtdbin: rtdb_get failed',51)
         
         nfcoor = nf
         nfscoo = nfs

       end

       subroutine mm_get_nfc(nf,nfs)
         implicit none
#include "md_common.fh"
         integer nf,nfs
         
         nf  = nfcoor
         nfs =nfscoo

       end

       subroutine mm_get_nfrest(nf)
         implicit none
#include "md_common.fh"
         integer nf
         
         nf  = nfrest

       end

       subroutine mm_set_keep_restart()
         implicit none
#include "md_common.fh"
         
         keepr = 1

       end

       subroutine mm_unset_keep_restart()
         implicit none
#include "md_common.fh"
         
         keepr = 0

       end

       subroutine mm_get_keep_restart_index(i)
         implicit none
#include "md_common.fh"
         integer i
         
         i = keepr-1

       end

       subroutine mm_get_msa(msa0)
         implicit none
#include "md_common.fh"
         integer msa0
         
         msa0=msa

       end

       subroutine mm_get_temp(temp0)
         implicit none
#include "md_common.fh"
         double precision temp0
         
         temp0 = temp 

       end

       subroutine mm_get_stime(stime0)
         implicit none
#include "md_common.fh"
         double precision stime0
         
         stime0 = stime

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

       subroutine mm_get_pwdt(pwdt)
         implicit none
#include "md_common.fh"
         integer pwdt
         
         pwdt = i_iw+(lwdyn-1)*mwm

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

       subroutine mm_get_psgm(psgm)
         implicit none
#include "md_common.fh"
         integer psgm
         
         psgm = i_is+(lssgm-1)*msa

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

       function mm_classical(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_classical
c        local variables
         integer psdt

         call mm_get_psdt(psdt)
         mm_classical = iand(int_mb(psdt+i-1),lquant).eq.0
         
       end

       function mm_ww(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         double precision  mm_ww
c        local variables

         mm_ww = dbl_mb(i_wws+i-1)
         
       end

       function mm_record()
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_record
c        local variables

         mm_record = lxw.or.lxs
         
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
c         mm_bqzone = mm_bqzone .and. iand(int_mb(psdt+i-1),lqlink).eq.0
         mm_bqzone = mm_bqzone .and. int_mb(psz+i-1).gt.0
         
       end

      function mm_active_bqs(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_bqs
c        local variables
         integer psdt,psz

         call mm_get_psdt(psdt)
         call mm_get_psz(psz)

         mm_active_bqs = (iand(int_mb(psdt+i-1),mfixed).ne.lfixed).and.
     >               (iand(int_mb(psdt+i-1),lquant).eq.0).and.
     >               (int_mb(psz+i-1).gt.0)
         
       end

      function mm_active_bqw(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_bqw
c        local variables
         integer pwdt,pwz

         call mm_get_pwdt(pwdt)
         call mm_get_pwz(pwz)

         mm_active_bqw = (iand(int_mb(pwdt+i-1),mfixed).ne.lfixed).and.
     >               (int_mb(pwz+i-1).gt.0)
         
       end

      function mm_fixed_w(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_fixed_w
c        local variables
         integer pwdt

         call mm_get_pwdt(pwdt)

         mm_fixed_w = iand(int_mb(pwdt+i-1),mfixed).eq.lfixed
         
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

      function mm_active(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active
c        local variables
         integer psdt

         call mm_get_psdt(psdt)

         mm_active = iand(int_mb(psdt+i-1),mfixed).ne.lfixed
         
       end

      function mm_active_qm(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_qm
c        local variables
         integer psdt

         call mm_get_psdt(psdt)

         mm_active_qm = (iand(int_mb(psdt+i-1),mfixed).ne.lfixed)
     >                   .and.
     >                  (iand(int_mb(psdt+i-1),lquant).ne.0) 
         
       end

      function mm_active_classical(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_classical
c        local variables
         integer psdt

         call mm_get_psdt(psdt)

         mm_active_classical = 
     >                (iand(int_mb(psdt+i-1),mfixed).ne.lfixed)
     >                .and.
     >               (iand(int_mb(psdt+i-1),lqlink).eq.0) 
     >                .and.
     >               (iand(int_mb(psdt+i-1),lquant).eq.0) 
         
       end

      function mm_active_link(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_link
c        local variables
         integer psdt

         call mm_get_psdt(psdt)

         mm_active_link = (iand(int_mb(psdt+i-1),mfixed).ne.lfixed)
     >                .and.
     >               (iand(int_mb(psdt+i-1),lqlink).ne.0) 
         
       end

      function mm_active_w(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_active_w
c        local variables
         integer pwdt

         call mm_get_pwdt(pwdt)

         mm_active_w = iand(int_mb(pwdt+i-1),mfixed).ne.lfixed
         
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

      function mm_qmbqzone(i)
         implicit none
#include "mafdecls.fh"
#include "md_common.fh"
         integer i
         logical mm_qmbqzone
c        local variables
         integer psdt,psz

         call mm_get_psdt(psdt)
         call mm_get_psz(psz)

         mm_qmbqzone = iand(int_mb(psdt+i-1),lquant).eq.0 .or.
     >                 int_mb(psz+i-1).gt.0
         
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


c $Id$
