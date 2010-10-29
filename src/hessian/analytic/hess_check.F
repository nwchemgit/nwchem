c
      logical function hess_check(rtdb)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer rtdb      ! [input] rtdb handle
c
      integer geom1, basis1
      integer ipol, mult
      logical cdfit, dmd_rel, dk_rel, zora_rel
      logical xc_chktau,xc_getxcfac,cam_exch
      external xc_chktau,xc_getxcfac
c
      hess_check = .false.
      cdfit = .false.
      ipol = 1
c
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_check: theory not specified',555,
     &       RTDB_ERR)
      if (theory.eq.'hyb') theory = 'dft'
c
c  Check for relativistic all-electron methods
c
      if (.not. rtdb_get(rtdb, 'dyall_mod_dir', MT_LOG, 1, dmd_rel))
     $   dmd_rel = .false.
      if (.not. rtdb_get(rtdb, 'doug_kroll', MT_LOG, 1, dk_rel))
     $   dk_rel = .false.
      if (.not. rtdb_get(rtdb, 'zora', MT_LOG, 1, zora_rel))
     $   zora_rel = .false.
c 
      if (dmd_rel.or.dk_rel.or.zora_rel) then
         if (ga_nodeid().eq.0)
     $    write(LuOut,*) 'Relativistic all-electron method detected'
         return
      endif
c
c  Check for ROHF in SCF
c
      if (theory.eq.'scf') then
        if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     *    nopen = 0
        if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
          if (nopen .eq. 0) then
            scftype = 'RHF'
          else if (scftype .ne. 'UHF') then
            scftype = 'ROHF'
          endif
        endif
        if (scftype.eq.'ROHF') then
           if (ga_nodeid().eq.0)
     *       write(LuOut,*) 'ROHF detected'
           return
        endif
      endif
c
c  Check for UHF in DFT
c
      if (theory.eq.'dft') then
        if (.not. rtdb_get(rtdb, 'dft:nopen', MT_INT, 1, nopen))
     $    nopen = 0
        if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1, mult))
     $    mult = 1
        if (.not. rtdb_get(rtdb, 'dft:ipol', MT_INT, 1, ipol))
     $    ipol = 1
        if (ipol.gt.1.or.nopen.gt.0.or.mult.gt.1) then
          if (ga_nodeid().eq.0)
     $      write(LuOut,*) 'Open shell DFT detected'
          return
        endif
      endif
c
c  Check for CD fit
c
      if (theory.eq.'dft') then
         if (.not. geom_create(geom1, 'geometry'))
     $        call errquit('hess_check: geom_create problem', 555,
     &       GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom1, 'geometry'))
     $        call errquit('hess_check: no geometry ', 555,
     &       RTDB_ERR)
         if (.not. bas_create(basis1, 'cd basis'))
     $        call errquit('hess_check: bas_create?',102, BASIS_ERR)
         cdfit=bas_rtdb_load(rtdb, geom1, basis1, 'cd basis')
         if (.not. geom_destroy(geom1))
     $        call errquit('hess_check:not able to destroy geometry:',
     $        geom1, GEOM_ERR)
         if (.not. bas_destroy(basis1))
     $        call errquit('hess_check:not able to destroy CD_bas:',
     $                   basis1, BASIS_ERR)
      endif
c
c     check for metaGGAs
c
      if(theory.eq.'dft') then
         if(xc_getxcfac(rtdb)) then
            if(xc_chktau()) then
               if (ga_nodeid().eq.0)
     $              write(LuOut,*) 'MetaGGA xc detected'
               return
            endif
         endif
      endif
c
c     check for attenuated exchange-correlation forms
c
      cam_exch = .false.
      if(theory.eq.'dft') then
        if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
         if (cam_exch) then
             if (ga_nodeid().eq.0) 
     &              write(LuOut,*) 'CAM/LC xc detected'
             return
         end if
      end if
c
c     check for cdfit
c      
      if (cdfit) then
        if (ga_nodeid().eq.0)
     $    write(LuOut,*) 'CD fit detected'
        return
      endif
c
      hess_check = .true.
c
      return
      end
