      subroutine get_vecF1(g_vecF1,    ! out:
     &                     g_vecF1_im, ! out:
     &                     g_rhs,      ! in :
     &                     g_rhs_im,   ! in :
     &                     vectors,    ! in : MO vectors
     &                     nbf,        ! in : nr. basis functions
     &                     nmo,        ! in : nr. MOs
     &                     ncomp,      ! in :
     &                     npol,       ! in : nr. polarizations
     &                     lifetime,   ! in : = (.true.,.false.) with/out damping
     &                     nocc,       ! in : nr. occupied MOs
     &                     nvirt,      ! in : nr. virtual  MOs
     &                     debug)      ! in : = .true. for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer ispin,npol,ncomp,nocc(npol),nvirt(npol)
c Note.- g_vecF1(npol,ncomp) does not work for two dimensions
      integer g_vecF1(2,2),g_vecF1_im(2,2),
     &        g_rhs(ncomp),g_rhs_im(ncomp),
     &        g_u
      integer vectors(npol)
      integer ipm,ifld,ndir,nbf,nmo,disp
      logical lifetime,debug
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        plo(3), phi(3),
     &        qlo(3), qhi(3)
      character*256 cstemp
      ndir=3 ! = nr directions (x,y,z)
c ---- Create (g_vecF1,g_vecF1_im) -------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
      do ispin=1,npol
       do ipm=1,ncomp
        write (cstemp,'(a,i1)') 'aor vecF1 ',ipm
        if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),
     &    alo,g_vecF1(ispin,ipm) ))
     &     call errquit('fiao_f1: nga_create failed vecF1',0,GA_ERR)
        call ga_zero(g_vecF1(ispin,ipm))
        if (lifetime) then
          write (cstemp,'(a,i1)') 'aor vecF1_Im ',ipm
          if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),
     &       alo,g_vecF1_im(ispin,ipm) ))
     &       call errquit('fiao_f1: nga_create failed E1_im',0,GA_ERR)
          call ga_zero(g_vecF1_im(ispin,ipm))
        end if
       enddo ! end-loop-ipm
      enddo ! end-loop-ispin
c ---- Create (g_vecF1,g_vecF1_im) -------- END

      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf  
       ahi(2) = nocc(ispin)
       ahi(3) = 1
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('giao_b1: nga_create failed g_u',0,GA_ERR)
c     C1 = C0 * U10
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)  
       alo(3) = 1
       ahi(3) = 1       
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nmo 
       blo(3) = 1
       bhi(3) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nocc(ispin)
       disp=nocc(1)*nvirt(1)*(ispin-1)
       plo(1) = disp+1
       phi(1) = disp+nocc(ispin)*nvirt(ispin)
       qlo(1) = nocc(ispin)+1
       qhi(1) = nmo
       qlo(2) = 1
       qhi(2) = nocc(ispin)
       qlo(3) = 1
       qhi(3) = 1
       do ifld = 1,ndir
        clo(3) = ifld
        chi(3) = ifld
c      Make C1           
        do ipm = 1,ncomp
c ======== Including g_u_ov (g_rhs --> g_u) ==== START
        plo(2) = ifld
        phi(2) = ifld
        call ga_zero(g_u)
        call nga_copy_patch('n',g_rhs(ipm),plo,phi,
     &                          g_u       ,qlo,qhi)
c ======== Including g_u_ov (g_rhs --> g_u) ==== END
         call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                         vectors(ispin)    ,blo,bhi,
     &                         g_u               ,alo,ahi,
     &                         g_vecF1(ispin,ipm),clo,chi)          
         if (lifetime) then
c ======== Including g_u_ov_im (g_rhs_im --> g_u_im) == START
         call ga_zero(g_u)
         call nga_copy_patch('n',g_rhs_im(ipm),plo,phi,
     &                           g_u          ,qlo,qhi)
c ======== Including g_u_ov_im (g_rhs_im --> g_u_im) == END
         call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                         vectors(ispin)       ,blo,bhi,
     &                         g_u                  ,alo,ahi,
     &                         g_vecF1_im(ispin,ipm),clo,chi)  
         end if ! end-if-lifetime (damping)     
        enddo ! end-loop-ipm
       enddo ! end-loop-ifld
      if (.not.ga_destroy(g_u)) call 
     &    errquit('fiao_b1: ga_destroy failed g_d1',0,GA_ERR)
      enddo ! end-loop-ispin
      return
      end
c
c     JEM: Updated to include imaginary part of rhs
      subroutine update_rhs_dipole(
     &                     g_rhs,    ! in/out: 
     &                     g_rhs_im, ! in/out:
     &                     vectors,  ! in : MO vectors
     &                     rtdb,     ! in : rtdb  handle
     &                     basis,    ! in : basis handle
     &                     lvelocity,! in : logical var
     &                     nat,      ! in : nr. atoms
     &                     npol,     ! in : nr of polarizations
     &                     nocc,     ! in : nr. occ  shells
     &                     nvirt,    ! in : nr. virt shells
     &                     nbf,      ! in : nr. basis functions
     &                     nmo,      ! in : nr. MOs   
     &                     ncomp,    ! in : nr components of ...
     &                     debug,    ! in : logical for debugging
     &                     omega,    ! in : freq value
     &                     lifetime) ! in : lifetime
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "prop.fh"
#include "util.fh"
      integer npol,ncomp
      integer g_dipole,g_dipole1,g_rhs(2),g_rhs_im(2)
      integer rtdb,basis
      integer vectors(npol)
      integer ifld,ndir,nbf,nmo,ipm,nat,
     &        nocc(npol),nvirt(npol)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      integer shift,disp,ispin
      integer alo(3), ahi(3),
     &        blo(3), bhi(3) 
      logical lvelocity,oskel,debug
c  DIM/QM stuffs
      logical use_dimqm, ldimqm, lrsp, lclfld, lifetime
      integer g_dipel_i, g_dipel_i_tmp
      double precision omega
c  DIM/QM END
      oskel = .false.
      ndir=3 ! nr directions (x,y,z)
c     Get dipole integrals in GA and transform to MO set (virt,occ)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
      if (.not.nga_create(MT_DBL,3,ahi,'dip matrix',alo,g_dipole1)) call 
     &    errquit('fiao_f1: nga_create failed g_dipole',0,GA_ERR)
      ahi(3) = ndir*npol
      if (.not.nga_create(MT_DBL,3,ahi,'dip matrix',alo,g_dipole)) call 
     &    errquit('fiao_f1: nga_create failed g_dipole',0,GA_ERR)
c jbecca START DIMQM
      use_dimqm=util_module_avail('dimqm')

      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlrsp(lrsp)
         call dimqm_getlclfld(lclfld)
      else
         ldimqm   =  .false.
         lrsp     =  .false.
         lclfld   =  .false.
      endif

      if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
       if (.not.nga_create(MT_DBL,3,ahi,'dip matrix',alo,g_dipel_i_tmp)) 
     $ call errquit('fiao_f1: nga_create failed g_dipole',0,GA_ERR)
       call ga_zero(g_dipel_i_tmp)
      endif

c jbecca END DIMQM

      if (debug) write (luout,*) 'fiao_f1: dipole matrix allocated'
c     Get H10 in GA, using g_dipole array
c     note: origin has been set to (0,0,0) for multipole integs.
      call ga_zero(g_dipole1)
      if (lvelocity) then 
        call int_giao_1ega(basis,basis,
     &                     g_dipole1,'velocity', ! out
     &                     origin,nat, oskel)
        call ga_scale (g_dipole1, -1d0)
      else
        call int_mpole_1ega(basis,basis,
     &                      g_dipole1,'dipole', ! out
     &                      origin,oskel)
      end if
      if (debug) write (luout,*) 'fiao_f1: AO integrals done'
c jbecca START ------------- DIM/QM -----------------
c For DIM/QM, we need to add in local field contributions here, if any
c so that they affect the linear response of the system. It is important
c to do here, because they are calculated in AO basis. Adding in here
c will allow for this routine to handle converting to MO basis. 
c NOTE: The following code assumes that the user will supply a damping
c        parameter if using a frequency dependent local field

      if (ldimqm .and. lrsp .and. lclfld) then 
         !write(luout,*)'calling dimqm_lclfld for resp.eq.'
         call dimqm_lclfld(g_dipole1,omega,lifetime,g_dipel_i)
      endif

c jbecca END   ------------- DIM/QM -----------------
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs.
c     (the rhs is NOT divided by (e_a - e_i -/+ omega), this
c     will be considered in the CPKS solver, in the precon-
c     ditioner and the 1e part of the "product" routine)
c --------- g_dipole1 --> g_dipole --------- START
      call ga_zero(g_dipole)       
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 1
       bhi(3) = ndir
      do ispin=1,npol  
       disp=ndir*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       alo(3) = disp+1
       ahi(3) = disp+ndir  
       call nga_copy_patch('n',g_dipole1,blo,bhi,
     &                         g_dipole ,alo,ahi) 
c  jbecca START: DIM/QM
       if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
         call nga_copy_patch('n',g_dipel_i,blo,bhi,
     &                         g_dipel_i_tmp ,alo,ahi) 
       endif
c  jbecca END
      enddo ! end-loop-ispin
c --------- g_dipole1 --> g_dipole --------- END
      call giao_aotomo(g_dipole,vectors,nocc,nvirt,npol,ndir,nbf)
c jbecca START: DIM/QM
      if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
        call giao_aotomo(g_dipel_i_tmp,vectors,nocc,nvirt,npol,ndir,nbf)
      endif
c jbecca END: DIM/QM
      
      do ispin=1,npol
       shift=ndir*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+ndir
       disp=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = ndir
       blo(3) = 1
       bhi(3) = 1
       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_dipole-nw-------- START'
        call ga_print(g_dipole)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_dipole-nw-------- END'
       endif ! end-if-debug
       do ipm = 1,ncomp
        call nga_add_patch(1.0d0,g_rhs(ipm),blo,bhi,
     &                     1.0d0,g_dipole  ,alo,ahi,
     &                           g_rhs(ipm),blo,bhi)
c  jbecca START DIM/QM
        if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
            call nga_add_patch(1.0d0,g_rhs_im(ipm),blo,bhi,
     &                         1.0d0,g_dipel_i_tmp  ,alo,ahi,
     &                           g_rhs_im(ipm),blo,bhi)
        endif
c  jbecca END DIM/QM
       enddo ! end-loop-ipm
      enddo ! end-loop-ispin
c jbecca start DIMQM
      if (ldimqm .and. lrsp .and. lclfld) call dimqm_tidy()
      if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
         if (.not.ga_destroy(g_dipel_i_tmp)) call errquit
     $      ('fiao_f1: ga_destroy failed g_dipel_i_tmp',0,GA_ERR)
         if (.not.ga_destroy(g_dipel_i)) call errquit
     $      ('fiao_f1: ga_destroy failed g_dipel_i',0,GA_ERR)
      endif
c jbecca END
      if (debug) write (luout,*) 'fiao_f1: dipole added to rhs'
c     Cleanup g_dipole as we do not need it right now
      if (.not.ga_destroy(g_dipole)) call 
     &    errquit('fiao_f1: ga_destroy failed g_dipole',0,GA_ERR)
      if (.not.ga_destroy(g_dipole1)) call 
     &    errquit('fiao_f1: ga_destroy failed g_dipole',0,GA_ERR)
      return
      end

      subroutine get_nocc(rtdb,   ! in : rtdb handle
     &                    nocc,   ! out: nr occupations
     &                    npol,   ! out: nr of polarization
     &                    nclosed,! in : nr closed shells
     &                    nopen,  ! in : nr open shells
     &                    nvirt,  ! in : nr virtual MOs
     &                    scftyp, ! in : string = UHF or RHF
     &                    ntot)   ! out: sum_{i,npol} nocc(i)*nvirt(i)
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c                

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "prop.fh"
      character*3 scftyp
      integer rtdb,ispin,
     &        npol,nocc(2),nclosed(2),
     &        nopen(2),nvirt(2),ntot
        if      (scftyp .eq. 'UHF') then
         npol=2
         nocc(1)=nopen(1)   
         nocc(2)=nopen(2)    
c ------ Store nopen in rtdb so that CPHF routine is happy ---- START
c In scf_init(): to avoid error message: 
c          ===>  scf: no. of closed-shell electrons is not even!
          if (.not. rtdb_put(rtdb, 'scf:nopen', 
     &         MT_INT, 1, nocc(1)-nocc(2)))
     *         call errquit('get_nocc:rtdbput nopen failed',
     &         nocc(1)-nocc(2),
     &       RTDB_ERR)
c ------ Store nopen in rtdb so that CPHF routine is happy ---- END
        else if (scftyp .eq. 'RHF') then
         npol=1
         nocc(1)=nclosed(1)
         nocc(2)=0
        endif
        ntot=0
        do ispin=1,npol
         ntot=ntot+nocc(ispin)*nvirt(ispin)
        enddo
      return
      end
c $Id$
