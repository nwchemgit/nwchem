czora...
czora...Scale zora eigenvalues and energy
czora...
          subroutine dft_zora_scale_so(
     &                   rtdb,
     &                   g_dens_ini,nexc,
     &                   geom, 
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   nbf_mo,
     &                   g_dens,
     &                   g_s,
     &                   g_moso,
     &                   g_zora_scale_sf,
     &                   g_zora_scale_so,
     &                   evals,
     &                   focc,
     &                   noc,
     &                   nocc, ! FA
     &                   ipol,
     &                   ener_scal)

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "zora.fh"
#include "rtdb.fh" 

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external get_rhoS_so         ! FA
      external hnd_efgmap_Z4_so    ! FA
      external print_EFGZ4_version ! FA
      integer g_dens(2)
      integer g_moso(2)
      integer g_orb(2)
      integer g_dens_so(2)
      integer g_scr(2)
      integer g_s
      integer g_zora_scale_sf(2)
      integer g_zora_scale_so(3)
      double precision numelecs, ener_kin
      integer l_vecsre, k_vecsre
      integer l_vecsim, k_vecsim
      integer iorb
      double precision eval_scal
      double precision ener_scal
      double precision ener_sf
      double precision ener_so_x
      double precision ener_so_y
      double precision ener_so_z
      double precision ener_tot
      integer noc
      integer ispin
      integer ipol

      integer geom
      integer ao_bas_han
      integer nbf
      integer nbf_ao
      integer nbf_mo

      double precision focc(*)  ! occupation no.
      double precision evals(*) ! eigenvalues
      integer do_replnxn,mem_replnxn,val_in
      integer l_evals_in, k_evals_in,
     &          l_zora_scale_sf, k_zora_scale_sf,
     &          l_zora_scale_so, k_zora_scale_so,
     &          l_dens_so, k_dens_so,
     &          l_vtmp, k_vtmp
      double precision zora_denom          ! Added by FA
      integer rtdb                         ! Added by FA
      integer nexc,g_dens_ini(2)           ! Added by FA
      integer j,j1,iorb1,nocc(2)           ! Added by FA
      logical do_zgc_old                   ! Added by FA
      integer ac_AB(2),iter                ! Added by FA
      character*1  soxyz(3)                ! Added by FA
      double precision ener_sf_1,ener_sf_2 ! Added by FA
      logical status                       ! Added by FA 
      integer noelfgz4                     ! Added by FA 
      data noelfgZ4/1/                     ! Added by FA 
      soxyz(1)='z'
      soxyz(2)='y'
      soxyz(3)='x'
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
      status=rtdb_get(rtdb,'prop:efieldgradZ4',MT_INT,1,noelfgZ4) ! FA
      do_zgc_old=do_zora_get_correction
      do_zora_get_correction= .true. !-1  ! FORCE-ZORA-CALC-INTS
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
c ----- execute this code ONLY when EFGZ4 is requested -- START
      if(noelfgZ4.eq.0) then
       if (ga_nodeid().eq.0) call print_EFGZ4_version()
       if (.not.ga_create(MT_DBL,2,nocc(1)+nocc(2),      ! FA
     &      'dft_zora_scale: g_Ci',0,0,g_Ci))            ! FA
     &  call errquit('dft_zora_scale_so: g_Ci',0,GA_ERR) ! FA
       call ga_zero(g_Ci)                                ! FA
      endif
c ----- execute this code ONLY when EFGZ4 is requested -- END
c     allocate memory

      if (.not.MA_Push_Get(MT_Dbl, 2*nbf_mo, 'real vec aux',
     &          l_vecsre, k_vecsre))
     & call errquit('dft_zora_scale_so: cannot allocate vec',0, MA_ERR)
c
#if 0
      if (.not.MA_Push_Get(MT_Dbl, nbf_mo, 'imag vec aux',
     &          l_vecsim, k_vecsim))
     & call errquit('dft_zora_scale_so: cannot allocate vec',0, MA_ERR)
#else
      k_vecsim=k_vecsre+nbf_mo
#endif
c 
c     spin-orbit vector - real
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'orbs re',0,0,g_orb(1)))
     & call errquit('dft_zora_scale_so: orb real error',0, GA_ERR)
      call ga_zero(g_orb(1))
c
c     spin-orbit vector - imag
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'orbs im',0,0,g_orb(2)))
     & call errquit('dft_zora_scale_so: orb imag 2 error',0, GA_ERR)
      call ga_zero(g_orb(2))
c
c     spin-orbit density matrix - real
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'denmxre',0,0,
     & g_dens_so(1)))
     & call errquit('dft_zora_scale_so: dens real error',0, GA_ERR)
      call ga_zero(g_dens_so(1))
c
c     spin-orbit density matrix - imag
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'denmxim',0,0,
     & g_dens_so(2)))
     & call errquit('dft_zora_scale_so: dens imag error',0, GA_ERR)
      call ga_zero(g_dens_so(2))
c
c     scratch array
       if(.not.ga_duplicate(g_dens(1),g_scr(1),'scratch 1'))
     &  call errquit('dft_zora_scale_so: ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_scr(1))
       if(.not.ga_duplicate(g_dens(2),g_scr(2),'scratch 2'))
     &  call errquit('dft_zora_scale_so: ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_scr(2))
      j=0 ! FA
      ac_AB(1)=1
      ac_AB(2)=1
      ener_scal = 0.d0
cwhy??      ener_tot  = 0.d0
c
c     check if enough memory is available to allocate a few nxn objects
c     
      do_replnxn=0
      mem_replnxn=
     C     nbf_mo*nbf_mo*2  +   ! dens_so
     C     nbf_ao*nbf_ao    +   ! vtmp
     C     nbf_ao*nbf_ao*2  +   ! zora_scale_sf
     C     nbf_ao*nbf_ao*3  +   ! zora_scale_so
     C     nbf_mo               ! evals_in
      if (MA_inquire_avail(mt_dbl)*8/10.gt.mem_replnxn) then
         do_replnxn=1
      else
         if(ga_nodeid().eq.0) then
            write(6,*) ' mem needed    in MB',mem_replnxn*
     / 8/1024/1024
            write(6,*) ' mem available in MB ',MA_inquire_avail(mt_dbl)*
     / 8/1024/1024
         endif
      endif
      if(noelfgZ4.eq.0) then
         write(6,*) ' need to compute elfgZ4, do_replnxn=0 '
         do_replnxn=0
      endif
      if(rtdb_get(rtdb,'zora:do_replnxn',MT_INT,1,val_in))
     c    do_replnxn=val_in
      call ga_igop(19481984, do_replnxn, 1, 'min')
      if(do_replnxn.eq.0) then
         if(ga_nodeid().eq.0)
     c        write(6,*) ' zoraso: NOT enough mem for repl'
      do iorb = 1, noc
       call ga_get(g_moso(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)
       call ga_zero(g_orb(1))
       call  ga_put(g_orb(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)
       call ga_get(g_moso(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1)
       call ga_zero(g_orb(2))
       call  ga_put(g_orb(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1)
       call dft_densm_so(g_dens_so, g_orb, nbf_ao, noc)
       call ga_zero(g_scr(1))
       call ga_zero(g_scr(2))
       call ga_dens_sf(g_scr, g_dens_so(1), nbf_ao)
        ener_sf_1=ga_ddot(g_scr(1),g_zora_scale_sf(1))
        ener_sf_2=ga_ddot(g_scr(2),g_zora_scale_sf(2))
        ener_sf=ener_sf_1+ener_sf_2
        ener_tot=ener_sf
        do iter=1,3
         call ga_zero(g_scr(1))
         call ga_dens_so(g_scr(1),g_dens_so,nbf_ao,soxyz(iter))
         ener_tot=ener_tot+
     &            ga_ddot(g_scr(1),g_zora_scale_so(iter)) 
        end do
       zora_denom=1.d0+ener_tot
       eval_scal = evals(iorb)
       eval_scal = eval_scal/zora_denom
       ener_scal = ener_scal - eval_scal*ener_tot 
       evals(iorb) = eval_scal
c ----- execute this code ONLY when EFGZ4 is requested -- START
        if(noelfgZ4.eq.0) then
         j1=j+1                          ! FA
         if (ac_AB(j1).le.nocc(j1)) then ! FA
          call ga_fill_patch(g_Ci,j1,j1,ac_AB(j1),ac_AB(j1),
     &                       zora_denom) ! FA
          ac_AB(j1)=ac_AB(j1)+1          ! FA
         end if                          ! FA
        j=1-j                            ! FA
        endif
c ----- execute this code ONLY when EFGZ4 is requested -- END
      enddo ! end-loop-iorb 
      else
         if(ga_nodeid().eq.0)
     c        write(6,*) ' zoraso: enough mem for repl'
      if (.not.ma_push_get(mt_dbl, nbf_mo, 'evals_in',
     &          l_evals_in, k_evals_in))
     & call errquit('dft_zora_scale_so: cannot allocate evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl, nbf_ao*nbf_ao*2, 'zora_scale_sf',
     &          l_zora_scale_sf, k_zora_scale_sf))
     & call errquit('dft_zora_scale_so:cant all zora_scale_sf',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl, nbf_ao*nbf_ao*3, 'zora_scale_so',
     &          l_zora_scale_so, k_zora_scale_so))
     & call errquit('dft_zora_scale_so:cant all zora_scale_so',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl, nbf_mo*nbf_mo*2, 'dens_so',
     &          l_dens_so, k_dens_so))
     & call errquit('dft_zora_scale_so:cant all dens_so',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl, nbf_ao*nbf_ao, 'vtmp',
     &          l_vtmp, k_vtmp))
     & call errquit('dft_zora_scale_so:cant all vtmp',0,MA_ERR)
c
c     replicate zora_scale_sx
c
      call util_mygabcast(g_zora_scale_sf(1), 
     M  nbf_ao,nbf_ao, dbl_mb(k_zora_scale_sf),nbf_ao)
      call util_mygabcast(g_zora_scale_sf(2), 
     M  nbf_ao,nbf_ao, dbl_mb(k_zora_scale_sf+nbf_ao*nbf_ao),nbf_ao)

      call util_mygabcast(g_zora_scale_so(1), 
     M  nbf_ao,nbf_ao, dbl_mb(k_zora_scale_so),nbf_ao)
      call util_mygabcast(g_zora_scale_so(2), 
     M  nbf_ao,nbf_ao, dbl_mb(k_zora_scale_so+nbf_ao*nbf_ao),nbf_ao)
      call util_mygabcast(g_zora_scale_so(3), 
     M  nbf_ao,nbf_ao, dbl_mb(k_zora_scale_so+2*nbf_ao*nbf_ao),nbf_ao)

c     copy input evals
      call dcopy(nbf_mo,evals,1,dbl_mb(k_evals_in),1)
c     zero evals [output] to get dgop working
      call dcopy(noc,0d0,0,evals,1)
      do iorb = ga_nodeid()+1, noc, ga_nnodes()
      call dft_zora_nogaga(iorb,nbf_ao,nbf_mo,noc,g_moso,
     C     soxyz,
     D     dbl_mb(k_zora_scale_sf),dbl_mb(k_zora_scale_so),
     E     ener_tot,
     D     dbl_mb(k_dens_so),
     O     dbl_mb(k_vtmp),dbl_mb(k_vecsre))
       zora_denom=1.d0+ener_tot
       eval_scal = dbl_mb(k_evals_in+iorb-1)/zora_denom
       ener_scal = ener_scal - eval_scal*ener_tot 
       evals(iorb) = eval_scal
#if 0
NOT WORKING YET
c ----- execute this code ONLY when EFGZ4 is requested -- START
        if(noelfgZ4.eq.0) then
         j1=j+1                          ! FA
         if (ac_AB(j1).le.nocc(j1)) then ! FA
          call ga_fill_patch(g_Ci,j1,j1,ac_AB(j1),ac_AB(j1),
     &                       zora_denom) ! FA
          ac_AB(j1)=ac_AB(j1)+1          ! FA
         end if                          ! FA
        j=1-j                            ! FA
        endif
#endif
c ----- execute this code ONLY when EFGZ4 is requested -- END
      enddo ! end-loop-iorb 
c dgop on evals
      call ga_dgop(20166102, evals, noc, '+')
      call ga_dgop(20166103, ener_scal, 1, '+')
      endif
c ----- execute this code ONLY when EFGZ4 is requested -- START
       if(noelfgZ4.eq.0) then
c       if (ga_nodeid().eq.0) write(*,*) "---g_Ci_so -- START"
c       call ga_print(g_Ci)
c       if (ga_nodeid().eq.0) write(*,*) "---g_Ci_so -- END"
c       write(luout,*) "ener_scal:",ener_scal
        call get_rhoS_so(rtdb,g_dens_ini,nexc,
     &                   geom,ao_bas_han,
     &                   nbf,nbf_ao,nbf_mo,
     &                   g_dens,g_moso,
     &                   noc,nocc,ipol)
c -------- computing EFGs --------- START
        call hnd_efgmap_Z4_so(rtdb,ao_bas_han,geom,
     &                        nbf,nbf_ao,nbf_mo,
     &                        g_dens,g_moso,
     &                        noc,nocc)
c -------- computing EFGs --------- END
       endif
c ----- execute this code ONLY when EFGZ4 is requested -- START

c     deallocate memory
#if 0
      if (.not. ma_chop_stack(l_vecsim)) 
     & call errquit('dft_zora_scale_so:l_vecsim', 0, MA_ERR)
#endif
      if (.not. ma_chop_stack(l_vecsre))
     & call errquit('dft_zora_scale_so:l_vecsre', 0, MA_ERR)

      if (.not. ga_destroy(g_orb(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_orb(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)

      if (.not. ga_destroy(g_dens_so(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_dens_so(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)

      if (.not. ga_destroy(g_scr(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_scr(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
c +++++++++++++++ RESTORE VALUE +++++++++++++START
      do_zora_get_correction=do_zgc_old
c +++++++++++++++ RESTORE VALUE +++++++++++++END  
      return
      end
c
czora...Inquire if the zora correction file is present
c
      logical function dft_zora_inquire_file_so(filename)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character*(*) filename
      logical found
c
      call ga_sync()
c
c     Inquire if file is present
      inquire(file=filename,exist=found)
      dft_zora_inquire_file_so = found
c
      call ga_sync()
c
      return
      end
c
czora...Read in the zora atomic corrections from disk
c
      logical function dft_zora_read_so(filename, nbf, nsets, nmo, 
     &               mult, g_zora_sf, g_zora_scale_sf, 
     &               g_zora_so, g_zora_scale_so) 
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character*(*) filename
      integer iset              ! restricted or unrestricted
c
      integer g_zora_sf(2)
      integer g_zora_scale_sf(2)
      integer g_zora_so(3)
      integer g_zora_scale_so(3)
c
      integer nsets             ! restricted or unrestricted
      integer nbf               ! No. of functions in basis
      integer nmo(nsets)       
      integer mult
      integer ok, jset, i, j
      integer l_zora, k_zora
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c     
      integer nsets_read
      integer nbf_read
      integer mult_read
c
      l_zora = -1               ! An invalid MA handle
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c
c     Print a message indicating the file being read
      write(6,22) filename(1:inp_strlen(filename))
 22   format(/' Read atomic ZORA corrections from ',a/)
      call util_flush(luout)
c
c     Open the file
      open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nsets_read
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) mult_read
c
c      Error checks
       if ((nsets_read .ne. nsets) 
     &  .or. (nbf_read .ne. nbf) 
     &  .or. (mult_read .ne. mult) ) goto 1003
c
c     Allocate the temporary buffer
      if (.not.ma_push_get(mt_dbl,nbf,'dft_zora_read_so',l_zora,k_zora))
     $        call errquit('dft_zora_read_so: ma failed', nbf, MA_ERR)
c
c     Read in g_zora_sf
      do iset = 1,2  ! 2 components
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_scale_sf
      do iset = 1,2  ! 2 components
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &   dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_so(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
       end do
      end do
c
c     Read in g_zora_scale_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1,nbf
        call sread(unitno, dbl_mb(k_zora), nbf)
        call ga_put(g_zora_scale_so(iset), 1, nbf, i, i, 
     &   dbl_mb(k_zora), 1)
       end do
      end do
c
c     Close the file
      close(unitno,err=1002)
      ok = 1
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora)) call errquit
     $      ('dft_zora_read_so: pop failed', l_zora, MA_ERR)
c
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      dft_zora_read_so = ok .eq. 1
      return
c
 1000 write(6,*) 'dft_zora_read_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_read_so: failed to read ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1003 write(6,*) 'dft_zora_read_so: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_read_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
czora...Write out the zora atomic corrections to disk
c
      logical function dft_zora_write_so(rtdb, basis, filename,
     &     nbf, nsets, nmo, mult, 
     &     g_zora_sf, g_zora_scale_sf,
     &     g_zora_so, g_zora_scale_so)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nsets             ! [input] No. of sets of matrices
      integer nmo(nsets)        ! [input] No. of mos in each set
      integer mult
c
      integer g_zora_sf(2)    
      integer g_zora_scale_sf(2)  
      integer g_zora_so(3)    
      integer g_zora_scale_so(3)  
c
      integer unitno
      parameter (unitno = 77)
c
      integer lentit
      integer lenbas
      integer l_zora, k_zora
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, key
      character*20 scftype20 
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize
c
      l_zora = -1               ! An invalid MA handle
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
      ok = 0
c
c     Read routines should be consistent with this
c
c     Write out the atomic zora corrections
c
      if (ga_nodeid() .eq. 0) then
c
c     Open the file
      open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c     Write out the number of sets and basis functions
      write(unitno, err=1001) nsets
      write(unitno, err=1001) nbf
      write(unitno, err=1001) mult
c
c     Allocate the temporary buffer
      if (.not.ma_push_get(mt_dbl,nbf,'dft_zora_write_so',
     & l_zora,k_zora))
     &   call errquit('dft_zora_write_so: ma failed', nbf, MA_ERR)
c
c     Write out g_zora_sf
      do iset = 1,2  ! 2 components
       do i = 1, nbf
        call ga_get(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_sf
      do iset = 1,2  ! 2 components
       do i = 1, nbf
        call ga_get(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &              dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1, nbf
        call ga_get(g_zora_so(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_so
      do iset = 1,3  ! 3 components: x,y,z
       do i = 1, nbf
        call ga_get(g_zora_scale_so(iset), 1, nbf, i, i, 
     &              dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora))
     $  call errquit('dft_zora_write_so: ma pop failed', l_zora, MA_ERR)
c
c     Close the file
      close(unitno,err=1002)
c
      ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      dft_zora_write_so = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote atomic ZORA corrections to ',a/)
         call util_flush(luout)
      endif
c
      return
c
 1000 write(6,*) 'dft_zora_write_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_write_so: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_write_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c $Id$
      subroutine dft_zora_nogaga(iorb,nbf_ao,nbf_mo,noc,g_moso,
     C     soxyz,
     D     zora_scale_sf,zora_scale_so,
     E     ener_tot,
     D     dens_so,
     O     vtmp,vecs)
      implicit none
      integer iorb
      integer nbf_mo,nbf_ao ! [in] nbf_mo=2*nbf_ao
      integer noc
      double precision zora_scale_sf(nbf_ao,nbf_ao,2)! replicated[in]
      double precision zora_scale_so(nbf_ao,nbf_ao,3)! replicated[in]
      integer g_moso(2) ! [in]
      double precision vecs(nbf_mo,2) ! [scratch]
      double precision dens_so(nbf_mo,nbf_mo,2) ! repl [scratch]
      double precision vtmp(nbf_ao,nbf_ao) ! repl [scratch]
      double precision ener_tot ! [out]
      character*1  soxyz(3)                

c
      integer xyz
      integer j,i
c
      external ddot
      double precision ddot
c

       call ga_get(g_moso(1),1,nbf_mo,iorb,iorb,vecs(1,1),1)
       call ga_get(g_moso(2),1,nbf_mo,iorb,iorb,vecs(1,2),1)
       call dftv_densm_so(dens_so,vecs,nbf_mo,noc)
c        ener_sf_1
c strided copy
        do j=1,nbf_ao
           call dcopy(nbf_ao,dens_so(1,j,1),1,vtmp(1,j),1)
        enddo
        ener_tot=ddot(nbf_ao*nbf_ao,
     D      vtmp,1,zora_scale_sf(1,1,1),1)
cstrided copy
        do j=1+nbf_ao,2*nbf_ao
           call dcopy(nbf_ao,dens_so(1+nbf_ao,j,1),1,vtmp(1,j-nbf_ao),1)
        enddo
c        ener_sf_2
        ener_tot=ener_tot+ddot(nbf_ao*nbf_ao,
     D      vtmp,1,zora_scale_sf(1,1,2),1)
        do xyz=1,3
         call vec_dens_so(vtmp,dens_so,nbf_ao,soxyz(xyz))
         ener_tot=ener_tot+ddot(nbf_ao*nbf_ao,
     &            vtmp,1,zora_scale_so(1,1,xyz),1) 
        end do
        return
        end
