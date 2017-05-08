c
czora...Scale zora eigenvalues and energy
          subroutine dft_zora_scale(
     &                   rtdb,g_dens_at,nexc, ! FA
     &                   geom, 
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   g_dens,
     &                   g_s,
     &                   g_movecs,
     &                   g_zora_scale_sf,
     &                   evals,
     &                   focc,
     &                   noc,
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
#include "geom.fh" 
      integer rtdb,geom,ao_bas_han ! handles
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external get_rhoS,print_EFGZ4_version,
     &         get_EPRg,print_EPRg_version,
     &         get_NMRCS_SRZORA,print_NMRCS_SRZORA_version,
     &         get_NMRHFine_ZORA,print_NMRHypFine_version,
     &         get_densZ4,
     &         dft_zorashield_read,dft_zorashield_write      
      integer dft_zorashield_read,dft_zorashield_write
      integer g_dens(2),g_movecs(2),g_zora_scale_sf(2)
      integer g_s,ga_Fji,g_orb,g_dens_sf,g_zora_scr(2),
     &        g_dens_at(2) 
      integer l_vecs,k_vecs
      integer noc(2),iorb,ispin,ipol,i,nexc,noc1
      integer nbf,nbf_ao,nat_slc  
      double precision eval_scal,ener_scal,scf_dbl
      double precision focc(nbf*ipol)  ! occupation no.
      double precision evals(ipol*nbf) ! eigenvalues                     
      double precision zora_eint,zora_denom,
     &                 scale,scale_MO,scale_MO1  
      logical status,done_Fji,do_zgc_old    
      integer nospin,nogshift,nogiao,noelfgz4,efgfile,switch_focc
      data nospin,nogshift,nogiao,noelfgZ4/1,1,1,1/                       
      data switch_focc/0/ ! =0 not using occ keyword =1 using occ keyword
      data efgfile/0/     ! =0 no NLMO analysis      =1 doing NLMO analysis
      character*255 zorafilename
      integer dims(3),chunk(3) ! for g_Cifull(i) i=1,ipol
      integer g_densZ4(3),ind_occ
      logical skip_csAOev,skip_gshiftAOev,
     &        skip_hypAOev,skip_efgz4AOev
c ---------- nlmo definitions ----------------- START
      integer g_munuV6,g_munu_rhoS,g_densx,       
     &        g_zora_scale_munu(2),ipolmunu 
      integer ndir
      logical dft_zoraEFGZ4_NLMOAnalysis_write
      external dft_zoraEFGZ4_NLMOAnalysis_write
     &         util_file_name
c ---------- nlmo definitions ----------------- END
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
      status=rtdb_get(rtdb,'focc:occ-switch',mt_int,1,switch_focc)  ! check-occupations-keyword
      status=rtdb_get(rtdb,'prop:efieldgradZ4',MT_INT,1,noelfgZ4)   
      status=rtdb_get(rtdb,'prop:giao'        ,MT_INT,1,nogiao)     
      status=rtdb_get(rtdb,'prop:gshift'      ,MT_INT,1,nogshift)   
      status=rtdb_get(rtdb,'prop:hyperfine'   ,MT_INT,1,nospin)     
      do_zgc_old=do_zora_get_correction                             
      do_zora_get_correction=.true. ! -1  ! FORCE-ZORA-CALC-INTS
c
      if (ipol.eq.1) then
        scf_dbl=2.0d00
        noc1=noc(1)
      else if (ipol.eq.2) then 
        scf_dbl=1.0d00   
        noc1=noc(1)+noc(2)    
      endif
c
      if (noelfgZ4.eq.0 .or. nospin.eq.0 .or. 
     &    nogshift.eq.0 .or. nogiao.eq.0) then
        if (ga_nodeid().eq.0) then
          write(LuOut,*)
          call util_print_centered(LuOut,
     $   'Commencing ZORA Property Calculations', 23, .true.)
          write(LuOut,*)
        end if
      end if
c
c      if (nospin.eq.0 .and. ga_nodeid().eq.0)
c     &   call print_NMRHypFine_version()
c      if (nogiao.eq.0 .and. ga_nodeid().eq.0)
c     &   call print_NMRCS_SRZORA_version()
c      if (nogshift.eq.0 .and. ga_nodeid().eq.0)
c     &   call print_EPRg_version()
c      if (noelfgZ4.eq.0 .and. ga_nodeid().eq.0) 
c     &   call print_EFGZ4_version()
c
      if (noelfgZ4.eq.0 .or. nospin.eq.0 .or. 
     &    nogshift.eq.0 .or. nogiao.eq.0) then
       if (.not.ga_create(MT_DBL,ipol,noc1,'dft_zora_scale: g_Ci',
     $       0,0,g_Ci)) call errquit('dft_zora_scale: g_Ci',0,
     &       GA_ERR)
       call ga_zero(g_Ci)
      endif
c
      if (nospin  .eq.0 .or. nogiao  .eq.0 .or.
     &    nogshift.eq.0 .or. noelfgZ4.eq.0) then
c ------- create g_Cifull for scaling MO vectors
       dims(1) =nbf
       chunk(1)=nbf 
       do ispin=1,ipol
        if (.not. nga_create(mt_dbl,1,dims,"g_Cifull",
     &       chunk,g_Cifull(ispin)))
     $       call errquit('dft_zora_scale: g_Cifull', 0,
     &                    GA_ERR)
        call ga_fill(g_Cifull(ispin),1.0d0)
       enddo ! end-loop-ispin
      endif   
c -- execute this code for EFGSRZ4+NLMO analysis is requested -- START
        if (noelfgZ4.eq.0 .and. efgfile.eq.1) then ! -- if-noelfgZ4.eq.0---- START
         ipolmunu=ipol
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'hnd_elfcon_symm: g_munu',0,0,g_zora_scale_munu(1)))
     $    call errquit('hnd_elfcon_symm:',0,GA_ERR)
         call ga_copy(g_zora_scale_sf(1),g_zora_scale_munu(1))
         if(ipol.gt.1) then
          if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'hnd_elfcon_symm: g_munu',0,0,g_zora_scale_munu(2)))
     $    call errquit('hnd_elfcon_symm:',0,GA_ERR)
          call ga_copy(g_zora_scale_sf(2),g_zora_scale_munu(2))
         endif
c ----- Store efgz4 data in a file ------- START
c       Note.- lbl_nlmo defined in zora.fh
        ndir=6
        call util_file_name(lbl_nlmo,.false.,.false.,zorafilename)
c ---------> Write NMLO analysis data: 1 of 3 ----- START
        if (.not. dft_zoraEFGZ4_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &              nlist, ! in: list of selected atoms 
     &                  1, ! in: writing order =1,2,3
     &           ipolmunu, ! in: write for ndada=1
     &  g_zora_scale_munu, ! in: write for ndada=1
     &        g_munu_rhoS, ! in: write for ndata=2
     &            g_densx, ! in: write for ndata=2
     &           g_munuV6))
     &     call errquit('get_rhoS: dft_zoraNLMO_write failed',
     &                  0,DISK_ERR)
c ---------> Write NMLO analysis data: 1 of 3 ----- END
        if (.not. ga_destroy(g_zora_scale_munu(1))) call errquit(
     &    'dft_zora_utils: ga_destroy failed ',0, GA_ERR)
         if (ipolmunu.gt.1) then
         if (.not. ga_destroy(g_zora_scale_munu(2))) call errquit(
     &    'dft_zora_utils: ga_destroy failed ',0, GA_ERR)
         endif
        endif ! -------- if-noelfgZ4.eq.0---- END
c -- execute this code for EFGSRZ4+NLMO analysis is requested -- END
      if (.not. ma_alloc_get(mt_dbl,nbf_ao,'vec aux',l_vecs, k_vecs))
     &    call errquit(
     &          'dft_zora_utils: cannot allocate vec aux',911,MA_ERR)
      g_orb     = ga_create_atom_blocked(geom, ao_bas_han, 'orbs')
      g_dens_sf = ga_create_atom_blocked(geom, ao_bas_han,'g_dens_sf')
      do i=1,3 ! for scaled density matrix used in NMR-shieldigns,hyperfine,gshift,EFGZ4
       g_densZ4(i)=ga_create_atom_blocked(geom,ao_bas_han,'densZ4')
       call ga_zero(g_densZ4(i))
      end do
c
c     scale the eigenvalues and energy
      ener_scal = 0.d0
      do ispin = 1,ipol
       do iorb = 1,nbf
        call ga_get(g_movecs(ispin), 1, nbf, iorb, iorb,
     &              dbl_mb(k_vecs), 1)
        call ga_zero(g_orb)
        call ga_put(g_orb, 1, nbf, iorb, iorb,
     &              dbl_mb(k_vecs), 1)
        call ga_zero(g_dens_sf)
        call ga_dgemm('n','t',nbf,nbf,noc(ispin),
     &                1.0d00,g_orb,g_orb,
     &                0.0d00,g_dens_sf)
        eval_scal = evals(iorb)
        if (ispin.gt.1) eval_scal = evals(iorb + nbf)
        zora_eint = ga_ddot(g_dens_sf,g_zora_scale_sf(ispin)) ! FA
        if ((do_NonRel) .or. (not_zora_scale)) then
c +++++++++ Create non-scaled density matrix +++++++ START
         zora_denom=1.0d0
         ener_scal = ener_scal-eval_scal ! for skipping scaling
         if (iorb .le. noc(ispin)) then        
          scale_MO=scf_dbl
          if (switch_focc.eq.1) then ! adding focc()
            ind_occ=nbf*(ispin-1)+iorb
            scale_MO=scf_dbl*focc(ind_occ)
          endif
          call ga_add(1.0d0   ,g_densZ4(ispin),
     &                scale_MO,g_dens_sf,
     &                g_densZ4(ispin))
         endif
c +++++++++ Create non-scaled density matrix +++++++ END
        else
c +++++++++ Create scaled density matrix +++++++ START
         zora_denom= 1.d0+zora_eint                            
         eval_scal = eval_scal/zora_denom                     
         ener_scal = ener_scal-eval_scal*zora_eint
         if (iorb .le. noc(ispin)) then
          scale=1.0d0/zora_denom
          call ga_scale(g_dens_sf,scale)
          scale_MO=scf_dbl
          if (switch_focc.eq.1) then ! adding focc()
            ind_occ=nbf*(ispin-1)+iorb
            scale_MO=scf_dbl*focc(ind_occ)
          endif
          call ga_add(1.0d0   ,g_densZ4(ispin),
     &                scale_MO,g_dens_sf,
     &                g_densZ4(ispin))
         endif
c +++++++++ Create scaled density matrix +++++++ END
        endif
        if (ispin.le.1) then
            evals(iorb)     = eval_scal
        else
            evals(iorb+nbf) = eval_scal
        end if
c ----- execute this code ONLY when EFGZ4/NMRCSZ4/EPRgshift is requested -- START
         if(noelfgZ4.eq.0 .or. nospin.eq.0 .or. 
     &      nogshift.eq.0 .or. nogiao.eq.0) then
          if ((ipol.eq.1 .and. iorb.le.noc(1)).or.
     &        (ipol.eq.2 .and.
     &      (ispin.eq.1 .and. iorb.le.noc(1)).or.
     &      (ispin.eq.2 .and. iorb.le.noc(2)))) then
           call ga_fill_patch(g_Ci,ispin,ispin,iorb,iorb,
     &                        zora_denom)
          endif
         endif
         if(nogiao  .eq.0 .or.
     &      nogshift.eq.0 .or. nospin.eq.0) then
          if ((ipol.eq.1 .and. iorb.le.noc(1)).or.
     &        (ipol.eq.2 .and.
     &      (ispin.eq.1 .and. iorb.le.noc(1)).or.
     &      (ispin.eq.2 .and. iorb.le.noc(2)))) then
            if (nogiao.eq.0 .or. nogshift.eq.0 .or.
     &          nospin .eq. 0) then
             scale_MO1=1.0d0/zora_denom
             if (switch_focc.eq.1) then ! adding focc()
              ind_occ=nbf*(ispin-1)+iorb
              scale_MO1=1.0d0/zora_denom*focc(ind_occ)
c              if (ga_nodeid().eq.0) then
c               write(*,5) 1.0d0/zora_denom,ind_occ,focc(ind_occ)
c5              format('In dft_zora_utils:: (scl_MO,ind_occ,focc)=(',
c    &                  f15.8,',',i4,',',f15.8,')') 
c              endif
             endif     
             call ga_fill_patch(g_Cifull(ispin),iorb,iorb,1,1,
     &                          scale_MO1)
            endif
          end if
         endif
c ----- execute this code ONLY when EFGZ4/NMRCSZ4/EPRgshift is requested -- END
       end do ! orbital loop
      end do ! polarization loop
c--- Create total scaled-density matrix ---- START
      call ga_zero(g_densZ4(3))
      call ga_copy(g_densZ4(1),g_densZ4(3))
      if (ipol.gt.1) then 
       call ga_add(1.0d00,g_densZ4(1),
     &             1.0d00,g_densZ4(2),
     &                    g_densZ4(3))      
      end if
c--- Create total scaled-density matrix ---- END
c     double the energy for closed-shell calculations
      if (ipol.eq.1) ener_scal = 2.d0*ener_scal     
c      write(luout,*) "scaled_energy:",ener_scal
c ----- execute this code ONLY when NMR-CS-Z4 is requested -- START
      done_Fji=.false.
      if(.not.rtdb_get(rtdb,'zora:skip_csAOev',
     &                 mt_log,1,skip_csAOev))        
     &  skip_csAOev = .false.       
      if (nogiao.eq.0 .and. .not.(skip_csAOev)) then
       call get_NMRCS_SRZORA(rtdb,g_dens_at,nexc,
     &                       geom,ao_bas_han,
     &                       nbf,noc,ipol,
     &                       ga_Fji,
     &                       g_densZ4)
       done_Fji=.true. ! I do not need to calculate again 
                       ! if calling get_EPRg
      endif
c ----- execute this code ONLY when NMR-CS-Z4 is requested -- END
c ----- execute this code ONLY when NMR-spinspin is requested -- START
      if(.not.rtdb_get(rtdb,'zora:skip_hypAOev',
     &                 mt_log,1,skip_hypAOev))        
     &  skip_hypAOev = .false.    
      if (nospin.eq.0 .and. .not.(skip_hypAOev)) then
        call get_NMRHFine_ZORA(rtdb,g_dens_at,nexc,
     &                         geom,ao_bas_han,
     &                         nbf,focc,noc,ipol,
     &                         g_densZ4) 
      endif
c ----- execute this code ONLY when NMR-spinspin is requested -- END
c ----- execute this code ONLY when EPR-gshift is requested -- START
      skip_gshiftAOev = .false. 
      if(.not.rtdb_get(rtdb,'zora:skip_gshiftAOev',
     &                 mt_log,1,skip_gshiftAOev))        
     &  skip_gshiftAOev = .false.  
      
      if (ga_nodeid().eq.0)
     &  write(*,22) nogshift,skip_gshiftAOev,done_Fji
 22   format('(nogshift,skip_gshiftAOev,done_Fji)=(',
     &        i2,',',L1,',',L1,')')

      if (nogshift.eq.0 .and. .not.(skip_gshiftAOev)) then
       call get_EPRg(rtdb,g_dens_at,nexc,
     &               geom,ao_bas_han,
     &               nbf,focc,noc,ipol,
     &               done_Fji,ga_Fji,
     &               g_densZ4)
      endif   
c ----- execute this code ONLY when EPR-gshift is requested -- END
      if ((nogiao  .eq.0 .and. .not.(skip_csAOev)).or. 
     &    (nogshift.eq.0 .and. .not.(skip_gshiftAOev))) then
        if (.not.ga_destroy(ga_Fji)) call 
     &    errquit('dft_zora_utils: ga_destroy failed ga_Fji',
     &            0,GA_ERR)
      endif
c ----- execute this code ONLY when EFGZ4 is requested -- START
      if(.not.rtdb_get(rtdb,'zora:skip_efgz4AOev',
     &                 mt_log,1,skip_efgz4AOev))        
     &  skip_efgz4AOev = .false.  
      if (noelfgZ4.eq.0 .and. .not.(skip_efgz4AOev)) then
       call get_rhoS(rtdb,g_dens_at,nexc,
     &               geom,ao_bas_han,nbf,nbf_ao,
     &               noc,ipol)     
      endif 
c ----- execute this code ONLY when EFGZ4 is requested -- END
c     deallocate memory
      if (.not.ma_free_heap(l_vecs)) call errquit
     &   ('dft_zora_utils, ma_free_heap of l_vecs failed',911,MA_ERR)
      if (.not. ga_destroy(g_orb)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_dens_sf)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)
      do i=1,3
        if (.not. ga_destroy(g_densZ4(i))) call errquit(
     &    'dft_zora_utils: ga_destroy failed ',0, GA_ERR)
      enddo
c
      do_zora_get_correction=do_zgc_old   ! restore value
c
      return
      end
c
czora...Inquire if the zora correction file is present
c
      logical function dft_zora_inquire_file(filename)
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
      dft_zora_inquire_file = found 
c
      call ga_sync()
c
      return
      end
c
czora...Read in the zora atomic corrections from disk
c
      logical function dft_zora_read(filename, nbf, nsets, nmo, 
     &               mult, g_zora_sf, g_zora_scale_sf)
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
      integer nsets             ! restricted or unrestricted
      character*(*) filename
      integer iset              ! restricted or unrestricted
      integer g_zora_sf(nsets)
      integer g_zora_scale_sf(nsets)
c
      integer nbf               ! No. of functions in basis
      integer nmo(nsets)       
      integer mult              ! multiplicity
      integer ok, jset, i, j
      integer l_zora, k_zora
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c
      integer nsets_read       ! No. of sets  from file
      integer nbf_read         ! No. of functions from file 
      integer mult_read
c     
c     Initialise to invalid MA handle
      l_zora = -1
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read atomic ZORA corrections from ',a/)
       call util_flush(luout)
c
c      Open the file
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
c      Allocate the temporary buffer
       if (.not. ma_push_get(mt_dbl,nbf,'dft_zora_read',l_zora,k_zora))
     $        call errquit('dft_zora_read: ma failed', nbf, MA_ERR)
c
c      Read in g_zora_sf
       do iset = 1, nsets
        do i = 1, nmo(iset)
         call sread(unitno, dbl_mb(k_zora), nbf)
         call ga_put(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
        end do
       end do
c
c      Read in g_zora_scale_sf
       do iset = 1, nsets
        do i = 1, nmo(iset)
         call sread(unitno, dbl_mb(k_zora), nbf)
         call ga_put(g_zora_scale_sf(iset), 1, nbf, i, i, 
     & dbl_mb(k_zora), 1)
        end do
       end do
c
c      Close the file
       close(unitno,err=1002)
       ok = 1
c
c      Deallocate the temporary buffer
       if (.not. ma_pop_stack(l_zora)) call errquit
     $      ('dft_zora_read: pop failed', l_zora, MA_ERR)
c
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
c       write(6,*)' g_zora_sf(1) from dft_scf'
c       call ga_print(g_zora_sf(1))
c       write(6,*)' g_zora_scale_sf(1) from dft_scf'
c       call ga_print(g_zora_scale_sf(1))
c
      dft_zora_read = ok .eq. 1
      return
c
 1000 write(6,*) 'dft_zora_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1003 write(6,*) 'dft_zora_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
czora...Write out the zora atomic corrections to disk
c
      logical function dft_zora_write(rtdb, basis, filename,
     &     nbf, nsets, nmo, mult, g_zora_sf, g_zora_scale_sf)
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
      integer g_zora_sf(nsets)    
      integer g_zora_scale_sf(nsets)  
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
      if (.not. ma_push_get(mt_dbl,nbf,'dft_zora_write',l_zora,k_zora))
     $        call errquit('dft_zora_write: ma failed', nbf, MA_ERR)
c
c     Write out g_zora_sf
      do iset = 1, nsets
       do i = 1, nmo(iset)
        call ga_get(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_sf
      do iset = 1, nsets
       do i = 1, nmo(iset)
        call ga_get(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &            dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora))
     $  call errquit('dft_zora_write: ma pop failed', l_zora, MA_ERR)
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
c       write(6,*)' g_zora_sf(1) from dft_scf'
c       call ga_print(g_zora_sf(1))
c       write(6,*)' g_zora_scale_sf(1) from dft_scf'
c       call ga_print(g_zora_scale_sf(1))
c
      dft_zora_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote atomic ZORA corrections to ',a/)
         call util_flush(luout)
      endif
c
      return
c
 1000 write(6,*) 'dft_zora_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
c     Convergence check for zora calculations
c
      subroutine dft_zora_scfcvg(rms, derr, etold, etnew, e_conv, 
     &                      d_conv, g_conv, ipol, iter, iterations, 
     &                      idone, rtdb, converged, diising)
c
c     $Id$
c
      implicit none
#include "errquit.fh"
c
      double precision rms(2)   ! [input]
      double precision derr(2)  ! [input]
      double precision etold    ! [input]
      double precision etnew    ! [input]
      double precision e_conv   ! [input]
      double precision d_conv   ! [input]
      double precision g_conv   ! [input]
      integer ipol              ! [input]
      integer iter              ! [input]
      integer iterations        ! [input]
      integer idone             ! [output]
      integer rtdb              ! [input]
      logical converged         ! [output]
      logical diising           ! [input]
c
#include "mafdecls.fh"
#include "rtdb.fh"
c     
      logical e_conv_logical, d_conv_logical, g_conv_logical
      logical ENERGY, DENSITY, GRADIENT
      double precision de, abde
c
      converged = .false.
c
      e_conv_logical = .false.
      d_conv_logical = .false.
      g_conv_logical = .false.
c
      ENERGY = e_conv.gt.0
      DENSITY = d_conv.gt.0
      GRADIENT = g_conv.gt.0
c
      idone = 0
c
c     Evaluate change in energy.
c     
      de = etnew - etold
      etold=etnew
      abde = dabs(de)
c     
c     Check to see if energy is converged.
c
      if (ENERGY)then     
         if (abde.lt.e_conv)e_conv_logical = .true.
      else
         e_conv_logical = .true.
      endif
c     
c     Check for density matrix convergence.
c     
      if (DENSITY)then
         if (dsqrt(rms(1)).le.d_conv)d_conv_logical = .true.
         if (ipol.eq.2)then
           if (dsqrt(rms(2)).le.d_conv) then
             d_conv_logical = d_conv_logical.and..true.
          else
             d_conv_logical = d_conv_logical.and..false.
          endif
         endif
      else
         d_conv_logical = .true.
      endif
c     
c     Check for gradient convergence.
c     
      if (GRADIENT.and.diising)then
         if (derr(1).le.g_conv)g_conv_logical = .true.
         if (ipol.eq.2)then
           if (derr(2).le.g_conv) then
              g_conv_logical = g_conv_logical.and..true.
             else
                g_conv_logical = g_conv_logical.and..false.
             endif
         endif
      else
         g_conv_logical = .true.
      endif
c
c     Check over-all convergence.
c
      converged = e_conv_logical
      if (converged)idone = 1
c     
c     Check iteration value.
c     
      if (iter.lt.1)then
         return
      elseif (iter.eq.iterations)then
         idone = 1
      endif
c
c     If all convergence criterion met or number of iterations has been
c     exceeded, write "converged" to RTDB.
c
      if (idone.eq.1)then
         if (.not.rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &      call errquit('dft_scfcvg: rtdb_put failed', 1, RTDB_ERR)
      endif
      return
      end
