      subroutine hnd_hyperfine_zora(rtdb,basis,geom)
c $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"
#include "zora.fh"  
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer nat,nat_slc,typeprop,
     &        ixy,ix,iy,iz,iatom,iocc,ifld,ioff
      integer alo(3),ahi(3),blo(3),bhi(3),clo(3),chi(3)
      integer dlo(3), dhi(3)
      integer l_occ,k_occ,l_eval,k_eval
      integer l_dia,k_dia,l_para,k_para,l_tot,k_tot
      integer l_xyz,k_xyz,l_zan,k_zan
      integer l_cfine,k_cfine
      integer l_AtNr,k_AtNr
      integer icalczora,type_nmrdata
      integer g_dens(3),type_NMR,
     &        g_d1, g_rhs, g_rhs0, g_u
      integer ga_dia_hfine,ga_para1,
     &        ga_h01_hfine,ga_Fji_hfine,g_AtNr1
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val
      double precision a(6),axs(3,3),eig(3),conv2MHz
      logical dft_zoraNMR_read
      character*255 zorafilename
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      integer ld(2),cbuf,ind
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
      external giao_aotomo,mat_transpose ! FA-09-16-10
      double precision val1,coeffpol,ac_occ(2),small,trace,
     &                 isotr,aniso,
     &                 span,skew,asym,rtemp
      data small /1.0d-10/
      logical     oskel, status
      data tol2e   /1.0d-10/
      integer ic,npol 
c ----- Definitions for NLMO analysis ---- START
      integer i1,j1,acc_vec,l_tvec,k_tvec,hypfile          
      integer g_tvec ! for nbo analysis
c ----- Definitions for NLMO analysis ---- END
      external get_par_hfine,get_dia_hfine,
     &         add_H10_hfine,skip_cphf,get_P10,
     &         get_convfactor_hfine,sort_asc,
     &         dft_zoraNMR_read,create_munu4nbo_hyp,
     &         get_densZ4,
     &         dft_zoraCPHF_write,dft_zoraCPHF_read

      integer ntot,ispin,indA,indB,nocc(2),nind_jk
      integer debug_giaox
      logical skip_cphf_ev_hyp

      debug_giaox=0 ! =1 for debugging print outs of matrices
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        call util_print_centered(LuOut,
     $   'Hyperfine Tensor (in au)', 23, .true.)
        write(LuOut,*)
      end if
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb   
c
c     Integral initialization
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
c     Find out from rtdb which atoms we need to calculate hyperfine for
c     Get number of atoms (all or number from rtdb)
c     Get which atoms (all or some read from rtdb)
c     Allocate arrays which will hold atomic information (k_zan and k_xyz)

      status = rtdb_parallel(.true.)
c ------- Read (nat,atmnr) --------- START
         status=geom_ncent(geom,nat)   
      if (.not.ma_alloc_get(
     &       mt_int,nat,'nmt tmp',l_AtNr,k_AtNr))
     &    call errquit('hnd_hyperfine_zora: ma failed',0,MA_ERR)
         typeprop=3 ! =1 EFG =2 Shieldings =3 Hyperfine  
         call get_slctd_atoms(nat_slc,       ! out: selected atoms
     &                        int_mb(k_AtNr),! out: list of selected atom nr.     
     &                        nat,           ! in : total nr atoms in molecule            
     &                        rtdb,          ! in : rdt  handle
     &                        typeprop)      ! in : =1,2,3=EFG,Shieldings,Hyperfine
      if (ga_nodeid().eq.0) then
       write(*,8) nat_slc
 8     format('nat_slc=',i4)
       do i=1,nat_slc
        write(*,7) i,int_mb(k_AtNr+i-1)
 7      format(' In hnd_hyperfine_zora:: atomnr(',i3,')=',i5)
       enddo
      endif
c ------- Read (nat,atmnr) --------- END
      if (.not. ma_push_get(mt_dbl,nat_slc,'nmr hfine',l_cfine,k_cfine)) 
     &    call errquit('hnd_hfine: ma_push_get failed k_cfine',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat_slc,'nmr at',l_xyz,k_xyz)) 
     &    call errquit('hnd_hfine: ma_push_get failed k_xyz',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat_slc,'nmr zan',l_zan,k_zan)) 
     &    call errquit('hnd_hfine: ma_push_get failed k_zan',0,MA_ERR)
      hypfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
      if (hypfile.eq.1) then ! ------- hypfile-if++++ START
c  Allocate memory for l_tvec,k_tvec --- start 
       if (.not.ma_alloc_get(mt_dbl,nat_slc*3*3,'tvec',l_tvec,k_tvec)) 
     &     call errquit('hnd_hfine: ma failed',0,MA_ERR) 
       call dcopy(nat_slc*3*3,0.0d0,0,dbl_mb(k_tvec),1) ! reset 
c  Allocate memory for l_tvec,k_tvec --- end
      endif ! ------------------------ hypfile-if++++ END

      do ixy = 0, nat_slc-1
         if (.not. geom_cent_get(geom, int_mb(k_AtNr+ixy), tag, 
     &                           dbl_mb(k_xyz+3*ixy),dbl_mb(k_zan+ixy)))
     &       call errquit('hnd_hfine: geom_cent_tag failed',0, GEOM_ERR)
      enddo 

      call get_convfactor_hfine(geom,           ! in  : geom handle
     &                          dbl_mb(k_cfine),! out : conversion factor for hyperfine calc.
     &                          dbl_mb(k_zan),
     &                          nat_slc,int_mb(k_AtNr)) 
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_hfine: could not get nbf',0, BASIS_ERR)
c ++++++ Reading hyperfine data from file ++++++ START
c       Note.- lbl_nmrhyp defined in zora.fh
        call util_file_name(lbl_nmrhyp,.false.,.false.,zorafilename)
      icalczora = 0  ! initialize the flag
      type_nmrdata=2 ! =1,2,3=shieldings,hyperfine,gshift       
      if (.not.dft_zoraNMR_read(zorafilename,
     &     type_nmrdata,
     &     nbf,nat_slc,
     &     g_AtNr1,
     &     ga_dia_hfine,
     &     ga_para1,
     &     ga_h01_hfine,
     &     ga_Fji_hfine)) icalczora=1 
c Note.- If I print the GAs here it gets freezed
c ++++++ Reading hyperfine data from file ++++++ END
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('hnd_hfine: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('hnd_hfine: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c ------ define npol ----- START
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, npol)) then
        if      (scftyp .eq. 'UHF') then
         npol=2
        else if (scftyp .eq. 'RHF') then
         npol=1
        endif
      endif
c ------ define npol ----- END
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'npol=',npol

      if      (npol.eq.1) then
       nocc(1)=nclosed(1)
       nocc(2)=0
      else if (npol.eq.2) then
       nocc(1)=nopen(1)   
       nocc(2)=nopen(2)     
      endif
      if (debug_giaox.eq.1) then
       write(*,10) nocc(1),nocc(2),
     &             nopen(1),nopen(2),
     &             nclosed(1),nclosed(2),
     &             nvirt(1),nvirt(2),scftyp
 10    format('nocc =(',i3,',',i3,') ',
     &        'nopen=(',i3,',',i3,') ',
     &        'nclos=(',i3,',',i3,') ',
     &        'nvirt=(',i3,',',i3,') scftyp=',a,')')
      endif
c ---- calculate total occupation alpha and beta --- START
       do i=1,npol
        val=0.0d0
        ind=nbf*(i-1)
        do j=1,nocc(i)         
         val=val+dbl_mb(k_occ+ind+j-1)
        enddo
        ac_occ(i)=val
       enddo
c       if (ga_nodeid().eq.0) then
c        write(*,2) ac_occ(1),ac_occ(2),nocc(1),nocc(2)
c 2      format('In hnd_hfine: ac_occ=(',
c     &         f15.8,',',f15.8,')  nocc=(',
c     &         i4,',',i4,')')
c       endif
c ---- calculate total occupation alpha and beta --- END
c ----- calc coeff_pol --------- START
      if (ac_occ(1) .ne. ac_occ(2)) then
       coeffpol=1.0d0/(ac_occ(1)-ac_occ(2))
      else
       write(*,1) nocc(1),nocc(2)
 1     format('Error in hnd_hyperfine_zora(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       stop
      endif 
c      if (ga_nodeid().eq.0) 
c     &  write(*,*) 'coeffpol=',coeffpol
c ----- calc coeff_pol --------- END

c ------ Store nopen in rtdb so that CPHF routine is happy ---- START
c WARNING: For restricted calc nocc(1)=9 nocc(2)=0 <=== PROBLEM
        if (npol.eq.2) then   
          if (.not. rtdb_put(rtdb, 'scf:nopen', 
     &         MT_INT, 1, nocc(1)-nocc(2)))
     *         call errquit('scfinit:rtdbput nopen failed',
     &         nocc(1)-nocc(2),
     &       RTDB_ERR)
        endif
c ------ Store nopen in rtdb so that CPHF routine is happy ---- END
c
c     Get Unperturbed Density Matrix
      call hnd_prp_get_dens(rtdb,geom,basis,
     &                       g_dens, ! out
     &                       ndens,scftyp,
     &                       nclosed,nopen,nvirt)
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
      enddo
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs))
     &   call errquit('hnd_gshift: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)

       if (debug_giaox.eq.1) then
        write(*,*) 'after-creating g_rhs ...'
        write(*,102)npol,nocc(1),nocc(2),
     &              nclosed(1),nclosed(2),
     &              nvirt(1),nvirt(2),scftyp,ndens
 102   format('BEF pre-fock::npol=',i3,' nocc=(',i3,',',i3,') ',
     &       'nclos=(',i3,',',i3,') ',
     &       'nvirt=(',i3,',',i3,') scftyp=',a,',ndens=',i3,')')
       endif
       call add_H10_hfine(g_rhs, ! out: ga_rhs(a,i)=ga_rhs(a,i)+H10(a,i)
     &                    ga_Fji_hfine,vectors, 
     &                    basis,nbf,nmo,npol,nocc,nvirt) 
       if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,*) 'Aft. add_H10()'
         write(*,*) '---- g_rhs-aft-add_H10 -------- START'
        endif
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-aft-add_H10 -------- END'
       endif

c ----> SKIPPING CPHF ------------------ START
c      if (ga_nodeid().eq.0)
c     &   write(*,*) 'WARNING : Using skip_cphf ...'
c            call skip_cphf(g_rhs,          ! IN/OUT
c     &                     dbl_mb(k_eval), ! IN: energy eigenvalues
c     &                     nocc,nvirt,nmo,npol)
c ----> SKIPPING CPHF ------------------ END

c -------free allocated memory -------------- START
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_occ',0,MA_ERR)
c -------free allocated memory -------------- END
      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- END'
      endif
c ------- creating dummy var -------- START
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nocc(ispin)
      enddo
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs0))
     &   call errquit('get_prelim_fock: ga_create failed g_rhs0',
     &                 0,GA_ERR)
      call ga_zero(g_rhs0)
c ------- creating dummy var -------- END
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'WARNING: SKIP cphf ...'
c        goto 123
      if(.not.rtdb_get(rtdb,'zora:skip_cphf_ev_hyp',
     &                 mt_log,1,skip_cphf_ev_hyp))        
     &  skip_cphf_ev_hyp = .false.       
      if (.not.(skip_cphf_ev_hyp)) then ! ----- if-skip-cphf--START
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'COMPUTE cphf hyp data ...'
c     Write ga_rhs to disk 
       call cphf_fname('cphf_rhs',cphf_rhs)
       call cphf_fname('cphf_sol',cphf_sol)
       if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_hfine: could not write cphf_rhs',0, DISK_ERR)
       call schwarz_tidy()
       call int_terminate()
c
c     Call the CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
       if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('hnd_hfine: failed to write skew ', 0, RTDB_ERR)
       if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_hfine: failure in cphf ',0, RTDB_ERR)
       if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('hnd_hfine: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
       call ga_zero(g_rhs)
       if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_hfine: could not read cphf_rhs',0, DISK_ERR)  
c ----- write CPHF data to file ----------------- START
c       Note.- lbl_cphfhyp defined in zora.fh
       call util_file_name(lbl_cphfhyp,.false.,.false.,zorafilename)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-hyp --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-hyp --------- END'
       call dft_zoraCPHF_write(
     &           zorafilename, ! in: filename
     &           npol,         ! in: nr polarization
     &           nocc,         ! in: nr occupied MOs
     &           nvirt,        ! in: nr virtual  MOs
     &           nbf,          ! in: nr basis functions
     &           vectors,      ! in: MOs
     &           g_rhs0,       ! out: (ntot,3)       GA matrix
     &           g_rhs)        ! in: (nocc*nvirt,3) GA matrix
c ----- write CPHF data to file ----------------- END
      else
       call ga_zero(g_rhs)
       call ga_zero(g_rhs0)
       do i=1,npol
        call ga_zero(vectors(i))
       enddo
       if (ga_nodeid().eq.0)
     &  write(*,*) 'READ cphf hyperfine data from file ...'
       call util_file_name(lbl_cphfhyp,.false.,.false.,zorafilename)
       call dft_zoraCPHF_read(
     &           zorafilename, !  in: filename
     &           npol,         !  in: nr polarization
     &           nocc,         !  in: nr occupied MOs
     &           nvirt,        !  in: nr virtual  MOs
     &           nbf,          !  in: nr basis functions
     &           vectors,      ! out: MOs
     &           g_rhs0,       ! out: (ntot,3)       GA matrix
     &           g_rhs)        ! out: (nocc*nvirt,3) GA matrix
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-hyp-read --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-hyp--read --------- END'
      endif ! ----- if-skip-cphf--END

 123  continue ! SKIPPING cphf

      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- START'
      call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- END'  
      endif
      type_NMR=2 ! =1,2,3=shieldings,hyperfine,gshift
      call get_P10(g_d1,     ! out: Perturbed density matrix occ-occ  contrib
     &             type_NMR, !  in: =1,2,3=shieldings,hyperfine,gshift
     &             g_rhs,    !  in: accumulated rhs expression
     &             g_rhs0,   !  in: from get_prelim_fock()
     &             vectors,g_CiFull, 
     &             nbf,nmo,npol,nocc,nvirt,
     &             do_zora,do_NonRel,not_zora_scale,rtdb,
     &             lbl_nlmohyp) 
       if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_hfine_zora: ga_destroy failed g_rhs',0,GA_ERR)
       if (.not.ga_destroy(g_rhs0)) call 
     &    errquit('hnd_hfine_zora: ga_destroy failed g_rhs0',0,GA_ERR)

c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the hyperfine tensor
      if (.not. ma_push_get(mt_dbl,9*nat_slc,'sh para',l_para,k_para)) 
     &    call errquit('hnd_hfine: ma_push_get failed k_para',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*nat_slc,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_hfine: ma_push_get failed k_dia',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*nat_slc,'sh para',l_tot,k_tot)) 
     &    call errquit('hnd_hfine: ma_push_get failed k_tot',0,MA_ERR)
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call get_par_hfine(dbl_mb(k_para), ! OUT: para-hfine tensor
     &                   ga_h01_hfine,   ! IN : h01 matrix
     &                   coeffpol,
     &                   g_d1,           ! IN : Perturbed density matrix
     &                   basis,nbf,nat_slc,rtdb)
      call get_dia_hfine(dbl_mb(k_dia),  ! OUT: dia-hfine tensor
     &                   ga_dia_hfine,   ! IN
     &                   coeffpol,
     &                   basis,nbf,nat_slc,rtdb)
c +++++++++ print-total-pardia-transferred +++ START
      if (ga_nodeid().eq.0)      
     &   write(*,105) nat_slc
 105     format('NATOMS=',i3) 
      ic=1
      do iatom = 1, nat_slc
       do ix = 1, 3
        do iy = 1, 3
         if (ga_nodeid().eq.0) then
           write(*,19) ix,iy,iatom,
     &                 dbl_mb(k_dia +ic-1),
     &                 dbl_mb(k_para+ic-1),
     &                 dbl_mb(k_dia +ic-1)+
     &                 dbl_mb(k_para+ic-1)
 19        format('NW:(dia,par,dia+par)(',
     &             i1,',',i1,',',i1,')=(',
     &             f15.8,' ',f15.8,' ',f15.8,' )')
         endif
         ic=ic+1
        enddo ! end-loop-iy
       enddo ! end-loop-ix
      enddo ! end-loop-iatom
c +++++++++ print-total-dia-transferred +++ END
       do ispin=1,ndens 
        if (.not.ga_destroy(g_dens(ispin))) call 
     &    errquit('hnd_hfine_zora: ga_destroy failed g_dens',
     &    0,GA_ERR)
       enddo
c
c     Print out tensor information, and write to Ecce file if necessary
      status = rtdb_parallel(.false.)   
      if (ga_nodeid().gt.0) goto 300
      acc_vec=0 ! For NMLO analysis
      call ecce_print_module_entry('nmr')
      do iatom = 1, nat_slc
         ioff = (iatom-1)*9
         if (.not. geom_cent_get(geom, int_mb(k_AtNr+iatom-1), tag, 
     &       dbl_mb(k_xyz), dbl_mb(k_zan))) call 
     &       errquit('hnd_hfine: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(tag, symbol, element, atn)) call
     &       errquit('hnd_hfine: geom_tag_to_element failed',0,GEOM_ERR)
c
c      Print tensor pieces and sum for total hyperfine tensor
         if (ga_nodeid().eq.0) then
c---- 2nd format: 2col hyperfine tensor (1scol: au 2ndcol: MHz) ---- START
          conv2MHz=dbl_mb(k_cfine+iatom-1)
          write(luout,9700) int_mb(k_AtNr+iatom-1),symbol ! Showing original atm nr.
          write(luout,97)
          write(luout,112) ' FC, SD terms:'
          trace=(dbl_mb(k_dia+ioff  ) +
     &           dbl_mb(k_dia+ioff+4) +
     &           dbl_mb(k_dia+ioff+8))/3.0d0
          write(luout,99) trace,trace*conv2MHz
  99      format('  isotropic',t16,'A(au)',
     &           f12.4,t60,'A(MHz)',f12.4)
          do iy=1,3
           iz=3*(iy-1)
           write (luout,98) (dbl_mb(k_dia+ioff+ix-1),ix=iz+1,iz+3), 
     &                      (dbl_mb(k_dia+ioff+ix-1)*conv2MHz,
     &                                               ix=iz+1,iz+3)
          enddo
          write(luout,112) ' PSO-Spin-Orbit terms:'
          trace=(dbl_mb(k_para+ioff  ) +
     &           dbl_mb(k_para+ioff+4) +
     &           dbl_mb(k_para+ioff+8))/3.0d0
          write(luout,99) trace,trace*conv2MHz
          do iy=1,3
           iz=3*(iy-1)
           write (luout,98) (dbl_mb(k_para+ioff+ix-1),ix=iz+1,iz+3), 
     &                      (dbl_mb(k_para+ioff+ix-1)*conv2MHz,
     &                                                ix=iz+1,iz+3)
          enddo
          do ix = 0, 8 
            dbl_mb(k_tot+ioff+ix) = dbl_mb(k_dia+ioff+ix) + 
     &                              dbl_mb(k_para+ioff+ix)
          enddo
          write(luout,112) ' Total hyperfine coupling tensor:'
          trace=(dbl_mb(k_tot+ioff  ) +
     &           dbl_mb(k_tot+ioff+4) +
     &           dbl_mb(k_tot+ioff+8))/3.0d0
          write(luout,99) trace,trace*conv2MHz
          do iy=1,3
           iz=3*(iy-1)
           write (luout,98) (dbl_mb(k_tot+ioff+ix-1),ix=iz+1,iz+3), 
     &                      (dbl_mb(k_tot+ioff+ix-1)*conv2MHz,
     &                                               ix=iz+1,iz+3)
          enddo
c---- 2nd format: 2col hyperfine tensor (1scol: au 2ndcol: MHz) ---- END
         endif

         do ix = 0, 8 
            dbl_mb(k_para+ioff+ix) = dbl_mb(k_dia+ioff+ix) + 
     &                               dbl_mb(k_para+ioff+ix)
         enddo
c
c     Print total hyperfine tensor

         if (ga_nodeid().eq.0) then
c            write(luout,9802) (dbl_mb(k_para+ioff+ix-1),ix=1,9)
c
c     Diagonalize total tensor
c     Order in a: xx xy yy xz yz zz 
            a(1) = dbl_mb(k_para+ioff)     
            a(2) = dbl_mb(k_para+ioff+1)
            a(3) = dbl_mb(k_para+ioff+4)
            a(4) = dbl_mb(k_para+ioff+2)
            a(5) = dbl_mb(k_para+ioff+5)
            a(6) = dbl_mb(k_para+ioff+8)
            ij = 0
            do 241 i = 1, 3
            do 241 j = 1, i
               ij = ij + 1
               axs(i,j) = a(ij)
               axs(j,i) = a(ij)
  241       continue
            call hnd_diag(axs,eig,3,.true.,.true.)
            call sort_asc(eig,axs) ! sort in ascending order
            isotr =(eig(1) + eig(2) + eig(3))/3.0d0
            aniso = eig(1) -(eig(2) + eig(3))/2.0d0
c ++++++++++ get eigenvectors for NMLO analysis +++++ START
            hypfile=0 ! not doing NLMO analysis by default
            status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
            if (hypfile.eq.1) then
             do i1=1,3
              do j1=1,3
               dbl_mb(k_tvec+acc_vec)=axs(i1,j1)
               acc_vec=acc_vec+1
              enddo 
             enddo 
            endif
c ++++++++++ get eigenvectors for NMLO analysis +++++ END
c ---- Calc: span,skew,asymm ----- START
             span  = eig(3) - eig(1)
             skew = 0.0d0
            if (abs(span)>small) 
     &       skew = 3.0d0*(isotr-eig(2))/span
            rtemp = eig(3) - isotr
             asym = 0.0d0
            if (abs(rtemp)>small) 
     &       asym = (eig(2)-eig(1))/rtemp
             write(luout,101)
 101         format(/1x,'In principal axis representation:',
     &                 ' total span,skew and asymm')
             write(luout,100) span,skew,asym,
     &                        span*conv2MHz,
     &                        skew*conv2MHz,
     &                        asym*conv2MHz
 100         format(t9,'span',t21,'skew',t32,'asymm',t52,
     &                  'span',t64,'skew',t75,'asymm'/3f12.4,
     &                   t44,3f12.4/)
c ---- Calc: span,skew,asymm ----- END
            write(luout,9986) (ix,ix=1,3),(ix,ix=1,3)
            write(luout,4489) (eig(ix),ix=1,3),
     &                        (eig(ix)*conv2MHz,ix=1,3)
 4489       format(3f12.4,t44,3f12.4,/)
            do iy=1,3
              write(luout,9983) iy,axs(iy,1),axs(iy,2),axs(iy,3)
            enddo
            write(luout,'(//)')
c
c     Print Ecce information

            call ecce_print1_char('atom name',symbol,1)
            call ecce_print2('hyperfine tensor',MT_DBL,
     &                       dbl_mb(k_para+ioff),3,3,3)
            call ecce_print1('hyperfine isotropic',MT_DBL,isotr,1)
            call ecce_print1('hyperfine anisotropy',MT_DBL,aniso,1)
            call ecce_print1('hyperfine eigenvalues',MT_DBL,eig,3)
            call ecce_print2('hyperfine eigenvectors',MT_DBL,axs,
     &                       3,3,3)
         endif
      enddo
      call ecce_print_module_exit('nmr','ok')
300   call ga_sync()
      status = rtdb_parallel(.true.)   

c ---- Destroy stored ga arrays ------ START
           if (.not. ga_destroy(g_d1)) call errquit(
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)  
           if (.not. ga_destroy(ga_dia_hfine)) call errquit(
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_hfine)) call errquit(
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)   
        if (.not. ga_destroy(ga_Fji_hfine)) call errquit(
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)             
c ---- Destroy stored ga arrays ------ END

      hypfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
      if (hypfile.eq.1) then ! ------- hypfile-if++++ START
         if (.not. ga_create(mt_dbl,1,nat_slc*9,
     &       'munu4nbo: g_tvec',0,0,g_tvec)) 
     $       call errquit('munu4nbo: g_tvec', 0,GA_ERR)
        call ga_dgop(msg_efgs_col,dbl_mb(k_tvec),nat_slc*9,'+')
        call ga_put(g_tvec,1,1,1,nat_slc*9,dbl_mb(k_tvec),1)      
        call create_munu4nbo_hyp(rtdb,g_tvec,
     &                           coeffpol,
     &                           nat_slc,int_mb(k_AtNr),
     &                           basis,npol,nocc,nvirt,nmo,
     &                           dbl_mb(k_cfine)) 
         if (.not. ga_destroy(g_tvec)) call errquit( ! destroy GA
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)   
       if (.not.ma_free_heap(l_tvec)) call
     &     errquit('hnd_hyperfine_zora: ma_free_heap l_tvec',0,MA_ERR)    
      endif ! ------------------------ hypfile-if++++ END
c ---- Destroy stored ga arrays ------ START
c          if (.not. ga_destroy(g_d1)) call errquit(
c     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)  
c           if (.not. ga_destroy(ga_dia_hfine)) call errquit(
c     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)  
c        if (.not. ga_destroy(ga_h01_hfine)) call errquit(
c     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)   
c        if (.not. ga_destroy(ga_Fji_hfine)) call errquit(
c     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)             
c ---- Destroy stored ga arrays ------ END
c
c     Clean up all remaining memory
      if (.not.ma_pop_stack(l_tot)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed k_tot',
     &            0,MA_ERR)
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed k_dia',
     &            0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed k_para',
     &            0,MA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('hnd_hyperfine_zora: ga_destroy failed vectors',
     &             0,GA_ERR)
      enddo
        if (.not. ga_destroy(g_AtNr1)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)   
      if (.not.ma_pop_stack(l_zan)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed k_zan',
     &             0,MA_ERR)
      if (.not.ma_pop_stack(l_xyz)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed k_xyz',
     &             0,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('hnd_hyperfine_zora: ma_free_heap l_tmp',
     &             0, MA_ERR)
      if (.not.ma_pop_stack(l_cfine)) call
     &    errquit('hnd_hyperfine_zora: ma_pop_stack failed l_cfine',
     &             0,MA_ERR)
      call schwarz_tidy()
      call int_terminate()
      return

 97   format(/' ----- A-tensor in atomic units ----',t44,
     1        '  ------------- A (MHz) ------------')
 98   format(3f12.4,t44,3f12.4) ! for showing 2col hyperfine tensor
 112  format(/a)
 7000 format(/,10x,'NMR hyperfine cannot be calculated for UHF',
     1      ' or ROHF wave functions: use RHF')
 9700 format(6x,'Atom: ',i4,2x,a2)
 9800 format(8x,'Diamagnetic',/,3(3F12.4,/))
 9801 format(8x,'Paramagnetic',/,3(3F12.4,/))
 9802 format(8x,'Total hyperfine Tensor',/,3(3F12.4,/))
 9983 format(i1,1x,f10.4,f12.4,f12.4)
c 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
c 9986 format(10x,'Principal Components and Axis System',/,10x,
c     1       3(7x,i1,4x))
 9986 format(1x,'Principal Components and Axis System',/,1x,
     1       3(7x,i1,4x),t44,3(7x,i1,4x))
 9987 format(10x,' isotropic = ',f12.4,/,
     1       10x,'anisotropy = ',f12.4,/)
 9999 format(
     1 /,10x,41(1h-),/,
     2 10x,'Hyperfine Tensors (in au)',/,
     3 10x,41(1h-),/)
      end
c
c --------------- For NMLO analysis ------------------ START
      subroutine get_tvec(dia,     ! in : dia tensor
     &                    tvec,    ! out: eigenvectors
     &                    acc_vec, ! in/out: counter for eigenvectors
     &                    iat,     ! in : atom nr.
     &                    hypfile) ! in : NMLO-hyp flag
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "bas.fh"
#include "stdio.fh"

         integer iat,acc_vec
         double precision dia(*),tvec(*)
         integer i,j,ic,hypfile
         double precision vec(3,3), eig(3)    
         external hnd_diag
         ic=9*(iat-1)+1
         do i = 1, 3
          do j = 1, 3
           vec(i,j)= dia(ic)
           ic=ic+1
          enddo ! end-loop-i
         enddo ! end-loop-j
         call hnd_diag(vec,eig,3,.true.,.false.)
c ------- copy eigenvectors, vec for NLMO analysis ----- START
          if (hypfile.eq.1) then
           do i=1,3
            do j=1,3
             tvec(acc_vec)=vec(i,j)
             acc_vec=acc_vec+1
            enddo 
           enddo 
          endif
      return
      end
c --------------- For NMLO analysis ------------------ END
      subroutine get_par_hfine(
     &                   par,         ! OUT : paramagnetic tensor
     &                   ga_h01_hfine,! IN :
     *                   coeffpol,    ! IN  : scaling factor 1/(nA-nB)
     &                   g_d1,        ! IN  : Perturbed density matrix 
     &                   basis,       ! IN  : basis handle
     &                   nbf,         ! IN  : nr. basis functions
     &                   natoms,      ! IN  : nr. atoms
     &                   rtdb)
c Purpose : Assemble NMR Hyperfine: paramagnetic tensor
c Author  : Fredy Aquino
c Date    : 03-09-11
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "bas.fh"
#include "stdio.fh"
#include "zora.fh"  
c      Global arrays variables defined in zora.fh
c      ga_dia_hfine,ga_h01_hfine,ga_Fji_hfine
      integer g_d1 ! input: Perturbed density matrices 
      integer ga_h01_hfine
      integer rtdb
      integer basis
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer ld(2),cbuf,cbuf1,cbuf2
      integer nbf,natoms
      integer iatom,ix,iy,ixy,ind
      integer debug_get_par
      double precision val,par(*),coeffpol
      logical oskel
       debug_get_par=0
       do ixy = 1, 9*natoms
         par(ixy) = 0.0d0  ! initialize
       enddo
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nbf
       ixy = 0
       blo(3) = 0
       bhi(3) = 0
       do iatom = 1, natoms
        do iy = 1, 3
         blo(3) = blo(3) + 1
         bhi(3) = bhi(3) + 1
         do ix = 1, 3
          alo(3) = ix
          ahi(3) = ix
          ixy = ixy + 1
          val = nga_ddot_patch(g_d1        ,'n',alo,ahi,
     &                         ga_h01_hfine,'n',blo,bhi)*
     *          coeffpol*
     &          (-2.0d0) ! factor from (+2) g_N B_N/(nA-nB)
                         !             (-1) from h_01 operator       
          cbuf=9*(iatom-1)+3*(ix-1)+iy ! transpose  
          par(cbuf)=val
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
c --------- symmmetrize para ------- START
       do iatom = 1, natoms
        do iy = 1,3
         do ix = iy+1,3
          cbuf1=9*(iatom-1)+3*(ix-1)+iy ! transpose  
          cbuf2=9*(iatom-1)+3*(iy-1)+ix   
          val=(par(cbuf1)+par(cbuf2))/2.0d0
          par(cbuf1)=val
          par(cbuf2)=val
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
c --------- symmmetrize para ------- END
      return
      end

      subroutine get_dia_hfine(dia,         ! OUT : dia-hfine tensor
     &                         ga_dia_hfine,! IN :
     &                         coeffpol,    ! IN  : scaling factor
     &                         basis,       ! IN  : basis handle
     &                         nbf,         ! IN  : nr. basis functions
     &                         natoms,      ! IN  : nr. atoms
     &                         rtdb)

c Purpose : Assemble NMR Hyperfine: diamagnetic tensor
c Author  : Fredy Aquino
c Date    : 03-09-11
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "zora.fh"  
c      Global arrays variables defined in zora.fh
c      --> ga_dia_hfine,ga_h01_hfine,ga_Fji_hfine
      integer rtdb
      integer basis     ! [input] basis handle
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer ld(2),cbuf,cbuf1,cbuf2
      integer nbf,natoms,ispin
      integer iatom,ix,iy,ixy
      integer ga_dia_hfine
      double precision dia(*),coeffpol
      logical oskel
      double precision val
      do ixy = 1, 9*natoms
         dia(ixy) = 0.0d0  ! initialize
      enddo
c ---- STORE: ga_dia --> dia
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=natoms
       ld(1)=3
       ld(2)=3 
       call ga_scale(ga_dia_hfine,coeffpol)
       call nga_get(ga_dia_hfine,alo,ahi,dia(1),ld)
c --------- symmmetrize dia ------- START
       do iatom = 1, natoms
        do iy = 1,3
         do ix = iy + 1,3
          cbuf1=9*(iatom-1)+3*(ix-1)+iy ! transpose  
          cbuf2=9*(iatom-1)+3*(iy-1)+ix   
          val=(dia(cbuf1)+dia(cbuf2))/2.0d0
          dia(cbuf1)=val
          dia(cbuf2)=val
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
c --------- symmmetrize dia ------- END
      return
      end

      subroutine add_H10_hfine( 
     &                    g_rhs, ! out: accumulated rhs expression
     &                   ga_Fji, !  in: Fock 1st-deriv from operator: p K x p
     &                  vectors, !  in: MO  coeffs
     &                    basis, !  in: basis handle
     &                      nbf, !  in: nr. basis functions
     &                      nmo, !  in: nr. MOs (occ+virt)
     &                     npol, !  in: nr. of polarizations
     &                     nocc, !  in: nr. occ     MOs
     &                    nvirt) !  in: nr. virtual MOs
c =============  Computing  ==============
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c ========================================
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol    ! nr. of polarizations
      logical do_zora ! logical to check if NOT doing zora calc
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,disp
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3)
      integer vectors(npol)
      integer g_rhs,ga_Fji
      integer g_s10 ! scratch ga arrays
      integer basis ! basis handle
      logical oskel
      external giao_aotomo
      oskel = .false.
c ----------- Create scratch ga-arrays ------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*npol
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                   alo,g_s10)) 
     &    call errquit('add_H10: nga_create failed g_s10',
     &                  0,GA_ERR)
       call ga_zero(g_s10)
c ----------- Create scratch ga-arrays ------- END
c     Transform H10 to MO and add to g_rhs
c     Copy: g_Fji --> g_s10
       blo(1) = 1
       bhi(1) = nmo
       blo(2) = 1
       bhi(2) = nmo
       blo(3) = 1
       bhi(3) = 3
      do ispin=1,npol   
       shift=3*(ispin-1)
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nmo
       alo(3) = shift+1
       ahi(3) = shift+3   
       call nga_copy_patch('n',ga_Fji,blo,bhi,
     &                          g_s10,alo,ahi) 
      enddo ! end-loop-ispin
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,3,nbf)
      do ispin=1,npol
       shift=3*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+3   
c ------ definitions for g_rhs -------- START
       disp=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c ------ definitions for g_rhs -------- END
       call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs,blo,bhi)
      enddo ! end-loop-ispin
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('add_H10: ga_destroy failed g_s10',0,GA_ERR)
      return
      end
      subroutine get_convfactor_hfine(geom,       ! in  : geom handle
     &                                const_hfine,! out : conversion factor for hyperfine calc.
     &                                Zan,        ! in  : atomic numbers
     &                                natoms,     ! in  : nr. of selected atoms
     &                                atmnr)      ! in  : list of atom nr. indices
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer geom
      integer iat,iat1,isonr
      character*16 at_tag 
      integer natoms,atmnr(natoms)
      double precision gnuc,hbar,ge,emf,vl,auev,evmhz,gmhz,
     &                pi,fac,betae,betan,conv,convf,con,gnu,
     &                Zan(natoms)
      double precision const_hfine(natoms)
      logical status
      logical atom_gfac
      external atom_gfac
      data gnuc   /5.05078343d-27/      ! Nuclear magneton
      data hbar   /1.05457168d-34/      ! Planck constant over 2 pi
      data ge     /2.002319304386d+00/  ! Electron g-factor
      data emf    /1836.152701d+00/     ! Proton-electron mass ratio
      data vl     /137.0359895d+00/     ! Speed of light in au
      data auev   /27.2113961d+00/      ! Conversion from au to eV
      data evmhz  /2.41798836d+08/      ! Conversion from eV to MHz
      data gmhz   /2.8025d+00/          ! Conversion from Gauss to MHz

c     --- calculate constants and conversion terms ---START
      pi   = acos(-1.0d0)
      fac  =(4.0d0*pi/3.0d0)
      betae=1.0d0/(2.0d0*vl)
      betan=betae/emf
      convf=auev*evmhz
      con  =ge*betae*betan*convf
c     --- calculate constants and conversion terms ---END
c----- Allocate memory - FA
       do iat1 = 1, natoms
         iat=atmnr(iat1)
        if (.not. atom_gfac(Zan(iat1),gnu,isonr)) call
     &      errquit('hnd_hfine: atom_gfac failed',0, UERR)
        const_hfine(iat1)=con*gnu
        if (ga_nodeid().eq.0)
     &   write(*,1) iat1,iat,con,gnu,isonr,const_hfine(iat1)
 1       format('(con,gnu,isonr,const_hfine)(',i3,',',i3,')=(',
     &          f15.8,',',f15.8,',',i3,',',f15.8,')')
       enddo ! end loop-iat1
      return
      end

      subroutine sort_asc(eig,vec)
c     Purpose: Sort in ascending order
c     Taken forn hnd_diag
c     FA-03-17-11
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
         integer i,j,jj
         double precision eig(3),vec(3,3),xx
         do 50 i=1,3
            jj=i
            do 30 j=i,3
               if( eig(j).le. eig(jj)) jj=j
   30       continue
            if(jj.eq.i) go to 50

            xx=eig(jj)
            eig(jj)=eig(i)
            eig(i)=xx
            do 40 j=1,3
               xx=vec(j,jj)
               vec(j,jj)=vec(j,i)
               vec(j,i)=xx
   40       continue
   50    continue
      return 
      end
c --------------------------------------------------------
c -------------- for NMLO analysis ----------------- START
      subroutine create_munu4nbo_hyp(
     &             rtdb,        ! in: rtdb handle
     &             g_tvec,      ! in: eigenvectors or T diagonalizing matrix
     &             coeffpol,    ! in: scaling factor for hyperfine constants
     &             nat,         ! in: nr. atoms
     &             atmnr,       ! in: list of selected atoms 
     &             basis,       ! in: basis handle
     &             npol,        ! in: nr. polarizaitons
     &             nocc,        ! in: nr. occ   nocc(i) i=1,npol
     &             nvirt,       ! in: nr. virt nvirt(i) i=1,npol
     &             nmo,         ! in: nr. MO
     &             const_hfine) ! in: conv2MHz

      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh" 
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "msgids.fh" 
#include "zora.fh"   
c FA: Revised on 06-22-11
c ------Main outputs -------- START
      integer g_munu_rot,         ! hyp-FCSD  dia  part
     &        g_munu_rot1,g_acc2  ! hyp-PSOSO para part
c ------Main outputs -------- END
      integer rtdb,basis
      integer g_sdens,g_tvec
      integer g_munuFCSD,g_munuPSOSO,
     &        g_tnp,g_acc,
     &        vectors(2),
     &        g_tnp1,g_acc1    
      integer ispin  
      integer nat  
      double precision ac_val,val1,sign,
     &                 const_hfine(nat)   

      integer iat1,iat,i,j,k,m,n,ndir,ndir1
      integer jlo,jhi,s,nbf,nmo,nsize,nsize1,nsize2
      integer npol,ntot,atmnr(nat)
      integer ind,nlst,count,nocc(npol),nvirt(npol)
      integer Natoms_munu,Ndir_munu,atmnr_munu(nat)
c
c     Ndir_munu, Nr. of directions stored
c                =3  xx yy zz
      double precision coeff,fact,para_rot(9),coeffpol
      double precision tmn(2),chcdata(3)
      integer jlo1,jhi1,jlo2,jhi2
      integer g_sdens1,g_c1,
     &        g_p10,g_munuPSOSO2d
      integer iind(2),jind(2),icalczora,type_NMR  
      integer alo(3),ahi(3),elo(3),ehi(3),flo(3),fhi(3)
      logical dft_zoraHYP_NLMOAnalysis_read ! for read-nlmo-mat
      character*255 zorafilename              ! for read-nlmo-mat
      integer arr_ind(6,2)
       data ((arr_ind(j,i),i=1,2),j=1,6)
     &  /1,1,2,2,3,3,1,2,1,3,2,3/
      external dft_zoraHYP_NLMOAnalysis_read,get_P10_rot,
     &         fill_munuPSOSO_1,get_par_hfine_rot,
     &         get2dmat,whypfile,
     &         get_par_hfine_rot1
c     --> To store ONLY munu principal components xx,yy,zz 
c     g_munuFCSD    is created in dft_zora_Hyperfine.F
c     size(g_munuFCSD)=nlst*ndir*nat (linear vector)
c     g_sdens, spin density matrix
c     nbf x nbf elements (bidimensional matrix)
c     Legend:
c     ndir=6 = xx, yy, zz, xy, xz, yz
c     nbf, Nr of basis functions
c     nlst=nbf*(nbf+1)/2
c     nat, nr. of selected atoms out of nat

      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('munu: bas_numbf failed',555, BASIS_ERR)
      Natoms_munu=nat
      do i=1,Natoms_munu
       atmnr_munu(i)=atmnr(i)
      enddo
      Ndir_munu=3
      nlst=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
c ++++++ Read NLMO matrices +++++++++ START
      ndir=6
      ndir1=3
      call util_file_name(lbl_nlmohyp,.false.,.false.,zorafilename)    
      icalczora = 0  ! initialize the flag
      if (.not.dft_zoraHYP_NLMOAnalysis_read(
     &       zorafilename, ! in : filename
     &                nbf, ! in : nr basis functions
     &               ndir, ! in : nr of directions: 6 = xx yy zz xy xz yz
     &              ndir1, ! in : nr of directions: 3 = x y z
     &                nat, ! in : list of selected atoms 
     &               nocc, ! in : nocc(i) i=1,2 nr. occupations in alpha and beta
     &               npol, ! in: nr polarizations
     &         g_munuFCSD, ! out: munu matrix of Fermi Contact + Spin Dipolar term
     &        g_munuPSOSO, ! out: munu matrix of PSOSO
     &            vectors, ! out: MOs
     &               g_c1, ! out: perturbed MO
     &            g_sdens)) icalczora=1 
      call ga_scale(g_munuFCSD ,       coeffpol) ! scaling dia
      call ga_scale(g_munuPSOSO,-2.0d0*coeffpol) ! scaling para     
      goto 10
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuFCSD -- START'
        call ga_print(g_munuFCSD)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuFCSD -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuPSOSO -- START'
        call ga_print(g_munuPSOSO)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuPSOSO -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_sdens ----- START'
        call ga_print(g_sdens)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_sdens ----- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- START'
        call ga_print(g_c1)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- END'
 10   continue
c ++++++ Read NLMO matrices +++++++++ END
      call get_unique_elmat(g_sdens,g_sdens1,nlst,nbf)   ! out: g_sdens1
      ndir=6 ! Nr. of directions: xx,yy,zz,xy,xz,yz
      nsize=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
      nsize1=nsize*ndir   ! size of whole munu per atom
      nsize2=nsize*ndir1  ! size of whole munu per atom
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_tnp',0,0,g_tnp)) 
     $    call errquit('munu4nbo: g_tnp', 0,GA_ERR)
        call ga_zero(g_tnp)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_tnp1',0,0,g_tnp1)) 
     $    call errquit('munu4nbo: g_tnp1', 0,GA_ERR)
        call ga_zero(g_tnp1)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_acc',0,0,g_acc)) 
     $    call errquit('munu4nbo: g_acc', 0,GA_ERR)
        call ga_zero(g_acc)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_acc1',0,0,g_acc1)) 
     $    call errquit('munu4nbo: g_acc1', 0,GA_ERR)
         call ga_zero(g_acc1)       
         alo(1) = nbf
         alo(2) = -1
         alo(3) = -1
         ahi(1) = nbf
         ntot=nocc(1)+nocc(2)
         ahi(2) = ntot
         ahi(3) = 3*nat
         if (.not.nga_create(MT_DBL,3,ahi,'g_acc2 matrix',
     &       alo,g_acc2)) call 
     &        errquit('g_acc2: nga_create failed g_acc2',0,GA_ERR)
         call ga_zero(g_acc2)
        if (.not. ga_create(mt_dbl,1,nlst*3*nat,
     &                       'munu4nbo: g_munu_rot',0,0,g_munu_rot)) 
     $    call errquit('munu4nbo: g_munu_rot', 0,GA_ERR)
        if (.not. ga_create(mt_dbl,1,nlst*3*nat,
     &                       'munu4nbo: g_munu_rot1',0,0,g_munu_rot1)) 
     $    call errquit('munu4nbo: g_munu_rot1', 0,GA_ERR)
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
      if (.not.nga_create(mt_dbl,3,ahi,'g_munuPSOSO2d matrix',
     &    alo,g_munuPSOSO2d)) call 
     &    errquit('munu4nbo: nga_create failed g_munuPSOSO2d',0,GA_ERR)
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &              ' NW-Hyperfine: Summary C+HC data [MHz] ',
     &              'oooooooooooooooo START'
      do iat1=1,nat
        call ga_zero(g_munuPSOSO2d)
        iat=atmnr(iat1)
        do n=1,3  ! xx,yy,zz
         m=n ! For principal components ONLY
         call ga_zero(g_acc)
c ----- Do: A'= T^t A T, calculate only [A']_pp --> (do n=1,3 m=n)
c       a_pp'=    \sum_i t_ip a_ii t_ip + 
c               2 \sum_{j=2}^n \sum_{i=1}^{j-1} t_jp a_ji t_ip
c       g_munu_rot = A'
c       WARNING: g_munu_rot, contains several rotated matrices
c                since the matrices are symmetric I store only
c                the main diagonal + lower (upper) triangular 
c                matrix in a format that looks like:
c                a_11 a_22 ... a_nn 
c                a_21
c                a_31 a_32
c                a_41 a_42 a_43
c                ...
c                a_n1 a_n2 ... a_{n(n-1)}
c      There are two additional transformations on g_munu_rot
c      before leaving this routine and entering wefgfile()
c      1. I make the diagonalized matrix traceless
c ===== Transform xx_munu to 2xx_munu-(yy_munu+zz_munu) = START
c       or                   3xx_munu-(xx_munu+yy_munu+zz_munu)
c      2. I need to do a reordering of elements so that it is
c         compatible with wefgfile()
c        call reorder_munu(g_munu_rot,nat,nlst,nbf,Ndir_munu)
c --------------------------------------------------------------
         do s=1,6
c ------- get coeff() --- START
          iind(1)=1
          iind(2)=1
          jind(1)=9*(iat1-1)+3*(arr_ind(s,1)-1)+m
          jind(2)=9*(iat1-1)+3*(arr_ind(s,2)-1)+n   
          call ga_gather(g_tvec,tmn,iind,jind,2)
          fact=1.0d0
          if (s.gt.3) fact=2.0d0
          coeff=fact*tmn(1)*tmn(2)      
c ------- get coeff() --- END
c Note.- g_munuFCSD will be the (hyp-diag)_uv matrix
           jlo=nsize1*(iat1-1)+nsize*(s-1)+1
           jhi=jlo+nsize-1
           call ga_copy_patch('n',g_munuFCSD,1,1,jlo,jhi,
     &                            g_tnp     ,1,1,1  ,nsize)
           call ga_add(1.0d0,g_acc,coeff,g_tnp,g_acc)    
         enddo ! end-loop-s
c Note.- g_acc = \widetilde{H}_{mu nu}^{(m,m)}
c        it is the rotated munu matrix using:  T^t H T
         call ga_zero(g_acc1)
         do s=1,3
c ------- get coeff() --- START
          iind(1)=1
          jind(1)=9*(iat1-1)+3*(s-1)+m
          call ga_gather(g_tvec,tmn,iind,jind,1)

c          if (ga_nodeid().eq.0) then
c           write(*,1) m,s,tmn(1)
c 1         format('(m,s,tvec)=(',i5,',',i5,',',f15.8,')')
c          endif

          coeff=tmn(1) 
c ------- get coeff() --- END
c Note.- g_munuPSOSO will be the (hyp-para)_uv matrix
          jlo=nsize2*(iat1-1)+nsize*(s-1)+1
          jhi=jlo+nsize-1
          call ga_copy_patch('n',g_munuPSOSO,1,1,jlo,jhi,
     &                           g_tnp1     ,1,1,1  ,nsize)
          call ga_add(1.0d0,g_acc1,coeff,g_tnp1,g_acc1)   
c ----- Calculate rotated perturbed MO: g_acc2 ----- START
c  \sum_{s=1,3} t_{sm} C_{ri}^{(s) sigma} --> g_acc2
          elo(1) = 1
          ehi(1) = nbf
          elo(2) = 1
          ehi(2) = ntot
          elo(3) = s
          ehi(3) = s
          flo(1) = 1
          fhi(1) = nbf
          flo(2) = 1
          fhi(2) = ntot
          flo(3) = 3*(iat1-1)+m
          fhi(3) = 3*(iat1-1)+m
          call nga_add_patch(1.0d0,g_acc2,flo,fhi,
     &                       coeff,g_c1  ,elo,ehi,
     &                             g_acc2,flo,fhi)
c ----- Calculate rotated perturbed MO: g_acc2 ----- END 
         enddo ! end-loop-s
c Note: g_acc1 = \widetilde{H}_{mu nu}^{(m)}  m=x,y,z
c       it is the rotated munu matrix using: T H
c ====== Store final munu matrices === START
         jlo2=nlst*Ndir_munu*(iat1-1)+
     &        nlst*(n-1)+1
         jhi2=jlo2+nlst-1
         call ga_scale(g_acc ,const_hfine(iat1)) ! scaling to conv. to MHz
         call ga_scale(g_acc1,const_hfine(iat1)) ! scaling to conv. to MHz
         call ga_copy_patch('n',g_acc     ,1,1,   1,nlst,
     &                          g_munu_rot,1,1,jlo2,jhi2)      
         call ga_copy_patch('n',g_acc1    ,1,1,   1,nlst,
     &                         g_munu_rot1,1,1,jlo2,jhi2)     
c ====== Store final munu matrices === END
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== START
c ==== sum (g_acc .* g_sdens1 + Nuclear CONTRIB) 
c      = TOTAL HYP diagonalized   
         jlo1=1+nbf
         jhi1=nsize
         call ga_scale_patch(g_acc,1,1,jlo1,jhi1,2.0d0)
         chcdata(n)=ga_ddot(g_acc,g_sdens1)
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
        enddo ! end-loop-n
        if (ga_nodeid().eq.0) then
          write(*,23) iat,
     &                chcdata(1),            ! dia-x
     &                chcdata(2),chcdata(3)  ! dia-y,z
 23       format(' CHC  FCSD(xx,yy,zz)(',i3,')=(', 
     &           f15.8,',',f15.8,',',f15.8,')') 
        endif
c ++++++++++ CHECK diagonalization in para hyperfine ++++++++ START
c Note.- Variables defined in zora.fh:
c        g_CiFull, zora scaling factors 
c                 filled out with values in dft_zora_scale
c        zora switches: do_zora,do_NonRel,not_zora_scale
        type_NMR=2 ! =1,2,3=shieldings,hyperfine,gshift
        call get_P10_rot(
     &          g_p10,            ! out: Perturbed density matrix (munu nbf x nbf x 3 square matrix)
     &          type_NMR,         !  in: =1,2,3=shieldings,hyperfine,gshift
     &          g_acc2,           !  in: rotated perturbed MO vector
     &          vectors,g_CiFull, !  in: to build zora scaled MO vector 
     &          iat1,             !  in: index for selected atom nr =1,nlist
     &          nbf,nmo,npol,nocc,nvirt,
     &          do_zora,do_NonRel,not_zora_scale,rtdb) 
      call fill_munuPSOSO_1(  ! g_munuPSOSO --> g_munuPSOSO2d
     &        g_munu_rot1   , ! in: array with unique elements
     &        g_munuPSOSO2d , !out: nbf x nbf x 3 munu matrix for ith atom
     &        iat1          , ! in: = 1,2,..., nlist
     &        2,              ! in: type_symm = 1 symm  = 2 antisymm
     &        nbf) 
c +++++++++ NOW: do ddot product to get diagonalized tensor +++ START
      call get_par_hfine_rot(
     &                   g_munuPSOSO2d, ! IN : h01 matrix
     &                   g_p10,         ! IN : Perturbed density matrix
     &                   basis,nbf,iat,rtdb)
c +++++++++ NOW: do ddot product to get diagonalized tensor +++ END
c ++++++++++ CHECK diagonalization in para hyperfine ++++++++ END
      if (.not. ga_destroy(g_p10)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      enddo ! end-loop-iat
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &   ' NW-Hyperfine: Summary C+HC data [MHz] ',
     &              'oooooooooooooooo END'
c --> Main outputs: g_acc2     ,rotated perturbed MO nbf*ntot*ndir*nlist
c                               ndir=1,2,3=x,y,z
c                               nbf, nr of basis functions
c                               ntot=nocc(1)+nocc(2)
c                               nlist, nr of selected atoms
c                   g_munu_rot1,rotated perturbed AO matrix
c                               storing only diag + off-diag elements
c                               Reminder: this comes from an antisymmetrix matrix
c                                         in case we want to pull back the 2d munu-matrix
c                   g_munu_rot, rotated AO matrix for dia part
c                               storing only diag + off-diag elements
c                               Reminder: this comes from a  symmetric matrix
c                                         in case we want to pull back the 2d munu-matrix
      call reorder_munu(g_munu_rot ,nat,nlst,nbf,Ndir_munu) ! reoder-munu matrix
      call reorder_munu(g_munu_rot1,nat,nlst,nbf,Ndir_munu) ! reoder-munu matrix
c ------ destroy unnecessary GAs 
      if (.not. ga_destroy(g_munuFCSD)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuPSOSO)) call errquit(
     &  'create_munu4nbo-2: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuPSOSO2d)) call errquit(
     &  'create_munu4nbo-7: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp)) call errquit(
     &  'create_munu4nbo-5: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc)) call errquit(
     &  'create_munu4nbo-6: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp1)) call errquit(
     &  'create_munu4nbo-5a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc1)) call errquit(
     &  'create_munu4nbo-6a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_sdens)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_sdens1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_c1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('create_munu4nbo: ga_destroy failed vectors',
     &             0,GA_ERR)
      enddo
c ------ Free allocated memory
      call  whypfile(
     &         rtdb,
     &         g_munu_rot,
     &         g_munu_rot1,
     &         g_acc2,
     &         nlst,       
     &         Ndir_munu,
     &         Natoms_munu,
     &         atmnr_munu)
      if (.not. ga_destroy(g_munu_rot)) call errquit( 
     &  'whypfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munu_rot1)) call errquit( 
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc2)) call errquit(     
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      return
      end  

      subroutine get_P10_rot(  
     &                     g_p1, ! out: Perturbed (spin)density matrix
     &                 type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift              
     &                     g_c1, ! in: rotated perturbed MO
     &                  vectors, ! in: MO                coeffs
     &                 g_CiFull, ! in: MO zora weighting coeffs
     &                     iat1, ! in: index of selected atoms =1,nlist
     &                      nbf, ! in: nr. basis functions
     &                      nmo, ! in: nr. MOs (occ+virt)
     &                     npol, ! in: nr. of polarizations
     &                     nocc, ! in: nr. occ     MOs
     &                    nvirt, ! in: nr. virtual MOs
     &                  do_zora, ! in: .true.  if doing zora calc
     &                do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &           not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &                     rtdb)
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
c Purpose: Compute rotated P10 (perturbed density matrix)
c          Using a rotated perturbed MO (g_c1)
c          This routine is modified from get_P10
c          and is intended to be used in hyperfine coupling constants
c          NLMO analysis.
c Author : Fredy W. Aquino
c Date   : 07-07-11
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol   ! nr. of polarizations
      integer g_p1   ! OUT: Perturbed density matrix
      integer g_u,g_d1 ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,type_NMR,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz,iat1
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      double precision coeff(2)
      logical status,do_zora,do_NonRel,not_zora_scale
c ------ for Hyperfine NMLO analysis --------- START
      integer ntot,g_c1,  ! g_c1 , collects perturbed MO coeffs C1
     &        plo(3),phi(3),qlo(3),qhi(3),
     &        nlist
c ------ for Hyperfine NMLO analysis --------- END          
       if      (type_NMR.eq.1) then ! Shieldings
         coeff(1)=  1.0d0
         coeff(2)=  1.0d0
       else if (type_NMR.eq.2 .or.  ! Hyperfine
     &          type_NMR.eq.3) then ! g-shifts
         coeff(1)=  1.0d0
         coeff(2)= -1.0d0
       else
        write(*,*) 'Error in get_P10_1:',
     &             ' Calc. should be giao, gshift or hyperfine.'
        stop
       endif
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'get_P10_rot: g_d1 matrix',
     &                    alo,g_d1)) call 
     &     errquit('g_d1: nga_create failed g_d1',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'get_P10_rot: g_p1 matrix',
     &                     alo,g_p1)) call 
     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
       call ga_zero(g_p1)
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('g_u: nga_create failed g_u',0,GA_ERR)
       call ga_zero(g_u)
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dlo(2) = 1
       dhi(1) = nbf
       dhi(2) = nocc(ispin)
c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then
        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
c +++++++++++ here g_c1 --> g_u +++++++++++ START
          shift=nocc(1)*(ispin-1)
          qlo(1) = 1
          qhi(1) = nbf
          qlo(2) = shift+1
          qhi(2) = shift+nocc(ispin)
          qlo(3) = 3*(iat1-1)+xyz
          qhi(3) = 3*(iat1-1)+xyz
          plo(1) = 1
          phi(1) = nbf
          plo(2) = 1
          phi(2) = nocc(ispin)
          plo(3) = xyz
          phi(3) = xyz
          call nga_copy_patch('n',g_c1,qlo,qhi,
     &                            g_u ,plo,phi)     
c +++++++++++ here g_c1 --> g_u +++++++++++ END
        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1      
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u,alo,ahi,
     &                                  g_d1,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1,clo,chi)
       enddo ! end-loop-xyz
       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
       call nga_add_patch(1.0d0       ,g_p1,elo,ehi,
     &                    coeff(ispin),g_d1,elo,ehi,
     &                                 g_p1,elo,ehi)
       if (.not.ga_destroy(g_u)) call 
     &    errquit('get_d1: ga_destroy failed g_u',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
      enddo ! end-loop-ispin
      if (npol.eq.1 .and. type_NMR.eq.1) then ! this happens ONLY for NMR-restricted calc.
        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1,elo,ehi,2.0d0)
      endif    
       if (.not.ga_destroy(g_d1)) call 
     &    errquit('get_d1: ga_destroy failed g_d1',0,GA_ERR)      
      return
      end 

      subroutine get_par_hfine_rot(   ! for checking rotated Pert. MO and rotated munu AO
     &                   ga_h01_hfine,! IN :
     &                   g_d1,        ! IN  : Perturbed density matrix 
     &                   basis,       ! IN  : basis handle
     &                   nbf,         ! IN  : nr. basis functions
     &                   iat,         ! IN  : atom nr
     &                   rtdb)        ! IN  : rtdb handle
c Purpose : Assemble NMR Hyperfine: paramagnetic tensor
c Author  : Fredy Aquino
c Date    : 07-08-11
c Note    : This routine is taken from get_par_hfine()
c           The only difference is that here natoms=1
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "bas.fh"
#include "stdio.fh"

      integer g_d1 ! input: Perturbed density matrices 
      integer ga_h01_hfine
      integer rtdb
      integer basis
      integer alo(3), ahi(3)
      integer ld(2),cbuf 
      integer nbf,iat
      integer iy
      double precision val,chcdata(3)
      logical oskel

       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
        do iy = 1, 3
          alo(3) = iy
          ahi(3) = iy
          val = nga_ddot_patch(g_d1        ,'n',alo,ahi,
     &                         ga_h01_hfine,'n',alo,ahi)  
          chcdata(iy)=val
        enddo ! end-loop-iy
        if (ga_nodeid().eq.0) then
          write(*,1) iat,
     &               chcdata(1),            ! par-x
     &               chcdata(2),chcdata(3)  ! par-y,z
 1        format(' CHC PSOSO(xx,yy,zz)(',i3,')=(', 
     &            f15.8,',',f15.8,',',f15.8,')') 
        endif
      return
      end
      subroutine get2dmat(a,n,b)
c     a : size(a)=n*(n+1)/2) symmetric matrix
      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
       integer a,b,n,n1,i,j,count,
     &         l_mat,k_mat

        if (.not. ga_create(mt_dbl,n,n,
     &     'get2dmat: b',0,0,b)) 
     $    call errquit('get2dmat: b', 0,GA_ERR)
       call ga_zero(b)

       n1=n*(n+1)/2
       if (.not.ma_alloc_get(mt_dbl,n1,
     &            'fcsd',l_mat,k_mat))
     &  call errquit('get2dmat: ma failed',0,MA_ERR)       
        call ga_get(a,1,1,1,n1,
     &              dbl_mb(k_mat),1)      
       count=0
       do i=1,n
        call ga_put(b,i,i,i,i,dbl_mb(k_mat+count),1)
        count=count+1
       enddo
       do i=2,n
        do j=1,i-1
         call ga_put(b,i,i,j,j,dbl_mb(k_mat+count),1)
         call ga_put(b,j,j,i,i,dbl_mb(k_mat+count),1)   
         count=count+1
        enddo ! end-loop-j
       enddo ! end-loop-i
      if (.not.ma_free_heap(l_mat)) call
     &    errquit('get2dmat:: ma_free_heap l_mat',0, MA_ERR)           
      return
      end
c -------------- for NMLO analysis ----------------- END
c --------------------------------------------------------
