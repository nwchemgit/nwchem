c Subject : Implementation of NMR Hyperfine Tensors
c           using ZORA.
c Author  : Fredy Aquino
c Date    : 03-10-11
      subroutine get_NMRHFine_ZORA(rtdb,
     &                             g_dens_at,
     &                             nexc,
     &                             geom, 
     &                             ao_bas_han,
     &                             nbf,
     &                             focc,
     &                             noc,
     &                             ipol, 
     &                             g_densZ4)
       implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      integer rtdb    
      integer g_densZ4(3)
      integer g_sdens,g_dens_at(2) 
      integer geom
      integer ao_bas_han    
      integer nbf,noc(2),ipol,nexc
      double precision focc(nbf*ipol) ! contains occupation values if modified
      integer nat_slc,typeprop,nat,l_AtNr,k_AtNr
      integer alo(3),ahi(3),ld(2)
      logical status  
      integer i,j,type_nmrdata
      logical dft_zoraNMR_write
      integer ga_dia_hfine,ga_para1,
     &        ga_h01_hfine,ga_Fji_hfine
      logical is_atom,ofinite,Knucl
      double precision atmass,AtNr_dbl
      character*255 zorafilename

      external get_Htensor_slow,get_HFine_F1ji,
     &         dft_zoraNMR_write,
     &         get_slctd_atoms
      logical do_prntNMRCS
      if(.not.rtdb_get(rtdb,'zora:do_prntNMRCS',        ! FA
     &                 mt_log,1,do_prntNMRCS))          ! FA
     &  do_prntNMRCS= .false.   
c ------ Read Knucl   for including ONLY nuclear part in K ZORA ----- START
c Note.- stored in rel_input.F(rel_input(rtdb))
         Knucl=.false.
         status=rtdb_get(rtdb,'zora:Knucl',mt_log,1,Knucl) ! Check if gaussian nucl model requested
         if (ga_nodeid().eq.0)
     &     write(*,*) 'In get_NMRHFine_ZORA:: zora:Knucl=',Knucl
c ------ Read Knucl   for including ONLY nuclear part in K ZORA ----- END    
c ------ Read ofinite to be used by HFine finite calc ---FA-03-21-11-- START
c Note.- stored in geom_input.F (geom_input(rtdb))
         ofinite=.false.
         status=rtdb_get(rtdb,'prop:ofinite',mt_log,1,ofinite) ! Check if gaussian nucl model requested
c ------ Read ofinite to be used by HFine finite calc ---FA-03-21-11-- END
       if (ga_nodeid().eq.0) then
        write(*,*) 'dft_zora_Hypefine: ofinite=',ofinite
       endif
c ---- get spin-densty: g_sdens -------- START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'NMRHFine: g_sdens',
     $                      0,0,g_sdens))
     $       call errquit('NMRHFine: g_sdens', 0,
     &                    GA_ERR)
       call ga_add( 1.0d0,g_densZ4(1),
     &             -1.0d0,g_densZ4(2),g_sdens)
c ---- get spin-densty: g_sdens -------- END
c === allocate arrays to store diamagnetic tensor === START
c ------- Read (nat,atmnr) --------- START
         status=geom_ncent(geom,nat)   
      if (.not.ma_alloc_get(
     &       mt_int,nat,'nmt tmp',l_AtNr,k_AtNr))
     &    call errquit('dft_zora_Hyperfine: ma failed',0,MA_ERR)
         typeprop=3 ! =1 EFG =2 Shieldings =3 Hyperfine  
         call get_slctd_atoms(nat_slc,       ! out: selected atoms
     &                        int_mb(k_AtNr),! out: list of selected atom nr.     
     &                        nat,           ! in : total nr atoms in molecule            
     &                        rtdb,          ! in : rtdb  handle
     &                        typeprop)      ! in : =1,2,3=EFG,Shieldings,Hyperfine
         if (.not. ga_create(mt_dbl,1,nat_slc,
     &                  'dft_zora_Hyperfine: g_AtNr',0,0,g_AtNr))
     $     call errquit('dft_zora_Hyperfine: g_AtNr', 0,GA_ERR)

       do i=1,nat_slc
        AtNr_dbl=int_mb(k_AtNr+i-1)
        call ga_put(g_AtNr,1,1,i,i,AtNr_dbl,1)
       enddo
      if (ga_nodeid().eq.0) then
       write(*,*) 'nat_slc=',nat_slc
       do i=1,nat_slc
        write(*,7) i,int_mb(k_AtNr+i-1)
 7      format('In dft_zora_Hyperfine:: atomnr(',i3,')=',i5)
       enddo
      endif
c ------- Read (nat,atmnr) --------- END
       call get_Htensor_fast(
     &          ga_dia_hfine,  ! OUT: dia hyperfine tensor
     &          ga_h01_hfine,  ! OUT: h01 hyperfine munu AO matrix
     &          g_sdens,       ! IN : spin density
     &          ofinite,       ! IN : = .true. if requesting Gaussian Nucl. Model for charge
     &          Knucl,         ! in : = .true. for including ONLY nuclear part in K ZORA
     &          nat,nat_slc,int_mb(k_AtNr),
     &          rtdb,g_dens_at,nexc,
     &          geom,ao_bas_han,nbf,
     &          noc,ipol)
c ---Destroying ga arrays ------- START
         if (.not. ga_destroy(g_sdens)) call errquit(
     &     'get_NMRHFine_ZORA: ga_destroy failed ',0, GA_ERR)
      call get_HFine_F1ji(
     &          ga_Fji_hfine, ! OUT: munu-mat-Fji
     &          rtdb,g_dens_at,
     &          nat,       ! in: nr. atoms 
     &          ofinite,   ! in: = .true. if Gaussian Nucl. Model of charges requested
     &          Knucl,     ! in: = .true. if K_ZORA(V=Nuclear pot. only)     
     &          nexc,
     &          geom,ao_bas_han,nbf)
c       == get filename for the zora data ==
        type_nmrdata=2 ! =1,2,3=shieldings,hyperfine,gshift
c       Note.- lbl_nmrhyp defined in zora.fh
        call util_file_name(lbl_nmrhyp,.false.,.false.,zorafilename)
        if (.not.dft_zoraNMR_write(
     &          zorafilename,
     &          type_nmrdata, ! =1,2,3=shieldings,hyperfine,gshift
     &          nbf,
     &          nat_slc,
     &          g_AtNr,
     &          ga_dia_hfine,
     &          ga_para1,
     &          ga_h01_hfine,
     &          ga_Fji_hfine))
     &     call errquit('get_NMRHFine_ZORA: dft_zoraNMR_write failed',
     &                  0,DISK_ERR)
c ---- Destroy stored ga arrays ------ START
           if (.not. ga_destroy(ga_dia_hfine)) call errquit(
     &    'get_NMRHFine_ZORA: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_hfine)) call errquit(
     &    'get_NMRHFine_ZORA: ga_destroy failed ',0, GA_ERR)   
        if (.not. ga_destroy(ga_Fji_hfine)) call errquit(
     &    'get_NMRHFine_ZORA: ga_destroy failed ',0, GA_ERR)   
c        if (.not. ga_destroy(g_AtNr)) call errquit(
c     &    'get_NMRHFine_ZORA: ga_destroy failed ',0, GA_ERR)              
c ---- Destroy stored ga arrays ------ END
       if (.not.ma_free_heap(l_AtNr)) call
     &     errquit('get_NMRHFine_ZORA: ma_free_heap l_AtNr',0,MA_ERR)    
      return
      end

      subroutine get_Htensor_fast(
     &                  gFCSD,      ! OUT :      dia hyperfine tensor
     &                  gPSOSO,     ! OUT : h01 para hyperfine munu AO matrix
     &                  g_sdens,    ! IN  : spin density
     &                  ofinite,    ! IN  : = .true. if requesting Gaussian Nucl. Model for charge
     &                  Knucl,      ! IN  : = .true. for including ONLY nuclear part in K ZORA
     &                  nat,        ! IN  : total nr of atoms in molecule
     &                  nat_slc,    ! IN  : nr of selected atoms
     &                  atmnr_slc,  ! IN  : selected atom numbers
     &                  rtdb,
     &                  g_dens_at,
     &                  nexc,
     &                  geom,       ! IN  : geom  handle
     &                  ao_bas_han, ! IN  : basis handle 
     &                  nbf,        ! IN  : nr. basis functions
     &                  noc,        ! IN  : nr. occupied MOs
     &                  ipol)       ! IN  : nr. of polarizations
c    Purpose : Compute
c    1. A_{u,v}^{FC+SD} = 2 g_N B_N /(nA-nB) 
c                         \sum_{r,s} P_{r,s}^{A-B} < chi_r| h^{(u,v)}| chi_s>
c       GA array gFCSD (dimension 3 x 3 x nlist (nlist, nr of atoms selected)
c    2. h_{u}^{PSOSO}= < chi_r | h_{u}^{PSOSO}|chi_s
c    h_{u}^{PSOSO}= 2c^2 i . h_{Au,rs}^{ZPSO}=
c    \int dr K/r_A^3 (\vec{r}_A x[chi_r^* \nabla chi_s-\nabla chi_r^* chi_s])_u
c     (Eq. 56 in JA's write-up on Sept 17,2007)
c
c    Author : Fredy Aquino
c    Date   : 03-10-11
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      integer g_sdens,vectors(2) 
      integer gFCSD, ! out
     &        gPSOSO ! out
      integer rtdb    
      integer g_dens_at(2),g_densZ4(3)
      integer noc(2),noc1,ipol
      integer geom,ao_bas_han                      
      integer ispin,nexc,iat,iat1,
     &        nat,nat_slc,typeprop
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt
      integer atmnr_slc(nat_slc)
      logical status,ofinite,Knucl             
      character*16 at_tag  
      integer stat_read
      integer alo(3),ahi(3),ld(2)
      integer u,v,i,j,k,t,a,nbf
      integer dims(3),chunk(3)
      double precision val,factor
      double precision xyz_NMRQcoords(3),atmass,
     &                 zetanuc_slc 
      character*255 zorafilename
      integer g_fcsd(3,3),g_scr,g_t1,g_t2,g_t3
      integer g_zpso(6)  
      integer l_buf,k_buf,cbuf,
     &        l_buf1,k_buf1,
     &        l_zetanuc,k_zetanuc   
c +++++++++ definitions for NLMO analysis ++++++++ START
      integer ndir,ndir1,n_munu,count,count1,hypfile,
     &        g_munuFCSD,g_munuPSOSO, ! GA contains matrices
     &        g_c1,ndata
      logical dft_zoraHYP_NLMOAnalysis_write
c +++++++++ definitions for NLMO analysis ++++++++ END
      external util_wallsec
      double precision util_wallsec
      if(.not.ma_alloc_get(mt_dbl,3*3,'HTensor:buf',
     &                    l_buf,k_buf))
     &    call errquit('HTensor: ma failed',911,MA_ERR)
      if(.not.ma_alloc_get(MT_DBL,nbf*nbf,'HTensor',
     &                    l_buf1,k_buf1))
     $ call errquit('HTensor: ma failed',911, MA_ERR)
c +++++++++++++++creating ga_arrays ++++++++START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_t1',0,0,g_t1))
     $    call errquit('HFine: g_t1',0,GA_ERR)
         call ga_zero(g_t1)  
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_t2',0,0,g_t2))
     $    call errquit('HFine: g_t2',0,GA_ERR)
         call ga_zero(g_t2)   
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_t3',0,0,g_t3))
     $    call errquit('HFine: g_t3',0,GA_ERR)
         call ga_zero(g_t3)     
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_scr',0,0,g_scr))
     $    call errquit('HFine: g_scr',0,GA_ERR)
         call ga_zero(g_scr)    
      do i=1,6 
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_zpso',0,0,g_zpso(i)))
     $    call errquit('HFine: g_zpso',0,GA_ERR)
         call ga_zero(g_zpso(i))
      enddo
      do i=1,3
         do j=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_fcsd',0,0,g_fcsd(i,j)))
     $    call errquit('HFine: g_fcsd',0,GA_ERR)
         call ga_zero(g_fcsd(i,j))    
         enddo ! end-loop-j
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
c +++++ Read Atom Nr for NMR calc ++
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat_slc,
     &             'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('HFine: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat_slc,
     &             'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('HFine: ma failed',911,MA_ERR)
c === allocate arrays to store diamagnetic tensor === START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*nat_slc
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,gPSOSO)) call 
     &    errquit('get_d2p1: nga_create failed gPSOSO_num',0,GA_ERR)
      call ga_zero(gPSOSO)
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  nat_slc ! Total nr. atoms requested
      if (.not.nga_create(MT_DBL,3,ahi,'gFCSD matrix',
     &                    alo,gFCSD)) 
     &  call errquit('HFine: nga_create failed gFCSD all',
     &               0,GA_ERR)
c === allocate arrays to store diamagnetic tensor === END
       if (ofinite) then
         if (.not.ma_alloc_get(mt_dbl,nat,
     &                  'zetanuc1',l_zetanuc,k_zetanuc))
     &   call errquit('HFine: ma failed',0,MA_ERR)
         call get_zetanuc_arr(geom,nat,dbl_mb(k_zetanuc)) ! zetanuc_arr(i) i=1,natoms
         do iat = 1,nat ! == loop over the atoms ==
          if (ga_nodeid().eq.0)
     &     write(*,1) iat,dbl_mb(k_zetanuc+iat-1) 
 1        format('In get_Htensor_fast:: zetanuc_arr(',i3,')=',f35.8)
          dbl_mb(k_zetanuc+iat-1)=dsqrt(dbl_mb(k_zetanuc+iat-1)) ! Calc sqrt(zetanuc)
         enddo ! end-loop-iat
       endif
c +++++++++ NLMO analysis : create diag matrix ++++++++ START
      hypfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
       if (hypfile.eq.1) then ! --- if-hypfile--START
        ndir=6 ! xx,yy,zz,xy,xz,yz
        n_munu=nbf*(nbf+1)/2*ndir*nat_slc
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_Htensor_fast: g_munu',0,0,g_munuFCSD))
     $    call errquit('get_Htensor_fast:',0,GA_ERR)
        call ga_zero(g_munuFCSD)
        ndir1=3 ! x,y,z
        n_munu=nbf*(nbf+1)/2*ndir1*nat_slc
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_Htensor_fast: g_munu',0,0,g_munuPSOSO))
     $    call errquit('get_Htensor_fast:',0,GA_ERR)
        call ga_zero(g_munuPSOSO)
       endif ! ------------if-hypfile---------END
       count=1
       count1=1 ! for storing g_munu_PSOSO
c +++++++++ NLMO analysis : create diag matrix ++++++++ END
      do iat1=1,nat_slc  ! nat_slc <= nat
       iat=atmnr_slc(iat1)
       status=geom_cent_get(geom,iat,at_tag,
     &                      dbl_mb(k_xyzpt+3*(iat1-1)),
     &                      dbl_mb(k_zanpt+iat1-1))
       if(.not.geom_mass_get(geom, iat, atmass)) call
     &    errquit(' mass_get  failed ',iat, GEOM_ERR)
       call get_znuc(atmass,zetanuc_slc) 
       xyz_NMRQcoords(1)= dbl_mb(k_xyzpt  +3*(iat1-1))
       xyz_NMRQcoords(2)= dbl_mb(k_xyzpt+1+3*(iat1-1))
       xyz_NMRQcoords(3)= dbl_mb(k_xyzpt+2+3*(iat1-1))

       if (ga_nodeid().eq.0) then
          write(luout,'(a,f10.2,a)') ' starting at time ',
     &         util_wallsec(),'s'
        write(*,153) iat,ofinite,atmass,zetanuc_slc
 153    format('CHECK:(atom,ofinite,atmass,zetanucl_slc)=(',
     &         i4,',',l1,',',f15.8,',',f35.8,')')
       endif

       call zora_getv_HFine_fast(rtdb,g_dens_at, 
     &                           ofinite,
     &                           dbl_mb(k_zetanuc),
     &                           zetanuc_slc,
     &                           Knucl,
     &                           xyz_NMRQcoords,
     &                           g_zpso, ! ZPSO
     &                           g_fcsd, ! FC+SD (v,u) term
     &                           nexc) 
c ---- prepare g_fcsd-final --------- START
       call ga_zero(g_scr)
       do u=1,3
         call ga_add(1.0d0,g_scr,1.0d0,g_fcsd(u,u),g_scr)
       enddo
       do u=1,3
         call ga_add(-1.0d0,g_scr,+1.0d0,g_fcsd(u,u),
     &               g_fcsd(u,u))
       enddo
c ---- prepare g_fcsd-final ---------  END
       do u=1,3
c ---- 2nd method (fastest) to calculate g_zpso ---- START
c ----- g_h01 = < chi_{mu} | h_t^{01}| chi_{nu} >
c       h_t^{01}=K/(2c) (\vec{r} x \vec{p})_t/r_Q^3 +
c                       (\vec{r} x \vec{p})_t/r_Q^3 K/(2c)
        call ga_add(1.0d0,g_zpso(u)  ,
     &             -1.0d0,g_zpso(u+3),g_t1) ! i.e. u=1  mn-nm=23-32, etc.
        call ga_transpose(g_t1,g_t2)
        call ga_zero(g_t3)
        call ga_add(1.0d0,g_t2,-1.0d0,g_t1,g_t3) ! g_t2=g_h01
c ---- 2nd method (fastest) to calculate g_zpso ---- END
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=nbf
        alo(3)=3*(iat1-1)+u
        ahi(3)=3*(iat1-1)+u
        ld(1)=nbf
        ld(2)=nbf
        call ga_scale(g_t3,-0.5d0) ! including (-1/2) from h^{u} operator
c ------------ NLMO analysis para term (PSOSO) ------------- START
        hypfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
        if (hypfile.eq.1) then
         call fill_munuPSOSO(g_munuPSOSO, !out   : array with matrices (g_t3--> g_munuPSOSO)
     &                       count1,      !in/out: counting data stored in g_munuPSOSO
     &                       g_t3,        ! in:  nbf
     &                       nbf)
        endif
c ------------ NLMO analysis para term (PSOSO) ------------- END
        call ga_get(g_t3,1,nbf,1,nbf,dbl_mb(k_buf1),nbf)
        call nga_put(gPSOSO,alo,ahi,dbl_mb(k_buf1),ld) ! store gPSOSO_u
        do v=1,3
         val=ga_ddot(g_sdens,g_fcsd(u,v))
         cbuf=k_buf+(u-1)*3+v-1
         dbl_mb(cbuf)=val ! missing :2 g_N beta_N /(n_a-n_b) * (1/2)
                          ! Note.- (1/2) is from h^{(u,v)} operator
          if (ga_nodeid().eq.0) then
           write(*,10) iat,u,v,val
10         format('gFCSD(',i3,',',i3,',',i3,')=',f15.8)
          endif
        enddo ! end-loop-v
       enddo ! end-loop-u
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=iat1
       ahi(3)=iat1
       ld(1)=3
       ld(2)=3
       call nga_put(gFCSD,alo,ahi,dbl_mb(k_buf),ld)
c +++++++++ NLMO analysis : store diag matrix ++++++++ START
        hypfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
        if (hypfile.eq.1) then
c ----- Symmetrize g_fcsd() ------ START
        do u=1,3
         do v=u+1,3
          call ga_add(0.5d0,g_fcsd(u,v),
     &                0.5d0,g_fcsd(v,u),g_t1)
          call ga_copy(g_t1,g_fcsd(u,v))
          call ga_copy(g_t1,g_fcsd(v,u))
         enddo
        enddo
c ----- Symmetrize g_fcsd() ------ END
         call fill_munuFCSD(g_munuFCSD, !out: array with matrices (g_fcsd--> g_munuFCSD)
     &                      count,      !in/out: counting data stored in g_munuFCSD
     &                      g_fcsd,     ! in:  nbf
     &                      nbf) 
        endif
c +++++++++ NLMO analysis : store diag matrix ++++++++ END
      end do ! iat loop
c ---Destroying ga arrays ------- START
        if (.not. ga_destroy(g_t1)) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
        if (.not. ga_destroy(g_t2)) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_t3)) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_scr)) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
       do i=1,6
        if (.not. ga_destroy(g_zpso(i))) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
       enddo
       do i=1,3
        do j=1,3
        if (.not. ga_destroy(g_fcsd(i,j))) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
        enddo
       enddo
        hypfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
        if (hypfile.eq.1) then   
         ndata=1 !  =1 write FCSD,PSOSO,sdens =2 write g_c1
         call util_file_name(lbl_nlmohyp,.false.,.false.,zorafilename)
         if (.not.dft_zoraHYP_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munuPSOSO
     &            nat_slc, ! in: list of selected atoms 
     &                noc, ! in: dummy not used yet here
     &              ndata, ! in: =1 write FCSD,PSOSO,sdens =2 write g_c1
     &         g_munuFCSD, ! in: munu dia or Fermi Contact + Spin Dipolar term
     &        g_munuPSOSO, ! in: munu para or PSOSO term
     &               ipol, ! in: nr. of polarizations
     &            vectors, ! in: dummy not used yet here
     &               g_c1, ! in: dummy not used yet here
     &            g_sdens))! in: spin density
     &   call errquit('get_Htensor_fast: dft_zoraHYPNLMO_write failed',
     &                0,DISK_ERR)
         if (.not. ga_destroy(g_munuFCSD)) call errquit(
     &     'HFine: ga_destroy failed ',0, GA_ERR) 
         if (.not. ga_destroy(g_munuPSOSO)) call errquit(
     &     'HFine: ga_destroy failed ',0, GA_ERR) 
        endif
c ---Destroying ga arrays ------- END
c----deallocate memory 
       if (.not.ma_free_heap(l_buf1)) call errquit
     &    ('HFine, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('HFine, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('HFine:, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('HFine:, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (ofinite) then
        if (.not.ma_free_heap(l_zetanuc)) call
     &     errquit('HFine:: ma_free_heap l_zetanuc',0, MA_ERR)
      endif
      return
      end

      subroutine get_Htensor_slow(
     &                       gFCSD,    ! OUTPUT : Fermi Contact - Spin Dipolar
     &                       gPSOSO,   ! OUTPUT : 
     &                       g_sdens,  ! IN  : spin density
     &                       ofinite,  ! IN  : = .true. if requesting Gaussian Nucl. Model for charge
     &                       rtdb,g_dens_at,nexc,
     &                       geom, 
     &                       ao_bas_han,
     &                       nbf,
     &                       noc,
     &                       ipol)
c    Purpose : Compute
c    1. A_{u,v}^{FC+SD} = 2 g_N B_N /(nA-nB) 
c                         \sum_{r,s} P_{r,s}^{A-B} < chi_r| h^{(u,v)}| chi_s>
c       GA array gFCSD (dimension 3 x 3 x nlist (nlist, nr of atoms selected)
c    2. h_{u}^{PSOSO}= < chi_r | h_{u}^{PSOSO}|chi_s
c    h_{u}^{PSOSO}= 2c^2 i . h_{Au,rs}^{ZPSO}=
c    \int dr K/r_A^3 (\vec{r}_A x[chi_r^* \nabla chi_s-\nabla chi_r^* chi_s])_u
c     (Eq. 56 in JA's write-up on Sept 17,2007)
c
c    Author : Fredy Aquino
c    Date   : 11-06-10
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      integer g_sdens
      integer gFCSD, ! out
     &        gPSOSO ! out
      integer rtdb    
      integer g_dens_at(2),g_densZ4(3)    
      integer noc(2),noc1,ipol
      integer geom,ao_bas_han                      
      integer ispin,nexc,iat,iat1,nat
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt    
      integer l_AtNr,k_AtNr 
      logical status,ofinite             
      character*16 at_tag  
      integer stat_read,read_SLCTD_HFine_Atoms
      integer alo(3),ahi(3),ld(2)
      integer u,v,i,j,k,t,a,nbf
      integer dims(3),chunk(3)
      double precision val,factor
      double precision xyz_NMRQcoords(3),atmass
      integer g_fcsd(3,3),g_scr,g_t1,g_t2,g_t3
      integer g_zpso(3) 
      integer l_buf,k_buf,cbuf,
     &        l_buf1,k_buf1     
      external zora_getv_HFine,    
     &         read_SLCTD_HFine_Atoms,
     &         zora_getv_HFine1
      if(.not.ma_alloc_get(mt_dbl,3*3,'HTensor:buf',
     &                    l_buf,k_buf))
     &    call errquit('HTensor: ma failed',911,MA_ERR)
      if(.not.ma_alloc_get(MT_DBL,nbf*nbf,'HTensor',
     &                    l_buf1,k_buf1))
     $ call errquit('HTensor: ma failed',911, MA_ERR)
      status=geom_ncent(geom,nat) ! Get nat, # of atoms
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('HFine: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('HFine: ma failed',911,MA_ERR)
c +++++++++++++++creating ga_arrays ++++++++START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_scr',0,0,g_scr))
     $    call errquit('HFine: g_scr',0,GA_ERR)
         call ga_zero(g_scr)    
      do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_zpso',0,0,g_zpso(i)))
     $    call errquit('HFine: g_zpso',0,GA_ERR)
         call ga_zero(g_zpso(i))
      enddo
      do i=1,3
         do j=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'HFine: g_fcsd',0,0,g_fcsd(i,j)))
     $    call errquit('HFine: g_fcsd',0,GA_ERR)
         call ga_zero(g_fcsd(i,j))    
         enddo ! end-loop-j
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
c +++++ Read Atom Nr for NMR calc ++
        if (.not. ga_create(mt_dbl,1,nat,
     &   'HFine: g_AtNr',0,0,g_AtNr)) 
     $   call errquit('HFine: g_AtNr',0,GA_ERR)
        call ga_zero(g_AtNr)
       stat_read=read_SLCTD_HFine_Atoms(rtdb,nat,nlist,g_AtNr)
c === allocate arrays to store diamagnetic tensor === START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*nlist
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,gPSOSO)) call 
     &    errquit('get_d2p1: nga_create failed gPSOSO_num',0,GA_ERR)
      call ga_zero(gPSOSO)
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  nlist ! Total nr. atoms requested
      if (.not.nga_create(MT_DBL,3,ahi,'gFCSD matrix',
     &                    alo,gFCSD)) 
     &  call errquit('HFine: nga_create failed gFCSD all',
     &               0,GA_ERR)
c === allocate arrays to store diamagnetic tensor === END
c  Allocate memory for l_AtNr,k_AtNr 
      if (.not.ma_alloc_get(mt_dbl,nat,
     &  'AtNr',l_AtNr,k_AtNr))
     &  call errquit('HFine: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
      do iat1=1,nlist  ! nlist <= nat
       iat=dbl_mb(k_AtNr+iat1-1)
       status=geom_cent_get(geom,iat,at_tag,
     &                      dbl_mb(k_xyzpt+3*(iat-1)),
     &                      dbl_mb(k_zanpt+iat-1))
         if(.not.geom_mass_get(geom, iat, atmass)) call
     &        errquit(' mass_get  failed ',iat, GEOM_ERR)
       xyz_NMRQcoords(1)= dbl_mb(k_xyzpt  +3*(iat-1))
       xyz_NMRQcoords(2)= dbl_mb(k_xyzpt+1+3*(iat-1))
       xyz_NMRQcoords(3)= dbl_mb(k_xyzpt+2+3*(iat-1))

       if (ga_nodeid().eq.0) then
        write(*,153) iat,ofinite,atmass
 153    format('CHECK:(atom,ofinite,atmass)=(',
     &         i4,',',l1,',',f15.8,')')
       endif

       call zora_getv_HFine_slow(rtdb,g_dens_at, 
     &                           ofinite,
     &                           atmass,
     &                           xyz_NMRQcoords,
     &                           g_zpso, ! ZPSO
     &                           g_fcsd, ! FC+SD (v,u) term
     &                           nexc) 
c ---- prepare g_fcsd-final --------- START
      call ga_zero(g_scr)
      do u=1,3
        call ga_add(1.0d0,g_scr,1.0d0,g_fcsd(u,u),g_scr)
      enddo
      do u=1,3
        call ga_add(-1.0d0,g_scr,+1.0d0,g_fcsd(u,u),
     &              g_fcsd(u,u))
      enddo
c ---- prepare g_fcsd-final ---------  END
       do u=1,3
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=nbf
        alo(3)=3*(iat1-1)+u
        ahi(3)=3*(iat1-1)+u
        ld(1)=nbf
        ld(2)=nbf
        call ga_scale(g_zpso(u),0.5d0)
        call ga_get(g_zpso(u),1,nbf,1,nbf,dbl_mb(k_buf1),nbf)
        call nga_put(gPSOSO,alo,ahi,dbl_mb(k_buf1),ld) ! store g_h01
        do v=1,3
         val=ga_ddot(g_sdens,g_fcsd(u,v))
         cbuf=k_buf+(u-1)*3+v-1
         dbl_mb(cbuf)=val
         if (ga_nodeid().eq.0) then
          write(*,10) iat,u,v,val
10        format('gFCSD(',i3,',',i3,',',i3,')=',f15.8)
         endif
        enddo ! end-loop-v
       enddo ! end-loop-u

       goto 17

       do j=1,3
       if (ga_nodeid().eq.0)
     &  write(*,1) iat,j
 1      format('-----ZPSO munu-mat(',i3,',',i3,')----- START')
        call ga_print(g_zpso(j))
       if (ga_nodeid().eq.0)
     &  write(*,2) iat,j
 2      format('-----ZPSO munu-mat(',i3,',',i3,')----- END')
       enddo ! end-loop-j (directions xyz)
       do j=1,3
        do i=1,3
       if (ga_nodeid().eq.0)
     &  write(*,3) iat,j,i
 3      format('-----FC+SD munu-mat(',i3,',',i3,',',i3,')----- START')
        call ga_print(g_fcsd(j,i))  
       if (ga_nodeid().eq.0)
     &  write(*,4) iat,j,i
 4      format('-----FC+SD munu-mat(',i3,',',i3,',',i3,')-----END')    
        enddo       
       enddo

 17    continue

      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=iat1
      ahi(3)=iat1
      ld(1)=3
      ld(2)=3
      call nga_put(gFCSD,alo,ahi,dbl_mb(k_buf),ld)

      end do ! iat loop
c      if (ga_nodeid().eq.0)
c    &  write(*,*) '------- gFC+SD tensor ------ START'
c      call ga_print(gFCSD)
c      if (ga_nodeid().eq.0)
c    &  write(*,*) '------- gFC+SD tensor ------ END'
c      if (ga_nodeid().eq.0)
c    &  write(*,*) '------- gPSOSO tensor ------ START'
c      call ga_print(gPSOSO)
c      if (ga_nodeid().eq.0)
c    &  write(*,*) '------- gPSOSO tensor ------ END'

c ---Destroying ga arrays ------- START
        if (.not. ga_destroy(g_scr)) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
       do i=1,3
        if (.not. ga_destroy(g_zpso(i))) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
       enddo
       do i=1,3
        do j=1,3
        if (.not. ga_destroy(g_fcsd(i,j))) call errquit(
     &    'HFine: ga_destroy failed ',0, GA_ERR) 
        enddo
       enddo
c ---Destroying ga arrays ------- END
c----deallocate memory 
       if (.not.ma_free_heap(l_buf1)) call errquit
     &    ('HFine, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('HFine, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('HFine:, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('HFine:, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('HFine:: ma_free_heap l_AtNr',0, MA_ERR)
      return
      end

      subroutine get_HFine_F1ji(
     &          g_Fji,        !out: munu-mat-Fji
     &          rtdb,g_dens_at,
     &          nat,          ! in: nr. atoms 
     &          ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &          Knucl,        ! in: = .true. if K_ZORA(V=Nuclear pot. only)     
     &          nexc,
     &          geom,ao_bas_han,nbf)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
      integer g_Fji ! OUTPUT
      logical status        
      integer rtdb,g_dens_at(2),geom,ao_bas_han       
      integer nexc            
      integer i,j,k   
      integer nbf
      integer g_hfine(3)
      integer l_buf,k_buf
      integer alo(3),ahi(3),ld(2)
      integer iat,nat
      integer l_zetanuc,k_zetanuc
      logical ofinite,Knucl
      double precision zetanuc_arr(nat)
      external zora_getv_NMRHFine_F1ji,
     &         get_zetanuc_arr

       if (ofinite) then
         if (.not.ma_alloc_get(mt_dbl,nat,
     &                  'zetanuc1',l_zetanuc,k_zetanuc))
     &   call errquit('HFine: ma failed',0,MA_ERR)
         call get_zetanuc_arr(geom,nat,dbl_mb(k_zetanuc)) ! zetanuc_arr(i) i=1,natoms
         do iat = 1,nat ! == loop over the atoms ==
          if (ga_nodeid().eq.0)
     &     write(*,1) iat,dbl_mb(k_zetanuc+iat-1) 
 1        format('In get_Htensor_fast:: zetanuc_arr(',i3,')=',f35.8)
          dbl_mb(k_zetanuc+iat-1)=dsqrt(dbl_mb(k_zetanuc+iat-1)) ! Calc sqrt(zetanuc)
         enddo ! end-loop-iat
       endif
c +++++++++++++++creating ga_arrays ++++++++START
      do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gFji: g_hfine',0,0,g_hfine(i)))
     $    call errquit('gFji: g_hfine',0,GA_ERR)
         call ga_zero(g_hfine(i))
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
       call zora_getv_NMRHFine_F1ji(rtdb,g_dens_at, 
     &                              g_hfine,
     &                              nat,
     &                              ofinite,
     &                              dbl_mb(k_zetanuc),
     &                              Knucl,
     &                              nexc) 
      if(.not.ma_alloc_get(mt_dbl,nbf*nbf,'gFji:buf',
     &                    l_buf,k_buf))
     $     call errquit('gFji: ma failed',911, MA_ERR)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'Fji matrix',alo,g_Fji)) 
     &    call 
     &    errquit('gFji: nga_create failed g_Fji',
     &            0,GA_ERR)
      call ga_zero(g_Fji)
      do k=1,3   
       alo(1)=1
       ahi(1)=nbf
       alo(2)=1
       ahi(2)=nbf
       alo(3)=k
       ahi(3)=k
       ld(1) =nbf
       ld(2) =nbf
       call ga_get(g_hfine(k),1,nbf,1,nbf,dbl_mb(k_buf),nbf)
       call nga_put(g_Fji,alo,ahi,dbl_mb(k_buf),ld)
      enddo ! end-loop-k
c ---Destroying ga arrays ----- START
       do i=1,3
        if (.not. ga_destroy(g_hfine(i))) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR) 
       enddo
c ---Destroying ga arrays ----- END
c----deallocate memory
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('gFij, ma_free_heap of l_buf failed',911,MA_ERR)
      if (ofinite) then
        if (.not.ma_free_heap(l_zetanuc)) call
     &     errquit('HFine:: ma_free_heap l_zetanuc',0, MA_ERR)
      endif
      return
      end
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++ NLMO analysis routines +++++++++++++++++ START
      subroutine fill_munuFCSD(g_munuFCSD, !out: array with matrices
     &                         count     , !in/out: counting data stored in g_munuFCSD
     &                         g_fcsd,     ! in:  nbf
     &                         nbf)        ! in: nr. basis functions
      ! Purpose: g_fcsd(u,v) --> g_munuFCSD
      ! Note: g_fcsd(u,v) corresponds to ith atom
      !       g_munuFCSD  contains matrices for all atoms
      ! g_munuFCSD contains unique munu
      ! in sequence: xx,yy,zz,xy,xz,yz
      ! each chunk is nbf*(nbf+1)/2 in size
      ! Author: Fredy W. Aquino 
      ! Date  : 06-29-11
      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"

       integer u,v,ii,i,j,k,indx,indx1
       integer nbf,count
       double precision val
       integer g_munuFCSD,g_fcsd(3,3)
       integer l_fcsd,k_fcsd
       integer ndir ! in: =6 = xx,yy,zz,xy,xz,yz
       integer nlist_uv(2,6)  ! ndir=6
       data nlist_uv /1,1,  ! xx
     &                2,2,  ! yy
     &                3,3,  ! zz
     &                1,2,  ! xy
     &                1,3,  ! xz
     &                2,3 / ! yz
       ndir=6 !  = xx,yy,zz,xy,xz,yz
       if (.not.ma_alloc_get(mt_dbl,nbf*nbf,
     &            'fcsd',l_fcsd,k_fcsd))
     &  call errquit('get_munuFCSD: ma failed',0,MA_ERR)
       do ii=1,ndir
         u=nlist_uv(1,ii)
         v=nlist_uv(2,ii)
         call ga_get(g_fcsd(u,v),1,nbf,1,nbf,
     &               dbl_mb(k_fcsd),nbf)
          do i=1,nbf
           indx=k_fcsd+nbf*(i-1)+i-1
           val=dbl_mb(indx)
           call ga_fill_patch(g_munuFCSD,1,1,count,count,val)
           count=count+1
          enddo ! end-loop-i
          do i=2,nbf
           do j=1,i-1
            indx=k_fcsd+nbf*(j-1)+i-1
            val=dbl_mb(indx)
            call ga_fill_patch(g_munuFCSD,1,1,count,count,val)
            count=count+1
           enddo ! end-loop-j
          enddo ! end-loop-i
       enddo ! end-loop-ii
      if (.not.ma_free_heap(l_fcsd)) call
     &    errquit('fill_munuFCSD:: ma_free_heap l_fcsd',0, MA_ERR)
      return
      end

      subroutine fill_munuPSOSO(g_munuPSOSO, !out: array with matrices
     &                          count      , !in/out: counting data stored in g_munuFCSD
     &                          g_psoso    , ! in:  nbf
     &                          nbf)         ! in: nr. basis functions
      ! Purpose: g_psoso --> g_munuPSOSO
      ! Note: g_psoso      corresponds to uth component of ith atom (g_psoso_u u=1,2,3)
      !       g_munuPSOSO  contains matrices for all atoms
      ! g_munuPSOSO contains unique munu elements
      ! main diag. elements + off diagonal elements
      ! Author: Fredy W. Aquino 
      ! Date  : 07-04-11
      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"

       integer i,j,indx
       integer nbf,count
       double precision val
       integer g_munuPSOSO,g_psoso
       integer l_psoso,k_psoso

       if (.not.ma_alloc_get(mt_dbl,nbf*nbf,
     &            'fcsd',l_psoso,k_psoso))
     &    call errquit('get_munuPSOSO: ma failed',0,MA_ERR)
         call ga_get(g_psoso,1,nbf,1,nbf,
     &               dbl_mb(k_psoso),nbf)
          do i=1,nbf
           indx=k_psoso+nbf*(i-1)+i-1
           val=dbl_mb(indx)
           call ga_fill_patch(g_munuPSOSO,1,1,count,count,val)
           count=count+1
          enddo ! end-loop-i
          do i=2,nbf
           do j=1,i-1
            indx=k_psoso+nbf*(j-1)+i-1
            val=dbl_mb(indx)
            call ga_fill_patch(g_munuPSOSO,1,1,count,count,val)
            count=count+1
           enddo ! end-loop-j
          enddo ! end-loop-i
      if (.not.ma_free_heap(l_psoso)) call
     &    errquit('fill_munuPSOSO:: ma_free_heap l_psoso',0, MA_ERR)
      return
      end

      subroutine fill_munuPSOSO_1(
     &             g_munuPSOSO   , ! in: array with matrices
     &             g_munuPSOSO2d , !out: nbf x nbf x 3 munu matrix for ith atom
     &             iat,            ! in: atom index = 1, 2, nlist
     &             type_symm     , ! in: =1 symmetric =2 antisymmetric
     &             nbf)            ! in: nr. basis functions
      ! Purpose: g_munuPSOSO --> g_munuPSOSO2d
      ! Note: g_munuPSOSO2d      corresponds to uth component of ith atom (g_munuPSOSO2d_u u=1,2,3)
      !       g_munuPSOSO  contains matrices for all atoms
      ! g_munuPSOSO contains unique munu elements
      ! main diag. elements + off diagonal elements
      ! Author: Fredy W. Aquino 
      ! Date  : 07-04-11
      ! from (upper/lower) triang matrix --> full symmetric matrix
      ! g_munuPSOSO = (11) (22) (33) ...(nn) (21) (31) (32) (41) (42) (43) ...(n (n-1))
      !       n_elements(g_munuPSOSO) = n (n-1)/2    n=nbf nr basis functions
      ! g_munuPSOSO2d is n x n antisymmetric matrix
      ! g_munuPSOSO ---> g_munuPSOSO2d
      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"

       integer i,j,indx,alo(3),ahi(3)
       integer nbf,nlst,iat,shift,ndir,ntot,xyz,type_symm
       double precision val,valneg
       integer g_munuPSOSO,g_munuPSOSO2d
       integer l_mat,k_mat

       if (.not.(type_symm.eq.1 .or. 
     &           type_symm.eq.2)) then
        if (ga_nodeid().eq.0)
     &    write(*,*) 'Error in fill_munuPSOSO_1: type_symm ne 1 or 2'
          stop
       endif

       nlst=nbf*(nbf+1)/2
       ndir=3 ! x, y, z
       ntot=nlst*ndir
       shift=ntot*(iat-1)
       if (.not.ma_alloc_get(mt_dbl,ntot,
     &            'fcsd',l_mat,k_mat))
     &    call errquit('get_munuPSOSO: ma failed',0,MA_ERR)
        call ga_get(g_munuPSOSO,1,1,shift+1,shift+ntot,
     &              dbl_mb(k_mat),1)
         indx=k_mat
         do xyz=1,ndir  
          do i=1,nbf
           val=dbl_mb(indx)
           alo(1)=i
           ahi(1)=i
           alo(2)=i
           ahi(2)=i
           alo(3)=xyz
           ahi(3)=xyz
           call nga_fill_patch(g_munuPSOSO2d,alo,ahi,val)
           indx=indx+1
          enddo ! end-loop-i
          do i=2,nbf
           do j=1,i-1
            val=dbl_mb(indx)
            alo(1)=i
            ahi(1)=i
            alo(2)=j
            ahi(2)=j
            alo(3)=xyz
            ahi(3)=xyz
            call nga_fill_patch(g_munuPSOSO2d,alo,ahi,val)
            if      (type_symm.eq.1) then
             valneg= val
            else if (type_symm.eq.2) then
             valneg=-val
            endif
            alo(1)=j
            ahi(1)=j
            alo(2)=i
            ahi(2)=i
            alo(3)=xyz
            ahi(3)=xyz
            call nga_fill_patch(g_munuPSOSO2d,alo,ahi,valneg)
            indx=indx+1
           enddo ! end-loop-j
          enddo ! end-loop-i
         enddo ! end-loop-xyz
      if (.not.ma_free_heap(l_mat)) call
     &    errquit('fill_munuPSOSO_1:: ma_free_heap l_mat',0, MA_ERR)
           
      return
      end

      logical function dft_zoraHYP_NLMOAnalysis_write(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr of directions: 3 = x y z for g_munuPSOSO
     &              nlist, ! in: list of selected atoms 
     &               nocc, ! in: nocc(i) i=1,2 nr. occupations
     &              ndata, ! in: writing order =1,2
     &         g_munuFCSD, ! in: munu dia
     &        g_munuPSOSO, ! in: munu para or PSOSO term
     &               npol, ! in: nr. of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: perturbed MO coeffs
     &            g_sdens) ! in: spin density
c Description: Collecting three matrices to be used
c              in wefgfile(rtdb) and wnbofile(rtdb)
c              Those routines are called in prop.F
c              after hnd_property(rtdb)
c              The info collected is:
c              g_munuFCSD, FCSD munu matrix 
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nlist,  ! = nr. slc atoms
     &        nmat,   ! = nlst*ndir*nlist  (ndir=6)
     &        nmat2,  ! = nlst*ndir1*nlist (ndir1=3)
     &        nmat1,  ! = nbf*nbf
     &        npol
      integer g_munuFCSD,g_munuPSOSO,g_sdens,
     &        vectors(npol),g_c1    
      integer ndir,ndir1,nbf,nlst,ndata,ntot,nocc(2)
      integer unitno
      parameter (unitno = 77)
      integer l_mat ,k_mat,
     &        l_mat1,k_mat1,
     &        l_mat2,k_mat2,
     &        l_c1,k_c1,
     &        l_mo,k_mo
      integer ok,iset,i,j,alo(3),ahi(3),ld(2)
      integer inntsize

      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Allocate the temporary buffer
      if (ndata.eq.1) then ! First time writing
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
       nlst=nbf*(nbf+1)/2
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writehyp',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_writehyp: k_mat failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1*nlist
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat2,'dft_zoraNLMO_writehyp',
     &        l_mat2,k_mat2))
     $  call errquit('dft_zoraNLMO_writehyp: k_mat2 failed', 
     &               nmat2, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writehyp',
     &        l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_writehyp: k_mat1 failed', 
     &               nmat1, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
c     Open the file - 1st time
        open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000)
c     Write out the number of sets and basis functions
        write(unitno, err=1001) nbf
        write(unitno, err=1001) nlst
        write(unitno, err=1001) ndir
        write(unitno, err=1001) ndir1
        write(unitno, err=1001) nlist
        call ga_get(g_munuFCSD,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_munuPSOSO,1,1,1,nmat2,
     &              dbl_mb(k_mat2),1)
        call swrite(unitno,dbl_mb(k_mat2),nmat2)
        call ga_get(g_sdens,1,nbf,1,nbf,
     &              dbl_mb(k_mat1),nbf)
        call swrite(unitno,dbl_mb(k_mat1),nmat1)
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
       if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_writehyp: l_mat free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat1))
     $  call errquit('dft_zoraNLMO_writehyp: l_mat1 free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat2))
     $  call errquit('dft_zoraNLMO_writehyp: l_mat2 free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
      else if (ndata.eq.2) then
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_writehyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_writehyp: l_mo failed', 
     &               nbf,MA_ERR)
        ndir = 3 ! x,y,z
        ntot = nocc(1)+nocc(2)    
        nmat = nbf*ndir*ntot
        if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writehyp',l_c1,k_c1))
     $   call errquit('dft_zoraNLMO_writehyp: ma failed', 
     &                nmat, MA_ERR)
c     Open the file - 2nd time
        open(unitno, status='unknown', form='unformatted',
     $       file=filename, err=1000,position='append')
        write(unitno, err=1001) nocc(1)
        write(unitno, err=1001) nocc(2)  
        write(unitno, err=1001) ntot
        write(unitno, err=1001) nmat
        write(unitno, err=1001) npol
        write(unitno, err=1001) nbf
c ----- Add MOs in file ----- START
        do i=1,npol
         do j=1,nbf
         call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
         call ga_get(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
         call swrite(unitno,dbl_mb(k_mo),nbf)         
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Add MOs in file ----- END
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=ntot
        alo(3)=1
        ahi(3)=3
        ld(1)=nbf
        ld(2)=ntot
        call nga_get(g_c1,alo,ahi,dbl_mb(k_c1),ld)
        call swrite(unitno,dbl_mb(k_c1),nmat)
        if (.not. ma_free_heap(l_mo))
     $   call errquit('dft_zoraNLMO_writehyp: ma free_heap failed', 
     &               911, MA_ERR)
        if (.not. ma_free_heap(l_c1))
     $   call errquit('dft_zoraNLMO_writehyp: ma free_heap failed', 
     &               911, MA_ERR)
      endif
c     Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraHYP_NLMOAnalysis_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA NLMO HYP data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraNLMO_writehyp: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNLMO_writehyp: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNLMO_writehyp: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraHYP_NLMOAnalysis_read(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr od directions: 3 = x y z for g_munuPSOSO
     &              nlist, ! in: list of selected atoms 
     &               nocc, ! in: nocc(i) i=1,2
     &               npol, ! in: nr polarizations
     &         g_munuFCSD, ! out: munu dia
     &        g_munuPSOSO, ! out: munu para
     &            vectors, ! out: MOs
     &               g_c1, ! out: perturbed MO coeffs
     &            g_sdens) ! out: spin density
      implicit none
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

      character*(*) filename    ! [input] File to write to
      integer nmat,nmat1,nmat2,npol
      integer g_munuFCSD,g_munuPSOSO,g_sdens,
     &        vectors(npol),g_c1
      integer nbf,nbf_read,nlst,nlst_read,
     &        ndir,ndir_read,
     &        ndir1,ndir1_read,   
     &        nlist,nlist_read,
     &        ntot,ntot_read,
     &        nocc(2),nocc_read(2),
     &        n_c1,n_c1_read,
     &        npol_read,
     &        alo(3),ahi(3),ld(2)
      integer unitno
      parameter (unitno = 77)
      integer l_mat ,k_mat,
     &        l_mat1,k_mat1,
     &        l_mat2,k_mat2,
     &        l_c1,k_c1,
     &        l_mo,k_mo
      integer ok,iset,i,j
      integer inntsize
c     Initialise to invalid MA handle
      nlst=nbf*(nbf+1)/2
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      nmat=nlst*ndir*nlist
       if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munuFCSD',0,0,g_munuFCSD)) 
     $   call errquit('dft_zoraNLMO_read: g_munuFCSD',0,GA_ERR)
        call ga_zero(g_munuFCSD)  
      nmat2=nlst*ndir1*nlist
       if (.not. ga_create(mt_dbl,1,nmat2,
     &   'dft_zoraNLMO_read: g_munuPSOSO',0,0,g_munuPSOSO)) 
     $   call errquit('dft_zoraNLMO_read: g_munuPSOSO',0,GA_ERR)
        call ga_zero(g_munuPSOSO)   
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_sdens',0,0,g_sdens)) 
     $   call errquit('dft_zoraNLMO_read: g_sdens',0,GA_ERR)
        call ga_zero(g_sdens)  
       do i=1,npol
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_sdens',0,0,vectors(i))) 
     $   call errquit('dft_zoraNLMO_read: vectors',0,GA_ERR)
        call ga_zero(vectors(i))  
       enddo
       ntot=nocc(1)+nocc(2)
       n_c1=nbf*3*ntot
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = ntot
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'c1 matrix',alo,g_c1)) call 
     &     errquit('dft_zoraNMR_read: nga_create failed g_c1',
     &             0,GA_ERR)
       call ga_zero(g_c1)
 
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA HYP NLMO data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nlst_read
       read(unitno, err=1001, end=1001) ndir_read
       read(unitno, err=1001, end=1001) ndir1_read
       read(unitno, err=1001, end=1001) nlist_read
c      Error checks
       if ((nbf_read      .ne. nbf)   .or.
     &     (nlst_read     .ne. nlst)  .or.
     &     (ndir_read     .ne. ndir)  .or.
     &     (ndir1_read    .ne. ndir1) .or.
     &     (nlist_read    .ne. nlist)) goto 1003
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat2, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat2,k_mat2))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat2, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(mt_dbl,nmat1, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat1, MA_ERR)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munuFCSD,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat2),nmat2)
       call ga_put(g_munuPSOSO,1,1,1,nmat2,dbl_mb(k_mat2),1)
       call sread(unitno,dbl_mb(k_mat1),nmat1)
       call ga_put(g_sdens,1,nbf,1,nbf,dbl_mb(k_mat1),nbf)
      if (.not. ma_free_heap(l_mat))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat1))      ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat2))      ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c -------- Read Perturbed MO coeffs ------------- START
       read(unitno, err=1001, end=1001) nocc_read(1)
       read(unitno, err=1001, end=1001) nocc_read(2)
       read(unitno, err=1001, end=1001) ntot_read
       read(unitno, err=1001, end=1001) n_c1_read
       read(unitno, err=1001, end=1001) npol_read
       read(unitno, err=1001, end=1001) nbf_read
c      Error checks
       if ((nocc_read(1) .ne. nocc(1)) .or.
     &     (nocc_read(2) .ne. nocc(2)) .or.
     &     (ntot_read    .ne. ntot)    .or.
     &     (nbf_read     .ne. nbf)     .or.
     &     (npol_read    .ne. npol)    .or.
     &     (n_c1_read    .ne. n_c1)) goto 1003
c ----- Read MOs ----- START
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_readhyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_readhyp: ma failed', 
     &               nbf,MA_ERR)
        do i=1,npol
         do j=1,nbf
          call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
          call sread(unitno,dbl_mb(k_mo),nbf)   
          call ga_put(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Read MOs ----- END
       if (.not. ma_alloc_get(mt_dbl,n_c1, ! allocate memory
     &    'dft_zoraNLMO_read',l_c1,k_c1))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               n_c1, MA_ERR)
! Note: n_c1 = nbf*3*ntot    [ ntot=nocc(1)+nocc(2) ]
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=ntot
        alo(3)=1
        ahi(3)=3
        ld(1) =nbf
        ld(2) =ntot
        call sread(unitno,dbl_mb(k_c1),n_c1)
        call nga_put(g_c1,alo,ahi,dbl_mb(k_c1),ld)
      if (.not. ma_free_heap(l_mo))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_c1))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c -------- Read Perturbed MO coeffs ------------- END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraHYP_NLMOAnalysis_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraNLMO_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNLMO_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 'dft_zoraNLMO_read: file inconsistent',
     &           ' with calculation',
     $           filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNLMO_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c ++++++++++++ NLMO analysis routines +++++++++++++++++ END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function read_SLCTD_HFine_Atoms
     &                 (rtdb,nat,nlist,g_AtNr)
c---- GA output: g_AtNr   

      implicit none
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
#include "rtdb.fh" 
#include "context.fh"

      integer rtdb,ii,nlist,nat,g_AtNr
      integer atomnr(nat)
      integer hfineatoms
      double precision AtNr_dbl

      if (.not. rtdb_get(rtdb, 'hfine:natoms',mt_int,
     $                     1,hfineatoms))
     &  hfineatoms=0 ! reset       
         if (hfineatoms.eq.0) then
          hfineatoms=nat
           nlist=hfineatoms
           do ii=1,hfineatoms
            AtNr_dbl=ii
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo     
         else
          if (.not. rtdb_get(rtdb, 'hfine:atom list',mt_int,
     $                     hfineatoms,atomnr))
     $      call errquit('prop_input-hfine: rtdb_get failed',
     $                   555, RTDB_ERR)
           nlist=hfineatoms
           do ii=1,hfineatoms
            AtNr_dbl=atomnr(ii)
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo
         endif 
       read_SLCTD_HFine_Atoms = 1
       return
       end   
c
      subroutine print_NMRHypFine_version()
c
      implicit none
#include "stdio.fh"
c
      write(LuOut,*)
      call util_print_centered(LuOut,
     $   'ZORA NMR Hyperfine Spin-Spin Coupling Constants', 23, .true.)
      write(LuOut,*)
c
      return
      end

c $Id$
