c=======================================================
c================ Fredy Aquino's routines ======== START
c=======================================================
      subroutine get_rhoS(rtdb,g_dens_ini,nexc,
     &                    geom, 
     &                    ao_bas_han,
     &                    nbf,nbf_ao,
     &                    noc,
     &                    ipol)
c -- Purpose: Calculation of small component density, \rho_S
c             for Electric Field Gradient calculation
c             Source: van Lenthe, et.al.,JCP, V112,N19,Y2000
c -- Author : Fredy Aquino 12-07-09

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 

      integer rtdb,geom,ao_bas_han   ! handles
      integer g_efgz4(2),g_dens_ini(2),g_densZ4(3),
     &        g_zora_scale_munu(2)  
      integer noc(2),noc1,ipol,ispin,pos_dat,
     &        nbf,nbf_ao,pos,ndir,
     &        do_zgc_old,nexc
      integer iat,iat1,nat,ipolmunu,zora_Qpq  
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt,    
     &        l_ilst,k_ilst,     
     &        l_jlst,k_jlst,  
     &        l_AtNr,k_AtNr,
     &        nlst,i,j,count
      integer g_rhoS,g_munu_rhoS,g_munuV6,
     &        efgfile
      logical dft_zoraEFGZ4_write
c     g_munu_rhoS, contains xx,yy,zz,xy,xz,yz
c                  munu half-triangle matrices
c                  including main diagonal
      logical status               
      character*16 element, at_tag  
      character*255 zorafilename
      double precision scf_dbl,zora_eint,toac,
     &                 xyz_EFGcoords(3)  
      integer stat_read,read_SLCTD_EFG_Atoms
      logical dft_zoraEFGZ4_NLMOAnalysis_write
      external get_densZ4,read_SLCTD_EFG_Atoms,
     &         zora_getv_EFGZ4_SR,dft_zoraEFGZ4_write,
     &         util_file_name,
     &         dft_zoraEFGZ4_NLMOAnalysis_write
c     Note.- g_munu_rhoS is created in dft_zora_rhos.F
c            g_munuV6    is created in hnd_elfcon_symm.F
c            g_munuV6 is not created/defined here yet     
c
c---- Input global variables (defined in zora.fh):
c     1. zora_calc_type  ! =3 -> ZORA4-EFG
c     2. so_term         ! =0 -> ZORA4-spin-free
c     3. xyz_EFGcoords(i) i=1,2,3
c     4. zora_Qpq=1,6 : xx,yy,zz,xy,xz,yz 
      if (ipol.eq.1) then
        scf_dbl=2.0d00
        noc1=noc(1)
      else if (ipol.eq.2) then 
        scf_dbl=1.0d00   
        noc1=noc(1)+noc(2)    
      endif
      status=geom_ncent(geom,nat)
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
c----- Allocate global arrays - FA
        do i=1,ipol
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'get_rhoS: g_efgz4',0,0,g_efgz4(i)))
     $       call errquit('get_rhoS: g_efgz4', 0,GA_ERR)
             call ga_zero(g_efgz4(i))
        enddo
c +++++ Read Atom Nr for EFG calc ++
         if (.not. ga_create(mt_dbl,1,nat,
     &                  'get_rhoS: g_AtNr',0,0,g_AtNr)) 
     $     call errquit('get_rhoS: g_AtNr', 0,GA_ERR)
        call ga_zero(g_AtNr)
       stat_read=read_SLCTD_EFG_Atoms
     &            (rtdb,nat,nlist,g_AtNr)
        if (ga_nodeid().eq.0)
     &      write(*,*) "==> nlist=",nlist
c  Allocate memory for l_AtNr,k_AtNr 
      if (.not.ma_alloc_get(mt_dbl,nat,'AtNr',l_AtNr,k_AtNr))
     &    call errquit('get_rhoS: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
c ------------ for NLMO analysis --------- START
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
      if (efgfile.eq.1) then ! ============= efgfile-if-START
c  Allocate memory for (l_ilst,k_ilst) (l_jlst,k_jlst)
      nlst=nbf*(nbf+1)/2  ! single triangle matrix
        if (.not.ma_alloc_get(mt_int,nlst,'ijlst',
     &                        l_ilst,k_ilst))
     &  call errquit('get_ijlst: ma failed',0,MA_ERR)
        if (.not.ma_alloc_get(mt_int,nlst,'ijlst',
     &                        l_jlst,k_jlst))
     &  call errquit('get_ijlst: ma failed',0,MA_ERR)
c ------- create ga_munu_rhoS --- START
c WARNING : ONLY for one atom
        ndir=6 ! Nr. directions: xx,yy,zz,xy,xz,yz
        if (.not. ga_create(mt_dbl,1,nlst*ndir*nlist,
     &     'get_munu_rhoS: g_munu_rhoS', 0,0,g_munu_rhoS)) 
     $     call errquit('get_rhoS: g_AtNr', 0,GA_ERR)
        call ga_zero(g_munu_rhoS)
c ------- create ga_munu_rhoS --- END
c    Define (ilst,jlst) indices to store munu-rhoS elem.
       count=0
       do i=1,nbf
        int_mb(k_ilst+count)=i
        int_mb(k_jlst+count)=i
        count=count+1      
       enddo
       do i=2,nbf
        do j=1,i-1
         int_mb(k_ilst+count)=i
         int_mb(k_jlst+count)=j
         count=count+1     
        enddo 
       enddo
      endif !=============================== efgfile-if-END
c ------------ for NLMO analysis --------- END
c ++++++++++++++++++++++++++++++++++
c--- About content of g_rhoS:
c--- 1st set of nat*6 elements corresponds to ZORA4-EFG
c--- 2nd set of nat*6 elements corresponds to NUM-EFG
         if (.not. ga_create(mt_dbl,1,nlist*6*2,
     &                      'get_rhoS: g_rhoS',
     $                       0,0,g_rhoS)) 
     &       call errquit('get_rhoS: g_rhoS', 0,
     &                     GA_ERR)
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis

      pos_dat=1 ! data-counter in g_rhoS
c     do zora_calc_type =3,4 !  ZORA4-EFG,NUM-EFG
      do zora_calc_type =4,3,-1 !  NUM-EFG,ZORA4-EFG
       do_NonRel=.false.
       if (zora_calc_type.eq.4) do_NonRel=.true.                   
       call get_densZ4(rtdb,ao_bas_han,geom,g_densZ4) 
       so_term=0           ! ZORA-spin-free
       pos=0  ! chunk-index counter for g_munu_rhoS
       do iat1=1,nlist  ! nlist <= nat
        iat=dbl_mb(k_AtNr+iat1-1)
        status=geom_cent_get(geom,iat,at_tag,
     &                       dbl_mb(k_xyzpt+3*(iat-1)),
     &                       dbl_mb(k_zanpt+iat-1))
        xyz_EFGcoords(1)= dbl_mb(k_xyzpt  +3*(iat-1))
        xyz_EFGcoords(2)= dbl_mb(k_xyzpt+1+3*(iat-1))
        xyz_EFGcoords(3)= dbl_mb(k_xyzpt+2+3*(iat-1))
        if (ga_nodeid().eq.0) then
        write(*,19) iat,xyz_EFGcoords(1),xyz_EFGcoords(2),
     &                  xyz_EFGcoords(3)
 19     format('xyz_EFG(',i2,')=(',f15.8,',',f15.8,',',
     &         f15.8,')')
        endif
        do zora_Qpq=1,6 ! xx,yy,zz,xy,xz,yz 
c------Generate munu A^{pq}_r ----- START
          do i=1,ipol
           call ga_zero(g_efgz4(i))
          enddo
          call zora_getv_EFGZ4_SR(rtdb,g_dens_ini, 
     &                            zora_calc_type,      
     &                            zora_Qpq,xyz_EFGcoords, 
     &                            g_efgz4, !  out: munu matrix
     &                            nexc)
         zora_eint=0.0d0
         do ispin=1,ipol
          toac=ga_ddot(g_densZ4(ispin),g_efgz4(ispin)) 
          zora_eint=zora_eint+toac
          if (ga_nodeid().eq.0) then
           write(*,15) zora_calc_type,iat,zora_Qpq,
     &                ispin,pos_dat,toac,zora_eint
 15        format('zora-efg(',i3,',',i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,')')  
          endif
         end do ! ispin-loop
c ++++++ for NLMO analysis ++++++++++++++++++ START
         if (zora_calc_type.eq.3 .and. efgfile.eq.1) then
          call get_munu_rhos_symm(g_efgz4,     ! input
     &                            ipol,
     &                            l_ilst,k_ilst,
     &                            l_jlst,k_jlst,
     &                            nlst,
     &                            g_munu_rhoS, ! output
     &                            pos)
         endif
         pos=pos+1
c ++++++ for NLMO analysis ++++++++++++++++++ END
         call ga_fill_patch(g_rhoS,1,1,pos_dat,pos_dat,
     &                      zora_eint) ! zora_eint --> g_rhoS(pos_dat)
c------Generate munu A^{pq}_r ----- END
         pos_dat=pos_dat+1   
        end do ! zora_Qpq loop
       end do ! iat loop
c ------------- destroy g_densZ4() ------------ START
        if (zora_calc_type.eq.4) then
         do i=1,3
          if (.not. ga_destroy(g_densZ4(i))) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
         enddo
        endif
c ------------- destroy g_densZ4() ------------ END
      end do ! zora_calc_type loop
c ----- Store efgz4 data in a file ------- START
c       Note.- lbl_efgz4 defined in zora.fh
        call util_file_name(lbl_efgz4,.false.,.false.,zorafilename)
        if (.not.dft_zoraEFGZ4_write(
     &           zorafilename,
     &           nlist,
     &           nat,
     &           g_AtNr,
     &           g_rhoS))
     &     call errquit('get_rhoS: dft_zoraNMR_write failed',
     &                  0,DISK_ERR)
c ----- Store efgz4 data in a file ------- END
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
      if (efgfile.eq.1) then ! ============ efgfile-if---START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '--------- g_munu_rhoS ---------- START'
c        call ga_print(g_munu_rhoS)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '--------- g_munu_rhoS ---------- END'
c ---------> Write NMLO analysis data: 2 of 3 ----- START
        call util_file_name(lbl_nlmo,.false.,.false.,zorafilename)
        if (.not.dft_zoraEFGZ4_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &              nlist, ! in: list of selected atoms 
     &                  2, ! in: writing order =1,2,3
     &           ipolmunu, ! in: write for ndata=1
     &  g_zora_scale_munu, ! in: write for ndata=1
     &        g_munu_rhoS, ! in: write for ndata=2
     &        g_densZ4(3), ! in: write for ndata=2
     &           g_munuV6))! in: write for ndata=3
     &     call errquit('get_rhoS: dft_zoraNLMO_write failed',
     &                  0,DISK_ERR)
        if (.not. ga_destroy(g_munu_rhoS)) call errquit(
     &    'dft_zora_utils: ga_destroy failed ',0, GA_ERR)
c ---------> Write NMLO analysis data: 2 of 3 ----- END
        if (.not.ma_free_heap(l_ilst)) call
     &      errquit('dft_zora_rhos: ma_free_heap l_ilst',0, MA_ERR)
        if (.not.ma_free_heap(l_jlst)) call
     &      errquit('dft_zora_rhos: ma_free_heap l_jlst',0, MA_ERR)
       endif! ============================= efgfile-if---END
c ------ destroy remaining from zora_calc_type=3 ------ START
        do i=1,3
         if (.not. ga_destroy(g_densZ4(i))) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
        enddo
c ------ destroy remaining from zora_calc_type=3 ------ END
c----deallocate memory
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('dft_zora_utils, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('dft_zora_utils, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('dft_zora_utils: ma_free_heap l_AtNr',0, MA_ERR)
        do i=1,ipol
         if (.not. ga_destroy(g_efgz4(i))) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
        enddo
         if (.not. ga_destroy(g_rhoS)) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
      return
      end

      subroutine get_gdensZ4_ith(basis,geom,iorb,nbf,ispin,
     &                           nocc,scftyp_dbl)
c     Purpose: Calculate Cmunu_iorb matrix
c              iorb, index   for MO
c              mu,nu indices for AO or Basis Set (BS)
c     Author : Fredy Aquino
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
      integer ispin,iorb,nbf
      integer nocc
      integer g_orb,prpvectors(2)      
      integer k_prpocc,l_prpocc    
      double precision scftyp_dbl    
      integer g_dens_sf
      integer  ga_create_atom_blocked   
      external ga_create_atom_blocked   
      common /Cmunu/g_dens_sf,prpvectors
c --- Cmunu --> g_dens_sf                         
      if (.not.ma_push_get(mt_dbl,nbf*2,'MO occ',l_prpocc,k_prpocc))
     &   call errquit('get_gdensZ4_ith:ma_push_get l_prpocc',
     &                0,MA_ERR)
      g_orb = ga_create_atom_blocked(geom,basis,'orbs')      
      call ga_get(prpvectors(ispin),1,nbf,iorb,iorb,
     &            dbl_mb(k_prpocc),1)
      call ga_zero(g_orb)
      call ga_put(g_orb,1,nbf,iorb,iorb,
     &            dbl_mb(k_prpocc),1)
      call ga_zero(g_dens_sf)
      call ga_dgemm('n','t',nbf,nbf,nocc,scftyp_dbl,
     &              g_orb,g_orb,0.d00,g_dens_sf)
      if (.not.ma_pop_stack(l_prpocc)) call
     &    errquit('get_gdensZ4_ith: ma_pop_stack l_occ',0, MA_ERR)
      if (.not. ga_destroy(g_orb)) call errquit(
     &  'get_gdensZ4_ith: ga_destroy failed ',0, GA_ERR)
      return
      end
c ------ calculate get_gdensZ4_jith -------------------- START
c --> To be used in evaluation of F_{ji}^{1k} NMRZ4
      subroutine get_gdensZ4_jith(basis,geom,
     &                            jorb,iorb,nbf,ispin,
     &                            nocc,scftyp_dbl)
c     Purpose: Calculate Cmunu_jiorb matrix
c              jorb, can be index of virtual  MO
c              iorb, can be index of occupied MO
c              mu,nu indices for AO or Basis Set (BS)
c     Output : g_dens_sf
c     Author : Fredy Aquino
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
      integer ispin,iorb,jorb,nbf
      integer nocc
      integer g_orbj,gorbi,g_dens_sf,prpvectors(2)      
      integer k_prpocc,l_prpocc    
      double precision scftyp_dbl    
      integer  ga_create_atom_blocked   
      external ga_create_atom_blocked   
      common /Cmunu/g_dens_sf,prpvectors
c --- Cmunu --> g_dens_sf                         
      if (.not.ma_push_get(mt_dbl,nbf*2,'MO occ',l_prpocc,k_prpocc))
     &   call errquit('get_gdensZ4_jith:ma_push_get l_prpocc',
     &                0,MA_ERR)
      g_orbj = ga_create_atom_blocked(geom,basis,'orbs')   
      g_orbi = ga_create_atom_blocked(geom,basis,'orbs')         
c ---- Getting Cj (g_orbj) MO coeffs
      call ga_get(prpvectors(ispin),1,nbf,jorb,jorb,
     &            dbl_mb(k_prpocc),1)
      call ga_zero(g_orbj)
      call ga_put(g_orbj,1,nbf,jorb,jorb,
     &            dbl_mb(k_prpocc),1)
c ---- Getting Ci (g_orbi) MO coeffs
      call ga_get(prpvectors(ispin),1,nbf,iorb,iorb,
     &            dbl_mb(k_prpocc),1)
      call ga_zero(g_orbi)
      call ga_put(g_orbi,1,nbf,iorb,iorb,
     &            dbl_mb(k_prpocc),1)
      call ga_zero(g_dens_sf)
      call ga_dgemm('n','t',nbf,nbf,nocc,scftyp_dbl,
     &              g_orbj,g_orbi,0.d00,g_dens_sf)
      if (.not.ma_pop_stack(l_prpocc)) call
     &    errquit('get_gdensZ4_jith: ma_pop_stack l_occ',0, MA_ERR)
      if (.not. ga_destroy(g_orbi)) call errquit(
     &  'get_gdensZ4_jith: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_orbj)) call errquit(
     &  'get_gdensZ4_jith: ga_destroy failed ',0, GA_ERR)
      return
      end
c ------ calculate get_gdensZ4_jith -------------------- END

      subroutine hnd_prp_get_vecs(rtdb,geom,basis,
     2                            prpvectors,
     2                            scftyp,
     2                            nclosed,nopen,nvirt)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "util.fh"
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      integer     rtdb          ! [input] database handle
      integer     geom          ! [input] geometry handle
      integer     basis         ! [input] handles to basis
      integer     ndens         ! [output] number of active density handles (RHF=1, UHF=3)
      character*3 scftyp        ! [output] type of wave function
      integer nclosed(2),nopen(2),nvirt(2) ! [output] occupation info
c
      integer nbf, nmo, k_prpocc, l_prpocc, k_prpeval, l_prpeval
      integer prpvectors(2)
c
c     Get vectors and other information
c     Arrays occ(nbf*2) and evals(nbf*2) are needed
c
      if (.not. bas_numbf(basis,nbf)) call 
     &    errquit('hnd_prp_get_dens: could not get nbf',0, BASIS_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*2,'MO eval',l_prpeval,k_prpeval)) 
     &   call errquit('hnd_prp_get_dens:ma_push_get l_prpeval',0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*2,'MO occ',l_prpocc,k_prpocc)) 
     &   call errquit('hnd_prp_get_dens:ma_push_get l_prpocc',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,nvirt,
     &                      scftyp,prpvectors,dbl_mb(k_prpocc),
     &                      dbl_mb(k_prpeval),nmo)
c
c     Make the density matrix
c
c     Cleanup of MA arrays that are not needed
c
      if (.not.ma_pop_stack(l_prpocc)) call
     &    errquit('hnd_prp_get_dens: ma_pop_stack l_occ',0, MA_ERR)
      if (.not.ma_pop_stack(l_prpeval)) call
     &    errquit('hnd_prp_get_dens: ma_pop_stack l_eval',0, MA_ERR)
c
      return
      end
      
      subroutine get_densZ4(rtdb,basis,geom,g_densZ4)
c
c     Purpose: Calculate Cmunu Z4 density matrix
c              mu,nu indices for AO or Basis Set (BS)
c     Author : Fredy Aquino

       implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "zora.fh" 
      character*3 scftyp  ! output from hnd_prp_get_vecs()
      integer rtdb
      integer basis       ! [Input] Basis set
      integer geom        ! [Input] Geometry
      integer ispin,ipol,iorb,porb,nbf
      integer nocc(2)
      integer prpvectors(2)      
      integer g_dens_sf,g_densZ4(3)
      integer l_Ci,k_Ci
      integer noc1 ! Total # of occupied MOs
      integer nclosed(2),nopen(2),nvirt(2) ! [output] occupation info
      double precision scftyp_dbl,scale 
      integer  ga_create_atom_blocked   
      external ga_create_atom_blocked   
      integer i,ii

      common /Cmunu/g_dens_sf,prpvectors

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('get_densZ4: could not get nbf',0, BASIS_ERR)      
      g_dens_sf=ga_create_atom_blocked(geom,basis,'sf-1')
      do ii=1,3
       g_densZ4(ii)=ga_create_atom_blocked(geom,basis,'densZ4-1')
      end do
      call hnd_prp_get_vecs(rtdb,geom,basis,prpvectors,scftyp,
     &                      nclosed,nopen,nvirt)
      if (scftyp.eq.'RHF') then
         ipol=1
         scftyp_dbl=2.0d00
         nocc(1)=nclosed(1)+nopen(1) ! is that right??
         noc1=nocc(1)
      else if (scftyp.eq.'UHF') then
         ipol=2
         scftyp_dbl=1.0d00
         nocc(1)=nopen(1)
         nocc(2)=nopen(2)    
         noc1=nocc(1)+nocc(2)
      endif
      if (.not.ma_alloc_get(mt_dbl,ipol*noc1,'Ci',l_Ci,k_Ci)) 
     &    call errquit('hnd_efgmap_Z4: ma failed',0,MA_ERR)
      call ga_get(g_Ci,1,ipol,1,noc1,dbl_mb(k_Ci),ipol)
      do ispin=1,ipol
      porb=ispin
       call ga_zero(g_densZ4(ispin))
       do iorb=1,nocc(ispin)
         call get_gdensZ4_ith(basis,geom,iorb,nbf,ispin,
     &                        nocc(ispin),scftyp_dbl)
         if ((do_NonRel) .or. (not_zora_scale)) then
             scale=1.0d0
         else
             scale=1.0d0/dbl_mb(k_Ci+porb-1)
         endif
         call ga_scale(g_dens_sf,scale)
         call ga_add(1.0d00,g_densZ4(ispin),
     &               1.0d00,g_dens_sf,g_densZ4(ispin))
         porb=porb+ipol
       end do     ! iorb 
      end do      ! ispin
      call ga_zero(g_densZ4(3))
      call ga_copy(g_densZ4(1),g_densZ4(3))
      if (ipol.gt.1) then 
       call ga_add(1.0d00,g_densZ4(1),
     &             1.0d00,g_densZ4(2),
     &                    g_densZ4(3))      
      end if
      if (.not. ga_destroy(g_dens_sf)) call errquit(
     &  'get_densZ4: ga_destroy failed ',0, GA_ERR)
      if (.not.ga_destroy(prpvectors(1))) call
     &    errquit('get_densZ4: ga_destroy vecs 1',0, GA_ERR)
      if (scftyp.eq.'UHF') then
         if (.not.ga_destroy(prpvectors(2))) call
     &       errquit('get_densZ4: ga_destroy vecs 2',0, GA_ERR)
      endif
      if (.not. MA_free_heap(l_Ci))
     &  call errquit('get_densZ4:cannot free heap',111, MA_ERR)
      return
      end

      integer function read_SLCTD_EFG_Atoms
     &            (rtdb,nat,nlist,g_AtNr)
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
      integer efgz4atoms
      double precision AtNr_dbl

       read_SLCTD_EFG_Atoms = 0
        if (.not. rtdb_get(rtdb, 'efgz4:natoms',mt_int,
     $                     1,efgz4atoms))
     &       efgz4atoms=0 ! reset
         if (efgz4atoms.eq.0) then
          efgz4atoms=nat
           nlist=efgz4atoms
           do ii=1,efgz4atoms
            AtNr_dbl=ii
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo     
         else
          if (.not. rtdb_get(rtdb, 'efgz4:atom list',mt_int,
     $                     efgz4atoms,atomnr))
     $      call errquit('prop_input-EFGZ4: rtdb_get failed',
     $                   555, RTDB_ERR)
           nlist=efgz4atoms
           do ii=1,efgz4atoms
            AtNr_dbl=atomnr(ii)
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo
         endif
       read_SLCTD_EFG_Atoms = 1
       return
       end    
       subroutine get_munu_rhos_symm(g_munu,        ! in  : munu matrix
     &                               ipol,          ! in  : nr. polarizations
     &                               l_ilst,k_ilst, ! in  : array of indices ilst      
     &                               l_jlst,k_jlst, ! in  : array of indices jlst                       
     &                               nlst,          ! in  : = nbf*(nbf+1)/2
     &                               g_munu_rhoS,   ! out : accumulate main-diag,off-diag munu elements
     &                               pos)
c    g_munu_rhoS accumulates xx,yy,zz,xy,xz,yz
c    matrices in the following format:
c    11 22 ... (nbf nbf)  -> main diagonal first
c    21 31 32 41 42 43    ->  off-diagonal later
c    nbf, number of basis functions
c    nlst=nbf*(nbf*+1)/2
c    1st chunk of nlst numbers corresponds to xx
c    2nd chunk --> yy,  3rd chunk --> zz,
c    4th chunk --> xy,  5th chunk --> xz,  6th chunk --> yz
c    Input: pos=0,1,..., chunk-index
c           --> the size of one chunk is nlst
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 

       integer g_munu_rhoS
       integer g_munu(2)
       integer i,ipol
       integer nlst,pos
       integer l_ilst,k_ilst,     
     &         l_jlst,k_jlst
       integer jlo,jhi
       double precision v(nlst)
       call ga_gather(g_munu(1),v,              ! g_munu ---> v
     &                int_mb(k_ilst),int_mb(k_jlst),
     &                nlst)
c      Now accumulates ith-chunk on g_munu_rhoS
       jlo=pos*nlst+1
       jhi=pos*nlst+nlst
       call ga_put(g_munu_rhoS,1,1,jlo,jhi,v,1) ! v --> g_munu_rhoS
       return
       end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++ READ/WRITE EFGZ4 data +++++++++++++ START
c Note.- Using modified versions of
c        dft_zora_read() and dft_zoraNMR_write()
c        --> located in dft_zora_utils.F
czora...Write out the zora NMR shieldings to disk

      logical function dft_zoraEFGZ4_write(
     &           filename, ! in: filename
     &              nlist, ! in: number of selected atoms
     &                nat, ! in: total number of atoms
     &            g_AtNr1, ! in: list of atoms to calc. shieldings
     &             g_rhoS) ! in: rhoS values
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
      integer nlist,   !       nr. slc atoms
     &        nat,     ! total nr. of  atoms
     &        ntot_efg ! = ntat*6*2 6=xx,yy,zz,xy,xz,yz 2=NonRel, EFGZ4
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_rhoS   ! rhoS values
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_rhoS,k_rhoS
      integer ok, iset, i, j
      integer inntsize
      integer alo(3),ahi(3),ld(2)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c     Write out the number of sets and basis functions
       write(unitno, err=1001) nlist
       write(unitno, err=1001) nat
c     Allocate the temporary buffer
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraEFGZ4_write',
     &                        l_AtNr,k_AtNr))
     $  call errquit('dft_zoraEFGZ4_write: ma failed', 
     &               nlist, MA_ERR)
       ntot_efg=nlist*6*2
       if (.not. ma_alloc_get(mt_dbl,ntot_efg,
     &                        'dft_zoraEFGZ4_write',
     &                        l_rhoS,k_rhoS))
     $  call errquit('dft_zoraEFGZ4_write: ma failed', 
     &               ntot_efg, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
       call ga_get(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
       call swrite(unitno,dbl_mb(k_AtNr),nlist)
       call ga_get(g_rhoS,1,1,1,ntot_efg,dbl_mb(k_rhoS),1)
       call swrite(unitno,dbl_mb(k_rhoS),ntot_efg)
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraEFGZ4_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_rhoS))
     $  call errquit('dft_zoraEFGZ4_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraEFGZ4_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA EFGZ4 data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraEFGZ4_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraEFGZ4_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraEFGZ4_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraEFGZ4_read(
     &           filename, ! in : filename
     &              nlist, ! in : number of selected atoms
     &                nat, ! in : total number of atoms
     &            g_AtNr1, ! out: list of atoms to calc. shieldings
     &             g_rhoS) ! out: rhoS values
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
      integer nlist,   !       nr. slc atoms
     &        nat,     ! total nr. of  atoms
     &        ntot_efg ! = ntat*6*2 6=xx,yy,zz,xy,xz,yz 2=NonRel, EFGZ4
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_rhoS
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_rhoS,k_rhoS
      integer ok,inntsize
      integer nat_read,nlist_read
c     Initialise to invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c---- Create GA arrays: g_AtNr1,g_rhoS --- START
         if (.not. ga_create(mt_dbl,1,nlist,
     &   'dft_zoraEFGZ4_read: g_AtNr1',0,0,g_AtNr1)) 
     $   call errquit('dft_zoraEFGZ4_read: g_AtNr1',0,GA_ERR)
        call ga_zero(g_AtNr1) 
        ntot_efg=nlist*6*2
         if (.not. ga_create(mt_dbl,1,ntot_efg,
     &   'dft_zoraEFGZ4_read: g_rhoS',0,0,g_rhoS)) 
     $   call errquit('dft_zoraEFGZ4_read: g_rhoS',0,GA_ERR)
        call ga_zero(g_rhoS)   
c---- Create GA arrays: g_AtNr1,g_rhoS --- END
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA EFGZ4 data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nlist_read 
       read(unitno, err=1001, end=1001) nat_read      
c      Error checks
       if ((nat_read   .ne. nat) .or.
     &     (nlist_read .ne. nlist) ) goto 1003
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraEFGZ4_read',
     &                     l_AtNr,k_AtNr))
     $  call errquit('dft_zoraEFGZ4_read: ma failed', 
     &               nlist, MA_ERR)
       if (.not. ma_alloc_get(mt_dbl,ntot_efg,'dft_zoraEFGZ4_read',
     &                     l_rhoS,k_rhoS))
     $  call errquit('dft_zoraEFGZ4_read: ma failed', 
     &               ntot_efg, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
       call sread(unitno,dbl_mb(k_AtNr),nlist)
       call ga_put(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
       call sread(unitno,dbl_mb(k_rhoS),ntot_efg)
       call ga_put(g_rhoS,1,1,1,ntot_efg,dbl_mb(k_rhoS),1)
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraEFGZ4_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_rhoS))
     $  call errquit('dft_zoraEFGZ4_read: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraEFGZ4_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraEFGZ4_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraEFGZ4_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 'dft_zoraEFGZ4_read: file inconsistent',
     &           ' with calculation',
     $           filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraEFGZ4_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraEFGZ4_NLMOAnalysis_write(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &              nlist, ! in: list of selected atoms 
     &              ndata, ! in: writing order =1,2,3
     &           ipolmunu, ! in: write for ndata=1
     &  g_zora_scale_munu, ! in: write for ndata=1
     &        g_munu_rhoS, ! in: write for ndata=2
     &             g_dens, ! in: write for ndata=2
     &           g_munuV6) ! in: write for ndata=3
c Description: Collecting three matrices to be used
c              in wefgfile(rtdb) and wnbofile(rtdb)
c              Those routines are called in prop.F
c              after hnd_property(rtdb)
c              The info collected is:
c              1. ipolmunu, nr. of polarizations 
c              2. g_zora_scale_munu, = g_zora_scale_sf
c                 --> Stored in dft_zora_scale() [dft_zora_utils.F]
c              3. g_munu_rhoS, small component density munu matrix
c                              about the format check info in header
c                              of get_munu_rhos_symm()
c                 --> Stored in get_rhoS() [dft_zora_rhos.F]
c              4. g_munuV6, electronic EFG munu matrix (without the
c                           small component density correction)
c                 --> Stored in get_munuV6() [hnd_elfcon_symm.F] 
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
     &        nlst,   ! = nbf*(nbf+1)/2
     &        nmat,   ! = nlst*ndir*nlist
     &        nmat1   ! = nbf*nbf
      integer ndata   ! = 1,2,3
      integer ipolmunu,g_zora_scale_munu(2) ! write for ndata=1
      integer g_munu_rhoS,g_dens            ! write for ndata=2
      integer g_munuV6                      ! write for ndata=3
      integer ndir,nbf
      integer unitno
      parameter (unitno = 77)
      integer l_mat,k_mat,l_mat1,k_mat1
      integer ok, iset, i, j
      integer inntsize
c WARNING: The sequence in which is called this routine is: 1->2->3.
c ---- check ndata=1,2,3 ----- START
      if (ndata.ne.1 .and. 
     &    ndata.ne.2 .and. ndata.ne.3) then
       write(*,*) 'Error in dft_zoraEFGZ4_NLMOAnalysis_write:',
     &            'ndata=',ndata,' It should be 1, 2 or 3'
       stop
      endif
c ---- check ndata=1,2,3 ----- END
      nlst=nbf*(nbf+1)/2
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Allocate the temporary buffer
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if      (ndata.eq.1) then
        nmat=nbf
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_write',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_write: ma failed', 
     &               nmat, MA_ERR)
       else if (ndata.eq.2) then
        nmat=nlst*ndir*nlist
        nmat1=nbf*nbf
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_write',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_write: ma failed', 
     &               nmat, MA_ERR)
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat1,'dft_zoraNLMO_write',
     &        l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_write: ma failed', 
     &               nmat1, MA_ERR)
       else if (ndata.eq.3) then
        nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_write',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_write: ma failed', 
     &               nmat, MA_ERR)
       endif
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
       if      (ndata.eq.1) then
c     Open the file - 1st time
        open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000)
c     Write out the number of sets and basis functions
        write(unitno, err=1001) nbf
        write(unitno, err=1001) nlst
        write(unitno, err=1001) ndir
        write(unitno, err=1001) ipolmunu
c     Write out g_zora_scale_sf
        do iset = 1, ipolmunu
         do i = 1, nbf
          call ga_get(g_zora_scale_munu(iset), 1, nbf, i, i, 
     &               dbl_mb(k_mat),1)
          call swrite(unitno, dbl_mb(k_mat), nbf)
         end do
        end do
       else if (ndata.eq.2) then
c     Open the file - 2nd time
        open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000,position='append')
        write(unitno, err=1001) nlist
        call ga_get(g_munu_rhoS,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_dens,1,nbf,1,nbf,
     &              dbl_mb(k_mat1),nbf)
        call swrite(unitno,dbl_mb(k_mat1),nmat1)
       else if (ndata.eq.3) then
c     Open the file - 3rd time
        open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000,position='append')
        call ga_get(g_munuV6,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
       endif
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if      (ndata.eq.1) then
      if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_write: ma free_heap failed', 
     &               911, MA_ERR)
      else if (ndata.eq.2) then
      if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat1))
     $  call errquit('dft_zoraNLMO_write: ma free_heap1 failed', 
     &               911, MA_ERR)
      else if (ndata.eq.3) then
      if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_write: ma free_heap failed', 
     &               911, MA_ERR)
      endif
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraEFGZ4_NLMOAnalysis_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA NLMO EFGZ4 data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraNLMO_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNLMO_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNLMO_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraEFGZ4_NLMOAnalysis_read(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &              nlist, ! in: list of selected atoms 
     &           ipolmunu, ! in: 
     &  g_zora_scale_munu, ! out: 
     &        g_munu_rhoS, ! out: 
     &             g_dens, ! out:
     &           g_munuV6) ! out:  
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
      integer nmat,nmat1
      integer g_zora_scale_munu(2),
     &        g_munu_rhoS,
     &        g_dens,
     &        g_munuV6
      integer nbf,nbf_read,nlst,nlst_read,
     &        ndir,ndir_read,nlist,nlist_read,
     &        ipolmunu,ipolmunu_read
      integer unitno
      parameter (unitno = 77)
      integer l_mat,k_mat,
     &        l_mat1,k_mat1
      integer ok, iset, i, j
      integer inntsize
c     Initialise to invalid MA handle
      nlst=nbf*(nbf+1)/2
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c---- Create GA arrays: g_AtNr1,g_rhoS --- START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_zora_scale_munu',
     &    0,0,g_zora_scale_munu(1))) 
     $   call errquit('dft_zoraNLMO_read: g_rhoS',0,GA_ERR)
        call ga_zero(g_zora_scale_munu(1))
       if (ipolmunu.gt.1) then
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_zora_scale_munu',
     &    0,0,g_zora_scale_munu(2))) 
     $   call errquit('dft_zoraNLMO_read: g_rhoS',0,GA_ERR)
        call ga_zero(g_zora_scale_munu(2))
       endif
       nmat=nlst*ndir*nlist
         if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_rhoS',0,0,g_munu_rhoS)) 
     $   call errquit('dft_zoraNLMO_read: g_rhoS',0,GA_ERR)
        call ga_zero(g_munu_rhoS) 
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_sdens',0,0,g_dens)) 
     $   call errquit('dft_zoraNLMO_read: g_dens',0,GA_ERR)
        call ga_zero(g_dens)    
         if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munuV6',0,0,g_munuV6)) 
     $   call errquit('dft_zoraNLMO_read: g_munuV6',0,GA_ERR)
        call ga_zero(g_munuV6)   
c---- Create GA arrays: g_AtNr1,g_rhoS --- END
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA EFGZ4 data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nlst_read
       read(unitno, err=1001, end=1001) ndir_read
       read(unitno, err=1001, end=1001) ipolmunu_read
c      Error checks
       if ((nbf_read      .ne. nbf)   .or.
     &     (nlst_read     .ne. nlst)  .or.
     &     (ndir_read     .ne. ndir)  .or.
     &     (ipolmunu_read .ne. ipolmunu)) goto 1003
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
c ------- Reading 1st matrix: g_zora_scale_munu ----- START
       nmat=nbf
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       do iset = 1, ipolmunu
        do i = 1, nbf
         call sread(unitno, dbl_mb(k_mat), nbf)
         call ga_put(g_zora_scale_munu(iset), 1, nbf, i, i, 
     &               dbl_mb(k_mat), 1)
        end do
       end do
      if (.not. ma_free_heap(l_mat))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c ------- Reading 1st matrix: g_zora_scale_munu ----- END
c ------- Reading 2nd matrix: g_munu_rhoS ----- START
       read(unitno, err=1001, end=1001) nlist_read  
       if (nlist_read    .ne. nlist) goto 1003
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(mt_dbl,nmat1, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat1, MA_ERR)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munu_rhoS,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat1),nmat1)
       call ga_put(g_dens,1,nbf,1,nbf,dbl_mb(k_mat1),nbf)
      if (.not. ma_free_heap(l_mat))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat1))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c ------- Reading 2nd matrix: g_munu_rhoS ----- END
c ------- Reading 3rd matrix: g_munuV6 ----- START
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munuV6,1,1,1,nmat,dbl_mb(k_mat),1)
      if (.not. ma_free_heap(l_mat))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c ------- Reading 3rd matrix: g_munuV6 ----- END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraEFGZ4_NLMOAnalysis_read = ok .eq. 1
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
c +++++++++++ READ/WRITE EFGZ4 data +++++++++++++ END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine print_EFGZ4_version()
c
      implicit none
c
#include "stdio.fh"
c
      write(LuOut,*)
      call util_print_centered(LuOut,
     $   'EPR Z4', 23, .true.)
      write(LuOut,*)
c
      return
      end
c $Id$
