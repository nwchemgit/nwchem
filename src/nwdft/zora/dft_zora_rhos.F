c=======================================================
c================ Fredy Aquino's routines ======== START
c=======================================================
      subroutine get_rhoS(
     &                    rtdb,g_dens_ini,nexc,
     &                    geom, 
     &                    ao_bas_han,
     &                    nbf,nbf_ao,
     &                    g_movecs,
     &                    g_zora_scale_sf,
     &                    noc,
     &                    ipol)
c -- Purpose: Calculation of small component density, \rho_S
c             for Electric Field Gradient calculation
c             Source: van Lenthe, et.al.,JCP, V112,N19,Y2000
c -- Author : Fredy Aquino 12-07-09
c
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 

      integer rtdb     
      integer g_orb,g_dens_sf
      integer g_zora_scr(2)        
      integer g_dens_ini(2)   
      integer g_zora_scale_sf(2)
      integer g_movecs(2)
      integer noc(2)
      integer ipol
      integer geom
      integer ao_bas_han                      
      integer ispin,pos_dat,
     &        do_zgc_old,nexc
      integer iat,nat,noc1   
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt    
      logical status               
      character*16 element, at_tag  
      double precision scf_dbl,zora_eint,toac 
c ++++++++++++++++++++++++++++
      character*100 fname_RSLC
      integer stat_read
      integer iat1,l_AtNr,k_AtNr    
c ++++++++++++++++++++++++++++
      integer g_densZ4(3)
      integer read_SLCTD_EFG_Atoms
      external get_densZ4,read_SLCTD_EFG_Atoms 
      integer nbf,nbf_ao
c
c---- Input global variables (defined in zora.fh):
c     1. zora_calc_type  ! =3 -> ZORA4-EFG
c     2. so_term         ! =0 -> ZORA4-spin-free
c     3. xyz_EFGcoords(i) i=1,2,3
c     4. zora_Qpq=1,6 : xx,yy,zz,xy,xz,yz 

c================================================   
c     if (ga_nodeid().eq.0) write(*,17) zoracutoff_EFG
c17   format('zoracutoff_EFG=',f15.8)
c================================================
      if (ipol.eq.1) then
        scf_dbl=2.0d00
        noc1=noc(1)
      else if (ipol.eq.2) then 
        scf_dbl=1.0d00   
        noc1=noc(1)+noc(2)    
      endif
c     if (ga_nodeid().eq.0)
c    &  write(*,18) ipol,scf_dbl,noc1
c18   format('(ipol,scf_dbl,noc1)=(',i3,',',
c    &       f15.8,',',i4,')')
      status=geom_ncent(geom,nat)
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
c----- Allocate global arrays - FA
       if(.not.ga_duplicate(g_zora_scale_sf(1),g_zora_scr(1),'sf 1'))
     &  call errquit('get_rhoS: ga_duplicate failed',1,GA_ERR)
        if(ipol.gt.1) then
         if(.not.ga_duplicate(g_zora_scale_sf(2),g_zora_scr(2),'sf 2'))
     &   call errquit('get_rhoS: ga_duplicate failed',1,GA_ERR)
        endif
c ++++++++++++++++++++++++++++++++++
c +++++ Read Atom Nr for EFG calc ++
         if (.not. ga_create(mt_dbl,1,nat,
     &                       'get_rhoS: g_AtNr',
     $                       0,0,g_AtNr)) 
     $       call errquit('get_rhoS: g_AtNr', 0,
     &                    GA_ERR)
        call ga_zero(g_AtNr)
       stat_read=read_SLCTD_EFG_Atoms
     &            (rtdb,nat,nlist,g_AtNr)
        write(*,*) "==> nlist=",nlist
        call ga_print(g_AtNr)
c  Allocate memory for l_AtNr,k_AtNr 
      if (.not.ma_alloc_get(mt_dbl,nat,'AtNr',
     &    l_AtNr,k_AtNr))
     &    call errquit('get_rhoS: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
c ++++++++++++++++++++++++++++++++++
c--- About content of g_rhoS:
c--- 1st set of nat*6 elements corresponds to ZORA4-EFG
c--- 2nd set of nat*6 elements corresponds to NUM-EFG
         if (.not. ga_create(mt_dbl,1,nat*6*2,
     &                      'get_rhoS: g_rhoS',
     $                       0,0,g_rhoS)) 
     &       call errquit('get_rhoS: g_rhoS', 0,
     &                     GA_ERR)
      call ga_zero(g_rhoS)
      pos_dat=1 ! data-counter in g_rhoS
c     do zora_calc_type =3,4 !  ZORA4-EFG,NUM-EFG
      do zora_calc_type =4,3,-1 !  NUM-EFG,ZORA4-EFG
      call get_densZ4(rtdb,ao_bas_han,geom,g_densZ4)
      so_term=0           ! ZORA-spin-free
      write(*,*) "zora_calc_type=",zora_calc_type
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
        do zora_Qpq=1,6 ! xx,yy,zz,xy,xz,yz - FA
c------Generate munu A^{pq}_r ----- START
        call ga_zero(g_zora_scr(1))
        if (ipol.gt.1) call ga_zero(g_zora_scr(2))
        call zora_getv(rtdb,g_dens_ini,g_zora_scr,nexc)
        zora_eint=0.0d0
         do ispin=1,ipol
c         write(*,1) iat,zora_Qpq,ispin,pos_dat
c1        format('--- g_densZ4(',i3,',',i3,',',
c    &           i3,',',i3,') ----- START')
c =======================================
          call ga_sync() ! Syncronizing 
c =======================================
c         call ga_print(g_densZ4(ispin))
c         call ga_print(g_zora_scr(ispin))
c         write(*,2) iat,zora_Qpq,ispin,pos_dat
c2        format('--- g_densZ4(',i3,',',i3,',',
c    &            i3,',',i3,') ----- END')
          toac=ga_ddot(g_densZ4(ispin),g_zora_scr(ispin)) 
          zora_eint=zora_eint+toac
          if (ga_nodeid().eq.0) then
           write(*,15) zora_calc_type,iat,zora_Qpq,
     &                ispin,pos_dat,toac,zora_eint
 15        format('zora-efg(',i3,',',i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,')')  
          endif
         end do ! ispin-loop
        call ga_fill_patch(g_rhoS,1,1,pos_dat,pos_dat,zora_eint) 
c------Generate munu A^{pq}_r ----- END
        pos_dat=pos_dat+1
        end do ! zora_Qpq loop
       end do ! iat loop
      end do ! zora_calc_type loop

c----deallocate memory - FA
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('dft_zora_utils, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('dft_zora_utils, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('dft_zora_utils: ma_free_heap l_AtNr',0, MA_ERR)
c ---------------- added by FA ---- end
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
c
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
c
       implicit none
c
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
      g_dens_sf=ga_create_atom_blocked(geom,basis,'orbs dens')
      do ii=1,3
       g_densZ4(ii)=ga_create_atom_blocked(geom,basis,'orbs dens')
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
c      do i=1,ipol*noc1
c       write(*,18) i,dbl_mb(k_Ci+i-1)
c  18   format('Ci(',i3,')=',f15.8)
c      end do

      do ispin=1,ipol
      porb=ispin
       call ga_zero(g_densZ4(ispin))
       do iorb=1,nocc(ispin)
         call get_gdensZ4_ith(basis,geom,iorb,nbf,ispin,
     &                        nocc(ispin),scftyp_dbl)
c         write(*,117) ispin,iorb,porb,dbl_mb(k_Ci+porb-1)
c 117     format('Ci(',i3,',',i3,',',i3,')=',f15.8)

         if (zora_calc_type.eq.3) scale=1.0d0/dbl_mb(k_Ci+porb-1)
         if (zora_calc_type.eq.4) scale=1.0d0  ! NUM-EFG

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
#include "rtdb.fh" 

#include "context.fh"

      integer rtdb,ii,nlist,nat,g_AtNr
      integer atomnr(nat)
      integer efgz4atoms
      double precision AtNr_dbl

       read_SLCTD_EFG_Atoms = 0
       efgz4atoms=0 ! reset
        if (.not. rtdb_get(rtdb, 'prop:atnr',mt_int,
     $                     1,efgz4atoms))
     $      call errquit('prop_input-EFGZ4: rtdb_get failed',
     $                   555, RTDB_ERR)
c          write(*,*) "in read-slc: efgz4atoms=",efgz4atoms
         if (efgz4atoms.eq.0) then
          efgz4atoms=nat
c          if (ga_nodeid().eq.0) 
c     $     write(*,*) "efgz4atoms=",efgz4atoms
           nlist=efgz4atoms
           do ii=1,efgz4atoms
            AtNr_dbl=ii
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
c            if (ga_nodeid().eq.0) 
c     $        write(*,482) ii,atomnr(ii)
c 482         format('atomnr(',i3,')=',i5)
           enddo     
         else
          if (.not. rtdb_get(rtdb, 'prop:atnr1',mt_int,
     $                     efgz4atoms,atomnr))
     $      call errquit('prop_input-EFGZ4: rtdb_get failed',
     $                   555, RTDB_ERR)
c          if (ga_nodeid().eq.0) 
c     $     write(*,*) "RHOS:efgz4atoms=",efgz4atoms
           nlist=efgz4atoms
           do ii=1,efgz4atoms
            AtNr_dbl=atomnr(ii)
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
c            if (ga_nodeid().eq.0) 
c     $        write(*,481) ii,atomnr(ii)
c 481         format('atomnr(',i3,')=',i5)
           enddo
         endif
       read_SLCTD_EFG_Atoms = 1
       return
       end    

       subroutine print_EFGZ4_version()
       implicit none
       write(*,*) "+++++++++++++++++++++++++++++++++"
       write(*,*) "+ EFG-Z4 Version 04-21-10 12:00 +"
       write(*,*) "+ Author: Fredy W. Aquino       +"
       write(*,*) "+++++++++++++++++++++++++++++++++"
       return
       end
c=======================================================
c================ Fredy Aquino's routines ======== END
c=======================================================
