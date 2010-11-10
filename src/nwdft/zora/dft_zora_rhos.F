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
      integer g_efgz4(2)        
      integer g_dens_ini(2)   
      integer noc(2)
      integer ipol
      integer geom
      integer ao_bas_han                      
      integer ispin,pos_dat,
     &        do_zgc_old,nexc
      integer iat,nat,noc1   
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt    
      integer l_ilst,k_ilst,     
     &        l_jlst,k_jlst,    
     &        nlst,i,j,count
      integer g_munu_rhoS,g_munuV6
c     g_munu_rhoS, contains xx,yy,zz,xy,xz,yz
c                  munu half-triangle matrices
c                  including main diagonal
      logical status               
      character*16 element, at_tag  
      double precision scf_dbl,zora_eint,toac 
c ++++++++++++++++++++++++++++
      integer stat_read
      integer iat1,l_AtNr,k_AtNr    
c ++++++++++++++++++++++++++++
      integer g_densZ4(3)
      integer read_SLCTD_EFG_Atoms
      external get_densZ4,read_SLCTD_EFG_Atoms,
     &         zora_getv_EFGZ4_SR 
      integer nbf,nbf_ao,pos,ndir
      integer g_zora_scale_munu(2)
      integer ipolmunu
      integer zora_Qpq
      double precision xyz_EFGcoords(3)   
      common/munu_EFG/g_munuV6,g_munu_rhoS,
     &                g_zora_scale_munu,ipolmunu ! munu-matrices
c     Note.- g_munu_rhoS is created in dft_zora_rhos.F
c            g_munuV6    is created in hnd_elfcon_symm.F
c            g_munuV6 is not created/defined here yet
      
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
      status=geom_ncent(geom,nat)
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('get_rhoS: ma failed',911,MA_ERR)
c----- Allocate global arrays - FA
        do i=1,ipol
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'get_rhoS: g_efgz4',
     $                       0,0,g_efgz4(i)))
     $       call errquit('get_rhoS: g_efgz4', 0,
     &                    GA_ERR)
             call ga_zero(g_efgz4(i))
        enddo
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
c  Allocat memory for (l_ilst,k_ilst) (l_jlst,k_jlst)
      nlst=nbf*(nbf+1)/2  ! single triangle matrix
      if (.not.ma_alloc_get(mt_int,nlst,'ijlst',
     &                      l_ilst,k_ilst))
     &    call errquit('get_ijlst: ma failed',0,MA_ERR)
      if (.not.ma_alloc_get(mt_int,nlst,'ijlst',
     &                      l_jlst,k_jlst))
     &    call errquit('get_ijlst: ma failed',0,MA_ERR)
c ------- create ga_munu_rhoS --- START
c WARNING : ONLY for one atom
         ndir=6 ! Nr. directions: xx,yy,zz,xy,xz,yz
         if (.not. ga_create(mt_dbl,1,nlst*ndir*nlist,
     &                       'get_munu_rhoS: g_munu_rhoS',
     $                       0,0,g_munu_rhoS)) 
     $       call errquit('get_rhoS: g_AtNr', 0,
     &                    GA_ERR)
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
       do_NonRel=.false.
       if (zora_calc_type.eq.4) do_NonRel=.true.                  
       call get_densZ4(rtdb,ao_bas_han,geom,g_densZ4)
       so_term=0           ! ZORA-spin-free
       write(*,*) "zora_calc_type=",zora_calc_type
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
        do zora_Qpq=1,6 ! xx,yy,zz,xy,xz,yz - FA
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
c ++++++++++++++++++++
         if (zora_calc_type.eq.3) then
          call get_munu_rhos_symm(g_efgz4,ipol,
     &                            l_ilst,k_ilst,
     &                            l_jlst,k_jlst,
     &                            nlst,
     &                            g_munu_rhoS,pos)
         endif
         pos=pos+1
c ++++++++++++++++++++
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
      if (.not.ma_free_heap(l_ilst)) call
     &    errquit('dft_zora_rhos: ma_free_heap l_ilst',0, MA_ERR)
      if (.not.ma_free_heap(l_jlst)) call
     &    errquit('dft_zora_rhos: ma_free_heap l_jlst',0, MA_ERR)
c ---------------- added by FA ---- end
        do i=1,ipol
         if (.not. ga_destroy(g_efgz4(i))) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
        enddo
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
      do ispin=1,ipol
      porb=ispin
       call ga_zero(g_densZ4(ispin))
       do iorb=1,nocc(ispin)
         call get_gdensZ4_ith(basis,geom,iorb,nbf,ispin,
     &                        nocc(ispin),scftyp_dbl)
         if (do_NonRel.eqv. .true. .or. 
     &       not_zora_scale .eqv. .true.) then
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
         if (efgz4atoms.eq.0) then
          efgz4atoms=nat
           nlist=efgz4atoms
           do ii=1,efgz4atoms
            AtNr_dbl=ii
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo     
         else
          if (.not. rtdb_get(rtdb, 'prop:atnr1',mt_int,
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

       subroutine print_EFGZ4_version()
       implicit none
       write(*,*) "+++++++++++++++++++++++++++++++++"
       write(*,*) "+ EFG-Z4 Version 10-29-10 12:00 +"
       write(*,*) "+ Author: Fredy W. Aquino       +"
       write(*,*) "+++++++++++++++++++++++++++++++++"
       return
       end

       subroutine get_munu_rhos_symm(g_munu,ipol,
     &                               l_ilst,k_ilst,     
     &                               l_jlst,k_jlst,                                        
     &                               nlst,
     &                               g_munu_rhoS,pos)
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
       integer nlst,count,pos
       integer l_ilst,k_ilst,     
     &         l_jlst,k_jlst
       integer jlo,jhi
       double precision v(nlst)

       if (ipol.gt.1) then
        call ga_add(1.0d00,g_munu(1),
     &              1.0d00,g_munu(2),g_munu(1))
       endif
       call ga_gather(g_munu(1),v,
     &                int_mb(k_ilst),int_mb(k_jlst),
     &                nlst)
c      Now accumulates ith-chunk on g_munu_rhoS
       jlo=pos*nlst+1
       jhi=pos*nlst+nlst
       call ga_put(g_munu_rhoS,1,1,jlo,jhi,v,1)
       return
       end
c=======================================================
c================ Fredy Aquino's routines ======== END
c=======================================================
