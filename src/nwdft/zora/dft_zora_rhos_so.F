c=======================================================
c================ Fredy Aquino's routines ======== START
c=======================================================
      subroutine get_rhoS_so(
     &                   rtdb,
     &                   g_dens_ini,nexc,
     &                   geom, 
     &                   ao_bas_han,
     &                   nbf,nbf_ao,nbf_mo,
     &                   g_dens,
     &                   g_moso,
     &                   g_zora_scale_sf,
     &                   noc,
     &                   nocc, 
     &                   ipol)
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

      integer g_dens(2)
      integer g_moso(2)
      integer g_zora_scr(2)        
      integer g_dens_ini(2)   
      integer g_zora_scale_sf(2)
      integer noc,nocc(2)
      integer ipol
      integer geom
      integer ao_bas_han                      
      integer i,pos_dat,
     &        do_zgc_old,nexc
      integer iat,nat,noc1,ac   
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt    
      logical status               
      character*16 element, at_tag  
      double precision zora_eint,toac(5) 
c ++++++++++++++++++++++++++++
      character*100 fname_RSLC
      integer stat_read
      integer iat1,l_AtNr,k_AtNr
c ++++++++++++++++++++++++++++
      integer rtdb
      integer g_densZ4so(5)
      integer g_acc  
      integer read_SLCTD_EFG_Atoms
      external get_densZ4_so,read_SLCTD_EFG_Atoms   
      integer nbf,nbf_ao,nbf_mo
c
c---- Input global variables (defined in zora.fh):
c     1. zora_calc_type  ! =3 -> ZORA4-EFG
c     2. so_term         ! =0 -> ZORA4-spin-free
c     3. xyz_EFGcoords(i) i=1,2,3
c     4. zora_Qpq=1,6 : xx,yy,zz,xy,xz,yz 
      if (ga_nodeid().eq.0)
     &  write(*,*) "nocc(1)=",nocc(1)," nocc(2)=",nocc(2)
c================================================   
      if (ga_nodeid().eq.0) write(*,17) zoracutoff_EFG
 17   format('zoracutoff_EFG=',f15.8)
c================================================
      if (ipol.eq.1) then
        noc1=nocc(1)
      else if (ipol.eq.2) then 
        noc1=nocc(1)+nocc(2)    
      endif
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
     &                     'get_rhoS: g_AtNr',
     $       0,0,g_AtNr))
     $       call errquit('get_rhoS: g_AtNr', 0,
     &       GA_ERR)
        call ga_zero(g_AtNr)
        stat_read=read_SLCTD_EFG_Atoms
     &             (rtdb,nat,nlist,g_AtNr)
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
     &                     'get_rhoS: g_rhoS',
     $       0,0,g_rhoS)) call errquit('get_rhoS: g_rhoS', 0,
     &       GA_ERR)
      call ga_zero(g_rhoS)
      pos_dat=1
c      do zora_calc_type =3,4 !  ZORA4-EFG,NUM-EFG
      do zora_calc_type =4,3,-1 !  NUM-EFG,ZORA4-EFG
       call get_densZ4_so(rtdb,geom,ao_bas_han,
     &                    nbf,nbf_ao,nbf_mo,
     &                    g_dens,g_moso,
     &                    noc,nocc,g_densZ4so)
       write(*,*) "zora_calc_type=",zora_calc_type
c      write(*,*) "---Z4-dft_scf-g_dens(1) array ------- START"
c      call ga_print(g_dens_ini(1))
c      write(*,*) "---Z4-dft_scf-g_dens(1) array ------- END"
       do iat1=1,nlist ! nlist<=nat
        iat=dbl_mb(k_AtNr+iat1-1)
        status=geom_cent_get(geom,iat,at_tag,
     &                       dbl_mb(k_xyzpt+3*(iat-1)),
     &                       dbl_mb(k_zanpt+iat-1))
        xyz_EFGcoords(1)= dbl_mb(k_xyzpt  +3*(iat-1))
        xyz_EFGcoords(2)= dbl_mb(k_xyzpt+1+3*(iat-1))
        xyz_EFGcoords(3)= dbl_mb(k_xyzpt+2+3*(iat-1))
        if (ga_nodeid().eq.0) then
         write(*,19) iat,xyz_EFGcoords(1),xyz_EFGcoords(2),
     &                   xyz_EFGcoords(3)
 19      format('xyz_EFG(',i2,')=(',f15.8,',',f15.8,',',
     &          f15.8,')')
        endif
        do zora_Qpq=1,6 ! xx,yy,zz,xy,xz,yz - FA
c------Generate munu A^{pq}_r ----- START
         so_term=0           ! ZORA-spin-free
         call ga_zero(g_zora_scr(1))
         call ga_zero(g_zora_scr(2))
         call zora_getv(rtdb,g_dens_ini,g_zora_scr,nexc)
c        write(*,*) "--------- efg-1 ------- START"
c        call ga_print(g_densZ4so(1))
c        call ga_print(g_zora_scr(1))
c        write(*,*) "--------- efg-1 ------- END"
c        write(*,*) "--------- efg-2 ------- START"
c        call ga_print(g_densZ4so(2))
c        call ga_print(g_zora_scr(2))
c        write(*,*) "--------- efg-2 ------- END"
         toac(1)=ga_ddot(g_densZ4so(1),g_zora_scr(1))         
         toac(2)=ga_ddot(g_densZ4so(2),g_zora_scr(2))   
         if (zora_calc_type.eq.3) then
          ac=3
          do so_term=1,3 ! z,y,x-component 
           call ga_zero(g_zora_scr(1))
           call zora_getv(rtdb,g_dens_ini,g_zora_scr,nexc)
c          write(*,*) "--------- efg-",ac," ------- START"
c          call ga_print(g_densZ4so(ac))
c          call ga_print(g_zora_scr(1))
c          write(*,*) "--------- efg-",ac," ------- END"
           toac(ac)=ga_ddot(g_densZ4so(ac),g_zora_scr(1)) 
           ac=ac+1
          end do ! so_term-loop
          zora_eint=toac(1)+toac(2)+         
     &              toac(3)+toac(4)+toac(5)
          if (ga_nodeid().eq.0) then
           write(*,155) zora_calc_type,iat,zora_Qpq,pos_dat,
     &                 toac(1),toac(2),toac(1)+toac(2),
     &                 toac(3),toac(4),toac(5),zora_eint
 155       format('zora-efg-so(',i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,',', f15.8,',',f15.8,',',
     &            f15.8,',', f15.8,',',f15.8,')')  
          endif
         else 
          zora_eint=toac(1)+toac(2)
          if (ga_nodeid().eq.0) then
           write(*,15) zora_calc_type,iat,zora_Qpq,pos_dat,
     &                 toac(1),toac(2),toac(1)+toac(2),
     &                 zora_eint
 15        format('zora-efg-so(',i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,',', f15.8,',',f15.8,')')  
          endif
         endif 
         call ga_fill_patch(g_rhoS,1,1,pos_dat,pos_dat,zora_eint) 
c------Generate munu A^{pq}_r ----- END
        pos_dat=pos_dat+1
        end do ! zora_Qpq loop
       end do ! iat-loop
      end do ! zora_calc_type loop
c----deallocate memory - FA
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('get_rhoS_so, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('get_rhoS_so, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('dft_zora_utils: ma_free_heap l_AtNr',0, MA_ERR)
      return
      end

          subroutine get_densZ4_so(
     &                   rtdb,
     &                   geom, 
     &                   ao_bas_han,
     &                   nbf,nbf_ao,nbf_mo,
     &                   g_dens,g_moso,
     &                   noc,nocc,
     &                   g_densZ4so)
c    Purpose : Calculation of density matrix-like
c              for EFG-ZORA-4
c    Output  : g_densZ4so
c    Author  : Fredy Aquino
c    Note    : Modified from dft_zora_scale_so()
c    Date    : 12-09-09
       implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "zora.fh"

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_dens(2)
      integer g_moso(2)
      integer g_orb(2)
      integer g_dens_so(2)
      integer g_scr(2)

      integer l_vecsre, k_vecsre
      integer l_vecsim, k_vecsim

      integer iorb
      integer noc
      integer ispin
      integer geom
      integer ao_bas_han
      integer nbf
      integer nbf_ao
      integer nbf_mo

      integer g_densZ4so(5)     ! FA - Main output
      integer rtdb              ! FA
      integer i,j,ac,noc1,count ! FA
      integer nocc(2)           ! FA
      integer l_Ci,k_Ci         ! FA
      character*1  soxyz(3)     ! FA
      character*7  lbls(6)      ! FA
      character*34 lbls1(6)     ! FA
      double precision scale    ! FA
c     Using global input : g_Ci
c     defined in zora.fh and calculated
c     in dft_zora_scale_so()       
      soxyz(1)='z'
      soxyz(2)='y'
      soxyz(3)='x'
      lbls(1)='orbs re'
      lbls(2)='denmxre'
      lbls(3)='scrch 1'
      lbls(4)='orbs im'
      lbls(5)='denmxim'
      lbls(6)='scrch 2'
      lbls1(1)='get_densZ4_so: orb  real     error'
      lbls1(2)='get_densZ4_so: dens real     error'
      lbls1(3)='get_densZ4_so: ga_duplicate failed'
      lbls1(4)='get_densZ4_so: orb  imag     error'
      lbls1(5)='get_densZ4_so: dens imag     error'
      lbls1(6)='get_densZ4_so: ga_duplicate failed'

      noc1=nocc(1)+nocc(2)
 
cc    allocate memory
       if (.not.MA_Push_Get(MT_Dbl, nbf_mo, 'real vec aux',
     &          l_vecsre, k_vecsre))
     & call errquit('get_densZ4_so: cannot allocate vec',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl, nbf_mo, 'imag vec aux',
     &          l_vecsim, k_vecsim))
     & call errquit('get_densZ4_so: cannot allocate vec',0, MA_ERR)

      if (.not.ma_alloc_get(MT_Dbl, 2*noc1, 'Ci',
     &          l_Ci, k_Ci))
     & call errquit('get_densZ4_so: ma failed',0,MA_ERR)
      
      call ga_get(g_Ci,1,2,1,noc1,dbl_mb(k_Ci),2)
    
c     do i=1,2*noc1
c      write(*,18) i,dbl_mb(k_Ci+i-1)
c18    format('Ci(',i3,')=',f15.8)
c     end do
c     stop
c 
      ac=1
      do i=1,2
c     spin-orbit vector -          i=1(real),2(imaginary)  
       if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,lbls(ac),0,0,
     &                   g_orb(i)))
     & call errquit(lbls1(ac),0, GA_ERR)
       call ga_zero(g_orb(i))
       ac=ac+1
c     spin-orbit density matrix -  i=1(real),2(imaginary)
       if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,lbls(ac),0,0,
     &                   g_dens_so(i)))
     & call errquit(lbls1(ac),0, GA_ERR)
       call ga_zero(g_dens_so(i))
       ac=ac+1
c     scratch array                i=1(real),2(imaginary)
       if(.not.ga_duplicate(g_dens(i),g_scr(i),lbls(ac)))
     &  call errquit(lbls1(ac),1, GA_ERR)
        call ga_zero(g_scr(i))
       ac=ac+1
      end do ! end i-loop

      do i=1,5  ! Added by FA   
       if(.not.ga_duplicate(g_dens(1),g_densZ4so(i),'densZ4so'))
     &  call errquit('get_densZ4_so: ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_densZ4so(i))
      end do 
c
      do ispin=1,2
       iorb=ispin
       do count=1,nocc(ispin) 
        call ga_get(g_moso(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)
        call ga_zero(g_orb(1))
        call ga_put(g_orb(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)

        call ga_get(g_moso(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1) 
        call ga_zero(g_orb(2))
        call ga_put(g_orb(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1)

        call dft_densm_so(g_dens_so,g_orb,nbf_ao,noc)

        call ga_zero(g_scr(1))
        call ga_zero(g_scr(2))
        call ga_dens_sf(g_scr, g_dens_so(1), nbf_ao)

c        scale=1.0d0/dbl_mb(k_Ci+iorb-1)
         if (zora_calc_type.eq.3) scale=1.0d0/dbl_mb(k_Ci+iorb-1)
         if (zora_calc_type.eq.4) scale=1.0d0  ! NUM-EFG

c        write(*,16) ispin,count,dbl_mb(k_Ci+iorb-1)
c 16     format('Ci(',i3,',',i3,')=',f15.8)

         do i=1,2
          call ga_scale(g_scr(i),scale)       
          call ga_add(1.0d0,g_densZ4so(i),1.0d0,g_scr(i),g_densZ4so(i))
         end do
        j=3
         do i=1,3
          call ga_zero(g_scr(1))
          call ga_dens_so(g_scr(1),g_dens_so,nbf_ao,soxyz(i))
          call ga_scale(g_scr(1),scale)   
          call ga_add(1.0d0,g_densZ4so(j),1.0d0,g_scr(1),g_densZ4so(j))
          j=j+1
         end do
        iorb=iorb+2
       end do  ! count-loop
      end do   ! ispin-loop

c     deallocate memory
      if (.not. ma_chop_stack(l_vecsim)) 
     & call errquit('get_densZ4_so:l_vecsim', 0, MA_ERR)
      if (.not. ma_chop_stack(l_vecsre))
     & call errquit('get_densZ4_so:l_vecsre', 0, MA_ERR)
      if (.not. MA_free_heap(l_Ci))
     &  call errquit('get_densZ4_so:cannot free heap',111, MA_ERR)
      do i=1,2
       if (.not. ga_destroy(g_orb(i))) 
     & call errquit('get_densZ4_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_dens_so(i))) 
     & call errquit('get_densZ4_so: ga_destroy failed',0, GA_ERR)
       if (.not. ga_destroy(g_scr(i))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      end do

      return
      end

      subroutine hnd_efgmap_Z4_so(rtdb,basis,geom,
     &                            nbf,nbf_ao,nbf_mo,
     &                            g_dens_4SO,g_moso,
     &                            noc,nocc)
c
c $Id: hnd_efgmap.F,v 1.6 2005-05-27 01:51:04 edo Exp $
c
c     This routine calculates the electric field gradient and 
c     the orientation of the EFG for a given density at the
c     atomic positions.
c
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
#include "zora.fh" ! Added by FA
c
      integer rtdb      ! [Input] rtdb  
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, j, ij
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_efgs, k_efgs
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision vec(3,3), eig(3), a(6)
      double precision pi, deg, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
      double precision rr, rr5, eta
      logical status
c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq

      integer nder                       ! Added by FA
      integer l_rhoS,k_rhoS              ! Added by FA
      integer ii,jj                      ! Added by FA   
      integer g_densZ4(3)                ! Added by FA
      double precision sum_efgs          ! Added by FA
      integer ipol,pos,indx,indy         ! Added by FA
      integer count_efgtyp               ! Added by FA
      integer g_densZ4so(5)              ! Added by FA
      integer g_dens_4SO(2)              ! Added by FA
      integer g_moso(2)                  ! Added by FA
      integer nbf,nbf_ao,nbf_mo,         ! Added by FA
     &        noc,nocc(*)                   ! Added by FA
      external get_densZ4_so             ! Added by FA
      integer iat1 
      integer l_AtNr,k_AtNr              ! Added by FA
 
c      if (ga_nodeid().eq.0) write(*,*) "---invoke-g_rhoS---START"
c      call ga_print(g_rhoS)
c      if (ga_nodeid().eq.0) write(*,*) "---invoke-g_rhoS---END"
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     ----- calculate electric field gradient -----
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)
c
      pi  = acos(-1.0d0)
      deg = 180.0d0/pi
c
      call ecce_print_module_entry('efg')
c
c     ----- define points for calculation -----
c           1. nuclei
c
      status=geom_ncent(geom,nat)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,6*nat,'efg pnt',l_efgs,k_efgs))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_efgmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
         status=geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &                        dbl_mb(k_zanpt+iat-1))
   30 continue
c
c  Allocate memory for l_rhoS,k_rhoS 
      if (.not.ma_alloc_get(mt_dbl,nat*6*2,'rhoS',
     &    l_rhoS,k_rhoS))
     &    call errquit('hnd_efgmap_Z4: ma failed',0,MA_ERR)
      call ga_get(g_rhoS,1,1,1,nat*6*2,dbl_mb(k_rhoS),1)
c  Allocate memory for l_AtNr,k_AtNr 
      if (.not.ma_alloc_get(mt_dbl,nat,'AtNr',
     &    l_AtNr,k_AtNr)) 
     &    call errquit('hnd_efgmap_Z4: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
c      count_efgtyp=1
      count_efgtyp=0 ! NUM-EFG, Z4-EFG
      do zora_calc_type=4,3,-1
       call get_densZ4_so(rtdb,geom,basis,
     &                    nbf,nbf_ao,nbf_mo,
     &                    g_dens_4SO,g_moso,
     &                    noc,nocc,g_densZ4so)
c ----- compute g_densZ4_all=g_densZ4so(1)+g_densZ4so(2)
       call ga_add(1.0d00,g_densZ4so(1),
     &             1.0d00,g_densZ4so(2),
     &                    g_densZ4so(1))   
      nder=2
      call hnd_elfcon(basis,geom,g_densZ4so(1),
     &                dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_efgs),nder)
       if (ga_nodeid().eq.0) then ! START-if-do-it-once
        goto 123
        if (zora_calc_type.eq.4) then
         write(*,*) "++++++++++++++++++"
         write(*,*) "+    EFG-ZORA    +"
         write(*,*) "++++++++++++++++++"
        endif
        if (zora_calc_type.eq.3) then
         write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
         write(*,*) "+  EFG-ZORA-4 (with small-component-density) +"
         write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
        endif
 123    continue ! SKIP printout
        write(*,112)
        do iat1=1,nlist 
          iat=dbl_mb(k_AtNr+iat1-1)
c ------- get Atom name: symbol ----------- START
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
c ------- get Atom name: symbol ----------- END
         indx=k_efgs+6*(iat-1)
         efgxx = dbl_mb(indx)
         efgyy = dbl_mb(indx+1)
         efgzz = dbl_mb(indx+2)
         efgxy = dbl_mb(indx+3)
         efgxz = dbl_mb(indx+4)
         efgyz = dbl_mb(indx+5)
         sum_efgs=(efgxx+efgyy+efgzz)/3.0d0
         efgxx=efgxx-sum_efgs
         efgyy=efgyy-sum_efgs
         efgzz=efgzz-sum_efgs
         indx=k_efgs+6*(iat-1)
c        indy=k_rhoS+6*(iat-1)+6*nat*count_efgtyp   
         indy=k_rhoS+6*(iat1-1)+6*nlist*count_efgtyp
         if (zora_calc_type.eq.3) then      
          dbl_mb(indx  )=efgxx+dbl_mb(indy  )
          dbl_mb(indx+1)=efgyy+dbl_mb(indy+1)
          dbl_mb(indx+2)=efgzz+dbl_mb(indy+2)
          dbl_mb(indx+3)=efgxy+dbl_mb(indy+3)
          dbl_mb(indx+4)=efgxz+dbl_mb(indy+4)
          dbl_mb(indx+5)=efgyz+dbl_mb(indy+5)
          write(*,113) symbol,efgxx,efgyy,efgzz,efgxy,
     &                 efgxz,efgyz
          write(*,114) symbol,dbl_mb(indy),dbl_mb(indy+1),
     &                 dbl_mb(indy+2),dbl_mb(indy+3),
     &                 dbl_mb(indy+4),dbl_mb(indy+5)  
          write(*,115) symbol,dbl_mb(indx),dbl_mb(indx+1),
     &                 dbl_mb(indx+2),dbl_mb(indx+3),
     &                 dbl_mb(indx+4),dbl_mb(indx+5)  
         endif
         if (zora_calc_type.eq.4) then
          dbl_mb(indx  )=efgxx
          dbl_mb(indx+1)=efgyy
          dbl_mb(indx+2)=efgzz
          dbl_mb(indx+3)=efgxy
          dbl_mb(indx+4)=efgxz
          dbl_mb(indx+5)=efgyz
          write(*,116) symbol,efgxx,efgyy,efgzz,efgxy,
     &                 efgxz,efgyz
          write(*,117) symbol,dbl_mb(indy),dbl_mb(indy+1),
     &                 dbl_mb(indy+2),dbl_mb(indy+3),
     &                 dbl_mb(indy+4),dbl_mb(indy+5)  
         endif        
        end do ! iat-loop
       end if ! END-if-do-it-once
c       count_efgtyp=count_efgtyp-1
       count_efgtyp=count_efgtyp+1 ! NUM-EFG, Z4-EFG
      end do ! zora_calc_type loop
c ------- All-FA-formats------------------------------------ START
 112     format('====> Electronic contribution to EFG',
     &          ' in molecular frame (a.u.)',/,
     & 21x,'XX',12x,'YY',12x,'ZZ',12x,'XY',12x,'XZ',12x,'YZ',/,
     & 16x,82(1h-))
 113      format('EFG-elec(',a2,')=(',f13.8,',',f13.8,',',
     &          f13.8,',',f13.8,',',f13.8,',',f13.8,')')
 114      format('EFG-rhoS(',a2,')=(',f13.8,',',f13.8,',',
     &          f13.8,',',f13.8,',',f13.8,',',f13.8,')')  
 115      format('EFG-tot (',a2,')=(',f13.8,',',f13.8,',',
     &          f13.8,',',f13.8,',',f13.8,',',f13.8,')')       
 116      format('  ANALYT(',a2,')=(',f13.8,',',f13.8,',',
     &          f13.8,',',f13.8,',',f13.8,',',f13.8,')')
 117      format('  NUMERI(',a2,')=(',f13.8,',',f13.8,',',
     &          f13.8,',',f13.8,',',f13.8,',',f13.8,')')  
c ------- All-FA-formats------------------------------------ END
c
c     get bq structures if any (MV)
c     -----------------------------
      dobq = .false.
      if(geom_extbq_on()) then
        dobq = .true.
        bq_ncent = geom_extbq_ncenter()
        i_cbq = geom_extbq_coord()
        i_qbq = geom_extbq_charge()
      end if
c
c     ----- collect and output results of all points -----
c
      if (ga_nodeid().gt.0) goto 300
c     if (ga_nodeid().eq.0) then  ! ---- node-0-START
c
      do 230  iat=1,nat
         xp = dbl_mb(k_xyzpt  +3*(iat-1))
         yp = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp = dbl_mb(k_xyzpt+2+3*(iat-1))
c     ----- add nuclear contribution -----
c        efgxx = dbl_mb(k_efgs  +6*(iat-1))/3.0d0
c        efgyy = dbl_mb(k_efgs+1+6*(iat-1))/3.0d0
c        efgzz = dbl_mb(k_efgs+2+6*(iat-1))/3.0d0
c        efgxy = dbl_mb(k_efgs+3+6*(iat-1))/3.0d0
c        efgxz = dbl_mb(k_efgs+4+6*(iat-1))/3.0d0
c        efgyz = dbl_mb(k_efgs+5+6*(iat-1))/3.0d0
         efgxx = 0.0d0 ! FA
         efgyy = 0.0d0 ! FA
         efgzz = 0.0d0 ! FA
         efgxy = 0.0d0 ! FA
         efgxz = 0.0d0 ! FA
         efgyz = 0.0d0 ! FA

         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr = sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            rr5=rr*rr*rr*rr*rr
            efgxx = efgxx - zan*xn*xn/rr5
            efgyy = efgyy - zan*yn*yn/rr5
            efgzz = efgzz - zan*zn*zn/rr5
            efgxy = efgxy - zan*xn*yn/rr5
            efgxz = efgxz - zan*xn*zn/rr5
            efgyz = efgyz - zan*yn*zn/rr5
  210    continue
c
c     ----- form -efc- contribution -----
c           from cosmo point charges !!!!
c
         if (cosmo_last) then
            if (.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc))
     &         call errquit('hnd_efgmap: rtdb get failed for nefc ',911,
     &         RTDB_ERR)
            if (.not.ma_push_get(mt_dbl,nefc*3,'efcc',l_efcc,k_efcc))
     &         call errquit('hnd_efgmap: malloc k_efcc fail',911,ma_err)
            if (.not.ma_push_get(mt_dbl,nefc,'efcz',l_efcz,k_efcz))
     &         call errquit('hnd_efgmap: malloc k_efcz fail',911,ma_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,
     &         dbl_mb(k_efcc))) call 
     &         errquit('hnd_efgmap: rtdb get failed efcc',912,rtdb_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,
     &         dbl_mb(k_efcz))) call 
     &         errquit('hnd_efgmap: rtdb get failed efcz',913,rtdb_err)
            do i = 1,nefc
               xn = dbl_mb(k_efcc+3*(i-1)  ) - xp
               yn = dbl_mb(k_efcc+3*(i-1)+1) - yp
               zn = dbl_mb(k_efcc+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if (rr.lt.1.0d-3) then
                  write(luout,9993) xp,yp,zp,i
               else
               rr5=rr*rr*rr*rr*rr
               efgxx = efgxx - dbl_mb(k_efcz+i-1)*xn*xn/rr5
               efgyy = efgyy - dbl_mb(k_efcz+i-1)*yn*yn/rr5
               efgzz = efgzz - dbl_mb(k_efcz+i-1)*zn*zn/rr5
               efgxy = efgxy - dbl_mb(k_efcz+i-1)*xn*yn/rr5
               efgxz = efgxz - dbl_mb(k_efcz+i-1)*xn*zn/rr5
               efgyz = efgyz - dbl_mb(k_efcz+i-1)*yn*zn/rr5
            endif
            enddo
 220        continue
            if (.not.ma_chop_stack(l_efcc)) call 
     &         errquit('hnd_efgmap: chop stack l_efcc',913,ma_err)
         endif
c
c        adding external bq contributions(MV)
c        ----------------------------------
         if (dobq) then
            do i = 1,bq_ncent
               xn = dbl_mb(i_cbq+3*(i-1)  ) - xp
               yn = dbl_mb(i_cbq+3*(i-1)+1) - yp
               zn = dbl_mb(i_cbq+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if (rr.lt.1.0d-3) then
                  write(luout,9993) xp,yp,zp,i
               else
               rr5=rr*rr*rr*rr*rr
               efgxx = efgxx - dbl_mb(i_qbq+i-1)*xn*xn/rr5
               efgyy = efgyy - dbl_mb(i_qbq+i-1)*yn*yn/rr5
               efgzz = efgzz - dbl_mb(i_qbq+i-1)*zn*zn/rr5
               efgxy = efgxy - dbl_mb(i_qbq+i-1)*xn*yn/rr5
               efgxz = efgxz - dbl_mb(i_qbq+i-1)*xn*zn/rr5
               efgyz = efgyz - dbl_mb(i_qbq+i-1)*yn*zn/rr5
               endif
            end do
         end if
c
c        dbl_mb(k_efgs  +6*(iat-1)) = 2.0d0*efgxx - efgyy - efgzz
c        dbl_mb(k_efgs+1+6*(iat-1)) = 2.0d0*efgyy - efgxx - efgzz
c        dbl_mb(k_efgs+2+6*(iat-1)) = 2.0d0*efgzz - efgxx - efgyy
c        dbl_mb(k_efgs+3+6*(iat-1)) = 3.0d0*efgxy
c        dbl_mb(k_efgs+4+6*(iat-1)) = 3.0d0*efgxz
c        dbl_mb(k_efgs+5+6*(iat-1)) = 3.0d0*efgyz
c ------- Adding modified electronic part + nuclear contribution
         indx=k_efgs+6*(iat-1)
         dbl_mb(indx  )=dbl_mb(indx  )+2.0d0*efgxx - efgyy - efgzz
         dbl_mb(indx+1)=dbl_mb(indx+1)+2.0d0*efgyy - efgxx - efgzz
         dbl_mb(indx+2)=dbl_mb(indx+2)+2.0d0*efgzz - efgxx - efgyy
         dbl_mb(indx+3)=dbl_mb(indx+3)+3.0d0*efgxy
         dbl_mb(indx+4)=dbl_mb(indx+4)+3.0d0*efgxz
         dbl_mb(indx+5)=dbl_mb(indx+5)+3.0d0*efgyz
c
c        ----- reorder into a as xx xy yy xz yz zz to form matrix -----
c
         a(1) = dbl_mb(k_efgs  +6*(iat-1))        
         a(2) = dbl_mb(k_efgs+3+6*(iat-1))        
         a(3) = dbl_mb(k_efgs+1+6*(iat-1))        
         a(4) = dbl_mb(k_efgs+4+6*(iat-1))        
         a(5) = dbl_mb(k_efgs+5+6*(iat-1))        
         a(6) = dbl_mb(k_efgs+2+6*(iat-1))        
         ij=0
         do 241 i = 1, 3
         do 241 j = 1, i
            ij = ij + 1
            vec(i,j) = a(ij) 
            vec(j,i) = a(ij)
  241    continue
c
c        ----- store ecce data -----
c
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
c        geom_tag_to_element returns false for Bq elements (MV)
c        -----------------------------------------------------
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
c
c         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
c     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)

c         call ecce_print1_char('atom name',symbol,1)
c         call ecce_print2('EFG tensor',MT_DBL,vec,3,3,3)

c
c        ----- print tensor components -----
c
         write(luout,9998) iat,symbol,xp,yp,zp
         write(luout,9997) 
         write(luout,9995) (dbl_mb(k_efgs+6*(iat-1)+i),i=0,5)
c
c        ----- diagonalize to get principal components and vectors -----
c
c        FA: I found that for some few cases it halts when calling hnd_diag()
c            I don't know why.
c            Example: AuI_SARC1-ZORA_BS_CNT_Z4_BLYP_4.out (04-16-10)
c        write(*,*) "BEF. hnd_diag ..."
         call hnd_diag(vec,eig,3,.true.,.false.)
c        write(*,*) "AFT. hnd_diag ..."
         eta  = abs( (eig(3)-eig(2)) / eig(1) )
c
c         call ecce_print1('EFG eigenvalues',MT_DBL,eig,3)
c         call ecce_print2('EFG eigenvectors',MT_DBL,vec,3,3,3)
c         call ecce_print1('EFG asymmetry',MT_DBL,eta,1)
c
         write(luout,9992)
         write(luout,9991) eig(1),eig(2),eig(3),eta
         write(luout,9988) ((vec(i,j),j=1,3),i=1,3)
         write(luout,*) ' '
c
  230 continue ! Assemblin and printing next atom
        call ecce_print_module_exit('EFG','ok')
        call util_flush(luout)
c     endif ! ga_nodeid().eq.0 ---- node-0-END

c      count_efgtyp=count_efgtyp-1  ! Added by FA 03-12-10
c     end do ! zora_calc_type loop    Added by FA 03-12-10
c
c     ----- release memory block -----
c
  300 call ga_sync()
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_efgs)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_efgs failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_efgmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
c ---------------- added by FA ---- start
c     Cleanup of MA arrays that are not needed
c     deallocate memory -------- Added by FA
      if (.not.ma_free_heap(l_rhoS)) call
     &    errquit('hnd_efgmap_Z4: ma_free_heap l_rhoS',0, MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('dft_zora_utils: ma_free_heap l_AtNr',0, MA_ERR)
c ---------------- added by FA ---- end
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,23(1h-),/,10x,'Electric field gradient',
     1       /,10x,23(1h-),/)
 9998 format(/,1x,60(1h-),/,3x,'Atom',6x,'X',9x,'Y',9x,'Z',/,1x,60(1h-),
     1       /,i5,1x,a2,3f10.5,/,1x,60(1h-),/)
 9997 format(1x,'Electric field gradient in molecular frame (a.u.)',/,
     2 9x,'XX',13x,'YY',13x,'ZZ',13x,'XY',13x,'XZ',13x,'YZ',/,
     3 1x,90(1h-))
 9996 format(' --- Warning - electric field gradient at ',
     1 3F10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(1x,6f15.6,/)
 9994 format(' 1 a.u. = 0.324123 10**(16) esu/cm**3 ',
     1       ' ( or statvolts/cm**2 )',' = 0.97174 10**(22) v/m**2 ',/)
 9993 format(' --- Warning - electric field gradient at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
 9992 format(1x,'Principal components (a.u.) and orientation ',
     1       /,' of principal axis w.r.t. absolute frame',
     2       22x,'Asymmetry parameter eta',/,1x,86(1h-))
 9991 format(1x,3f15.6,16x,f15.6,/)
 9988 format(1X,3F15.6)
      end
c=======================================================
c================ Fredy Aquino's routines ======== END
c=======================================================
