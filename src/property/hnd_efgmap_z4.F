      subroutine hnd_efgmap_z4(rtdb,basis,geom)
c
c $Id$
c
c     This routine calculates the electric field gradient and 
c     the orientation of the EFG for a given density at the
c     atomic positions.
c Description : Modified from hnd_efgmap()
c               Adding quasi-relativistic corrections through
c               small-component density \rhoS
c Author      : Fredy W. Aquino
c Date        : 06-22-11
c Date        : 08-13-11 (Revised and improved)
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

      integer rtdb      ! [Input] rtdb  
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
      character*255 zorafilename
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, j, ij
      integer l_xyzpt , k_xyzpt, 
     &        l_xyzpt1, k_xyzpt1, 
     &        l_zanpt , k_zanpt, 
     &        l_efgs  , k_efgs
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer l_AtNr,k_AtNr ! to store indices of selected atoms
      integer g_rhoS,g_Atnr1 
      logical dft_zoraEFGz4_read
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision vec(3,3), eig(3), a(6)
      double precision pi, deg, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
      double precision rr, rr5, eta
      logical status
c     bq variables (MV)
      logical dobq
      integer bq_ncent,i_cbq,i_qbq
      integer nat_slc,typeprop
      double precision elpotbq
c
      logical docosmo
      integer ncosbq
c
      integer nder,l_rhoS,k_rhoS
      integer ii,jj,count_efgtyp   
      integer g_densz4(3)               
      double precision sum_efgs           
      integer ipol,pos,indx,indy,indx1,indx2         
      external get_densz4 ,               
     &         dft_zoraEFGz4_read,util_file_name,
     &         hnd_elfcon,get_slctd_atoms
      integer iat1,icalczora ! Added by FA
      logical skip_efgz4AOev
c ----- Definitions for NLMO analysis ---- START
      external hnd_elfcon_symm
      integer acc_vec,l_tvec,k_tvec,efgfile             
      integer g_Nuc,g_tvec ! for nbo analysis
      integer l_EFGnuc,k_EFGnuc,indxEFGnuc ! for create_munu4nbo()
c ----- Definitions for NLMO analysis ---- END
c
c     Initialize integrals

      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     ----- calculate electric field gradient -----

      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)

      pi  = acos(-1.0d0)
      deg = 180.0d0/pi
      call ecce_print_module_entry('efg')
c
c     ----- define points for calculation -----
c           1. nuclei
c ------- Read (nat,atmnr) --------- START
         status=geom_ncent(geom,nat)   
      if (.not.ma_alloc_get(
     &       mt_int,nat,'nmt tmp',l_AtNr,k_AtNr))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
         typeprop=1 ! =1 EFG =2 Shieldings =3 Hyperfine  
         call get_slctd_atoms(nat_slc,       ! out: selected atoms
     &                        int_mb(k_AtNr), ! out: list of selected atom nr.     
     &                        nat,           ! in : total nr atoms in molecule            
     &                        rtdb,          ! in : rdt  handle
     &                        typeprop)      ! in : =1,2,3=EFG,Shieldings,Hyperfine
      if (ga_nodeid().eq.0) then
       write(*,*) 'nat_slc=',nat_slc
       do i=1,nat_slc
        write(*,7) i,int_mb(k_AtNr+i-1)
 7      format('In hnd_efgmap_z4:: atomnr(',i3,')=',i5)
       enddo
      endif
c ------- Read (nat,atmnr) --------- END
      if (.not.ma_alloc_get(
     &       mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
      if (.not.ma_alloc_get(
     &       mt_dbl,6*nat,'efg pnt',l_efgs,k_efgs))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
      if (.not.ma_alloc_get(
     &       mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
      do 30 iat=1,nat
         status=geom_cent_get(geom,iat,at_tag,
     &                        dbl_mb(k_xyzpt+3*(iat-1)),
     &                        dbl_mb(k_zanpt+iat-1))
   30 continue
c ======== ga arrays for munu_nbo (NLMO analysis) === START
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
c ======== ga arrays for munu_nbo (NLMO analysis) === END
c
c ++++++ Reading EFGz4 data from file ++++++ START
c       Note.- lbl_efgz4 defined in zora.fh
        call util_file_name(lbl_efgz4,.false.,.false.,zorafilename)
        icalczora = 0  ! initialize the flag
        if (.not.dft_zoraEFGz4_read(
     &              zorafilename,
     &              nat_slc,
     &              nat,
     &              g_AtNr1,
     &              g_rhoS)) icalczora=1 
c ++++++ Reading EFGz4 data from file ++++++ END
c      if (ga_nodeid().eq.0) 
c     &  write(*,*) '-------hnd_efgmat_z4: g_rhoS ---------- START'
c      call ga_print(g_rhoS)
c      if (ga_nodeid().eq.0) 
c     &  write(*,*) '-------hnd_efgmat_z4: g_rhoS ---------- END'
c  Allocate memory for l_rhoS,k_rhoS 
       if (.not.ma_alloc_get(
     &       mt_dbl,nat_slc*6*2,'rhoS',l_rhoS,k_rhoS))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
      call ga_get(g_rhoS,1,1,1,nat_slc*6*2,dbl_mb(k_rhoS),1)
c --------- for NLMO analysis --------------- START
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
      if (efgfile.eq.1) then ! ------- efgfile-if++++ START
c =========== allocate memory for EFGnuc == START
       if (.not.ma_alloc_get(
     &       mt_dbl,nat_slc*6,'EFGnuc',l_EFGnuc,k_EFGnuc))
     &    call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR)
       call dcopy(nat_slc*6,0.0d0,0,dbl_mb(k_EFGnuc),1) ! reset 
c =========== allocate memory for EFGnuc == END
c  Allocate memory for l_tvec,k_tvec --- start 
       if (.not.ma_alloc_get(
     &       mt_dbl,nat_slc*3*3,'tvec',l_tvec,k_tvec)) 
     &     call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR) 
       call dcopy(nat_slc*3*3,0.0d0,0,dbl_mb(k_tvec),1) ! reset 
c  Allocate memory for l_tvec,k_tvec --- end
      endif ! ------------------------ efgfile-if++++ END
c --------- for NLMO analysis --------------- END
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
c      count_efgtyp=1
      count_efgtyp=0 ! NUM-EFG, z4-EFG
      do zora_calc_type=4,3,-1
       do_NonRel=.false.
       if (zora_calc_type.eq.4) do_NonRel=.true.       
       call get_densz4(rtdb,basis,geom,g_densz4)     
       nder=2
c --------- for NLMO analysis ---------------- START
c         if (ga_nodeid().eq.0)
c     &    write(*,*) "Using elfcon with munu-symmetry:"
c ---- extract selected atoms coordinates ----- START
       if (.not.ma_alloc_get(
     &       mt_dbl,3*nat_slc,'xyz pnt1',l_xyzpt1,k_xyzpt1)) 
     &     call errquit('hnd_efgmap_z4: ma failed',0,MA_ERR) 
         do iat1=1,nat_slc
          iat=int_mb(k_AtNr+iat1-1)
          indx1=k_xyzpt1+3*(iat1-1)
          indx2=k_xyzpt +3*(iat -1)     
          dbl_mb(indx1  )= dbl_mb(indx2  )
          dbl_mb(indx1+1)= dbl_mb(indx2+1)
          dbl_mb(indx1+2)= dbl_mb(indx2+2)
          if (ga_nodeid().eq.0) then
           write(*,12) iat1,iat,
     &                 dbl_mb(indx1),
     &                 dbl_mb(indx1+1),
     &                 dbl_mb(indx1+2)
 12        format('Atom(',i3,',',i3,')=(',
     &            f15.8,',',f15.8,',',f15.8,')')
          endif
         enddo
c ---- extract selected atoms coordinates ----- END
         call hnd_elfcon_symm(basis,            ! in: basis handle
     &                        geom,             ! in: geom  handle
     &                        g_densz4(3),      ! in: electron density
     &                        dbl_mb(k_xyzpt1), ! in: (x,y,z) centers
     &                        nat_slc,          ! in: number of centers
     &                        dbl_mb(k_efgs),   !out: EFG values at centers
     &                        nder,             ! in: =2 for second derivative
     &                        efgfile)          ! in: efgfile=0,1= NO,YES NLMONBO analysis
       if (.not.ma_free_heap(l_xyzpt1)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_xyzpt1',0, MA_ERR)

       if (ga_nodeid().eq.0) then ! START-if-do-it-once
        write(*,112)
        do iat1=1,nat_slc 
          iat=int_mb(k_AtNr+iat1-1)
c ------- get Atom name: symbol ----------- START
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_efgmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) then 
            if(symbol.ne."bq") call
     &      errquit('hnd_efgmap: geom_tag_to_element failed',0,GEOM_ERR)
         end if
c ------- get Atom name: symbol ----------- END
         indx=k_efgs+6*(iat1-1)
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
         indx=k_efgs+6*(iat1-1)
         indy=k_rhoS+6*(iat1-1)+6*nat_slc*count_efgtyp   
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
       count_efgtyp=count_efgtyp+1 ! NUM-EFG, z4-EFG
        do i=1,3
         if (.not. ga_destroy(g_densz4(i))) call errquit(
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)
        enddo
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
      docosmo = .false.
      ncosbq = 0
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq).and.(ncosbq.gt.0))
     &   docosmo = .true.
c
c     ----- collect and output results of all points -----
c
      if (docosmo) then
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
      end if ! docosmo
c
      acc_vec=0 ! reset tvec accumulator <-- For NLMO analysis
c
      do 230  iat=1,nat_slc
         iat1=int_mb(k_AtNr+iat-1)
         xp = dbl_mb(k_xyzpt  +3*(iat1-1))
         yp = dbl_mb(k_xyzpt+1+3*(iat1-1))
         zp = dbl_mb(k_xyzpt+2+3*(iat1-1))
c
c     ----- add nuclear contribution -----
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

         if (docosmo) then
            do i = 1,nefc
               xn = dbl_mb(k_efcc+3*(i-1)  ) - xp
               yn = dbl_mb(k_efcc+3*(i-1)+1) - yp
               zn = dbl_mb(k_efcc+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if (rr.lt.1.0d-3) then
                  if (ga_nodeid().eq.0) write(luout,9993) xp,yp,zp,i
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
         end if ! docosmo
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
                  if (ga_nodeid().eq.0) write(luout,9993) xp,yp,zp,i
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
c ------- Adding modified electronic part + nuclear contribution
         indx=k_efgs+6*(iat-1)
         dbl_mb(indx  )=dbl_mb(indx  )+2.0d0*efgxx - efgyy - efgzz
         dbl_mb(indx+1)=dbl_mb(indx+1)+2.0d0*efgyy - efgxx - efgzz
         dbl_mb(indx+2)=dbl_mb(indx+2)+2.0d0*efgzz - efgxx - efgyy
         dbl_mb(indx+3)=dbl_mb(indx+3)+3.0d0*efgxy
         dbl_mb(indx+4)=dbl_mb(indx+4)+3.0d0*efgxz
         dbl_mb(indx+5)=dbl_mb(indx+5)+3.0d0*efgyz
c============= Store EFG nuclear for create_munu4nbo()== START
         if (efgfile.eq.1) then
          indxEFGnuc=k_EFGnuc+6*(iat-1)
          dbl_mb(indxEFGnuc  )=3.0d0*efgxx
          dbl_mb(indxEFGnuc+1)=3.0d0*efgyy
          dbl_mb(indxEFGnuc+2)=3.0d0*efgzz
          dbl_mb(indxEFGnuc+3)=3.0d0*efgxy
          dbl_mb(indxEFGnuc+4)=3.0d0*efgxz
          dbl_mb(indxEFGnuc+5)=3.0d0*efgyz
          write(*,52) iat1,
     &        dbl_mb(indxEFGnuc  ),dbl_mb(indxEFGnuc+1), 
     &        dbl_mb(indxEFGnuc+2),dbl_mb(indxEFGnuc+3), 
     &        dbl_mb(indxEFGnuc+4),dbl_mb(indxEFGnuc+5)
 52       format('OUT:EFGnuc(',i5,')=(',f15.8,',',f15.8,',',
     &           f15.8,',',f15.8,',',f15.8,',',f15.8,')' )
         endif
c============= Store EFG nuclear for create_munu4nbo()== END    
c
c        ----- reorder into a as xx xy yy xz yz zz to form matrix -----
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
         if (.not. geom_cent_tag(geom,iat1,at_tag)) call
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
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print2('EFG tensor',MT_DBL,vec,3,3,3)
c
c        ----- print tensor components -----
         if (ga_nodeid().eq.0) then
           write(luout,9998) iat1,symbol,xp,yp,zp
           write(luout,9997) 
           write(luout,9995) (dbl_mb(k_efgs+6*(iat-1)+i),i=0,5)
         end if 
c
c        ----- diagonalize to get principal components and vectors -----

         call hnd_diag(vec,eig,3,.true.,.false.)
c ------- copy eigenvectors, vec for NLMO analysis ----- START
          if (efgfile.eq.1) then
           do ii=1,3
            do jj=1,3
             dbl_mb(k_tvec+acc_vec)=vec(ii,jj)
             acc_vec=acc_vec+1
            enddo 
           enddo 
          endif
c ------- copy eigenvectors, vec for NLMO analysis----- END
         eta  = abs( (eig(3)-eig(2)) / eig(1) )
         call ecce_print1('EFG eigenvalues',MT_DBL,eig,3)
         call ecce_print2('EFG eigenvectors',MT_DBL,vec,3,3,3)
         call ecce_print1('EFG asymmetry',MT_DBL,eta,1)
         if (ga_nodeid().eq.0) then
           write(luout,9992)
           write(luout,9991) eig(1),eig(2),eig(3),eta
           write(luout,9988) ((vec(i,j),j=1,3),i=1,3)
           write(luout,*) ' '
         end if
c
  230 continue ! Assembling and printing next atom
c
      if (docosmo) then
        if (.not.ma_chop_stack(l_efcc)) call 
     &         errquit('hnd_efgmap: chop stack l_efcc',913,ma_err)
      endif
c
      call ecce_print_module_exit('EFG','ok')
      call util_flush(luout)
c
c     ----- release memory block -----

  300 call ga_sync()
c -------- For NLMO analysis ----------------------------- START
      status = rtdb_parallel(.true.)   ! FA-04-23-10
      efgfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
      if (efgfile.eq.1) then ! ++++++ efgfile-if++++ START
         if (.not. ga_create(mt_dbl,1,nat_slc*6,
     &       'munu4nbo: g_Nuc',0,0,g_Nuc)) 
     $       call errquit('munu4nbo: g_Nuc', 0,GA_ERR)
         if (.not. ga_create(mt_dbl,1,nat_slc*9,
     &       'munu4nbo: g_tvec',0,0,g_tvec)) 
     $       call errquit('munu4nbo: g_tvec', 0,GA_ERR)
c ===== store (EFGnuc,tvec) in global arrays === START
       call ga_dgop(msg_efgs_col,dbl_mb(k_EFGnuc),6*nat_slc,'+')
       call ga_dgop(msg_efgs_col,dbl_mb(k_tvec)  ,9*nat_slc,'+')
       call ga_put(g_Nuc ,1,1,1,nat_slc*6,dbl_mb(k_EFGnuc),1)
       call ga_put(g_tvec,1,1,1,nat_slc*9,dbl_mb(k_tvec),1)
c ===== store (EFGnuc,tvec) in global arrays === END
       call create_munu4nbo(rtdb,          ! in: rtd handle
     &                      g_tvec,        ! in: eigenvectors
     &                      nat_slc,       ! in: nr of selected atoms
     &                      int_mb(k_AtNr),! in: list of selected atoms
     &                      basis,         ! in: basis handle
     &                      g_Nuc)         ! in: EFG at nuclei
         if (.not. ga_destroy(g_Nuc)) call errquit(  ! destroy GA
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)   
         if (.not. ga_destroy(g_tvec)) call errquit( ! destroy GA
     &    'dft_zora_rhos: ga_destroy failed ',0, GA_ERR)         
c     ------- Deallocate MA memory ------
       if (.not.ma_free_heap(l_EFGnuc)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_EFGnuc',0, MA_ERR)
       if (.not.ma_free_heap(l_tvec)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_tvec',0, MA_ERR)
      endif ! +++++++++++++++++++++++ efgfile-if++++ END
c -------- For NLMO analysis ----------------------------- END
c
c ------- Destroy ga arrays -----------
      if (.not. ga_destroy(g_AtNr1)) call errquit(
     &  'hnd_efgmap_z4: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_rhoS)) call errquit(
     &  'hnd_efgmap_z4: ga_destroy failed ',0, GA_ERR)
c     efgfile=0 ! not doing NLMO analysis by default
c     status=rtdb_get(rtdb,'prop:efgfile',mt_int,1,efgfile) ! for NLMO analysis
c     if (efgfile.eq.0) then ! ONLY destroy if not doing NLMO analysis
c        if(.not.rtdb_get(rtdb,'zora:skip_efgz4AOev',
c    &             mt_log,1,skip_efgz4AOev)) then      
c           if (.not. ga_destroy(g_AtNr)) call errquit(
c    &         'hnd_efgmap_z4: ga_destroy failed ',0, GA_ERR)
c        endif
c     endif
c ------- Deallocate MA memory ------
       if (.not.ma_free_heap(l_rhoS)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_rhoS',0, MA_ERR)
       if (.not.ma_free_heap(l_zanpt)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_zanpt',0, MA_ERR)
       if (.not.ma_free_heap(l_efgs)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_efgs',0, MA_ERR)
       if (.not.ma_free_heap(l_xyzpt)) call
     &     errquit('hnd_efgmap_z4: ma_free_heap l_xyzpt',0, MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('hnd_efgmap_z4: ma_free_heap l_AtNr',0, MA_ERR)
c
c     Terminate integrals

      call schwarz_tidy()
      call int_terminate()
      return
 9999 format(/,10x,26(1h-),/,10x,'Z4-Electric field gradient',
     1       /,10x,26(1h-),/)
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
c +++++++++++++++=++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++ nbo utilities +++++++++++++++++++++++++++++++++++++++START
      subroutine create_munu4nbo(rtdb,   ! in: rtdb handle
     &                           g_tvec, ! in: eigenvectors
     &                           nat,    ! in: nr of selected atoms
     &                           atomnr, ! in: list of selected atoms
     &                           basis,  ! in: basis handle
     &                           g_Nuc)
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
      integer rtdb,basis
      integer g_dens,g_tvec,g_Nuc
      integer npol,ipolmunu,g_munuV6,g_munu_rhoS,
     &        g_zora_scale_munu(2),g_zsr
      integer g_tnp,g_acc,g_acc1
      integer iat1,iat,nat,atomnr(nat),
     &        i,j,k,m,n,ndir
      integer jlo,jhi,s,nbf,nsize,nsize1
      integer ind,nlst,count,Nel
c     Nel, Nr. of electrons
      integer Natoms_munu,Ndir_munu,atmnr_munu(nat)
c     Ndir_munu, Nr. of directions stored
c                =3  xx yy zz
      double precision coeff,fact,tm,tn
      double precision EFGmn,EFGmn_1(6),ac_val,scaleNuc
      double precision tmn(2),valNuc(2),chcdata(3)
      integer g_tnp1,jlo1,jhi1,jlo2,jhi2
      integer g_ov,g_ov1,g_dens1,g_munu_rot
      integer iind(2),jind(2),icalczora   
      logical dft_zoraEFGz4_NLMOAnalysis_read ! for read-nlmo-mat
      character*255 zorafilename              ! for read-nlmo-mat
      integer arr_ind(6,2)
       data ((arr_ind(j,i),i=1,2),j=1,6)
     &  /1,1,2,2,3,3,1,2,1,3,2,3/
      external dft_zoraEFGz4_NLMOAnalysis_read,wefgfile
c     --> To store ONLY munu principal components xx,yy,zz 
c     Note.- g_munu_rhoS is created in dft_zora_rhos.F
c            g_munuV6    is created in hnd_elfcon_symm.F
c     Legend:
c     nbf, Nr of basis functions
c     nlst=nbf*(nbf+1)/2
c     1. g_munuV6   , ANALYTIC munu (unique elements) size=nlst
c     2. g_munu_rhoS, EFG picture-change corrections munu
c                     munu (unique elements) size=nlst
c     3. g_zora_scale_munu(i) , scaling munu  
c                       i=1,ipolmunu (1, closed shell)
c                                    (2, open shell -not tested)
c                     munu (symmetric matrix) size=nbf*nbf
c     4. g_ov1      , overlap matrix (unique elements) size=nlst
c     5. Nel        , Number of electrons
c     6. EFGnuc(ndir,nat), Nuclear EFGs 
c                     ndir=6  xx yy zz xy xz yz
c                     nat, Nr. of atoms in molecular system
      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('munu: bas_numbf failed',555, BASIS_ERR)
      Natoms_munu=nat
      do i=1,Natoms_munu
       atmnr_munu(i)=atomnr(i)
      enddo
      Ndir_munu=3
      nlst=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
c ++++++ Read NLMO matrices +++++++++ START
      ndir=6
      ipolmunu=1
         if (.not. rtdb_get(rtdb, 'dft:ipol',mt_int,1,ipolmunu))
     $      call errquit('dft:ipol: rtdb_put failed',
     $                   555, RTDB_ERR)
      call util_file_name(lbl_nlmo,.false.,.false.,zorafilename)    
      icalczora = 0  ! initialize the flag
      if (.not.dft_zoraEFGz4_NLMOAnalysis_read(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &                nat, ! in: list of selected atoms 
     &           ipolmunu, ! in: 
     &  g_zora_scale_munu, ! out: 
     &        g_munu_rhoS, ! out: 
     &             g_dens, ! out:
     &           g_munuV6)) icalczora=1 
c ++++++ Read NLMO matrices +++++++++ END
c +++++++++ FA-04-29-10 +++++ START
c ===> Writing Number of electrons on rtdb
c      to be used in create_munu4nbo() defined
c      in hnd_efgmap_z4.F
         if (.not. rtdb_get(rtdb, 'prop:Nel',mt_int,1,Nel))
     $      call errquit('prop_input-EFGz4-nel: rtdb_put failed',
     $                   555, RTDB_ERR)
c +++++++++ FA-04-29-10 +++++ END
       call get_unique_elmat(g_dens,g_dens1,nlst,nbf)             ! out: g_dens1
       call get_unique_elmat(g_zora_scale_munu(1),g_zsr,nlst,nbf) ! out: g_zsr
c ================= compute overlap matrix ===== START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'munu4nbo: g_ov',
     $                       0,0,g_ov)) 
     $       call errquit('munu4nbo: g_ov', 0,
     &                    GA_ERR)
       call ga_sync()
       call ga_zero(g_ov)
       call int_1e_ga(basis,basis,g_ov,'overlap',.false.)
       call get_unique_elmat(g_ov,g_ov1,nlst,nbf) ! out: g_ov1
      if (.not. ga_destroy(g_ov)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
c Where: Z_uv = <u| p K^2/(4 m^2c^2) p |v>
       call ga_add(1.0d0,g_ov1,1.0d0,g_zsr,g_ov1) ! S_uv + Z_uv
c ================= compute overlap matrix ===== END
      ndir=6 ! Nr. of directions: xx,yy,zz,xy,xz,yz
      nsize=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
      nsize1=nsize*ndir   ! size of whole munu per atom
         if (.not. ga_create(mt_dbl,1,nsize,
     &                       'munu4nbo: g_tnp',
     $                       0,0,g_tnp)) 
     $       call errquit('munu4nbo: g_tnp', 0,
     &                    GA_ERR)
        call ga_zero(g_tnp)
         if (.not. ga_create(mt_dbl,1,nsize,
     &                       'munu4nbo: g_tnp1',
     $                       0,0,g_tnp1)) 
     $       call errquit('munu4nbo: g_tnp1', 0,
     &                    GA_ERR)
        call ga_zero(g_tnp1)
         if (.not. ga_create(mt_dbl,1,nsize,
     &                       'munu4nbo: g_acc',
     $                       0,0,g_acc)) 
     $       call errquit('munu4nbo: g_acc', 0,
     &                    GA_ERR)
        call ga_zero(g_acc)
         if (.not. ga_create(mt_dbl,1,nsize,
     &                       'munu4nbo: g_acc1',
     $                       0,0,g_acc1)) 
     $       call errquit('munu4nbo: g_acc1', 0,
     &                    GA_ERR)
      call ga_zero(g_acc1)
      call ga_add(1.0d0,g_munuV6,1.0d0,g_munu_rhoS,g_munuV6)
      if (.not. ga_destroy(g_munu_rhoS)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
         if (.not. ga_create(mt_dbl,1,nat*3*nlst,
     &                       'munu4nbo: g_munu_rot',
     $                       0,0,g_munu_rot)) 
     $       call errquit('munu4nbo: g_munu_rot', 0,
     &                    GA_ERR)
      do iat1=1,nat
        iat=atomnr(iat1)
        call ga_zero(g_acc1)
        do n=1,3  ! xx,yy,zz
         m=n ! For principal components ONLY
         call ga_zero(g_acc)
         ac_val=0.0d0
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
c ------- get coeffNuc --- START
c            ind=6*(iat-1)+s
            ind=6*(iat1-1)+s
            call ga_get(g_Nuc,1,1,ind,ind,valNuc,1)
            scaleNuc=valNuc(1)/(Nel*1.0d0)
c ------- get coeffNuc --- END
            jlo=nsize1*(iat1-1)+nsize*(s-1)+1
            jhi=jlo+nsize-1
            call ga_copy_patch('n',g_munuV6,1,1,jlo,jhi,
     &                             g_tnp   ,1,1,1  ,nsize)
            call ga_copy(g_tnp,g_tnp1)
c -------- in the future make this optional --- START
c Doing: (V_ee)_uv + (V_nuc)_uv
c        (V_ee)_uv  = (non-rel electronic efg)_uv + (picture-change corrections to efg)_uv 
c        (V_nuc)_uv = V_nuc / N_el [ S_uv + Z_uv ]
c Where: (V_ee)_uv  = g_tnp1
c        (V_nuc)_uv = scaleNuc * g_ov1   
           call ga_add(1.0d0,g_tnp1,scaleNuc,g_ov1,g_tnp1)   
c -------- in the future make this optional --- END
           call ga_sync()
           call ga_add(1.0d0,g_acc,coeff,g_tnp1,g_acc)    
c +++++++++++++++++++++++ CHECK-1 ++ START
           jlo1=nbf+1
           jhi1=nsize
           call ga_scale_patch(g_tnp1,1,1,jlo1,jhi1,2.0d0)
           EFGmn_1(s)=ga_ddot(g_tnp1,g_dens1)         
           ac_val=ac_val+EFGmn_1(s)*coeff
c +++++++++++++++++++++++ CHECK-1 ++ END
         enddo ! end-loop-s
          call ga_add(1.0d0,g_acc1,1.0d0,g_acc,g_acc1)   
c ====== Store final munu matrices === START
         jlo2=nlst*Ndir_munu*(iat1-1)+
     &        nlst*(n-1)+1
         jhi2=jlo2+nlst-1
         call ga_copy_patch('n',g_acc     ,1,1,   1,nlst,
     &                          g_munu_rot,1,1,jlo2,jhi2)        
c ====== Store final munu matrices === END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== START
c ==== sum (g_acc .* g_dens1 + Nuclear CONTRIB) 
c      = TOTAL EFG diagonalized   
         jlo1=1+nbf
         jhi1=nsize
         call ga_scale_patch(g_acc,1,1,jlo1,jhi1,2.0d0)
         EFGmn=ga_ddot(g_acc,g_dens1)
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
        enddo ! end-loop-n
c ===== Transform xx_munu to 2xx_munu-(yy_munu+zz_munu) = START
c       or                   3xx_munu-(xx_munu+yy_munu+zz_munu)
         do i=1,3
          jlo1=Ndir_munu*nlst*(iat1-1)+nlst*(i-1)+1
          jhi1=jlo1+nlst-1     
          call ga_add_patch(1.0d0            ,g_munu_rot,1,1,jlo1,jhi1,
     &                      -0.333333333333d0,g_acc1    ,1,1,1   ,nlst,
     &                                        g_munu_rot,1,1,jlo1,jhi1)  
         enddo ! end-loop-i
c ===== Transform xx_munu to 2xx_munu-(yy_munu+zz_munu) = END
      enddo ! end-loop-iat
c ===== CHECK EFG reproducing EFG eigenvalues from g_munu_rot== START
          if (ga_nodeid().eq.0)
     &     write(*,*) 'CHCooooooooooooo',
     &                ' NW-EFG: Summary C+HC data [au] ',
     &                'ooooooooooooooooooooooo START'
         do iat1=1,nat
          iat=atomnr(iat1)
          do i=1,3
           jlo=nlst*Ndir_munu*(iat1-1)+
     &         nlst*(i-1)+1
           jhi=jlo+nlst-1
           call ga_copy_patch('n',g_munu_rot,1,1,jlo,jhi,
     &                            g_tnp     ,1,1,  1,nlst) ! out: g_tnp
           jlo1=nbf+1
           jhi1=nlst
           call ga_scale_patch(g_tnp,1,1,jlo1,jhi1,2.0d0)
           chcdata(i)=ga_ddot(g_tnp,g_dens1)
          enddo ! end-loo-i    
          if (ga_nodeid().eq.0) then
           write(*,23) iat,
     &                 chcdata(1),            ! dia-x
     &                 chcdata(2),chcdata(3)  ! dia-y,z
 23        format(' CHC   EFG(xx,yy,zz)(',i3,')=(', 
     &            f15.8,',',f15.8,',',f15.8,')') 
          endif         
         enddo ! end-loop-iat1
          if (ga_nodeid().eq.0)
     &     write(*,*) 'CHCooooooooooooo',
     &                ' NW-EFG: Summary C+HC data [au] ',
     &                'ooooooooooooooooooooooo END'
c ===== CHECK EFG reproducing EFG eigenvalues from g_munu_rot== END
c ===== Reorder matrix before-output ==== START
        call reorder_munu(g_munu_rot,nat,nlst,nbf,Ndir_munu)
c ===== Reorder matrix before-output ==== END
c ------ destroy unnecessary GAs 
      if (.not. ga_destroy(g_munuV6)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_zsr)) call errquit(
     &  'create_munu4nbo-2: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp1)) call errquit(
     &  'create_munu4nbo-4: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp)) call errquit(
     &  'create_munu4nbo-5: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc)) call errquit(
     &  'create_munu4nbo-6: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc1)) call errquit(
     &  'create_munu4nbo-7: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_ov1)) call errquit(
     &  'create_munu4nbo-3: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_dens)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_dens1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)

      if (.not. ga_destroy(g_zora_scale_munu(1))) call errquit(
     &  'create_munu4nbo-8: ga_destroy failed ',0, GA_ERR)
       if (ipolmunu.gt.1) then
        if (.not. ga_destroy(g_zora_scale_munu(2))) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
       endif
c ------ Free allocated memory
       call wefgfile(rtdb,
     &               g_munu_rot,
     &               nlst,         
     &               Ndir_munu,
     &               Natoms_munu,
     &               atmnr_munu,
     &               Nel)
      if (.not. ga_destroy(g_munu_rot)) 
     &          call errquit( ! destroy GA - FA
     &  'wefgfile: ga_destroy failed ',0, GA_ERR)
      return
      end    
      
      subroutine get_unique_elmat(g_a,g_b,nlst,nbf)
c
c     Purpose : g_a  is a 2-d array that contains a symmetric matrix
c               dim(g_a) = nbf * nbf
c               g_b are unique elements of g_a
c               dim(g_b)= nbf*(nbf+1)/2
c               Structure of g_b: 
c                    11 22 ... nbf         (main diagonal elem)
c                    21 31 32 41 42 43 ... (off  diagonal elem)  
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh" 
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "msgids.fh" 

      integer nlst    ! INPUT : Nr. elements to extract   
      integer nbf     ! INPUT
      integer g_a     ! INPUT
      integer g_b     ! OUTPUT
      integer l_v,k_v
      integer i,j,count
      double precision val

      if (.not.ma_alloc_get(mt_dbl,nlst,'varr',
     &                      l_v,k_v))
     &    call errquit('get_ijlst: ma failed',0,MA_ERR)
       count=1
       do i=1,nbf
         call ga_get(g_a,i,i,i,i,val,1)
         dbl_mb(k_v+count-1)=val
        count=count+1
       enddo
       do i=2,nbf
        do j=1,i-1
         call ga_get(g_a,i,i,j,j,val,1)
         dbl_mb(k_v+count-1)=val
         count=count+1
        enddo
       enddo
c      Now copy unique elements to g_b
         if (.not. ga_create(mt_dbl,1,nlst,
     &                       'munu4nbo: g_b',
     $                       0,0,g_b)) 
     $       call errquit('munu4nbo: g_b', 0,
     &                    GA_ERR)
       call ga_put(g_b,1,1,1,nlst,dbl_mb(k_v),1) ! dbl_mb(k_v) --> g_b
       if (.not.ma_free_heap(l_v)) call
     &    errquit('unique_elmat: ma_free_heap l_v',0, MA_ERR)
      return
      end

      subroutine reorder_munu(g_a,nlist,nlst,nbf,ndir)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh" 
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "rtdb.fh"

      integer g_a,g_b
      integer nlist,nlst,nbf,ndir
      integer i,j,k,iat,indx,count
      integer jlo,jhi
      integer l_v,k_v
      integer me,nproc
      double precision arr(nbf,nbf)

      me = ga_nodeid()
      nproc = ga_nnodes()
      if (.not.ma_alloc_get(mt_dbl,nlst,'varr',
     &                      l_v,k_v))
     &    call errquit('reorder_munu: ma failed',0,MA_ERR)
      do iat=1,nlist
       do k=1,ndir    
         call dcopy(nlst,0.0d0,0,dbl_mb(k_v),1) ! reset
          jlo=nlst*ndir*(iat-1)+
     &        nlst*(k-1)+1
          jhi=jlo+nlst-1
         if (ga_nodeid().eq.0) then
          call ga_get(g_a,1,1,jlo,jhi,dbl_mb(k_v),1)
          count=0
          do i=1,nbf
              arr(i,i)=dbl_mb(k_v+count)
              count=count+1
          enddo
          do i=2,nbf
           do j=1,i-1
             arr(i,j)=dbl_mb(k_v+count)
             arr(j,i)=dbl_mb(k_v+count)
             count=count+1
           enddo
          enddo
c ------- Storing in correct order ('pack')
          count=0
          do i=1,nbf
           do j=1,i
            dbl_mb(k_v+count)=arr(i,j)
            count=count+1
           enddo
          enddo
         endif ! reordering done my master
c ------ Propagate dbl_mb(k_v) to other CPUs
         call ga_sync()
         call ga_dgop(msg_efgs_col,dbl_mb(k_v),nlst,'+')
         call ga_put(g_a,1,1,jlo,jhi,dbl_mb(k_v),1)        
       enddo ! end-loop-ndir
      enddo ! end-loop-iat
      if (.not.ma_free_heap(l_v)) call
     &    errquit('munu4nbo: ma_free_heap l_v',0, MA_ERR)
      return 
      end    
c +++++++++++ nbo utilities +++++++++++++++++++++++++++++++++++++++++  END
c +++++++++++++++=++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_slctd_atoms(nat_slc, ! out: selected atoms
     &                           atomnr,  ! out: list of selected atom nr.     
     &                           nat,     ! in : total nr atoms in molecule            
     &                           rtdb,    ! in : rdt  handle
     &                           typeprop)! in : =1,2,3=EFG,Shieldings,Hyperfine
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
        integer geom,rtdb
        integer i,nat,nat_slc,typeprop
        integer atomnr(nat)   
        character*12 lblprop1(3)
        character*15 lblprop2(3)
        data lblprop1/'efgz4:natoms',
     &                'giao:natoms ',
     &                'hfine:natoms'/
        data lblprop2/'efgz4:atom list',
     &                'giao:atom list ',
     &                'hfine:atom list'/
c ------- Read (nat,atmnr) --------- START
c         status=geom_ncent(geom,nat0)       
        if (.not. rtdb_get(rtdb,lblprop1(typeprop),
     &                     mt_int,1,nat_slc))
     &       nat_slc=0 ! reset
         if (nat_slc.eq.0) then
           nat_slc=nat
           do i=1,nat
            atomnr(i)= i      
           enddo     
         else
          if (.not. rtdb_get(rtdb,lblprop2(typeprop),
     &                       mt_int,nat_slc,atomnr))
     $      call errquit('get_slctd_atoms: rtdb_get failed',
     $                   555, RTDB_ERR)
         endif
c ------- Read (nat,atmnr) --------- END
      return
      end
