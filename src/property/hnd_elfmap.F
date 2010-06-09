      subroutine hnd_elfmap(rtdb,basis,geom)
c
c
c     This routine calculates the electric field
c     for a given cc density at the atomic positions.
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "bas.fh"
c
      integer rtdb      ! [Input] rtdb
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_efld, k_efld
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision elf(3), elfld
      double precision rr, rr3
c
c     property grid points variables (MV)
      integer h_prp_c,i_prp_c
      integer ma_prp_type          
      integer nprp
      character*26 prp_date
      logical do_points
      integer l_tefld,k_tefld
      logical rtdb_par
c->d3p975
      character*30 theory
      integer nbf
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical ao_1prdm_read
      external ao_1prdm_read
c<-d3p975
c
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)
c
c
c->d3p975
c
c    reading g_dens from a file   
c    1. allocate the memory
      if(theory.eq.'tce') then
         ndens = 1
c
        do i = 1, ndens
         g_dens(i) = ga_create_atom_blocked(geom,basis,'density matrix')
         call ga_zero(g_dens(i))
        enddo
c
c    2. read the data 
c
        if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_elfmap: could not get nbf',0, BASIS_ERR)
c
        if(.not.ao_1prdm_read(nbf,g_dens(ndens)))
     1  call errquit('hnd_elfmap: ao_1prdm_read failed',0,0)

      else 
         call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
      endif
c 
c<-d3p975
c
c
c     ----- calculate electric field -----
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)
c
       call ecce_print_module_entry('EField')
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elfmap: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elfmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c     
c     define points for the calculation now 
c     either custom grid or (default) nuclei positions (M.V.)
c     -------------------------------------------------
      if(rtdb_get_info(rtdb, "prop:xyz", ma_prp_type, 
     >                 nprp, prp_date)) then
        nprp = nprp/3
        if (.not. ma_push_get(mt_dbl,3*nprp,'prop:xyz',h_prp_c,i_prp_c))
     &    call errquit('hnd_elfmap: prop:xyz',911,MA_ERR)
        if (.not. rtdb_get(rtdb,'prop:xyz',mt_dbl,
     >                      3*nprp,dbl_mb(i_prp_c)))
     &    call errquit('hnd_elfmap: prop:xyz failed',911,RTDB_ERR)
        do_points = .true.
      else
        nprp = nat
        if (.not. ma_push_get(mt_dbl,3*nat,'prop:xyz',h_prp_c,i_prp_c))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
        call dcopy(3*nat,dbl_mb(k_xyzpt),1,dbl_mb(i_prp_c),1)
        do_points = .false.
      end if
c
      if (.not. ma_push_get(mt_dbl,6*nprp,'efld pnt',l_efld,k_efld))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
c
c     total electric field array (M.V.)
c     --------------------------------
      if (.not. ma_push_get(mt_dbl,3*nprp,'tot field',l_tefld,k_tefld))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp,
     &                dbl_mb(k_efld),1)
c
c     ----- collect and output results of all points -----
c
      if (ga_nodeid().gt.0) goto 300
c
      write(luout,9997)
      do 230  iat=1,nprp
         xp = dbl_mb(i_prp_c  +3*(iat-1))
         yp = dbl_mb(i_prp_c+1+3*(iat-1))
         zp = dbl_mb(i_prp_c+2+3*(iat-1))
c
c     ----- add nuclear contribution -----
c
         elf(1) = dbl_mb(k_efld  +3*(iat-1))
         elf(2) = dbl_mb(k_efld+1+3*(iat-1))
         elf(3) = dbl_mb(k_efld+2+3*(iat-1))
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.1.0d-3) go to 210
            rr3=rr*rr*rr
            elf(1) = elf(1) - zan*xn/rr3
            elf(2) = elf(2) - zan*yn/rr3
            elf(3) = elf(3) - zan*zn/rr3
  210    continue
c
c     ----- form -efc- contribution -----
c           from cosmo point charges !!!!
c
         if (cosmo_last) then
c
c           Switch the run time data base temporarily to sequential
c           as access as this code is executed on process 0 only.
c
            rtdb_par = rtdb_parallel(.false.)
            if (.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc))
     &         call errquit('hnd_elfmap: rtdb get failed for nefc ',911,
     &         RTDB_ERR)
            if (.not.ma_push_get(mt_dbl,nefc*3,'efcc',l_efcc,k_efcc))
     &         call errquit('hnd_elfmap: malloc k_efcc fail',911,ma_err)
            if (.not.ma_push_get(mt_dbl,nefc,'efcz',l_efcz,k_efcz))
     &         call errquit('hnd_elfmap: malloc k_efcz fail',911,ma_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,
     &         dbl_mb(k_efcc))) call
     &         errquit('hnd_elfmap: rtdb get failed efcc',912,rtdb_err)
            if (.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,
     &         dbl_mb(k_efcz))) call
     &         errquit('hnd_elfmap: rtdb get failed efcz',913,rtdb_err)
c
c           Switch run time data base access back to original setting.
c
            rtdb_par = rtdb_parallel(rtdb_par)
            do 220 i = 1,nefc
               xn = dbl_mb(k_efcc+3*(i-1)  ) - xp
               yn = dbl_mb(k_efcc+3*(i-1)+1) - yp
               zn = dbl_mb(k_efcc+3*(i-1)+2) - zp
               rr =  sqrt(xn*xn + yn*yn + zn*zn)
               if(rr.lt.1.0d-3) then
                  write(luout,9993) xp,yp,zp,i
                  go to 220
               else
                  rr3=rr*rr*rr
                  elf(1) = elf(1) - dbl_mb(k_efcz+i-1)*xn/rr3
                  elf(2) = elf(2) - dbl_mb(k_efcz+i-1)*yn/rr3
                  elf(3) = elf(3) - dbl_mb(k_efcz+i-1)*zn/rr3
               endif
  220       continue
            if (.not.ma_chop_stack(l_efcc)) call
     &         errquit('hnd_elfmap: chop stack l_efcc',913,ma_err)
         endif
c
         dbl_mb(k_tefld  +3*(iat-1)) = elf(1)
         dbl_mb(k_tefld+1+3*(iat-1)) = elf(2)
         dbl_mb(k_tefld+2+3*(iat-1)) = elf(3)

         elfld   = sqrt(elf(1)*elf(1) + elf(2)*elf(2) + elf(3)*elf(3))
         if(do_points) then
c
c        skipping printout due to potentially large number of points
c        -----------------------------------------------------------
c           write(luout,9995) iat,"point",xp,yp,zp,elf(1),elf(2),elf(3),
c     &                       elfld
         else
           if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &        errquit('hnd_elfmap: geom_cent_tag failed',0,GEOM_ERR)
           if (.not.geom_tag_to_element(at_tag,symbol,element,atn)) then
             if(symbol.ne."bq") call
     &         errquit('hnd_elfmap: geom_tag_to_element failed',
     &                 0,GEOM_ERR)
           end if
           write(luout,9995) iat,symbol,xp,yp,zp,elf(1),elf(2),elf(3),
     &                       elfld
c
c          ----- store ecce data -----
c
           call ecce_print1_char('atom name',symbol,1)
           call ecce_print1('EField components',MT_DBL,elf,3)
           call ecce_print1('EField',MT_DBL,elfld,1)
         
         end if
c
  230 continue ! Assembling and printing next atom
c
      call ecce_print_module_exit('EField','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  300 call ga_sync()

c
c     if custom grid is requested save final electric 
c     into rtdb(M.V.)
c     -----------------------------------------------
      if(do_points) then
        if (.not. rtdb_put(rtdb,'prop:efield_xyz',mt_dbl,
     >                        3*nprp,dbl_mb(k_tefld)))
     &      call errquit('hnd_elfmap: efield_xyz failed',0,RTDB_ERR)
      end if
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_tefld)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_efld failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_efld)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_efld failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_prp_c)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('elfmap: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,14(1H-),/,10x,'Electric field',
     1       /,10x,14(1H-),/)
 9997 format(3x,'Atom ',6x,'X',9x,'Y',9x,'Z',8x,1x,
     1 15x,'Electric field (a.u.)',/,
     2 39x,7x,'X',14x,'Y',14x,'Z',11x,'Field',/,
     3 2x,96(1h-))
 9995 format(i5,1x,a2,3f10.5,4f15.6)
 9996 format(' --- Warning - electric field at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9994 format(' 1 a.u. = 0.171524 10**(-8) dyn/esu ',/)
 9993 format(' --- warning - electric field at ',
     1 3f10.5,' . contribution from  -EFC-  ',i3,' ignored')
      end
