      subroutine hnd_spinspin(rtdb,basis,geom)
c $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "apiP.fh"
c
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
c
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer ixy
      integer icom_atom, ipairs, i_pert, i_resp
      integer l_occ, k_occ, l_eval, k_eval
      integer l_xyz, k_xyz, l_zan, k_zan, l_tr, k_tr
      integer l_pr, k_pr, l_un, k_un, l_data, k_data
      integer g_dens(3), vectors(2), geomnew, i, j, i_off, ii, jj
      integer isonra, isonrb
      double precision atn, isotr, aniso, znuc
      character*3 scftyp
      character*16 tag, tagi, tagj
      character*32 element
      character*2 symbola, symbolb
c
      logical  cphf2, file_write_ga, file_read_ga, cphf, atom_gfac
      external cphf2, file_write_ga, file_read_ga, cphf, atom_gfac
c
      logical     oskel, status, notthere
      double precision gfaci, gfacj, value
c
      double precision pi, clight, mubohr, munuc, munul, planck
      double precision hbar, gfactor, masselec, echarge, bohrrad
      double precision uamu, amu, massprot, pramu, light, autohz
      double precision factor, nfc, nsd, nsdfc, npso, ndso
c
      parameter(pi=3.14159265358979323846264338327950288419d0)
      parameter(clight=2.99782458d8)                                    ! Speed of light in m/s
      parameter(mubohr=9.27400949d-24)                                  ! Bohr magneton in A m^2
      parameter(munuc=5.050783431d-27)                                  ! Nuclear magneton in A m^2
      parameter(munul=4.0d0*pi*1.0d-7)                                  ! Permeability of vacuum in m kg s^-2 A^-2
      parameter(planck=6.6260693d-34)                                   ! Planck constant in J s
      parameter(hbar=planck/(2.0d0*pi))                                 ! Plank / 2Pi
      parameter(gfactor=2.0023193043718d0)                              ! Electronic g-factor ge
      parameter(masselec=9.1093826d-31)                                 ! Rest mass electron in kg
      parameter(echarge=1.60217653d-19)                                 ! Elementary charge in A s
      parameter(bohrrad=0.5291772108d-10)                               ! Meters to Bohr
      parameter(uamu=1.66053886d-27)                                    ! Unified atomic mass unit in kg
      parameter(amu=uamu/masselec)                                      ! Atomic mass unit
      parameter(massprot=1.67262171d-27)                                ! Proton mass in kg
      parameter(pramu=massprot/uamu)                                    ! Rest mass proton in amu
      parameter(light=clight*bohrrad*masselec/hbar)                     ! Speed of light in au
      parameter(autohz=6.5796838999d15)                                 ! au to hz conversion factor
      parameter(factor=1.0d0/(light**4*4.0d0*amu**2*pramu**2)*autohz)   ! Spin-spin prefactor
      parameter(nfc=(4.0d0/3.0d0*gfactor*pi)**2*factor)                 ! FC-FC term prefactor
      parameter(nsd=(gfactor/2.0d0)**2*factor)                          ! SD-SD term prefactor
      parameter(nsdfc=4.0d0/6.0d0*(gfactor**2)*pi*factor)               ! SD-FC term prefactor
      parameter(npso=factor, ndso=factor)                               ! PSO-PSO and DSO-DSO term prefactor
c
      if (ga_nodeid().eq.0) write(luout,9999)
c
c     Current CPHF does not handle symmetry 
c     
      oskel = .false.
c
c     Find out from rtdb for which atom combinations we need to calculate the spin-spin couplings
c     Get number of atoms (all or number from rtdb)
c     Get which atoms (all or some read from rtdb)
c     Allocate arrays which will hold atomic information (k_zan and k_xyz)
c
      icom_atom = 0
      status = rtdb_parallel(.true.)
      if (.not.rtdb_get(rtdb, 'spinspin:pairs', MT_INT, 1,ipairs)) then 
         if (.not.geom_ncent(geom, icom_atom)) call                   
     &       errquit('hnd_spinspin: geom_ncent failed ',0, GEOM_ERR)
         ipairs = icom_atom*(icom_atom-1)
      endif
      if (.not. ma_push_get(mt_int,ipairs*2,'pairs',l_pr,k_pr)) 
     &   call errquit('hnd_spinspin: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_int,ipairs*2,'unique',l_un,k_un)) 
     &   call errquit('hnd_spinspin: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_int,ipairs*2,'transl',l_tr,k_tr)) 
     &   call errquit('hnd_spinspin: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,3*2*ipairs,'at coord',l_xyz,k_xyz)) 
     &   call errquit('hnd_spinspin: ma_push_get failed k_xyz',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*ipairs,'zan',l_zan,k_zan)) 
     &   call errquit('hnd_spinspin: ma_push_get failed k_zan',0,MA_ERR)
c
c     Try to read the atom list from rtdb when pair list is defined
c     Stored with all perturbing atoms first and all responding atoms second
c
      if (icom_atom.gt.0) then
         ixy = 0
         do i = 1, icom_atom
            do j = 1, icom_atom
               if (i.ne.j) then
                  int_mb(k_pr+ixy)        = i
                  int_mb(k_pr+ipairs+ixy) = j
                  ixy = ixy + 1
               endif
            enddo
         enddo
      else
         if (.not.rtdb_get(rtdb,'spinspin:pair list',mt_int,2*ipairs,
     &       int_mb(k_pr))) call errquit('hnd_spinspin: rtdb_get failed'
     &                                   ,0,RTDB_ERR)
      endif
c
c     Get unique list of perturbing and responding atoms
c     and setup lookup table between unique list and full pair list
c
      i_pert = 1
      int_mb(k_un) = int_mb(k_pr)
      int_mb(k_tr) = 1
      do i = 2, ipairs
         notthere = .true.
         do j = 1, i_pert
            if (int_mb(k_un+j-1).eq.int_mb(k_pr+i-1)) then
               int_mb(k_tr+i-1) = j
               notthere = .false.
            endif
         enddo
         if (notthere) then
            i_pert = i_pert + 1
            int_mb(k_tr+i-1) = i_pert
            int_mb(k_un+i_pert-1) = int_mb(k_pr+i-1)
         endif
      enddo
      i_resp = 1
      int_mb(k_un+i_pert) = int_mb(k_pr+ipairs)
      int_mb(k_tr+ipairs) = 1
      do i = ipairs+2, 2*ipairs
         notthere = .true.
         do j = i_pert+1, i_pert+i_resp
            if (int_mb(k_un+j-1).eq.int_mb(k_pr+i-1)) then
               int_mb(k_tr+i-1) = j - i_pert
               notthere = .false.
            endif
         enddo
         if (notthere) then
            i_resp = i_resp + 1
            int_mb(k_tr+i-1) = i_resp
            int_mb(k_un+i_pert+i_resp-1) = int_mb(k_pr+i-1)
         endif
      enddo
c
c     Get info for the all unique atoms
c
      do ixy = 0, i_pert+i_resp-1
         if (.not. geom_cent_get(geom, int_mb(k_un+ixy), tag, 
     &                           dbl_mb(k_xyz+3*ixy),dbl_mb(k_zan+ixy)))
     &    call errquit('hnd_spinspin: geom_cent_tag failed',0, GEOM_ERR)
      enddo 
c
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues
c
      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_spinspin: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('hnd_spinspin: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('hnd_spinspin: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c
c     Return if scftyp equals RHF
c
      if (scftyp.eq.'RHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          goto 911
      endif
c
c     UHF or UDFT but still closed shell
c
      if (nopen(1).ne.nopen(2)) then
          if (ga_nodeid().eq.0) write(luout,7001)
          goto 911
      endif
c
c     Make sure the CPHF knows that a DFT run is also UHF
c
      if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, scftyp))
     $    call errquit('hnd_spinspin: put of scftyp failed',0, RTDB_ERR)
c
c     Allocate array for tensor data (5 x # of pairs x 9, one for each term)
c
      if (.not. ma_push_get(mt_dbl,9*5*ipairs,'eigenval',l_data,k_data))
     &   call errquit('hnd_spinspin: ma_push_get fails k_data',0,MA_ERR)
      call dcopy(45*ipairs,0.0d0,0,dbl_mb(k_data),1)
c
c     Call routine to do FC and SD CPHF (FC-FC, SD-SD, and FC-SD)
c
      if (ga_nodeid().eq.0) write(luout,9998)
      call hnd_sdfc(rtdb,geom,basis,vectors,nclosed,nopen,nvirt,nbf,nmo,
     &              int_mb(k_pr),int_mb(k_tr),ipairs,int_mb(k_un),
     &              i_pert,i_resp,dbl_mb(k_data),dbl_mb(k_xyz),
     &              nfc,nsd,nsdfc)
c
c     Call routine to do PSO-PSO CPHF
c
      if (ga_nodeid().eq.0) write(luout,9997)
      call hnd_psoten(rtdb,geom,basis,vectors,nclosed,nopen,nvirt,nbf,
     &                nmo,int_mb(k_pr),int_mb(k_tr),ipairs,int_mb(k_un),
     &                i_pert,i_resp,dbl_mb(k_data),
     &                dbl_mb(k_xyz),npso)
c
c     Get Unperturbed Density Matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Calculate DSO term
c
      call hnd_dsoten(rtdb,geom,basis,g_dens,ndens,int_mb(k_pr),
     &                int_mb(k_tr),ipairs,int_mb(k_un),i_pert,i_resp,
     &                dbl_mb(k_data),dbl_mb(k_xyz),nbf,ndso)
c
c     Get the nuclear g-factors and scale the results with those.
c
      if (ga_nodeid().eq.0) then
      write(LuOut,9996)
c
c     Print out tensor information, and write to Ecce file if necessary
c
      call ecce_print_module_entry('spinspin')
c
      do ixy = 1, ipairs
         i = int_mb(k_pr+ixy-1)
         j = int_mb(k_pr+ixy+ipairs-1)
         i_off = k_data+45*(ixy-1)
         if (.not. geom_cent_get(geom, i, tagi, dbl_mb(k_xyz), znuc))
     &       call errquit('hnd_spinspin: geom_cent_get',0, GEOM_ERR)
         if (.not. geom_tag_to_element(tagi,symbola,element,atn)) call
     &       errquit('hnd_spinspin: geom_tag_to_element failed',
     &       0,GEOM_ERR)
         if (.not. atom_gfac(znuc,gfaci,isonra)) call
     &       errquit('hnd_spinspin: atom_gfac failed',0, UERR)
         if (.not. geom_cent_get(geom, j, tagj, dbl_mb(k_xyz), znuc))
     &       call errquit('hnd_spinspin: geom_cent_get',0, GEOM_ERR)
         if (.not. geom_tag_to_element(tagj,symbolb,element,atn)) call
     &       errquit('hnd_spinspin: geom_tag_to_element failed',
     &       0,GEOM_ERR)
         if (.not. atom_gfac(znuc,gfacj,isonrb)) call
     &       errquit('hnd_spinspin: atom_gfac failed',0, UERR)
         call dscal(45,gfaci*gfacj,dbl_mb(i_off),1)
c
c     Print tensor components and total tensor
c
         write(LuOut,9000) i,isonra,symbola,j,isonrb,symbolb,gfaci,gfacj
         isotr=(dbl_mb(i_off)+dbl_mb(i_off+4)+dbl_mb(i_off+8))
         isotr=isotr/3.0d0
         write(LuOut,9001) (dbl_mb(i_off   +ii),ii=0,8),isotr
         isotr=(dbl_mb(i_off+9)+dbl_mb(i_off+13)+dbl_mb(i_off+17))
         isotr=isotr/3.0d0
         write(LuOut,9002) (dbl_mb(i_off+ 9+ii),ii=0,8),isotr
         isotr=(dbl_mb(i_off+18)+dbl_mb(i_off+22)+dbl_mb(i_off+26))
         isotr=isotr/3.0d0
         write(LuOut,9003) (dbl_mb(i_off+18+ii),ii=0,8),isotr
         isotr=(dbl_mb(i_off+27)+dbl_mb(i_off+31)+dbl_mb(i_off+35))
         isotr=isotr/3.0d0
         write(LuOut,9004) (dbl_mb(i_off+27+ii),ii=0,8),isotr
         isotr=(dbl_mb(i_off+36)+dbl_mb(i_off+40)+dbl_mb(i_off+44))
         isotr=isotr/3.0d0
         write(LuOut,9005) (dbl_mb(i_off+36+ii),ii=0,8),isotr
c
c     Sum components into FC section
c
         do ii = 0, 8
            value = 0.0d0
            do jj = 0, 4
               value = value + dbl_mb(i_off+9*jj+ii)
            enddo
            dbl_mb(i_off+ii) = value
         enddo
         isotr = (dbl_mb(i_off)+dbl_mb(i_off+4)+dbl_mb(i_off+8))/3.0d0
         write(LuOut,9006) (dbl_mb(i_off+ii),ii=0,8)
         write(LuOut,9803) isotr
         write(LuOut,'(/)')
c
c     Print Ecce information
c
         call ecce_print1('index atom one',MT_INT,i,1)
         call ecce_print1_char('symbol atom one',symbola,1)
         call ecce_print1('isotope atom one',MT_INT,isonra,1)
         call ecce_print1('index atom two',MT_INT,j,1)
         call ecce_print1_char('symbol atom two',symbolb,1)
         call ecce_print1('isotope atom two',MT_INT,isonrb,1)
         call ecce_print2('spinspin tensor',MT_DBL,
     &                     dbl_mb(i_off),3,3,3)
         call ecce_print1('spinspin isotropic',MT_DBL,isotr,1)
c
      enddo
      endif
c
      call ecce_print_module_exit('spinspin','ok')
      call ga_sync()
c
c     Clean up all remaining memory
c
      do i = 1, ndens
          if (.not.ga_destroy(g_dens(i))) call 
     &       errquit('hnd_spinspin: ga_destroy failed gdens',i,GA_ERR)
      enddo
      if (.not.ga_destroy(vectors(2))) call 
     &    errquit('hnd_spinspin: ga_destroy failed vectors',0,GA_ERR)
 911  if (.not.ga_destroy(vectors(1))) call 
     &    errquit('hnd_spinspin: ga_destroy failed vectors',0,GA_ERR)
c
      if (.not.ma_chop_stack(l_pr)) call
     &    errquit('hnd_spinspin: ma_pop_stack failed l_tmp',0,MA_ERR)
c
      return
 7000 format(/,6x,'Indirect spin-spin coupling not calculated for RHF',
     1      /,6x,'wave functions: needs UHF wave function')
 7001 format(/,6x,'Indirect spin-spin coupling not calculated for open',
     1      /,6x,'shell wave functions: needs closed shell UHF!')
 9000 format(6x,49(1h-),/,
     &       6x,'Atom ',i4,': ',i3,'-',a2,' and Atom ',i4,': ',
     &       i3,'-',a2,//,
     &       6x,'Respective Nuclear g-factors: ',f8.6,' & ',f8.6/,
     &       6x,49(1h-),/)
 9001 format(6x,'Fermi Contact (FC) Term (Hz):',/,3(3F12.4,/),
     &       8x,'Isotropic = ',f12.4/)
 9002 format(6x,'Spin-Dipole   (SD) Term (Hz):',/,3(3F12.4,/),
     &       8x,'Isotropic = ',f12.4/)
 9003 format(6x,'FC - SD  Cross Term (Hz):',/,3(3F12.4,/),
     &       8x,'Isotropic = ',f12.4/)
 9004 format(6x,'Paramagnetic Spin-Orbit Term (Hz):',/,3(3F12.4,/),
     &       8x,'Isotropic = ',f12.4/)
 9005 format(6x,'Diamagnetic Spin-Orbit Term (Hz):',/,3(3F12.4,/),
     &       8x,'Isotropic = ',f12.4/)
 9006 format(6x,'Total Spin-Spin Coupling Tensor (Hz):',/,3(3F12.4,/),/)
 9700 format(6x,'Atom: ',i4,2x,a2)
 9803 format(6x,'Isotropic Spin-Spin Coupling = ',f12.4,' Hz',/)
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9986 format(10x,'Principal Components and Axis System',/,10x,
     1       3(7x,i1,4x))
 9987 format(10x,' isotropic = ',f12.4,/,
     1       10x,'anisotropy = ',f12.4,/)
 9999 format(
     1 /,6x,35(1h-),/,
     2 6x,'Indirect Spin-Spin Coupling Tensors',/,
     3 6x,35(1h-),/)
 9998 format(
     1 /,6x,51(1h-),/,
     2 6x,'Calculating Fermi-Contact and Spin-Dipole Responses',/,
     3 6x,51(1h-),/)
 9997 format(
     1 /,6x,46(1h-),/,
     2 6x,'Calculating Paramagnetic Spin-Orbit Responses',/,
     3 6x,46(1h-),/)
 9996 format(
     1 /,6x,34(1h-),/,
     2 6x,'Indirect Spin-Spin Tensors (in Hz)',/,
     3 6x,34(1h-),/)
      end
