      subroutine dim_grid_quadv0(rtdb, g_vdim)
c
      implicit none
      integer rtdb
      integer g_vdim(*)

      call dim_grid_quadv0_gen(rtdb, g_vdim, 1, 1, .false., 0, 1)

      return
      end
c
      subroutine dim_grid_quadv0_gen(rtdb, g_vdim,
     &      nmat, calc_type, l3d, ipm, imag)

      implicit none

#include "errquit.fh"
c
      integer rtdb
      integer g_vdim

      Logical oprint,oprintw
      Logical xc_chktau 
      external xc_rep_close, xc_chktau
      logical xc_rep_close
      integer nmat              ! [input] number of DIM matrices to make;
c                               !         only used for lhs cpks
      integer calc_type         ! what to calculate [input]
c                               !   1 = XC matrices for SCF
c                               !   2 = XC matrices for CPKS LHS
      logical l3d               ! true if 3d GA is used
      integer ipm               ! Control over dipoles
      integer imag              ! Control over dipoles
c
#include "cgridfile.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "dft2drv.fh"
#include "steric.fh"
#include "dimqm.fh"
c
      logical rtdb_status       
      logical grid_reopen
      external grid_reopen
      logical lsa ! Qin
c
      integer me, nproc, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2,  
     &        ncoef_max,
     &        npol
c
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat,
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags,
     &        ltau, itau,  lmmat, immat, lttau, ittau
c
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat

      double precision nelec
c
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, l_Mmat2, k_Mmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2, len_Mmat2
      integer len_tau, len_mmat, len_ttau
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision job_grid_acc,rho_err
      integer ii,isp
c
      logical fon ! fractional occupancy
      double precision pstrace  ! a measure of the total charge
      integer i3, idir
      integer l_dimxyz, k_dimxyz
      integer l_muind,  k_muind
      character*(1) direction(3)
      character*(1) dpm(2)
      character*(1) dri(2)
      data direction /'x', 'y', 'z'/
      data dpm /'+', '-'/
      data dri /'r', 'i'/

c
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = ldebug 
c      if(oprint) then
c        write(luout,*) "dim_grid_quadv0 start"
c        write(*,*) "GEOM:", geom
c      end if
c      oprint = .true.
      oprintw = util_print('quadrature warning', print_medium)
c
      irdelchi_atom = 0
c
      npol = (ipol*(ipol+1))/2
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
c
c     Open grid pts file
c     
cng
      if (grid_written .and. n_rec_in_file .le. 0 ) then
c
c       Somebody (fock_xc?) closed the grid file, better reopen it
c       before continuing. Note, fock_xc tends to open the grid file
c       at the beginning and close it at the end.
c
        if (.not.grid_reopen(geom)) call errquit(
     &    'dim_grid_quadv0: could not open grid file',0,DISK_ERR)
      endif
cng
      if (.not. grid_written .or. n_per_rec .le. 0 ) then
         call grid_packing_info
      end if
c
      if (oprint) then      
        write(*,*) " *** dim_grid_quadv0: calc_type:",calc_type 
        write(*,*) " *** dim_grid_quadv0: grid_written:",grid_written 
      end if
c
c     Set up logicals denoting calculation type
c
      do_scf = calc_type .eq. 1
      do_cpks_l = calc_type .eq. 2
      if (.not. (do_scf .or. do_cpks_l)) then
         call errquit('dim_grid_quadv0_gen: illegal calculation type',
     &    1,INPUT_ERR)
      endif
c
      if (oprint) then
        write(*,*) " *** dim_grid_quadv0: do_scf:",do_scf 
        write(*,*) " *** dim_grid_quadv0: n_per_rec:",n_per_rec
      end if
c
      len_tau = 1
      len_mmat = 1
      len_ttau = 1
      len_Amat2 = 1
      len_Cmat2 = 1
      len_Mmat2 = 1
c
      if (do_scf) then
         len_rho = npol*n_per_rec
      elseif (do_cpks_l) then
         len_rho = n_per_rec*(nmat+1)+n_per_rec
      endif
      if (.not. geom_ncent(geom, nctrs)) 
     &   call errquit('grid_den: geom_ncent failed',geom, GEOM_ERR)
      nxyz = 3*nctrs
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('grid_den: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'charge',lcharge,icharge))
     &   call errquit('grid_den: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nctrs*16,'tags',ltags,itags))
     &   call errquit('grid_den: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, nctrs, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('grid_den: geom_cart_get failed',74, GEOM_ERR)
c
c      write(*,*)' Coordinates from grid_den.F '
c      call output(dbl_mb(ixyz), 1, 3, 1, nctrs, 3, nctrs, 1)
c     &                    
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
c
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c      if (do_nmrks_r) nbf_ao_mxnbf_ce2=nbf_ao_mxnbf_ce2*3
c
c     Divide the multi-center integrand into a set of single center integrands.
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_den: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('grid_den: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit('grid_den: cannot allocate atom list',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit('grid_den: cannot allocate atom nbf',0, MA_ERR)
      nexpo=nbf_ao_mxprim
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_den: cannot allocate expo',0, MA_ERR)
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit('grid_den: cannot allocate coef',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit('grid_den: cannot allocate DENSF',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit('grid_den: cannot allocate DENS',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit('grid_den: cannot allocate rho',0, MA_ERR)
c
c
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit('grid_den: cannot allocate Bmat',0,
     &       MA_ERR)
c
c     Assume this numerical density routine can be called from anywhere.
c     -compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_ao)
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('grid_den: cannot pop stack',1, MA_ERR)
c
c     allocate scratch array for grid compression
c
      nibf=nbf_ao
      if (.not.MA_Push_Get(mt_int, nibf, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit('grid_quadv0: cannot allocate ibf_ao',2,
     &       MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3,
     &       MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0, MA_ERR)
      irqlen=nctrs*n_per_rec
      if(.not.grid_written) irqlen=max(irqlen,nctrs*nqmax)
      if (.not.MA_Push_get(mt_dbl,irqlen,'rq',lrq,irq))
     &   call errquit('grid_den: cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit('grid_quadv0: cannot allocate bp',4, MA_ERR)
c
c     compute variables for scratch space in function evaluation
c     assumes cartesian; spherical set would only require 2*MXang + 1
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
c     
c     Allocate scratch space.
c     
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 
c
c     Allocate and read in DIM coordinates
      i3 = 3*nDIM
      if(.not.ma_push_get(mt_dbl,i3,'dimqm xyz',l_dimxyz,k_dimxyz))
     $   call errquit('dimqm_main malloc k_xyz failed',1,MA_ERR)
c   Read in coords from RTDB
      if(.not. rtdb_get(rtdb,'dimpar:coords',mt_dbl,i3,
     $                  dbl_mb(k_dimxyz)))
     $   call errquit('dimqm:xyz rtdb get failed',1,RTDB_ERR)
c
c     Allocate and read in induced dipoles      
      if(.not.ma_push_get(mt_dbl, i3*nmat, 'dimqm:muind',
     $                     l_muind, k_muind))
     $  call errquit('quadv0 malloc muind failed',2112, MA_ERR)

      if(nmat > 1) then ! only 1 perturbation for ground state, 3 for response
        do idir = 1, nmat
          if(ipm .gt. 0) then ! FD Response (real and imaginary, +/-)
            if(.not.rtdb_get(dimqm_rtdb,'dimqm:muind_'//direction(idir)
     $                       //'_'//dri(imag)//dpm(ipm), mt_dbl, 3*nDIM,
     $                       dbl_mb(k_muind+i3*(idir-1))))
     $       call errquit('get perturbed +/- dipoles failed',6,RTDB_ERR)
          else ! Static response
            if(.not.rtdb_get(dimqm_rtdb,'dimqm:muind_'//direction(idir),
     $                       mt_dbl, 3*nDIM,
     $                       dbl_mb(k_muind+i3*(idir-1))))
     $       call errquit('get perturbed xyz dipoles failed',6,RTDB_ERR)
          end if
        end do
      else ! Ground state
        if(.not.rtdb_get(dimqm_rtdb, 'dimqm:muind', mt_dbl, 3*nDIM,
     $                   dbl_mb(k_muind)))
     $    call errquit('get perturbed dipoles failed', 2112, RTDB_ERR)
      endif

      call ga_sync()
c      
c     Grid quadrature
c
         if (iAOacc .le. 0) then
          if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))then
             iAOacc=-nint(log(e_conv))
          else
             iAOacc=max(iAOacc,-nint(log(e_conv)))
          endif
         end if
c
         call dim_grid_quadv0a(ncoef_max,nctrs, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho),  
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        int_mb(icetobfr),
     ,        int_mb(iibf_ao), dbl_mb(irqbuf),
     ,        dbl_mb(irchi_atom), g_vdim,
     ,        nmat, do_scf, do_cpks_l,
     ,        l3d, ipm, imag, dbl_mb(k_dimxyz),
     $        dbl_mb(k_muind))
c     
      
      if (.not.ma_chop_stack(lqxyz))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      call ga_sync()
      if(do_scf.and.(.not.l3d)) then
          call ga_symmetrize(g_vdim)
      endif
c
      if (.not.MA_Free_Heap(lxyz))
     &   call errquit('grid_den: cannot free heap',0, MA_ERR)
c
 123  format(
     &   1x,'Grid integrated density: ',f20.12,/
     &   1x,'Requested integration accuracy: ',e10.2)
c
C      if(oprint) then
C        write(luout,*) "dim_grid_quadv0 end"
C      end if

      return
      end
