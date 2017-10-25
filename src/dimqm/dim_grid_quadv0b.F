      subroutine dim_grid_quadv0b(ictr_buf, qxyz, qwght, nqpts, rad,
     $                            ncoef_max, natoms, iniz, ifin, rho,
     $                            xyz, expo, Bmat, Fmat, Pmat, rq,
     $                            cetobfr, ibf_ao, bas_cent_info,
     $                            bas_cset_info, docset, rchi_atom,
     $                            g_vdim, nmat, do_scf, do_cpks_l, l3d,
     $                            ipm, imag, dimxyz, muind)
c
c$Id: dim_grid_quadv0b.F 19821 2010-12-14 07:46:49Z d3y133 $
c
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_sic.fh"
#include "dftps.fh"
#include "geom.fh"
#include "dimqm.fh"
c
c   Input Variables
      integer ictr_buf ! ctr of grid
      double precision qxyz(3, n_per_rec) ! Quadrature point coordinates
      double precision qwght(n_per_rec)   ! Quadrature point weights
      integer nqpts                       ! Number of quadrature points
      double precision rad
      integer ncoef_max
      integer natoms    ! Number of atoms
      integer iniz(natoms)
      integer ifin(natoms)
      double precision rho(*)
      double precision xyz(3, natoms)   ! Atom coordinates
      double precision expo(nbf_ao_mxprim)
      double precision Bmat(nbf_ao_mxnbf_ce*n_per_rec)
      double precision Fmat(*)
      double precision Pmat(*)
      double precision rq(n_per_rec, natoms)
      integer cetobfr(2, natoms)
      integer ibf_ao(nbf_ao)
      integer bas_cent_info(3,natoms)
      integer bas_cset_info(6,*)
      logical docset(*)
      double precision rchi_atom(natoms)
      integer g_vdim ! DIM potential for this block of the quadrature
      integer nmat ! Number of density matrices (number of perturbing directions)
      logical do_scf
      logical do_cpks_l
      logical l3d
      integer ipm
      integer imag
      double precision dimxyz(3, nDIM) ! DIM Coordinates
      double precision muind(3, nDIM, nmat) ! DIM Coordinates
c      
c   Local Variables 
      integer nbf_ao_mxnbf_ce2, mbf_ao, npol, ncontrset, maxdbas
      integer lchi_ao, ichi_ao, ldelchi_ao, idelchi_ao
      integer ldmat,idmat,i, dim_grid_nbfm
      double precision dabsmax
      double precision acc_ao_gauss
      external dabsmax
      external dim_grid_nbfm
      integer iscf_rho,iscf_delrho
      integer iscf_tau,iscf_ttau
      logical grid_sicinit,grid_sicend
      external grid_sicinit,grid_sicend
      double precision acc_xc_gauss
      integer lemat,iemat,lfmat,ifmat,k_scr,l_scr
      logical do_2nd
      logical stat
      logical dprint
c
      npol = (ipol*(ipol+1))/2
      acc_ao_gauss = dble(iaoacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(nqpts.eq.0) return
      dprint = ldebug
c      if(dprint) then
c        write(luout,*) "dim_grid_quadv0b start"
c      end if
       maxdbas = 0
       idelchi_ao = 0
c
c     Evaluate the AO basis set at each of the quad. points.
c     allocate arrays for exponents and contraction coefficients
c     within int_eval_gbsets
c     Evaluate the distances (squared) between each center and the points
c
      mbf_ao=nbf_ao
      if(natoms.gt.1) then
         call icopy(mbf_ao, 0,0, ibf_ao, 1)  

         mbf_ao=dim_grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao, docset, iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss,iatype_pt_chg)
         if (mbf_ao.eq.0) return
      else
        if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
        iniz(1)=1
        ifin(1)=nbf_ao
        do i=1,ncontrset
          docset(i)=.true.
        enddo
        do i=1,nbf_ao
          ibf_ao(i)=i
        enddo
      endif

      if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &     lchi_ao, ichi_ao))
     &     call errquit('grid_quadv0b: cannot allocate chi_ao',0,
     &       MA_ERR)

      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)

      call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &     0d0, 0d0, 0d0, rq, qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)
c
c     get reduced Xi(r) and dXi(r) over atoms
c
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(ichi_ao),rchi_atom)

      if(nqpts.eq.0) goto 2010

c
c     Evaluate the DIM potential
c     Set up pointer to the SCF density for the CPKS LHS case.
c
      call dim_eval_fnl(rho, nqpts, qwght,
     &     dimxyz, qxyz, 
     $     nmat, ipm, imag, muind)
c       rho(1:nqpts*nmat) = 0.0d0
c
c   Transform the potential from over the quadrature points to the shape of the Fock Matrix
      if(do_scf) then
c
        call dim_tabcd(0, l3d, Fmat, Pmat, rho, Bmat,
     $                   dbl_mb(ichi_ao), 1, nqpts, mbf_ao,
     $                   nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     $                   AO_bas_han, natoms, iniz, ifin, g_vdim,
     $                   ibf_ao, cetobfr)
c
      else if (do_cpks_l) then
c
        call dim_tabcd(1, l3d, Fmat, Pmat, rho, Bmat,
     $                   dbl_mb(ichi_ao), nmat, nqpts, mbf_ao,
     $                   nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     $                   AO_bas_han, natoms, iniz, ifin, g_vdim,
     $                   ibf_ao, cetobfr)
c
      endif
c         write(luout,*) "g_dim AFTER dim_tabcd"
c         call ga_print(g_vxc)
c
c      endif

 2010 continue

      if (sic_orb_index.eq.1) then
         if(.not.grid_sicend(l_vect1,ldelrho_sig))
     ,        call errquit(' grid_quadv0b: grid_sicend failed',0,
     &       CALC_ERR)
      endif

      if (.not.ma_pop_stack(lchi_ao))
     &     call errquit('grid_quadv0b: cannot pop stack', 3, MA_ERR)

c      if(dprint) then
c        write(luout,*) "end dim_grid_quadv0b"
c      end if
      return
      end subroutine dim_grid_quadv0b
c
c     function dim_grid_nbfm
c
      integer function dim_grid_nbfm(basis_hand,   mcenters, 
     &     ctr_quad,  rad_quad, xyz, qxyz,nq,
     &                  ibf, docset, iniz, ifin, zprim,
     .     minexp,ldiff,acc_gauss,iatype_pt_chg)
c
C$Id: grid_quadv0b.F 19821 2010-12-14 07:46:49Z d3y133 $
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand
      integer mcenters ! [input]
      double precision acc_gauss ! [input]
      double precision xyz(3,*) ! [input]
      double precision rad_quad ! something about radius of this integration shell [input]
      integer nq ! [in] # grid pts
      double precision qxyz(3,nq) ! [input] coord grid pts
      integer ctr_quad ! grid center  [input]
      integer ibf(*) ! [output]
      logical docset(*) ! [output]
      integer ldiff(*) ! [in]
      double precision minexp(*) ! [in]
      logical iatype_pt_chg(*) ! [in]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision r_q0 ! min dist
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(*)
      integer n1, icset, ictr,  nprim, ncontr, 
     &        isphere,  l, iprimo,npt,nshbf
      double precision zmin,acc_loc,x,y,z
      integer n,ic1,ic2,m,ibfhi,ibflo
      double precision alpha,logeps
      double precision gaussian_range,r2,r_arg
      logical veryfar
c
      gaussian_range(n,logeps,alpha) = 
     $    (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $    sqrt(-16.0d0*alpha*logeps)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
      if(acc_gauss.gt.25d0.or.acc_gauss.lt.3d0) call errquit(
     ' ' grid_nbfm: silly accgauss ',nint(acc_gauss), UNKNOWN_ERR)
      acc_loc=-acc_gauss
c
      n1 = 0
      npt=0
      do 400 ictr=1,mcenters ! Loop over all atoms
        if(iatype_pt_chg(ictr)) goto 400 ! Cycle if this atom is a point charge
        if (.not.bas_ce2cnr(basis_hand,ictr,ic1,ic2))
     &      call errquit('Exiting in xc_signf.',11, BASIS_ERR)

        r2=0d0
        if (ictr.ne.ctr_quad) then ! Calculate r-squared if this atom isn't the center of the quadrature
          x = xyz(1,ctr_quad) - xyz(1,ictr)
          y = xyz(2,ctr_quad) - xyz(2,ictr)
          z = xyz(3,ctr_quad) - xyz(3,ictr)
          r2 = sqrt(x*x + y*y + z*z)
        endif
        r_arg=0d0
        if (rad_quad.lt.r2) r_arg = (r2-rad_quad)
c
c     check on most diffuse fn
c
        veryfar=r_arg.gt.gaussian_range(ldiff(ictr),acc_loc,
     .                                  minexp(ictr)) 
        if(veryfar) then
          if(.not.bas_ce2bfr(basis_hand, ictr, ibflo, ibfhi))
     &        call errquit('Exiting in grid_nbfm',4, BASIS_ERR)
          nshbf=ibfhi-ibflo+1
          npt=npt+nshbf
        else
          r_q0=1d10
          do n=1,nq
            x = qxyz(1,n) - xyz(1,ictr)
            y = qxyz(2,n) - xyz(2,ictr)
            z = qxyz(3,n) - xyz(3,ictr)
            r_q0 = min(r_q0,sqrt(x*x + y*y + z*z))
          enddo
          do icset = ic1,ic2
            docset(icset) = .false.
c         
c       get info about current contraction set
            if (.not. bas_continfo(basis_hand, icset,  l ,nprim,
     &                             ncontr, isphere))
     &          call errquit('Exiting in grid_nbfm.',5, BASIS_ERR)
c       
c       get exponents and contraction coefficients for this contraction set
            if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &          call errquit('Exiting in grid_nbfm.',7, BASIS_ERR)
c       
c       Determine the minimum Gaussian exponent.
            zmin = 1.D+06
            do iprimo = 1,nprim
              zmin = min(zprim(iprimo),zmin)
            enddo
c       
c       Only include those basis functions that are "non-zero" for at least 
c       one  point in the sampling set.
            nshbf=ncontr*(((l+1)*(l+2))/2)
            if(isphere.eq.1) then
              nshbf=ncontr*(2*l+1)
            endif
c
c     pre-screening for worst case (max radius)
c
            if (r_q0.lt.gaussian_range(l,acc_loc,zmin)) then
              do m=1,nshbf
                ibf(n1+m) = npt+m
              enddo
              docset(icset) = .true.
              if (iniz(ictr).eq.0) iniz(ictr)=n1+1
              n1=n1+nshbf
            endif
            npt=npt+nshbf
          enddo ! icset = ic1, ic2
          ifin(ictr)= n1
        endif ! If veryfar
  400 continue ! End loop over atoms
c
      dim_grid_nbfm = n1
c
      return
      end
c
c     Precompute relevant basis info for XC
c
c     BGJ - 9/00
c
      Subroutine dim_make_basis_info(basis_hand, bas_cent_info,
     &     bas_cset_info, mcenters)
c
C$Id: grid_quadv0b.F 19821 2010-12-14 07:46:49Z d3y133 $
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand        ! [input]
      integer mcenters          ! [input]
      integer bas_cent_info(3,mcenters) ! [output]
      integer bas_cset_info(6,*) ! [output]
c
#include "bas.fh"
c
      integer ictr, icset, ifirst, ilast, nprim, ncontr, l, isphere
      integer ic1, ic2
c
      do ictr = 1,mcenters
         bas_cent_info(1,ictr) = 0 ! max angular momentum
         if (.not.bas_ce2cnr(basis_hand,ictr,ic1,ic2))
     &        call errquit('Exiting in xc_make_basis_info',1, BASIS_ERR)
         bas_cent_info(2,ictr) = ic1
         bas_cent_info(3,ictr) = ic2
c
         do icset = ic1,ic2
c         
c     Get info about current contraction set - first and last basis function,
c     angular momentum, number of primitives, number of contractions and
c     whether spherical harmomic
c       
            if (.not. bas_cn2bfr(basis_hand, icset, ifirst, ilast))
     &           call errquit('Exiting in xc_make_basis_info',1,
     &       BASIS_ERR)
            if (.not. bas_continfo(basis_hand, icset, l, nprim,
     &           ncontr, isphere))
     &           call errquit('Exiting in xc_make_basis_info',2,
     &       BASIS_ERR)
c
            if (l.lt.0) then
               call errquit('L code < 0 in xc_make_basis_info',1,
     &       BASIS_ERR)
            endif
c
            bas_cent_info(1,ictr) = max(bas_cent_info(1,ictr),l)
c
            bas_cset_info(1,icset) = ifirst
            bas_cset_info(2,icset) = ilast
            bas_cset_info(3,icset) = l
            bas_cset_info(4,icset) = nprim
            bas_cset_info(5,icset) = ncontr
            bas_cset_info(6,icset) = isphere
         enddo
      enddo
c
      return
      end
