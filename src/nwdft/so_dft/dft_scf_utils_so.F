c
c     Determine the HOMO-LUMO gap
c
      subroutine calc_homolumogap(k_eval,nelec,rlshift,homo,lumo,gap)
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer k_eval(2)      ! Pointers to the eigenvalue arrays. Only (1) is used
      integer nelec
      double precision rlshift
      double precision homo
      double precision lumo
      double precision gap
c
c     Determine the HOMO, LUMO and gap
c
      homo = Dbl_MB(k_eval(1)+nelec-1)   ! Extract the contents of the eigenvalue array
      lumo = Dbl_MB(k_eval(1)+nelec)
      gap = min(gap, (lumo-homo-rlshift))
c
      return
      end
c
c     Transform Fock To Orthonormal 
c
      subroutine trans_fock_to_ortho(g_tmp2,nbf_mo,g_sm12,g_fockso)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "consts.fh"
c
      integer g_tmp2          ! temp scratch array
      integer nbf_mo          ! number of molecular orbitals
      integer g_sm12          ! S^(-1/2)
      integer g_fockso(2)     ! 1: real, 2: imag
c
c     Real part of the Fock matrix
c
      call ga_zero(g_tmp2)
      call ga_dgemm('T', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_sm12, g_fockso(1), zero, g_tmp2)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_tmp2, g_sm12, zero, g_fockso(1))
c
c     Imag part of the Fock matrix
c
      call ga_zero(g_tmp2)
      call ga_dgemm('T', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_sm12, g_fockso(2), zero, g_tmp2)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_tmp2, g_sm12, zero, g_fockso(2))
c
      return
      end
c
c     Diagonalize complex Fock Matrix
c
      subroutine diag_fock(nbf_mo,ia,g_fockso,ibuff,g_moso,iwork,irwork,
     &     k_eval,trace,llwork,info)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer nbf_mo
      integer ia
      integer g_fockso(2)
      integer ibuff
      integer g_moso(2)
      integer iwork
      integer irwork
      integer k_eval(2)
      double precision trace
      integer llwork
      integer info
c
      integer i,j,i1
c
c     Prepare arrays for diagonalization
c
      trace = 0.d0
      do i = 1, nbf_mo
         do j = 1, nbf_mo
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=dcmplx(0.0, 0.0)
         enddo
      enddo
      do i = 1, nbf_mo
         call ga_get(g_fockso(1), 1,i, i,i, dbl_mb(ibuff),1)
         do j=1,i
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=
     =           dcmplx(dbl_mb(ibuff+j-1),0d0)
         enddo
         call ga_get(g_fockso(2), 1,i, i,i, dbl_mb(ibuff),1)
         do j=1,i
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=
     $               DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))
     $           +dcmplx(0d0,dbl_mb(ibuff+j-1))
         enddo
      enddo
      call ga_zero(g_moso(1))
      call ga_zero(g_moso(2))
c
c     Call the diagonalizer (complex diagonalizer)
c
      call zheev( 'V', 'U', nbf_mo, DCpl_mb(ia), nbf_mo,
     $            Dbl_mb(k_eval(1)),
     $            DCpl_mb(iwork), llwork, Dbl_mb(irwork), info )
      do i = 1, nbf_mo
         do j = 1, nbf_mo
            dbl_mb(ibuff+j-1)=0.0d0
            dbl_mb(ibuff+j-1)=dble(DCpl_mb(ia+nbf_mo*(i-1)+(j-1)))
         enddo
         i1=i
         call ga_put(g_moso(1),1,nbf_mo,i1,i1,dbl_mb(ibuff),1)
         trace = ddot(nbf_mo,dbl_mb(ibuff),1,dbl_mb(ibuff),1)
         do j = 1, nbf_mo
            dbl_mb(ibuff+j-1)=0.0d0
            dbl_mb(ibuff+j-1)=
     $             dimag(dcmplx(DCpl_mb(ia+nbf_mo*(i-1)+(j-1))))
         enddo
         i1=i
         call ga_put(g_moso(2),1,nbf_mo,i1,i1,dbl_mb(ibuff),1)
         trace = ddot(nbf_mo,dbl_mb(ibuff),1,dbl_mb(ibuff),1)
      enddo
c
      return
      end
c
c     Back-transform eigenvectors with S^-1/2.
c
      subroutine trans_vec_to_ao(nbf_mo,g_sm12,g_moso,g_fockso)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "consts.fh"
c
      integer nbf_mo
      integer g_sm12      ! contains S^-1/2
      integer g_moso(2)   ! MO vecs 1: real, 2: imag
      integer g_fockso(2) ! being used as scratch
c
c     Back-transform eigenvectors with S^-1/2.
c
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &     g_sm12, g_moso(1), zero, g_fockso(1))
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &     g_sm12, g_moso(2), zero, g_fockso(2))
c
c     Transfer into MO vec arrays
c
      call ga_zero(g_moso(1))
      call ga_zero(g_moso(2))
      call ga_copy(g_fockso(1), g_moso(1))
      call ga_copy(g_fockso(2), g_moso(2))
c
      return
      end
c
c     Calculate S powers: S,S^(-1/2),S^(1/2),S^(-1) 
c
      subroutine calc_s_powers(g_scr,g_tmp,nbf_ao,toll_s,svals,g_svecs,
     &   g_sp1,g_sm12,g_sp12,g_sm1)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "consts.fh"
c
      integer iw
      integer g_scr
      integer g_tmp
      integer nbf_ao
      double precision toll_s
      double precision svals(*)
      integer g_svecs
      integer g_sp1
      integer g_sm12
      integer g_sp12
      integer g_sm1
c
c     Calculate S^1: g_sp1
c
      iw = 1
      call ga_zero(g_scr)
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr,
     &     g_tmp, nbf_ao, iw)
      call ga_zero(g_sp1)
      call ga_fock_sf(g_scr, g_sp1, nbf_ao) ! Map onto larger array
c
c     Calculate S^-1/2: g_sm12
c
      iw = 2
      call ga_zero(g_scr)
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr,
     &     g_tmp, nbf_ao, iw)
      call ga_zero(g_sm12)
      call ga_fock_sf(g_scr, g_sm12, nbf_ao) ! Map onto larger array
c
c     Calculate S^+1/2: g_sp12
c
      iw = 3
      call ga_zero(g_scr)
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr,
     &     g_tmp, nbf_ao, iw)
      call ga_zero(g_sp12)
      call ga_fock_sf(g_scr, g_sp12, nbf_ao) ! Map onto larger array
c
c     Calculate S^-1: g_sm1
c
      iw = 4
      call ga_zero(g_scr)
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr,
     &     g_tmp, nbf_ao, iw)
      call ga_zero(g_sm1)
      call ga_fock_sf(g_scr, g_sm1, nbf_ao) ! Map onto larger array
c
      return
      end
c
c     Print energies 
c
      subroutine print_energies(etnew,enuc,ecore,ecoul,exc,
     &   nexc,rho_n,dft_time)
c
      implicit none
c
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
c
      double precision etnew
      double precision enuc
      double precision ecore
      double precision ecoul
      double precision exc(2)
      double precision nexc
      double precision rho_n
      double precision dft_time
c
c     Print  
c
      if (nexc.le.1) then
       write(luout,222)etnew+enuc,ecore,ecoul,exc(1),enuc
      else
       write(luout,223)etnew+enuc,ecore,ecoul,exc(1),exc(2),enuc
      end if
      write(luout,2222) rho_n
      write(luout,2223) dft_time
c
 222  format(//
     &     '      Total SO-DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '    Exchange-Corr. energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 223  format(//
     &     '      Total SO-DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '          Exchange energy =', f22.12/
     &     '       Correlation energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 2222 format(' Numeric. integr. density =', f22.12/)
 2223 format('     Total iterative time =', f9.1,'s'//)
c
      return
      end
c
c     Level shifting is implemented here (similarity 
c     transformation before standard eigensolver).  Note,
c     levelshifting is appropriate once a transformation
c     is available which makes the resulting Fock matrix 
c     diagonally dominant, e.g., in an approximate MO basis.  
c     Also note, there are many matrix multiplies with S^+-1/2 
c     which are redundant if one is sure that the former basis
c     is orthonormal.
c
      subroutine levelshift_fock(nbf_mo,ntotocc,g_tmp2,g_sp12,g_moso,
     &   g_scr2,g_fockso,rlshift)
c
      implicit none
c
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "consts.fh"
c
         integer nbf_mo
         integer ntotocc
         integer g_tmp2
         integer g_sp12
         integer g_moso(2)
         integer g_scr2
         integer g_fockso(2)
         double precision rlshift
c
         integer j
         integer me
         integer nproc
c
c        Preliminaries
c
         me = ga_nodeid()
         nproc = ga_nnodes()
c
c        Build a matrix which is diagonal in the "MO" rep,
c        back-transform, and shift the current Fock matrix
c     
c        Use S^+1/2 (g_sp12)  * old movecs (as a transform).
c     
c        Real part
         call ga_zero(g_tmp2)
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_sp12, g_moso(1), zero, g_tmp2)
         call ga_copy(g_tmp2,  g_moso(1))

c        Imag part
         call ga_zero(g_tmp2)
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_sp12, g_moso(2), zero, g_tmp2)
         call ga_copy(g_tmp2,  g_moso(2))
c     
c        Build diagonal matrix for the shift
         call ga_zero(g_tmp2)
         do j = nTotOcc+1+me, nbf_mo, nproc
            call ga_put(g_tmp2, j, j, j, j, rlshift, 1)
         enddo
c     
c        Transform this into "AO" basis and add to current 
c        Fock matrix
c     
c        Real part
         call ga_zero(g_scr2)  ! used as a work area
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_moso(1), g_tmp2, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_scr2, g_moso(1), one, g_fockso(1))
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_moso(2), g_tmp2, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_scr2, g_moso(2), one, g_fockso(1))
c
c        Imag part
         call ga_zero(g_scr2)  ! used as a work area
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_moso(1), g_tmp2, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, mone,
     &                 g_scr2, g_moso(2), one, g_fockso(2))
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_moso(2), g_tmp2, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one,
     &                 g_scr2, g_moso(1), one, g_fockso(2))
c
      return
      end
c $Id$
