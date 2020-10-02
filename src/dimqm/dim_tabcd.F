c
C$Id: xc_tabcd.F 20247 2011-04-28 18:58:49Z d3y133 $
c
c     Modified to handle CPKS while reusing existing code
c
c     BGJ - 8/98
c
      Subroutine dim_tabcd(what, l3d, Tmat, TTmat, Amat, Bmat, chi,
     $                     nmat, nq, nbf, max_at_bf, max_at_bf2, basis,
     $                     natoms, iniz, ifin, g_vdim, ibf, cetobfr)
c
c     We're using xc_tabcd for CPKS purposes
c
      implicit none
#include "mafdecls.fh"
#include "dftpara.fh"
#include "dft2drv.fh"
#include "dist.fh"
#include "dft_fdist.fh"
#include "stdio.fh"
#include "nwc_const.fh"
c
c     Input Variables
      integer what ! Control for SCF vs CPKS 
c                  ! 0 = SCF 
c                  ! 1 = CPKS_LHS
      logical l3d  ! Using 3d arrays?
      integer nq   ! Number of quadrature points
      integer nbf  ! Number of basis functions
      integer max_at_bf         ! Unsure
      double precision Tmat(*)  ! Work array
      double precision TTmat(*) ! Work array
      double precision Amat(nq,*) ! DIM potential calculated at each quadrature point
      double precision Bmat(nq,max_at_bf) ! Work array
      double precision chi(nq,nbf) ! Basis functions evaluated at each quadrature point
      integer nmat ! Number of perturbations
                   ! 1 for SCF
                   ! 3 for CPKS
      integer max_at_bf2 ! Unsure
      integer basis ! Basis set handle
      integer natoms ! Number of atoms in QM system
      integer iniz(natoms) ! Unsure
      integer ifin(natoms) ! Unsure
      integer g_vdim(*)   ! DIM potentials in Fock matrix form [IN/OUT]
      integer ibf(nbf) ! Unsure
      integer cetobfr(2,natoms) ! Unsure
c
c     Local Variables
      integer imat 
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast
      integer mu, mu1
      integer n, lastjat
      integer jrsh, jrsh2, idir, jdir
      integer sizeblk
      integer nonzero(nw_max_atom), natleft, iat0, jat0
      logical dofull
c
c       scf: l3d=.f.    & n3d=1
c       lhs: l3d=.true. & n3d=1
c
      natleft=0
      do  iat = 1, natoms
         if(iniz(iat).ne.0) then
            natleft=natleft+1
            nonzero(natleft)=iat
         endif
      enddo
      dofull=what.ne.0.or.l3d
c
      do 500 imat = 1,nmat ! Loop over perturbations
c
c     Compute the matrix product for the DIM potential:
c
c              T = Transpose(C * A) * C
c
c       repl stuff
c
        do 430 iat0=1,natleft ! Loop over iatoms
          iat=nonzero(iat0)
          inizia = iniz(iat)
          ifinia = ifin(iat)
          ifirst = cetobfr(1,iat)
          ilast = cetobfr(2,iat)
          nnia = ifinia - inizia + 1
          nbfia = ilast - ifirst + 1
          do mu = 1, nnia ! Loop over nnia
            mu1 = mu+inizia-1
            do n = 1, nq ! Loop over quadrature points
              Bmat(n,mu) = chi(n,mu1)*Amat(n,imat)
            enddo ! End loop over quadrature points
          enddo ! End loop over nnia
          lastjat=iat0
          do 168 jat0=1,lastjat ! Loop over jatoms
            jat=nonzero(jat0)
            inizja = iniz(jat)
            ifinja = ifin(jat)
            jfirst = cetobfr(1,jat)
            jlast = cetobfr(2,jat)
            nbfja = jlast - jfirst + 1
            nnja = ifinja - inizja + 1
            sizeblk=nbfia*nbfja
            call dgemm('T', 'N', nnia, nnja, nq, 1.0d0, Bmat,
     &                 nq, chi(1,inizja), nq, 0.d0, Tmat, nnia)
            call dfill(max_at_bf2, 0.0d0, TTmat, 1)
            if(dofull) then
              call scat_mat(TTmat, Tmat, nbfia, nbfja, nnia,
     &                    nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
            else
              call scat_matup(TTmat, Tmat, nbfia, nbfja, nnia,
     &                    nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
            endif
            if(what.eq.0) then
              jrsh=1
            elseif(what.eq.1) then
              jrsh=imat
            endif
            jrsh2=jrsh
            call updist(monitor_size_ga_acc1, sizeblk)
            if(l3d) then
              call dft_3dacc(g_vdim, TTmat, jrsh, jrsh2, ifirst,
     $                        ilast, jfirst, jlast, nbfia)
            else
              call upd_atom_block(g_vdim, basis, iat,
     $                            jat, TTmat)
            endif
            if(dofull) then ! If dofull
              if (iat .ne. jat) then ! iat =/= jat
                call dfill(max_at_bf2, 0.0d0, Tmat, 1)
                call transp_mat(TTmat, Tmat, nbfia, nbfja)
                if(l3d) then
                  call dft_3dacc(g_vdim, Tmat, jrsh, jrsh2,
     $                   jfirst, jlast, ifirst, ilast, nbfja)
                else
                  call upd_atom_block(g_vdim, basis, jat, iat, Tmat)
                endif
              endif ! End iat =/= jat
            endif ! If dofull
  168     continue ! End loop over jatoms
  430   continue ! End loop over natleft
  500 continue ! End loop over perturbations
      return
      end subroutine dim_tabcd
