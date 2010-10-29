      subroutine multipole_density(basis, center, lmax, 
     &     g_dens, moments, len_moments)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "cscfps.fh"
c     
c     Add into in moments() the cartesian multipole moments computed
c     as an expectation value using the given density.  A minus sign
c     corresponding to the sign of the electron charge is inserted.
c     
c     moments(q) -= <u|q|v> Duv , q = 0, ..., lmax
c     
c     L=0  charge
c     L=1  dipole        x, y, z
c     L=2  quadrupole    xx xy xz yy yz zz  
c     L=3  octopole      xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
c     ...
c     Run NWints/int/showxyz for more details
c     
c     There are (L+1)(L+2)/2 cartesian moments of order L
c     Sum(L=0,M) (L+1)(L+2)/2  =  (M+1)(M+2)(M+3)/6
c     
      integer basis             ! [input] basis 
      double precision center(3) ! [input] the expansion center
      integer g_dens            ! [input] 
      integer lmax              ! [input] up to 2^{lmax} multipoles
      integer len_moments       ! [input] size of moments
      double precision moments(len_moments) ! [output] multipoles
c     
      double precision one, two
      parameter (one=1.d0, two=2.d0)
c     
      integer geom, moff
      integer nshell, noperators, maxang, maxscr, me, nproc
      integer nbf_max, lmpmax, maxd, ishell, jshell, ijshell
      integer ilo, ihi, jlo, jhi, idim, jdim, ind, i, j, l, m, ioff
      integer l_scr, k_scr, l_mp, k_mp, l_dens, k_dens, l_m, k_m
c
      if (oscfps) call pstat_on(ps_multip)
c     
      noperators = (lmax+1)*(lmax+2)*(lmax+3)/6
      if (noperators .gt. len_moments)
     $     call errquit('multiple: moments is too short', len_moments,
     &       UNKNOWN_ERR)
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('multiplole: bad basis', 0, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('multipole: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_high_angular(basis,maxang)) call errquit(
     &     'multipole: bas_high_angular call failed',15, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit(
     &     'multipole: bas_nbf_cn_max failed',20, BASIS_ERR)
c     
c     
c     length of int_mpole integral output for full square list
c     includes l_pole = 0,...,lmax, where l_pole = 0 is simply
c     the 2-c overlap matrix.  (cartesian or spherical components).
c     
      maxd = nbf_max**2
      call int_mem_dipole(lmpmax,maxscr,basis,basis,lmax)
      maxscr = max(100000,maxscr)
c     
c     allocate necessary local temporary arrays on the stack
c     
      if(.not. ma_push_get(mt_dbl, lmpmax, 'mult:mp', l_mp, k_mp))
     &     call errquit('multipole: cannot allocate mp', lmpmax, MA_ERR)
      if(.not. ma_push_get(mt_dbl, maxscr, 'mult:scr', l_scr, k_scr))
     &     call errquit('multipole: cannot allocate scratch', maxscr,
     &       MA_ERR)
      if(.not. ma_push_get(mt_dbl, maxd, 'mult:dens', l_dens, k_dens))
     &     call errquit('multipole: cannot allocate dens', maxd, MA_ERR)
      if(.not. ma_push_get(mt_dbl, noperators, 'mult:mom', l_m, k_m))
     &     call errquit('multipole: cannot allocate mom', noperators,
     &       MA_ERR)
c     
      call dfill(noperators, 0.0d0, dbl_mb(k_m), 1)
      ijshell = -1
      me = ga_nodeid()
      nproc = ga_nnodes()
      do ishell = 1, nshell
c     
c     get basis info
c     
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &        ('multipole: bas_cn2bfr failed for basis', basis,
     &       BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c     
c     get basis info
c     
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('multipole: bas_cn2bfr', basis,
     &       BASIS_ERR)
               jdim = jhi - jlo + 1
c     
c     Get the density patch, make the integrals and contract
c     (note sign change when adding in since electrons are negative)
c     
               call ga_get(g_dens, jlo, jhi, ilo, ihi,
     $              dbl_mb(k_dens), jdim)
c     
               call int_mpole(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              lmpmax, dbl_mb(k_mp))
c     
c     output from int_mpole is: overlap, dipole, q-pole, ...
c     within a multipole block, the order is <i|m|j>  j fastest, 
c     then m, then i
c     
               ind = k_mp
               do l = 0, lmax
                  do i = 1, idim
                     ioff = (i-1)*jdim + k_dens - 1
                     do m = 1, ((l+1)*(l+2))/2
                        moff = (l*(l+1)*(l+2))/6 + k_m - 1 + m
                        do j = 1, jdim
                           dbl_mb(moff) = dbl_mb(moff) -
     $                          dbl_mb(j+ioff)*dbl_mb(ind)
                           ind = ind + 1
                        end do
                     enddo
                  end do
               end do
            end if
         end do
      end do
c     
c     Global sum of locate array and then add into result
c     
      call ga_sync()
      call ga_dgop(331, dbl_mb(k_m), noperators, '+')
      do i = 1, noperators
         moments(i) = moments(i) + dbl_mb(k_m+i-1)
      end do
c     
c     clean up stack
c     
      if (.not. ma_pop_stack(l_m)) call errquit('multipole: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_dens)) call errquit('multipole: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_scr)) call errquit('multipole: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_mp)) call errquit('multipole: ma?',0,
     &       MA_ERR)
c
      if (oscfps) call pstat_off(ps_multip)
c     
      end
      subroutine geom_nuc_mpole(geom, center, lmax,
     &     moments, len_moments)
      implicit none
#include "errquit.fh"
#include "geom.fh"
c     
c     Add the nuclear contribution to cartesian multipole moments
c     into the array moments().
c     
c     arguments
c     
      integer geom              ! [input] geometry handle
      double precision center(3) ! [input] Multipole center
      integer lmax              ! [input] Max order of multipole
      integer len_moments
      double precision moments(len_moments)
c     
c     locals
c     
      character*16 tag
      double precision coord(3), q
      integer ncenter, icent, ltotal, k, l, m, ind, atn
      double precision xk, xl, xm
      logical bqbq, isatom
      character*2 symbol
      character*20 element
      double precision pow_multi
      external pow_multi
c     
      if (.not. geom_ncent(geom, ncenter)) 
     $     call errquit('geom_nuc_mpole: bad geometry',0, GEOM_ERR)
c
      bqbq = geom_include_bqbq(geom)
c     
c     loop over centers, summing multipole moments
c     
      do icent = 1,ncenter
         if (.not. geom_cent_get(geom,icent,tag,coord,q))
     &        call errquit('geom_nuc_mpole:  geom_cent_get failed',10,
     &       GEOM_ERR)
         isatom = geom_tag_to_element(tag, symbol, element, atn)
         if (isatom .or. bqbq) then
c     
c     x^k y^l z^m, k+l+m=ltotal
c     
            ind = 0
            do ltotal = 0, lmax
               do k = ltotal, 0, -1
                  do l = (ltotal-k), 0, -1
                     m = ltotal - k - l
                     ind = ind + 1
c
c     On Fujitsu 0**0 is not 1
c
                     moments(ind) = moments(ind) + q*
     $                    pow_multi((coord(1)-center(1)),k) *
     $                    pow_multi((coord(2)-center(2)),l) * 
     $                    pow_multi((coord(3)-center(3)),m)
                  end do
               end do
            end do
         endif
      end do
c
      end 
      double precision function pow_multi(a,b)
      double precision a
      integer b
c
c     On Fujitsu 0**0 is not 1
c
      if (b .eq. 0) then
         pow_multi = 1.0d0
      else
         pow_multi = a**b
      endif
      end
