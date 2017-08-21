      subroutine tddft_dipole(basis, center, lmax, 
     &     g_dens, moments, len_moments, k_evl, ipolar,
     1     g_movecs)
*
*  Return GA with dipole integrals
*  Adapted multipole.F
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
cjb
      logical nodezero
      integer idx, ipos, size, unitn
      integer irec, icount
      integer nbf,l_dip,k_dip
      integer k_evl, ipolar

      integer g_movecs,reclen

      nodezero=(ga_nodeid().eq.0)
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
      if (.not. bas_numbf(basis,nbf)) call errquit(
     &     'multipole: bas_nbf_cn_max failed',20, BASIS_ERR)
c      if(nodezero) write(6,*)'NBF:',nbf
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
               call int_mpole(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              lmpmax, dbl_mb(k_mp))
c     
c     output from int_mpole is: overlap, dipole, q-pole, ...
c     within a multipole block, the order is <i|m|j>  j fastest, 
c     then m, then i
c    

      if(.not. ma_push_get(mt_dbl, idim*jdim, 'mult:mp', l_dip, k_dip))
     &     call errquit('multipole: cannot allocate mp', lmpmax, MA_ERR)

      call dfill(idim*jdim, 0.0d0, dbl_mb(k_dip), 1)

               ind = k_mp
               do l = 0, lmax
                  do i = 1, idim
                     ioff = (i-1)*jdim + k_dens - 1
                     do m = 1, ((l+1)*(l+2))/2
                        moff = (l*(l+1)*(l+2))/6 + k_m - 1 + m
                        do j = 1, jdim
c                           dbl_mb(moff) = dbl_mb(moff) -
c     $                          dbl_mb(j+ioff)*dbl_mb(ind)
c          if(nodezero.and.(l.eq.1)) then
c        write(6,"('<i|m|j>:',I5,I2,I5,'=',F18.12)")i+ilo-1,m,j+jlo-1,
c     1    dbl_mb(ind)
c          endif
          if((l.eq.1).and.(m.eq.ipolar)) then
                 dbl_mb(k_dip+j-1+(i-1)*jdim) = 
     2 dbl_mb(k_dip+j-1+(i-1)*jdim) + dbl_mb(ind)
          endif
                           ind = ind + 1
                        end do
                     enddo
                  end do
               end do

      call ga_put(g_dens, jlo,jhi, ilo, ihi, dbl_mb(k_dip), jdim)
    
      if (.not. ma_pop_stack(l_dip)) call errquit('multipole: ma?',0,
     &       MA_ERR)

            end if
         end do
      end do

      call ga_sync()

c ================================================================
cjbcjb
cjb      if(nodezero) then
cjb
cjb      if(.not. ma_push_get(mt_dbl,nbf*nbf, 'mult:mp', l_dip, k_dip))
cjb     &     call errquit('multipole: cannot allocate mp', lmpmax, MA_ERR)
cjb
cjb       inquire(iolength=reclen)dbl_mb(k_dip)
cjb
cjb       write(6,"('Record length is: ',I8)")reclen
cjb
cjb       unitn = 76
cjb       open(unit=unitn,file='DIP_AO.bin',
cjb     2 status='replace',access='direct',
cjb     2 form='unformatted',recl=reclen*nbf*nbf)
cjb
cjb      call ga_get(g_dens, 1, nbf, 1, nbf, dbl_mb(k_dip), nbf)
cjb
cjb       write(unitn,rec=1) dbl_mb(k_dip:k_dip+(nbf*nbf)-1)
cjb
cjb       do i =1,nbf*nbf
cjb         write(6,"('DIP_AO:',F18.12)")dbl_mb(k_dip+i-1)
cjb       enddo
cjb
cjb      if (.not. ma_pop_stack(l_dip)) call errquit('multipole: ma?',0,
cjb     &       MA_ERR)
cjb
cjb       close(unitn)
cjb
cjbc ----------------- MO ------------------------------------------
cjb
cjb       open(unit=unitn,file='LCAO.bin',
cjb     2 status='replace',access='direct',form='unformatted',
cjb     3 recl=reclen*nbf*nbf)
cjb
cjb      if(.not. ma_push_get(mt_dbl,nbf*nbf, 'mult:mp', l_dip, k_dip))
cjb     &     call errquit('multipole: cannot allocate mp', lmpmax, MA_ERR)
cjb
cjb      call ga_get(g_movecs, 1, nbf, 1, nbf, dbl_mb(k_dip), nbf)
cjb
cjb       write(unitn,rec=1) dbl_mb(k_dip:k_dip+(nbf*nbf)-1)
cjb
cjb       do i =1,nbf*nbf
cjb         write(6,"('LCAO:',F18.12)")dbl_mb(k_dip+i-1)
cjb       enddo
cjb
cjb      if (.not. ma_pop_stack(l_dip)) call errquit('multipole: ma?',0,
cjb     &       MA_ERR)
cjb
cjb       close(unitn)
cjb
cjb
cjb      endif
cjb
cjb      call ga_sync()
cjb
c ================================================================
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
