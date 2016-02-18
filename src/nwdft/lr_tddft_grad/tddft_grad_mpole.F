c
      subroutine tddft_grad_mpole(rtdb, ao_bas_han, ipol, g_adens, 
     &        g_bdens, g_ddens)
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb              ! [input]
      integer ao_bas_han        ! [input]
      integer ipol              ! [input] > 1 for open shell
      integer g_adens, g_bdens  ! [input] Alpha and beta densities (open)
c                                         or total density and junk (closed)
      integer g_ddens           ! [input] Difference density
c
c     Calculate the multipoles using the input density matrices
c
      integer geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments) 
      double precision alphamoments(len_moments)
      double precision betamoments(len_moments)
      double precision diffmoments(len_moments)
      double precision nuclmoments(len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      data center/3*0.0d0/
c
c     Preliminaries
c
      if (.not. bas_cando_mpoles(ao_bas_han)) return
      if (.not. bas_geom(ao_bas_han, geom)) call errquit
     $     ('dft_mpole: bas_geom failed',0, BASIS_ERR)
c
c     Initialize     
c
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, alphamoments, 1)
      call dfill(len_moments, 0.0d0, betamoments, 1)
      call dfill(len_moments, 0.0d0, diffmoments, 1)
c
c     Nuclear part
c
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
c
c     Electronic part from the alpha and beta densities
c
      call multipole_density(ao_bas_han, center, lmax, g_adens, 
     $     alphamoments, len_moments)
      if (ipol .gt. 1) then
         call multipole_density(ao_bas_han, center, lmax, g_bdens, 
     $     betamoments, len_moments)
      else
         do k = 1, len_moments
            alphamoments(k) = 0.5d0 * alphamoments(k)
            betamoments(k)  = alphamoments(k)
         enddo
      endif
c
c     Electronic part from difference density 
c
      call multipole_density(ao_bas_han, center, lmax, g_ddens, 
     $     diffmoments, len_moments)
c
      do k = 1, len_moments
         totalmoments(k) = alphamoments(k) + betamoments(k) +
     $        diffmoments(k) +  nuclmoments(k)
      enddo
c
c     Print multipoles
c
      if (ga_nodeid().eq.0) then
         write(LuOut,*)
         call util_print_centered(LuOut,
     $        'Multipole analysis of the density', 20, .true.)
         write(LuOut,*)
         write(LuOut,*) '    L   x y z        total'
         write(LuOut,*) '    -   - - -        -----'
         ind = 0
         do ltotal = 0, lmax
            do k = ltotal, 0, -1
               do l = ltotal-k, 0, -1
                  m = ltotal - k - l
                  ind = ind + 1
                  write(LuOut,12) ltotal, k, l, m, totalmoments(ind)
 12               format(4x,i2,2x,3i2,f14.6)
               enddo
            enddo
            write(LuOut,*)
         enddo
      endif ! ga_nodeid()
c
c     Store the dipole moment
c
      if (.not. rtdb_put(rtdb, 'tddft:dipole', mt_dbl, 
     $     3, totalmoments(2))) 
     &  call errquit('tddft_grad_mpole: rtdb?',0, RTDB_ERR)
c     
      end
