      subroutine grad_dens_so2 
     $     (geom, basis, g_dens,  
     $     evals, occ, ndens, nbf, nbf2, nopen, nclosed, scftype, 
     $     movecs, omp2, pdm2, coeff)

c     $Id: grad_dens_so2.F,v 1.1 2000-05-01 05:36:24 zzhang Exp $

C     get densities for exact-exchange sodft: Im(Daa), Im(Dbb), 
c     Re(Dab), Im(Dab)
C    

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"

      integer ndens,geom, basis, g_dens(ndens), 
     $     nbf, nbf2, nopen, nclosed, g_densso(2)
      double precision evals(*), occ(*), pdm2(*), coeff(*)
      character*(*) scftype
      character*(*) movecs
      double precision offset   ! EAF file offset
      character*255 pdmfilename ! EAF file name
      integer pdmfile           ! EAF file handle
      integer l_buf, k_buf

      logical omp2

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_so, file_read_ga
      external movecs_read_so, file_read_ga

      character*255 fname       ! filename for lagrangian (ROHF) 
                                ! or second-order contributions (MP2)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs(2), i
      double precision zero, one, mone, two, mtwo
      logical status

      parameter ( mtwo = -2.D0 )
      parameter ( mone = -1.D0 )
      parameter ( zero =  0.D0 )
      parameter ( one  =  1.D0 )
      parameter ( two  =  2.D0 )

C     get MO vectors from file (taken from Roberts scf_movecs_read)

      status = ga_create(mt_dbl, nbf2, nbf2,'Movecs Re', 0, 0, 
     &                   g_vecs(1))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'Movecs Re'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0)
      endif

      status = ga_create(mt_dbl, nbf2, nbf2,'Movecs Im', 0, 0, 
     &         g_vecs(2))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'Movecs Im'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0)
      endif

      if (.not. movecs_read_so(movecs, occ, evals, g_vecs)) 
     $     call errquit('grad_dens_so: could not read mo vectors', 110)
c
      if (scftype .eq. 'UHF') then
         nocc = 2*nclosed + nopen 
         status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &                   g_densso(1))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' grad_dens_so: ', 'DenMx Re'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('grad_dens_so: ga_create ', 0)
         endif
         status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &        g_densso(2))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' grad_dens_so: ', 'DenMx Im'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('grad_dens_so: ga_create ', 0)
         endif
C     density matrix
         call dft_densm_so(g_densso,g_vecs,nbf,nocc)
         call ga_dens_aabbi(g_dens(3), g_densso, nbf)
         call ga_dens_abri(g_dens(5), g_densso, nbf)
         call ga_sync()
      end if
         
c     

      status = ga_destroy ( g_vecs(1) )
      if ( scftype .eq. 'UHF' ) then
         status = ga_destroy ( g_vecs(2) )
      end if
      status = ga_destroy ( g_densso(1) )
      status = ga_destroy ( g_densso(2) )

      return
      end

