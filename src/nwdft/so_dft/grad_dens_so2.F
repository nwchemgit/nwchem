      subroutine grad_dens_so2 
     $     (g_dens,  
     $     evals, occ, ndens, nbf, nbf2, nopen, nclosed,
     $     movecs)

c     $Id$

C     get densities for exact-exchange sodft: Im(Daa), Im(Dbb), 
c     Re(Dab), Im(Dab)
C    

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"

      integer ndens, g_dens(ndens), 
     $     nbf, nbf2, nopen, nclosed, g_densso(2)
      double precision evals(*), occ(*)
      character*(*) movecs


      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_so, file_read_ga
      external movecs_read_so, file_read_ga

                                ! or second-order contributions (MP2)

C     local variables
      integer nocc, g_vecs(2)
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
         call errquit('grad_dens_so: ga_create ', 0, GA_ERR)
      endif

      status = ga_create(mt_dbl, nbf2, nbf2,'Movecs Im', 0, 0, 
     &         g_vecs(2))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'Movecs Im'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0, GA_ERR)
      endif

      if (.not. movecs_read_so(movecs, occ, evals, g_vecs)) 
     $     call errquit('grad_dens_so2:could not read mo vectors',110,
     &       DISK_ERR)
c
      nocc = 2*nclosed + nopen 
      status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &     g_densso(1))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'DenMx Re'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0, GA_ERR)
      endif
      status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &     g_densso(2))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'DenMx Im'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0, GA_ERR)
      endif
C     density matrix
c
      call dft_densm_so(g_densso,g_vecs,nbf,nocc)
c
      call ga_dens_aabbi(g_dens(3), g_densso, nbf)
      call ga_dens_abbar(g_dens(5), g_densso, nbf)
      call ga_dens_abbai(g_dens(7), g_densso, nbf)
c
      call ga_add(one, g_dens(3), one, g_dens(4), g_dens(3))
      call ga_scale(g_dens(4), 2.0d0)
      call ga_add(one, g_dens(3), mone, g_dens(4), g_dens(4))
c
      call ga_add(one, g_dens(5), one, g_dens(6), g_dens(5))
      call ga_scale(g_dens(6), 2.0d0)
      call ga_add(one, g_dens(5), mone, g_dens(6), g_dens(6))
c
      call ga_add(one, g_dens(7), one, g_dens(8), g_dens(7))
      call ga_scale(g_dens(8), 2.0d0)
      call ga_add(one, g_dens(7), mone, g_dens(8), g_dens(8))
      call ga_sync()
c     
      status = ga_destroy ( g_vecs(1) )
      status = ga_destroy ( g_vecs(2) )
      status = ga_destroy ( g_densso(1) )
      status = ga_destroy ( g_densso(2) )
      return
      end

