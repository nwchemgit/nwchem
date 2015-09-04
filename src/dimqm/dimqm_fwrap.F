      subroutine dimqm_fwrap(g_Ax)
c
c     Determines the induced dipoles for the given density matrix 
c    
c     Called from aoresponse.F
c
      implicit none
#include "cscf.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "constants.fh"
#include "dimqm.fh"
c
c     Input Variables
      integer g_Ax                 ! A matrix handle
c
c     Local variables
      double precision dx, dy, dz
      double precision dsum
      external dsum
      integer i3, n
      integer icmplx
      integer l_fld, k_fld
      integer nNuc
      integer l_qmxyz, k_qmxyz
c      integer l_pot, k_pot
      integer g_dim
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer dims(3), chunk(3)
      integer nvir, voff, xoff
      integer g_tmp1, g_tmp2, g_dcv, g_tmp3
      double precision Edim(2)
      integer nExc
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      data direction /'x', 'y', 'z'/

c
c     Common variables used
c
c     logical lcmplx
c     integer l_dimxyz, k_dimxyz
c     integer l_muind, k_muind
c     integer l_qind, k_qind
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'nfockdim:3ddim',chunk,
     &                        g_dim))
     &  call errquit('hessv_2e2: could not allocate g_dim',555,
     &                 GA_ERR)
      call ga_zero(g_dim)
      if(.not. ga_duplicate(g_Ax, g_tmp3, 'precond temp'))
     $  call errquit('fwrap: could not duplicate g_ax', 1, GA_ERR)
      Edim(1) = 0.0d0
      Edim(2) = 0.0d0
      nExc = 1
c      write(luout,*) "Calling dim_fock_xc"
      call fock_dim(geom, nbf, basis, 3, g_fdens, g_dim, .true.)
      xoff = 1
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'indD4: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'indD4: tmp2')
      g_dcv = ga_create_atom_blocked(geom, basis, 'indD4: dcv')
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      do n = 1, 3
        alo(1) = n
        ahi(1) = n
        xoff = 1
        call nga_copy_patch('N',g_dim,alo,ahi,g_dcv,blo,bhi)
        call ga_sync()
c        call ga_dscal(g_dcv, two)
        call ga_zero(g_tmp1)
        call ga_zero(g_tmp2)
        call ga_matmul_patch('n', 'n', four, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
        call ga_matmul_patch('t', 'n', four, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
        call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax,
     $                         xoff, n, 1.0d0, '+')
c        xoff = nvir*nclosed
c        call rohf_hessv_precond(g_tmp3, 0.0d0)
c        write(luout,*) "Done precon"
c        call ga_add_patch(1.0d0, g_Ax, 1, xoff, n, n,
c     $                    1.0d0, g_tmp3, 1, xoff, n, n,
c     $                    g_Ax, 1, xoff, n, n)
        call ga_sync()
      end do
      if(.not. ga_destroy(g_dcv))
     $  call errquit('dimqm_fwrap: destroy g_dcv failed', 1, GA_ERR)
      if(.not. ga_destroy(g_tmp2))
     $  call errquit('dimqm_fwrap: destroy g_tmp2 failed', 1, GA_ERR)
      if(.not. ga_destroy(g_tmp1))
     $  call errquit('dimqm_fwrap: destroy g_tmp1 failed', 1, GA_ERR)
      if(.not. ga_destroy(g_tmp3))
     $  call errquit('dimqm_fwrap: destory g_tmp3 failed', 1, GA_ERR)
      if(.not. ga_destroy(g_dim))
     $  call errquit('dimqm_fwrap: destroy g_dim failed', 1, GA_ERR)
      if(.not. ga_destroy(g_fdens))
     $  call errquit('dimqm_fwrap: destroy g_fdens failed', 1, GA_ERR)
c      write(luout,*) "Done fwrap"
      end subroutine dimqm_fwrap
