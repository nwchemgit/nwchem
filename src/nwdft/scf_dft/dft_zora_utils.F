c
czora...Scale zora eigenvalues and energy
c
          subroutine dft_zora_scale(geom, 
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   g_dens,
     &                   g_s,
     &                   g_movecs,
     &                   g_zora_scale_sf,
     &                   evals,
     &                   focc,
     &                   noc,
     &                   ipol,
     &                   ener_scal)

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_dens(2)
      integer g_movecs(2)
      integer g_s
      integer g_zora_scale_sf(2)

      double precision numelecs, ener_kin
      integer g_orb
      integer g_dens_sf
      integer l_vecs,k_vecs
      integer iorb
      double precision eval_scal
      double precision ener_scal
      integer noc(2)
      integer ispin
      integer ipol

      integer geom
      integer ao_bas_han
      integer nbf,nbf_ao

      double precision focc(nbf*ipol)  ! occupation no.
      double precision evals(ipol*nbf) ! eigenvalues
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     debugging
c     separate the density
c      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1))
c
c     electron number test
c      numelecs = ga_ddot(g_dens(1),g_s)
c      write(luout,*) "Number of spin up electrons:",numelecs
c      if (ipol.gt.1) then
c       numelecs = ga_ddot(g_dens(2),g_s)
c       write(luout,*) "Number of spin down electrons:",numelecs
c      end if
c
c     kinetic energy test
c      ener_kin = ga_ddot(g_dens(1),g_zora_scale_sf(1))   
c      write(luout,*) "Kinetic energy spin up:",ener_kin
c      if (ipol.gt.1) then
c        ener_kin = ga_ddot(g_dens(2),g_zora_scale_sf(2))
c        write(luout,*) "Kinetic energy spin down:",ener_kin
c      end if
c
c     consolidate the density before leaving
c     call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     allocate memory
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'vec aux',
     &            l_vecs, k_vecs))
     &     call errquit('xc_sic: cannot allocate vec aux',0, MA_ERR)
      g_orb = ga_create_atom_blocked(geom, ao_bas_han, 'orbs')
      g_dens_sf = ga_create_atom_blocked(geom, ao_bas_han,'orbs dens')
c
c     scale the eigenvalues and energy
      ener_scal = 0.d0
      do ispin=1,ipol
       do iorb = 1,nbf
c
        call ga_get(g_movecs(ispin), 1, nbf, iorb, iorb,
     &                      dbl_mb(k_vecs), 1)
        call ga_zero(g_orb)
        call ga_put(g_orb, 1, nbf, iorb, iorb,
     &                    dbl_mb(k_vecs), 1)
c
        call ga_zero(g_dens_sf)
        call ga_dgemm('n', 't', nbf, nbf, noc(ispin), 1.0d00,
     &                      g_orb, g_orb, 0.d00, g_dens_sf)
c
        eval_scal = evals(iorb)
        if (ispin.gt.1) eval_scal = evals(iorb + nbf)
c
        eval_scal 
     &    = eval_scal/(1.d0 + ga_ddot(g_dens_sf,g_zora_scale_sf(ispin)))
c
        ener_scal = ener_scal
     &     - eval_scal*ga_ddot(g_dens_sf,g_zora_scale_sf(ispin))
c
c       ener_scal = ener_scal + ga_ddot(g_orb_dens,g_zora_scale_sf(ispin))
c
        if (ispin.le.1) then
            evals(iorb) = eval_scal
c            write(luout,*) iorb,evals(iorb),eval_scal
        else
            evals(iorb+nbf) = eval_scal
c            write(luout,*) iorb,evals(iorb+nbf),eval_scal
        end if
c
       end do ! orbital loop
      end do ! polarization loop
c
c     double the energy for closed-shell calculations
      if (ipol.eq.1) ener_scal = 2.d0*ener_scal
      
c      write(luout,*) "scaled_energy:",ener_scal

c     deallocate emory
      if (.not.ma_pop_stack(l_vecs))
     &   call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      if (.not. ga_destroy(g_orb)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_dens_sf)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)

c      stop

      return
      end
c
czora...Inquire if the zora correction file is present
c
      logical function dft_zora_inquire_file(filename)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character*(*) filename
      logical found
c
      call ga_sync()
c
c     Inquire if file is present
      inquire(file=filename,exist=found)
      dft_zora_inquire_file = found 
c
      call ga_sync()
c
      return
      end
c
czora...Read in the zora atomic corrections from disk
c
      logical function dft_zora_read(filename, nbf, nsets, nmo, 
     &               mult, g_zora_sf, g_zora_scale_sf)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer nsets             ! restricted or unrestricted
      character*(*) filename
      integer iset              ! restricted or unrestricted
      integer g_zora_sf(nsets)
      integer g_zora_scale_sf(nsets)
c
      integer nbf               ! No. of functions in basis
      integer nmo(nsets)       
      integer mult              ! multiplicity
      integer ok, jset, i, j
      integer l_zora, k_zora
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c
      integer nsets_read       ! No. of sets  from file
      integer nbf_read         ! No. of functions from file 
      integer mult_read
c     
c     Initialise to invalid MA handle
      l_zora = -1
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read atomic ZORA corrections from ',a/)
       call util_flush(luout)
c
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nsets_read
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) mult_read
c
c      Error checks
       if ((nsets_read .ne. nsets)
     &  .or. (nbf_read .ne. nbf)
     &  .or. (mult_read .ne. mult) ) goto 1003
c
c      Allocate the temporary buffer
       if (.not. ma_push_get(mt_dbl,nbf,'dft_zora_read',l_zora,k_zora))
     $        call errquit('dft_zora_read: ma failed', nbf, MA_ERR)
c
c      Read in g_zora_sf
       do iset = 1, nsets
        do i = 1, nmo(iset)
         call sread(unitno, dbl_mb(k_zora), nbf)
         call ga_put(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora), 1)
        end do
       end do
c
c      Read in g_zora_scale_sf
       do iset = 1, nsets
        do i = 1, nmo(iset)
         call sread(unitno, dbl_mb(k_zora), nbf)
         call ga_put(g_zora_scale_sf(iset), 1, nbf, i, i, 
     & dbl_mb(k_zora), 1)
        end do
       end do
c
c      Close the file
       close(unitno,err=1002)
       ok = 1
c
c      Deallocate the temporary buffer
       if (.not. ma_pop_stack(l_zora)) call errquit
     $      ('dft_zora_read: pop failed', l_zora, MA_ERR)
c
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
c       write(6,*)' g_zora_sf(1) from dft_scf'
c       call ga_print(g_zora_sf(1))
c       write(6,*)' g_zora_scale_sf(1) from dft_scf'
c       call ga_print(g_zora_scale_sf(1))
c
      dft_zora_read = ok .eq. 1
      return
c
 1000 write(6,*) 'dft_zora_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1003 write(6,*) 'dft_zora_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
czora...Write out the zora atomic corrections to disk
c
      logical function dft_zora_write(rtdb, basis, filename,
     &     nbf, nsets, nmo, mult, g_zora_sf, g_zora_scale_sf)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nsets             ! [input] No. of sets of matrices
      integer nmo(nsets)        ! [input] No. of mos in each set
      integer mult
c
      integer g_zora_sf(nsets)    
      integer g_zora_scale_sf(nsets)  
c
      integer unitno
      parameter (unitno = 77)
c
      integer lentit
      integer lenbas
      integer l_zora, k_zora
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, key
      character*20 scftype20 
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize
c
      l_zora = -1               ! An invalid MA handle
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
      ok = 0
c
c     Read routines should be consistent with this
c
c     Write out the atomic zora corrections
c
      if (ga_nodeid() .eq. 0) then
c
c     Open the file
      open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c     Write out the number of sets and basis functions
      write(unitno, err=1001) nsets
      write(unitno, err=1001) nbf
      write(unitno, err=1001) mult
c
c     Allocate the temporary buffer
      if (.not. ma_push_get(mt_dbl,nbf,'dft_zora_write',l_zora,k_zora))
     $        call errquit('dft_zora_write: ma failed', nbf, MA_ERR)
c
c     Write out g_zora_sf
      do iset = 1, nsets
       do i = 1, nmo(iset)
        call ga_get(g_zora_sf(iset), 1, nbf, i, i, dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Write out g_zora_scale_sf
      do iset = 1, nsets
       do i = 1, nmo(iset)
        call ga_get(g_zora_scale_sf(iset), 1, nbf, i, i, 
     &            dbl_mb(k_zora),1)
        call swrite(unitno, dbl_mb(k_zora), nbf)
       end do
      end do
c
c     Deallocate the temporary buffer
      if (.not. ma_pop_stack(l_zora))
     $  call errquit('dft_zora_write: ma pop failed', l_zora, MA_ERR)
c
c     Close the file
      close(unitno,err=1002)
c
      ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
c       write(6,*)' g_zora_sf(1) from dft_scf'
c       call ga_print(g_zora_sf(1))
c       write(6,*)' g_zora_scale_sf(1) from dft_scf'
c       call ga_print(g_zora_scale_sf(1))
c
      dft_zora_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote atomic ZORA corrections to ',a/)
         call util_flush(luout)
      endif
c
      return
c
 1000 write(6,*) 'dft_zora_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) 'dft_zora_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) 'dft_zora_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c
c     Convergence check for zora calculations
c
      subroutine dft_zora_scfcvg(rms, derr, etold, etnew, e_conv, 
     &                      d_conv, g_conv, ipol, iter, iterations, 
     &                      idone, rtdb, converged, diising)
c
c     $Id: dft_scfcvg.F,v 1.14 2003-10-17 22:55:44 carlfahl Exp $
c
      implicit none
#include "errquit.fh"
c
      double precision rms(2)   ! [input]
      double precision derr(2)  ! [input]
      double precision etold    ! [input]
      double precision etnew    ! [input]
      double precision e_conv   ! [input]
      double precision d_conv   ! [input]
      double precision g_conv   ! [input]
      integer ipol              ! [input]
      integer iter              ! [input]
      integer iterations        ! [input]
      integer idone             ! [output]
      integer rtdb              ! [input]
      logical converged         ! [output]
      logical diising           ! [input]
c
#include "mafdecls.fh"
#include "rtdb.fh"
c     
      logical e_conv_logical, d_conv_logical, g_conv_logical
      logical ENERGY, DENSITY, GRADIENT
      double precision de, abde
c
      converged = .false.
c
      e_conv_logical = .false.
      d_conv_logical = .false.
      g_conv_logical = .false.
c
      ENERGY = e_conv.gt.0
      DENSITY = d_conv.gt.0
      GRADIENT = g_conv.gt.0
c
      idone = 0
c
c     Evaluate change in energy.
c     
      de = etnew - etold
      etold=etnew
      abde = dabs(de)
c     
c     Check to see if energy is converged.
c
      if (ENERGY)then     
         if (abde.lt.e_conv)e_conv_logical = .true.
      else
         e_conv_logical = .true.
      endif
c     
c     Check for density matrix convergence.
c     
      if (DENSITY)then
         if (dsqrt(rms(1)).le.d_conv)d_conv_logical = .true.
         if (ipol.eq.2)then
           if (dsqrt(rms(2)).le.d_conv) then
             d_conv_logical = d_conv_logical.and..true.
          else
             d_conv_logical = d_conv_logical.and..false.
          endif
         endif
      else
         d_conv_logical = .true.
      endif
c     
c     Check for gradient convergence.
c     
      if (GRADIENT.and.diising)then
         if (derr(1).le.g_conv)g_conv_logical = .true.
         if (ipol.eq.2)then
           if (derr(2).le.g_conv) then
              g_conv_logical = g_conv_logical.and..true.
             else
                g_conv_logical = g_conv_logical.and..false.
             endif
         endif
      else
         g_conv_logical = .true.
      endif
c
c     Check over-all convergence.
c
      converged = e_conv_logical
      if (converged)idone = 1
c     
c     Check iteration value.
c     
      if (iter.lt.1)then
         return
      elseif (iter.eq.iterations)then
         idone = 1
      endif
c
c     If all convergence criterion met or number of iterations has been
c     exceeded, write "converged" to RTDB.
c
      if (idone.eq.1)then
         if (.not.rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &      call errquit('dft_scfcvg: rtdb_put failed', 1, RTDB_ERR)
      endif
      return
      end
