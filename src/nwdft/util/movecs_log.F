      logical function movecs_to_logfile( rtdb, n_homo, n_lumo,
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "errquit.fh"
C$Id$
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
c     
c**** output the coefficients of a specified number of HOMOs and LUMOs 
c**** to standard output.  (For parsing into ECCE').   Hacked from 
c**** movecs_write
c**** 
c**** arguments
c**** 
      integer rtdb ! runtime database handle
      integer n_homo ! no. of highest occupied orbs to output
      integer n_lumo ! no. of lowest unoccupied orbs to output
      integer nbf ! No. of functions in orbital basis
      integer nsets ! No. of sets of vectors
      integer nmo(nsets) ! No. of vectors in each set
      integer ldocc
      integer ldevals
      double precision occ(ldocc, nsets) ! Occupation numbers of each set
      double precision evals(ldevals, nsets) ! Eigenvalues for each set
      integer g_vecs(nsets) ! Global array with eigen-vectors
c**** 
c**** Local variables
c**** 
      double precision zero
      parameter ( zero=0.d0 )
      double precision occ_number
      integer noc(2) ! no. of occupied orbitals
      integer iorb, i_homo, i_lumo
      integer l_vecs, k_vecs
c****      integer ok 
      integer iset, j
      integer i_below_gap
      character*6 alpha_or_beta(2)
c----------------------------------------------------------------------- 
      alpha_or_beta(1) = 'alpha '
      alpha_or_beta(2) = 'beta  '

      if ( nsets .eq. 1) then ! spin-restricted
        alpha_or_beta(1) = '      '
      endif 
c**** number of occupied orbitals
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &     call errquit(
     &     'movecs_to_logfile: rtdb_get of noc failed', 0, RTDB_ERR)

c**** vectors.F has a ga_sync at this point...
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
        if (nsets .gt. 2) then
          write(LuOut,*) ' movecs_write_logfile: nsets > 2 ', nsets
          call errquit('movecs_write_logfile: ',0, DISK_ERR)
        endif

c****   number of HOMOs and LUMOs to output, if any
c****   if (.not. rtdb_get(rtdb, 'dft:n_HOMO_to_logfile',
c****   &       mt_int, 1, n_homo)) n_homo = 0
c****   if (.not. rtdb_get(rtdb, 'dft:n_LUMO_to_logfile',
c****   &       mt_int, 1, n_lumo)) n_lumo = 0
c****   n_homo = min(n_homo,nbf)
c****   n_lumo = min(n_lumo,nbf)

c****   array to store one orbital

        if(.not. ma_push_get(mt_dbl, nbf, 'movecs', l_vecs, k_vecs))
     &    call errquit('movecs_log: cannot allocate movecs', 0, MA_ERR)

c****   write header to logfile
        write(LuOut,*)'<MO coefficients> start'

        write (LuOut,*) 'n_basis_functions: ',nbf,
     &       '      n_HOMO: ',
     &       n_homo, '      n_LUMO: ',n_lumo

c****   output n_homo highest occupied orbitals and n_lumo lowest
c****   unoccupied orbitals
        occ_number = 3-nsets
        do iset = 1, nsets ! alpha, then beta 

c****     output HOMOs
          occ_number = 3-nsets
          do i_homo = 1, n_homo
            iorb = noc(iset) - n_homo + i_homo
            i_below_gap = n_homo - i_homo + 1
            write (LuOut,'(a,a,a,i2)')
     &           'Orbital_Label: ',alpha_or_beta(iset),
     &           'HOMO ',i_below_gap
            write (LuOut,'(a,f4.1)')
     &           'Occupation_Number=', occ_number ! label
            write (LuOut,'(a,e14.7)')
     &           'Orbital_energy=', evals(iorb,iset)
            call ga_get(g_vecs(iset), 1, nbf, iorb, iorb,
     &           dbl_mb(k_vecs),1)
            write(LuOut,'(4(e18.7))') (dbl_mb(k_vecs+j), j=0,nbf-1)
          enddo

c****     output LUMOs
          occ_number = zero
          do i_lumo = 1, n_lumo
            iorb = noc(iset) + i_lumo 
            write (LuOut,'(a,a,a,i2)') 'Orbital_Label: ',
     &           alpha_or_beta(iset),
     &           'LUMO ',i_lumo
            write (LuOut,'(a,f4.1)')
     &           'Occupation_number=', occ_number ! label
            write (LuOut,'(a,e14.7)')
     &           'Orbital_energy=', evals(iorb,iset)
            call ga_get(g_vecs(iset), 1, nbf, iorb, iorb,
     &           dbl_mb(k_vecs),1)
            write(LuOut,'(4(e18.7))') (dbl_mb(k_vecs+j), j=0,nbf-1)
          enddo

        enddo

c****   write footer to logfile
        write(luout,*)'<MO coefficients> end'

        if(.not. ma_pop_stack(l_vecs))
     &    call errquit('movecs_log: cannot pop stack',0, MA_ERR)

c****   ok = 1
      endif

      movecs_to_logfile = .true.

      return 
c     
      end







