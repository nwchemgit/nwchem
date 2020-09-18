      subroutine fcidump_header(rtdb,occ_frozen,vir_frozen,iunit)
      implicit none
      integer rtdb
      integer occ_frozen,vir_frozen
      integer iunit
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "rtdb.fh"
c     FCIDUMP references      
c     https://doi.org/10.1016/0010-4655(89)90033-7
c     https://bitbucket.org/pjknowles/fci/src/master/
c     https://bitbucket.org/pjknowles/fcidump      
c FCIDUMP header      
C Variable Default value Description
C NORB         _         Number of orbitals                                                   
C NELEC        _         Number of electrons                                                  
C MS2          0         2S, where S is the spin quantum number                          
C ISYM         1         Spatial symmetry of wavefunction                                     
C ORBSYM-      1, 1,     Spatial symmetries of orbitals                                       
C (l:NORB)     
C IPRTIM       -1        If 0, print additional CPU timing analysis                           
C INT          5         Fortran stream from which integrals will be read                     
C MEMORY       100000    Size of workspace array in floating point words                      
C CORE         0.0       Core energy (may also be given in integral file) (floatingpoint)     
C MAXIT        25        Maximum number of iterations in Davidson diagonalisation             
C THR          10-5      Convergence threshold for Davidson diagonalisation (floating point)  
C THRRES       10-1      Threshold for printing final CI coefficients (floating point)        
C NROOT        1         Number of eigenvalues of Hamiltonian to be found
C IUHF         0         0 for Restricted HF, 1 for Unrestricted
C
C HF tz example      
C $FCI NORB=11,NELEC=8,
C  ORBSYM=1,1,1,1,1,1,1,2,2,3,3,
C  MEMORY=150000,
C $END
C   0.884276869077      1  1  1  1
C     ....
      character*16 pname,theory
      character*4 scftype
      parameter(pname='fcidump_header: ')
      integer geom,basis,mult
      double precision nuclear_charge,charge
      integer nelec  ! NELEC
      integer norb   ! NORB
      integer iuhf   ! IUHF
      integer ms2    ! MS2
      integer isym   ! ISYM
      integer, dimension(:), allocatable :: orbsym
      integer mymemory,iorb
      parameter (mymemory=1000000) ! 8MB
C
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     &  call errquit(pname//'theory not specified',0, RTDB_ERR)
c
c     extract NORB<-nbf
c     
c     
c     Basis and geometry
c     
      if (.not. geom_create(geom,'geometry'))
     $     call errquit(pname//'geom failed',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $     call errquit(pname//'geometry not found',0,
     &       RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit(pname//'basis failed',0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit(pname//'ao basis not found',0, RTDB_ERR)
      if (.not. bas_numbf(basis, norb)) 
     $     call errquit(pname//'failed to get nbf', 0, BASIS_ERR)
c     if orbitals are frozen, subract no. frozen orb. from norb
      norb=norb-occ_frozen-vir_frozen

      if (.not. bas_destroy(basis))
     $     call errquit(pname//'basis corrupt',0, BASIS_ERR)

c
c     extract NELEC
c      
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit(pname//'geom_nuc_charge failed', 0, GEOM_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit(pname//'geom corrupt',0, GEOM_ERR)
      nelec = nint(nuclear_charge - charge)
c     if orbitals are frozen, subract no. 2x frozen occ. orb. from elec
      nelec = nelec - 2*occ_frozen

c
c     extract IHUF<-ipol
c
      iuhf=0
      ms2=0
      if(theory(1:inp_strlen(theory)) .eq. 'dft')then
         if(.not. rtdb_cget(rtdb, 'dft:scftype', 1,scftype))
     $        call errquit(pname//'cget scftype failed',0, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1,mult))
     *        mult = 1
         if(mult.ne.1) then
            ms2=mult-1
         endif
         if(scftype(1:inp_strlen(scftype)).eq.'UHF') iuhf=1
      else
         write(luout,*) ' fcidump for theory ',
     &        theory(1:inp_strlen(theory)), ' not implemented yet'
      endif
c      
c      no symmetry for now
c
      isym=1
      allocate(orbsym(norb))
      do iorb=1,norb
         orbsym(iorb)=1
      enddo
C
      if(ga_nodeid().eq.0) then
         write(iunit,1) norb,nelec
         if(iuhf.ne.0) then
            write(iunit,3) iuhf,ms2
         else
            write(iunit,30) ms2
         endif
      endif
 11   format(1x,'$FCI NORB=',i7,',NELEC=',i6,',',/
     A 1x,'IUHF=',i1,',MS2=',i3,',',/
     A 1x,'MEMORY=',i10,/,
     A     1x'$END')
c      
 1    format(1x,'&FCI NORB=',i7,',NELEC=',i6,',')
c      
 3    format(1x,'IUHF=',i1,',MS2=',i3,',')
 30   format(1x,'MS2=',i3,',')
      return
      end
      subroutine fcidump_symm(iunit,norb,isym,orbsym)
      implicit none
#include "global.fh"
      integer iunit
      integer isym,norb
      integer orbsym(*)
      integer mymemory,iorb
      parameter (mymemory=6000000) ! 8MB
c     ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,
      if(ga_nodeid().eq.0) then
         if(isym.gt.1) then
            write(iunit,2) isym,
     c           (orbsym(iorb),iorb=1,norb)
         else
            write(iunit,2) isym,
     c           (1,iorb=1,norb)
         endif
         write(iunit,4)  mymemory
      endif
 2    format(1x,'ISYM='i2,',ORBSYM=',9999(I2,','))
 4    format(1x,'MEMORY=',i10,/,
     A     1x'&END')
      return
      end
