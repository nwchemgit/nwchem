      subroutine hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                            nvirt,scftyp,vectors,occ,evals,nmo)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer     rtdb               ! [input] rtdb handle
      integer     geom               ! [input] geometry handle
      integer     basis              ! [input] basis handle
      integer     nbf                ! [input] number of basis functions
      integer     nclosed(2)         ! [output] # of closed shell electrons
      integer     nopen(2)           ! [output] # of open shell electrons
      integer     nvirt(2)           ! [output] # of virtuals
      character*3 scftyp             ! [output] wave function type
      integer     vectors(2)         ! [output] handles for vectors
      double precision  occ(nbf,2)   ! [output] array with occupation numbers
      double precision  evals(nbf,2) ! [output] array with eigen values
      integer     nmo                ! [output] number of orbitals
c       
      character*255 title, movecs, title_molec, basis_name_molec
      character*20  scftype_molec
      integer i, j, ndens
      integer nbf_molec, nsets_molec, nmo_molec(2)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
c
c     Extract high level info from the data-base setting defaults
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'prop:input vectors',  1,
     $     movecs )) then
           call util_file_name('movecs', .false.,.false.,movecs )
      endif
      call util_file_name_resolve(movecs , .false.)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2)) call 
     &    errquit('hnd_vec_read: failed to read header', 911, DISK_ERR)
c
      scftyp = 'RHF'
      nmo    = nmo_molec(1)
      ndens  = 1
      if (nsets_molec.gt.1) then
         scftyp = 'UHF'
         nmo    = max(nmo_molec(1),nmo_molec(2))
         ndens  = 3
      endif
c
c     Create handles for vectors
c
      vectors(1) = ga_create_atom_blocked(geom,basis,'mo vectors' )
      if (.not.movecs_read(movecs,1,occ(1,1),evals(1,1),vectors(1)))
     &    call errquit('prop_vec_read: error reading mo vec',0,DISK_ERR)
c
      if (scftyp.eq.'UHF') then
         vectors(2) = ga_create_atom_blocked(geom,basis,'mo vectors 2')
         if (.not.movecs_read(movecs,2,occ(1,2),evals(1,2),vectors(2)))
     &       call errquit('prop_vec_read: error readin mo vec 2',0,
     &       DISK_ERR)
      endif
c
c     Determine number of closed, open and virtuals
c
      do i = 1, nsets_molec
         nclosed(i) = 0
         nopen(i)   = 0
         nvirt(i)   = 0
         do j = 1, nmo_molec(i)
            if (occ(j,i).eq.2.0d0) then
                nclosed(i) = nclosed(i)+1
            elseif(occ(j,i).ne.0.0d0) then
                nopen(i) = nopen(i)+1
            else
                nvirt(i) = nvirt(i)+1
            endif
         enddo
      enddo
c
      return
      end
c $Id$
