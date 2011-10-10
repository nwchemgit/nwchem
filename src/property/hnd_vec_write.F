      subroutine hnd_vec_write(rtdb,geom,basis,nbf,nclosed,nopen,
     &                            nvirt,scftyp,vectors,occ,evals,nmo,
     &   filename)
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
      integer     nclosed(2)         ! [iput] # of closed shell electrons
      integer     nopen(2)           ! [imput] # of open shell electrons
      integer     nvirt(2)           ! [input] # of virtuals
      character*3 scftyp             ! [input] wave function type
      integer     vectors(2)         ! [input] handles for vectors
      double precision  occ(nbf,2)   ! [input] array with occupation numbers
      double precision  evals(nbf,2) ! [input] array with eigen values
      integer     nmo                ! [input] number of orbitals
c       
      character*255 title, movecs, title_molec, basis_name_molec
      character*(*) filename
      character*20  scftype_molec
      integer i, j, ndens
      integer nbf_molec, nsets_molec, nmo_molec(2)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical  movecs_read_header, movecs_write
      external movecs_read_header, movecs_write
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
c     read MO vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2)) call 
     &    errquit('hnd_vec_write: failed to read header', 911, DISK_ERR)
c

c     write MO vectors to movecs file or a file with a user-specified
c     name

      if (filename.eq.' ' .or. filename.eq.'') then
        filename = movecs
      else
        continue ! let's hope the file name is OK ...
      end if

c      write (luout,*) 'writing localized MOs to file ',trim(movecs)

      if (.not. movecs_write(rtdb, basis, filename,
     $     scftype_molec, title_molec,
     $     nbf_molec, nsets_molec, nmo_molec,
     $     occ(1,1), nbf, evals(1,1),
     $     nbf, vectors)) call errquit
     $     ('hnd_vec_write: movecs_write failed', 666, DISK_ERR)

c
      return
      end
c $Id$
