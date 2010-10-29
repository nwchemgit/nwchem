C $Id$
c
c  Sets up NWChem and the DNTMC geometry information from the input
c  information.  When come out of here should have a geometry in the rtdb
c  corresponding to the cluster, and the initial coordinates (overlapping for
c  the same species), masses, total number of molecules, number of species,
c  number of atoms in a species, and list of unique species numbers will be
c  populated.
c  
C     version 5.0  Dec. 20, 2007 full release version (lcrosby)
C     June 10, 2008 modified to support for dynamic arrays MA (lcrosby)


      subroutine dntmc_getspecinfo(rtdb, nmols, l_natom, k_natom, 
     & l_rinit, k_rinit, l_uspc, k_uspc, nspecies, l_m, k_m, l_spc,
     & k_spc, l_tag, k_tag, nattotal, natmax)

      implicit none
#include "gibbsP.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
      integer rtdb     ! [input] rtdb handle
      integer nmols    ! [output] number of molecules
      integer nspecies ! [output] number of species
      integer l_spec, k_spec, i, j, k
      integer nattotal, natmax

      integer l_natom, k_natom
      integer l_uspc, k_uspc
      integer l_nmolunique, k_nmolunique
      integer l_spc, k_spc

      integer l_rinit, k_rinit
      integer l_m, k_m
      integer l_tag, k_tag

      external dntmc_get_iarray_element
      integer dntmc_get_iarray_element

      external dntmc_put_iarray_element

c
c  Get the number of species and molecules off the rtdb
c
      if (.not.rtdb_get(rtdb, 'dntmc:nspecies', MT_INT, 1, nspecies))
     &    call errquit('dntmc_getspecinfo: rtdb error for nspecies', 
     &                  0, RTDB_ERR)
c
c  Get the number of molecules for each species
c
      if (.not. ma_push_get(MT_INT,nspecies,'dntmc:natom',
     &    l_natom, k_natom)) call errquit
     &    ('entmc_getspecinfo: cannot allocate natom',0,MA_ERR)

      if (.not. ma_push_get(MT_INT,nspecies,'dntmc:uspc',
     &    l_uspc, k_uspc)) call errquit
     &    ('entmc_getspecinfo: cannot allocate natom',0,MA_ERR)


      if (.not. ma_push_get(MT_INT,nspecies,'dntmc:nmolunique',
     &    l_nmolunique, k_nmolunique)) call errquit
     &    ('dntmc_getspecinfo: cannot allocate nmolunique',0, MA_ERR)

      if (.not.rtdb_get(rtdb, 'dntmc:nmol', MT_INT, nspecies, 
     &    int_mb(k_nmolunique)))
     &    call errquit('dntmc_getspecinfo: rtdb error for nmolunique', 
     &                  0, RTDB_ERR)
c
c Sum up each of the seperate molecule numbers to get the total, nmols
c
      nmols = 0
      do i = 1, nspecies
         nmols = nmols + 
     &   dntmc_get_iarray_element(int_mb(k_nmolunique),nspecies,i)
      enddo

C Setup list of species numbers for each molecule (done in blocks)
      if (.not. ma_push_get(MT_INT,nmols,'dntmc:spc',l_spc,k_spc))
     & call errquit
     &('dntmc_getspecinfo: cannot allocate spc',0,MA_ERR)

      K = 1
      DO I = 1, nspecies
        DO J = K, K+
     & dntmc_get_iarray_element(int_mb(k_nmolunique),nspecies,i)-1
          call dntmc_put_iarray_element(int_mb(k_spc),nmols,J,I)
        END DO
        K = J
      END DO

c     write(6,*) 'in getspecinfo nmols is ', nmols
c
c  Get the geometry names out of the rtdb and load the geometries
c
      if (.not.ma_push_get(MT_BYTE, nspecies*255, 'dntmc: species',
     &    l_spec, k_spec)) call errquit
     &    ('dntmc_getspecgeom: cannot allocate species', 0, MA_ERR)
      call dntmc_getspecgeomnames(rtdb, nspecies, byte_mb(k_spec))
      call dntmc_getspecgeom(rtdb, nspecies, byte_mb(k_spec),
     & int_mb(k_uspc), int_mb(k_natom), l_rinit, k_rinit, l_m, k_m, 
     & int_mb(k_nmolunique), l_tag, k_tag, nattotal, natmax)
C Dont pop getspecgeom allocates additional entries on the stack
C Will be chopped later
C      if (.not.ma_pop_stack(l_spec)) call errquit
C     &    ('dntmc_getspecgeom: ma pop failed for species', 1, MA_ERR)
C       write(6,*) 'end of routine'
C       call flush(6)


      return
      end
c
C Little subroutine to help with the character-byte translation
c
      subroutine dntmc_getspecgeomnames(rtdb, nspecies, species)
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb   ! [input]
      integer nspecies ! [input]
      character*255 species(nspecies) ! [output] list of geometry names
c
c
      if (.not.rtdb_cget(rtdb, 'dntmc:species', nspecies, species))
     &   call errquit('dntmc_getspecgeomnames: rtdb failure for species'
     &     , 0, RTDB_ERR)

      return
      end
C
C Get information from the geometries that are stored on the rtdb
C This is really the routine that does most of the heavy lifting
C
      subroutine dntmc_getspecgeom(rtdb, nspecies, species, uspc, natom,
     &           l_rinit, k_rinit, l_ms, k_ms, nmolunique, l_tag, k_tag,
     &           nattotal, natmax)

      implicit none
#include "gibbsP.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
      integer rtdb     ! [input] rtdb handle
      integer nspecies ! [input] number of species
      integer nmolunique(nspecies)  ! [input] number of molecules in each species
      character*255 species(nspecies) ! [input] list of geometry names
      integer natom(nspecies)    ! [output] number of atoms in a particular species
      integer uspc(nspecies) ! [output] list of unique species numbers

      integer l_rinit, k_rinit
      integer l_ms, k_ms
      integer l_tag, k_tag
c
      integer i, nattotal, istart, geom, natuspc, natmax, nmols
      integer l_ct, k_ct, l_qt, k_qt, l_mt, k_mt, l_tt, k_tt ! cluster indices
      integer l_c, k_c, l_q, k_q, l_m, k_m, l_t, k_t ! species indices
c
c  Find out the total number of atoms overall
c  I do this in a seperate loop so that I can set up some temporary
c  space for the full cluster geometry.
c  There certainly a better way to do this
c
c     write (6,*) 'the number of species is ',nspecies
c     write(6,*) 'species(1) ', species(1)
c     write(6,*) 'nmolunique1 and nmolunique2 ', nmolunique(1),
c    &   nmolunique(2)
      nattotal = 0
      natuspc = 0
      natmax = 0
      nmols = 0
      do i = 1, nspecies
         if (.not. geom_create(geom, species(i)))
     &     call errquit('dntmc_getspecgeom: failed creating geometry?',
     &                   0, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom, species(i)))
     &     call errquit('dntmc_getspecgeom: no geometry ', 0,
     &                   RTDB_ERR)
         if (.not. geom_ncent(geom, natom(i)))
     &     call errquit('dntmc_getspecgeom: ncent error', 0,
     &                   GEOM_ERR)
         nattotal = nattotal + (natom(i)*nmolunique(i))
         natuspc = natuspc + natom(i)
         if (natmax .lt. natom(i)) natmax = natom(i)
         nmols = nmols + nmolunique(i)
         if (.not. geom_destroy(geom))
     &     call errquit('dntmc_getspecgeom: destroying geometry',0,
     &       GEOM_ERR)
      enddo

C Allocate some other arrays
      if (.not. ma_push_get(mt_dbl, 3*natmax*nmols,'dntmc:rinit',
     & l_rinit, k_rinit)) call errquit
     &('dntmc_getspecgeom: rinit not allocated',0,MA_ERR)

      if (.not. ma_push_get(mt_dbl, nspecies*natmax,'dntmc:ms',
     & l_ms, k_ms)) call errquit
     &('dntmc_getspecgeom: ms not allocated',0,MA_ERR)

      if (.not. ma_push_get(mt_byte, 16*nattotal,'dntmc:tag',
     & l_tag, k_tag)) call errquit
     &('dntmc_getspecgeom: tag not allocated',0,MA_ERR)

c
c Allocate some buffers for the full cluster geometry
c
c     write(6,*) 'the number of atoms is ',nattotal
      if (.not.ma_push_get(mt_dbl, 3*nattotal,'coordinate', l_ct, k_ct))
     &   call errquit('dntmc_getspecgeom: coordinate not allocated',
     &                 3*nattotal, MA_ERR)
      if (.not.ma_push_get(mt_dbl, nattotal, 'charge', l_qt, k_qt))
     &   call errquit('dntmc_getspecgeom: charge not allocated',
     &                 nattotal, MA_ERR)
      if (.not.ma_push_get(mt_dbl, nattotal, 'mass', l_mt, k_mt))
     &   call errquit('dntmc_getspecgeom: mass not allocated',
     &                 nattotal, MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nattotal,'tags',l_tt,k_tt))
     &   call errquit('dntmc_getspecgeom: tags not allocated',
     &                 nattotal, MA_ERR)
c
c  Loop over the species geometries, putting info in the appropriate arrays
c  Put info in the global arrays
c
      istart = 1
      do i = 1, nspecies
         uspc(i) = i
         if (.not. geom_create(geom, species(i)))
     &     call errquit('dntmc_getspecgeom: failed creating geometry?',
     &                   0, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom, species(i)))
     &     call errquit('dntmc_getspecgeom: no geometry ', 0, 
     &                   RTDB_ERR)
C Shouldn't need the following since we did it up above
C         if (.not. geom_ncent(geom, natom(i)))
C     &     call errquit('dntmc_getspecgeom: ncent error', 0,
C     &                   GEOM_ERR)
c
c get some local memory to handle temporary data
c
         if (.not.ma_push_get(mt_dbl, 3*natom(i), 'c', l_c, k_c))
     &     call errquit('dntmc_getspecgeom: c allocation error',
     &       3*natom(i), MA_ERR)
         if (.not.ma_push_get(mt_dbl, natom(i), 'q', l_q, k_q))
     &     call errquit('dntmc_getspecgeom: q allocation error',
     &       natom(i), MA_ERR)
         if (.not.ma_push_get(mt_dbl, natom(i), 'm', l_m, k_m))
     &     call errquit('dntmc_getspecgeom: m allocation error',
     &       natom(i), MA_ERR)
         if(.not.ma_push_get(mt_byte,16*natom(i),'t',l_t,k_t))
     &     call errquit('dntmc_getspecgeom: t not allocated',
     &       natom(i), MA_ERR)
c
         if (.not. geom_cart_get(geom, natom(i), byte_mb(k_t), 
     &     dbl_mb(k_c), dbl_mb(k_q)))
     &     call errquit('dntmc_getspecgeom: coord error',0, GEOM_ERR) 
         if (.not. geom_masses_get(geom, natom(i), dbl_mb(k_m)))
     &     call errquit('dntmc_getspecgeom: mass error',0, GEOM_ERR)
c
c put information into dntmc variables
c
         call dntmc_filldntmcinit(i, natom(i), dbl_mb(k_c), 
     &     dbl_mb(k_rinit), dbl_mb(k_m), dbl_mb(k_ms),
     &     nmols, nspecies, natmax)
c
c expand to the full list to be used for the generation of
c the cluster geometry object
c NOTE: istart get changed after the next call!!!
c
c        write(6,*) 'before the dntmc_speciestocluster call'
         call dntmc_speciestocluster(istart, nmolunique(i), natom(i),
     &         nattotal,
     &         dbl_mb(k_c), dbl_mb(k_ct), dbl_mb(k_q), dbl_mb(k_qt),
     &         dbl_mb(k_m), dbl_mb(k_mt), byte_mb(k_t), byte_mb(k_tt),
     &         byte_mb(k_tag))
c        write(6,*) 'after the dntmc_speciestocluster call'
c
c get rid of local memory for next round
c
         if (.not.ma_chop_stack(l_c))
     &     call errquit('dntmc_getspecgeom: chop of l_c did not work',
     &       0, MA_ERR)
c
         if (.not. geom_destroy(geom))
     &     call errquit('dntmc_getspecgeom: destroying geometry',0, 
     &       GEOM_ERR)

      enddo
c
c Create a geometry on the rtdb with the correct dimensions and other
c starting information
c
      call dntmc_setgeom(rtdb, nattotal, dbl_mb(k_ct), dbl_mb(k_qt),
     &      dbl_mb(k_mt), byte_mb(k_tt))
c
c get rid of cluster memory
c
         if (.not.ma_chop_stack(l_ct))
     &     call errquit('dntmc_getspecgeom: chop of l_ct did not work',
     &       0, MA_ERR)
c
      return
      end
C
C Copy the species geometry information to the cluster information
C
      subroutine dntmc_speciestocluster(iloc, nmols, natom, nattotal, 
     &     clocal, ccluster, qlocal, qcluster, mlocal, mcluster,
     &     tlocal, tcluster, tag)

      implicit none

      integer iloc  ! [input/output] index in the cluster coordinates
      integer nmols ! [input] number of molecules to be copied
      integer natom ! [input] number of atoms in the molecule
      integer nattotal ! [input] number of atoms in the cluster
      double precision clocal(3, natom) ! [input] coordinates of species
      double precision ccluster(3, nattotal) ! [input/output] coordinate of cluster
      double precision qlocal(natom) ! [input] charges for the species
      double precision qcluster(nattotal) ! [input/output] charges for the cluster
      double precision mlocal(natom) ! [input] masses of the species
      double precision mcluster(nattotal) ! [input/output] masses of the cluster
      character*16 tlocal(natom) ! [input] tags for the species
      character*16 tcluster(nattotal) ! [input/output] tags for the cluster
      character*16 tag(nattotal) ! [output] tags for cluster
c
      integer i, j
c
      do i = 1, nmols
         do j = 1, natom
            ccluster(1, iloc) = clocal(1, j)
            ccluster(2, iloc) = clocal(2, j)
            ccluster(3, iloc) = clocal(3, j)
            qcluster(iloc) = qlocal(j)
            mcluster(iloc) = mlocal(j)
            tcluster(iloc) = tlocal(j)
            tag(iloc) = tlocal(j)
            iloc = iloc + 1
         enddo
      enddo
c
      return
      end
C
C  Copy the species info to the DNTMC data structures
C
      subroutine dntmc_filldntmcinit(iloc, natom, clocal, 
     &     rinit, mlocal, m, nmols, nspecies, natmax)

      implicit none
#include "gibbsP.fh"

      integer iloc ! [input] index in the DNTMC species
      integer natom ! [input] number of atoms in the molecule
      integer nmols, nspecies, natmax ! to dimension rinit and m
      double precision clocal(3, natom) ! [input] coordinates of species
      double precision rinit(nmols, natmax, 3) ! [input/output] initial dntmc coordinates
      double precision mlocal(natom) ! [input] masses of the species
      double precision m(nspecies, natmax) ! [input/output] masses for dntmc
C Convert from au to Angstroms
      double precision scale
c
      integer i
c
      scale = 0.52917720832D+00
      do i = 1, natom
         rinit(iloc, i, 1) = clocal(1, i)*scale
         rinit(iloc, i, 2) = clocal(2, i)*scale
         rinit(iloc, i, 3) = clocal(3, i)*scale
         m(iloc, i) = mlocal(i)
      enddo
c
      return
      end
c
c set up the geometry the first time, takes the coordinates, labels, charges,
c and masses and sets up a geometry in the rtdb
c
      subroutine dntmc_setgeom(rtdb, natoms, coords, charges, masses,
     &   tags)

      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
      integer rtdb     ! [input] rtdb handle
      integer natoms   ! [input] number of atoms in the full cluster
      double precision coords (3, natoms)   ! [input] coordinates
      double precision charges(natoms)   ! [input] charges
      double precision masses(natoms)   ! [input] masses
      character*16 tags(natoms)        ! [input] atomic labels
c
      integer geom
c
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('dntmc_setgeom: geom_create failed',911,
     &       GEOM_ERR)
      if (.not.geom_set_user_units(geom,'angstroms'))
     &      call errquit
     &       ('dntmc_setgeom: geom_set_user_units failed',911,
     &       GEOM_ERR)
      if (.not.geom_cart_set(geom,natoms,tags,coords,charges))
     &      call errquit('dntmc_setgeom: geom_cart_set failed',911,
     &       GEOM_ERR)
      if (.not.geom_masses_set(geom,natoms,masses))
     &      call errquit('dntmc_setgeom:geom_masses_set failed',911,
     &       GEOM_ERR)
      call geom_compute_values(geom)
c
c store as default geometry
c
      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &      call errquit('dntmc_setgeom: geom_rtdb_store failed',911,
     &       RTDB_ERR)
      if (.not.geom_destroy(geom))
     &      call errquit('dntmc_setgeom: geom_destroy failed',
     &                    911, GEOM_ERR)

      return
      end
c
c set the coordinates assuming that a geometry already is set up in the rtdb
c
      subroutine dntmc_setcoords(rtdb, r, nmols, natom, spc,nspc,natmax)

      implicit none

#include "gibbsP.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"

      integer rtdb ! input
      integer nspc, natmax ! to dimension r and natom
      integer nmols ! input - number of molecules
      double precision r(nmols, natmax, 3) ! input - coordinates in DNTMC ordering
      integer natom(nspc) ! input - number of atoms for each species
      integer spc(nmols) ! input - species numbers
c
      integer geom, ncent
      integer l_ct, k_ct
c
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('dntmc_setcoords: geom_create failed',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &  call errquit('dntmc_setcoords: no geometry ', 555, RTDB_ERR)
      if (.not. geom_ncent(geom, ncent)) 
     &  call errquit('dntmc_setcoords: ncent failed', 555, GEOM_ERR)
c
c reorder the coordinates to that of internal NWChem
c
      if (.not.ma_push_get(mt_dbl, 3*ncent,'coordinate', l_ct, k_ct))
     &   call errquit('dntmc_setcoords: coordinate not allocated',
     &                 3*ncent, MA_ERR)
      call dntmc_dntcoordtonwchem(r, nmols, natom, spc, dbl_mb(k_ct),
     &   ncent, nspc, natmax)
c
c set the coordinates and store
c
      if (.not.geom_cart_coords_set(geom,dbl_mb(k_ct)))
     &      call errquit('dntmc_setcoords: geom_cart_coords_set',911,
     &       GEOM_ERR)
c     if (.not.geom_print(geom)) 
c    &     call errquit('dntmc_setcoords: geom_print failed', 555,
c    &       GEOM_ERR)
      if (.not.ma_chop_stack(l_ct))
     &     call errquit('dntmc_setcoords: chop of l_ct did not work',
     &       0, MA_ERR)
      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &      call errquit('dntmc_setcoords: geom_rtdb_store failed',
     &       911, RTDB_ERR)
      if (.not.geom_destroy(geom))
     &      call errquit('dntmc_setcoords: geom_destroy failed',
     &                    911, GEOM_ERR) 
     
      return
      end
c
c actually set the coordinates
c
      subroutine dntmc_dntcoordtonwchem(r, nmols, natom, spc, 
     &    coords, ncent, nspc, natmax)
c
      implicit none

#include "gibbsP.fh"
#include "errquit.fh"
c
      integer nmols ! input - number of molecules
      integer nspc, natmax !to dimension r and natom
      double precision r(nmols, natmax, 3) ! input - coordinates in DNTMC ordering
      integer natom(nspc) ! input - number of atoms for each species
      integer spc(nmols) ! input - species numbers
      integer ncent ! input - total number of atoms
      double precision coords(3, ncent) ! output - NWChem coordinates
C Convert from Angstroms to au
      double precision scale
c
      integer i, j, nat
c
c      write(6,*) 'inside dntcoordtonwchem and coords are'
      scale = 0.52917720832D+00
      nat = 0
      do i= 1, nmols
         do j= 1, natom(spc(i))
c           write(6,*) 'natom, spc(j), j ', natom(spc(j)), spc(j), j
            nat = nat + 1
            coords(1, nat) = r(i, j, 1)/scale
            coords(2, nat) = r(i, j, 2)/scale
            coords(3, nat) = r(i, j, 3)/scale
c           write(6,*) r(i, j, 1),r(i, j, 2),r(i, j, 3)
         enddo
      enddo
c     write(6,*) 'nmols is ', nmols
c
c little internal check on the number of atoms
c
c     write(6,*) 'nat and ncent ',nat,ncent
      if (nat.ne.ncent)
     &   call errquit('dntmc_dntcoordtonwchem: nat check failed', 
     &                 555, UNKNOWN_ERR)
c
      return
      end
