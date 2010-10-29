      subroutine hnd_cos_get_dens_dft(rtdb,geom_handle,basis_handle,
     1                               da,db,d1nw,d2nw,
     2                               scftyp,g_movecs)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
#include "cdft.fh"
      integer       rtdb
      integer       geom_handle
      integer       basis_handle
c
      character*8   scftyp
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_movecs
      dimension g_movecs(2)

      integer g_dens      ! density matrices (up to 3)
      dimension g_dens(3) ! max # of density matrices
 
      integer ndens
      integer nmo(2)
      integer noc(2)

c
      integer ir, iw, ip
      common/hnd_iofile/ir,iw,ip

      integer nap, iap
      common/hnd_lcapid/nap,iap
c
      double precision da(*),db(*),d1nw(*),d2nw(*)

      integer i
c
      logical dbug
      dbug=.false.
      dbug=dbug.and.(iap.eq.1)
c
c     ----- ... scftyp ... -----
c
      if (ipol.eq.2) then
         scftyp  = 'UHF     '
      else
         scftyp  = 'RHF     '
      endif
c
      ndens = 3
c
c     handles for density matrices
c
      do i=1, ndens
        g_dens(i) = ga_create_atom_blocked ( geom_handle, 
     $              basis_handle,'density matrix' )
      enddo
c

      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     $   call errquit('hnd_cos_get_dens_dft: rtdb_get of noc failed', 0,
     &       RTDB_ERR)

      nbf   =nbf_ao
      nmo(1)=noc(1)
      nmo(2)=noc(2)
c
      if(dbug) then
         write(iw,*) 'in hnd_cos_dens_read_dft'
         write(iw,*) 'scftyp  = ',scftyp 
         write(iw,*) 'nbf     = ',nbf
         write(iw,*) 'nmo     = ',nmo
         write(iw,*) 'ndens   = ',ndens
      endif
c
c     ----- make density matrices -----
c
      call hnd_cos_dens_make_dft(rtdb, g_movecs, g_dens, ndens, 
     $                           da, db, d1nw, d2nw )
c
      return
      end
      subroutine hnd_cos_dens_make_dft(rtdb, g_movecs, g_dens, ndens,
     $                                 da, db, d1nw, d2nw )
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cdft.fh"

      integer rtdb
      logical status

      integer g_dens, ndens
      dimension g_dens(ndens)

      integer  g_movecs(2)

      integer  noc(2)
      integer  nmo(2)
 
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer ir, iw, ip
      common/hnd_iofile/ir,iw,ip

      integer nap, iap
      common/hnd_lcapid/nap,iap

      double precision da, db, d1nw, d2nw
      dimension da(*), db(*), d1nw(*), d2nw(*)

c     local variables

      integer nclosed, nopen
      integer nalpha,  nbeta

      integer i, j, ij

c     calculate the density matrices from the g_movecs

      logical dbug
      dbug=.false.
      dbug=dbug.and.(iap.eq.1)

      if(dbug) then
         write(iw,*) 'in hnd_cos_dens_make_dft'
      endif
c
c     ----- basis set info -----
c
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     $   call errquit('hnd_cos_dens_make_dft: rtdb_get of noc failed',0,
     &       RTDB_ERR)

      nbf   =nbf_ao
      nmo(1)=noc(1)
      nmo(2)=noc(2)

c       density matrices ...

      if (ipol.eq.1) then
c
c     ----- rhf -----
c

        nclosed=noc(1)
        nopen  =0
        if(dbug) then
           write(iw,*) ' nbf, nclosed, nopen = ',nbf,nclosed,nopen
        endif

        if(nclosed.eq.0) then
           call ga_zero(g_dens(1))
        else
           call ga_matmul_patch('n', 't', 2.0D0, 0.0d0  ,
     $          g_movecs(1), 1, nbf,            1, nclosed,
     $          g_movecs(1), 1, nclosed,        1, nbf    ,
     $          g_dens(1)  , 1, nbf,            1, nbf    )
           if(nopen.gt.0) then
              call ga_matmul_patch('n', 't', 1.0D0, 0.0d0       ,
     $             g_movecs(1), 1, nbf, nclosed+1, nclosed+nopen,
     $             g_movecs(1), nclosed+1, nclosed+nopen, 1, nbf,
     $             g_dens(2)  , 1, nbf,            1, nbf       )
              call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2),
     $                                              g_dens(3) )
              call ga_dadd ( 1.D0, g_dens(3), 0.D0, g_dens(2),
     $                                              g_dens(1) )
           endif
        endif

        call ga_symmetrize(g_dens(1))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(da,nbf)
        endif

      else  ! UHF
c
c     ----- uhf -----
c

        nalpha=noc(1)
        nbeta =noc(2)
        if(dbug) then
           write(iw,*) 'nbf, nalpha, nbeta = ',nbf,nalpha,nbeta
        endif

        call ga_zero(g_dens(1))
        call ga_zero(g_dens(2))

        call ga_dgemm('n', 't', nbf, nbf, nalpha, 1.D0, 
     $                g_movecs(1), g_movecs(1), 0.D0, g_dens(1))
        if (nbeta .gt. 0) then
           call ga_dgemm('n', 't', nbf, nbf, nbeta, 1.D0, 
     $                   g_movecs(2), g_movecs(2), 0.D0, g_dens(2))
        endif

C       density total

        call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2), g_dens(3) )
        call ga_symmetrize(g_dens(1))
        call ga_symmetrize(g_dens(2))
        call ga_symmetrize(g_dens(3))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(da,nbf)
        endif

        call ga_get(g_dens(2),1,nbf,1,nbf,d1nw,nbf)
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              db(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(db,nbf)
        endif

c
      endif

c
      do i=1,3
         status = ga_destroy ( g_dens(i) )
      enddo
      return
      end
