      subroutine intd_1eh1(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lH1a,H1a)
C $Id: intd_1eh1.F,v 1.2 1994-06-09 09:08:42 d3e129 Exp $
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
c layer routine to compute the derivative 1 electron hamiltonian integrals
c for shells/contractions ish,jsh
c
c Order is...   nint*3*nat (3=> xyz, nat=number of atoms)
c
c  /                   |   
c | nint,d <ij>        |
c |      --------------|
c  \     d[idatom(1),x]|
c                          |
c       nint,d <ij>        |
c            --------------|
c            d[idatom(1),y]|
c                              |
c           nint,d <ij>        |
c                --------------|
c                d[idatom(1),z]|
c                                  |
c               nint,d <ij>        |
c                    --------------|
c                    d[idatom(2),x]|
c                                      |
c                   nint,d <ij>        |
c                        --------------|
c                        d[idatom(2),y]|
c                                           |
c                       nint,d <ij>         |
c                            -------------- |
c                            d[idatom(2),z] |
c
c                                  . . . 
c                                                            |
c                                         nint,d <ij>        |
c                                              --------------|
c                                            d[idatom(nat),x]|
c                                                                |
c                                             nint,d <ij>        |
c                                                  --------------|
c                                                d[idatom(nat),y]|
c                                                                    \
c                                                 nint,d <ij>         |
c                                                      -------------- |
c                                                    d[idatom(nat),z]/
c      
c::functions
      integer int_nint
      external int_nint
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer lH1a      ! [input] number of h1 integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(*)    ! [output] derivative integrals
c
c::local
      integer nint, offset, scrsize, nat
c
      nat = ncenter(ibs_geom((i_basis + BASIS_HANDLE_OFFSET)))
c
      nint = int_nint(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*nat.gt.lH1a) then
        write(6,*) 'nint*3*nat = ',nint*3*nat
        write(6,*) 'lH1a       = ',lH1a
        call errquit('intd_1eh1: nint>lH1a error',911)
      endif
c
      offset = nint*3*2       ! scratch for Ta array in intd_1eh1P routine
      scrsize = lscr - offset ! new scratch array size
      offset = offset + 1     ! increment for passing to intd_1eh1P
c
      call intd_1eh1P(i_basis,ish,j_basis,jsh,
     &       scrsize,scr(offset),nint,H1a,scr)
c
      end
      subroutine intd_1eh1P(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,H1a,Ta)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::external subroutines used
c... errquit
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer nint      ! [input] number of integrals in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(nint,3,*)    ! [output] derivative integrals (nint,3,n_atoms)
      double precision Ta(nint,3,2)     ! [scratch] space for kinetic integrals
c::local
      logical doT
      integer ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      ucont = ibs_cn2ucn(ish,ibas)
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP,ucont,ibas)
      itype = infbs_cont(CONT_TYPE,ucont,ibas) 
      igeom = ibs_geom(ibas)
      iatom = ibs_cn2ce(ish,ibas)
c
      ucont = ibs_cn2ucn(jsh,jbas)
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP,ucont,jbas)
      jtype = infbs_cont(CONT_TYPE,ucont,jbas) 
      jgeom = ibs_geom(jbas)
      jatom = ibs_cn2ce(jsh,jbas)
c
      if (igeom.ne.jgeom) then
        write(6,*)'intd_1eh1P.F: two different geometries for',
     &         ' derivatives??'
        call errquit('intd_1eh1P: geom error ',911)
      endif
c
      if (iatom.eq.jatom) then
        doT = .false.
      else
        doT = .true.
      endif
        
      call hf1d(
     &       coords(1,iatom,igeom),
     &       exndcf(iexp,ibas),
     &       exndcf(icf,ibas),
     &       inp,igen,itype,iatom,
c
     &       coords(1,jatom,jgeom),
     &       exndcf(jexp,jbas),
     &       exndcf(jcf,jbas),
     &       jnp,jgen,jtype,jatom,
c
     &       coords(1,1,igeom),charge(1,igeom),ncenter(igeom),
     &       scr,Ta,H1a,nint,
c............overlap, k-e,  pot-e,  canab,   dryrun
     &       .false., doT, .true., .false., .false.,
     &       scr,lscr)
c
c if needed add in Ta derivative integrals
      if (doT) then
        call daxpy(nint*3,1.0d00,Ta(1,1,1),1,H1a(1,1,iatom),1)
        call daxpy(nint*3,1.0d00,Ta(1,1,2),1,H1a(1,1,jatom),1)
      endif
      end
      
