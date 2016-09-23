c
      subroutine get_modelpotential_params(rtdb, natoms, 
     &         gelem, gexpo, gcoef)
c
      implicit none
c
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
c
      integer rtdb ! [input] handle to database
      integer natoms
      character*2 gelem(natoms)
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c    
      character*32 pname
c
      integer i,j
c
      pname = "get_modelpotential_params"
c
c     retrieve parameters from rtdb
      if (.not. rtdb_cget(rtdb,'modelpotential:gelem',natoms,gelem))
     &  call errquit(pname//'unable to access gelem',0,RTDB_ERR)
      if(.not. rtdb_get(rtdb,'modelpotential:gexpo',mt_dbl,natoms*50,
     & gexpo)) call errquit(pname//'unable to access gexpo',0,RTDB_ERR)
      if(.not. rtdb_get(rtdb,'modelpotential:gcoef',mt_dbl,natoms*50,
     & gcoef)) call errquit(pname//'unable to access gcoef',0,RTDB_ERR)
c
c     test retrieve zora parameters
c      do i = 1,natoms
c        write(luout,*) "gelem:",gelem(i)
c        do j = 1,50
c           write(luout,*) gexpo(i,j),gcoef(i,j)
c        end do
c      end do
c
      return
      end
