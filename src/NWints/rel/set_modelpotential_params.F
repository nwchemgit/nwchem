c
      subroutine set_modelpotential_params(rtdb, geom, natoms)
c
      implicit none
c
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "modelpotential.fh"
c
      integer rtdb ! [input] handle to database
      integer geom
      integer natoms
c
      character*2 gelem(natoms)
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c    
      character*32 pname
c
      integer i,j,k
      character*16 tagi      ! tag of atom i
      character*2 symi
      character*16 elemi
      double precision ci(3),chgi
      integer atni
      logical status
c
      logical geom_cent_get
      external geom_cent_get
c
      logical geom_tag_to_element
      external geom_tag_to_element
c
      pname = "set_modelpotential_params"
c
c     map the model potential onto the geometry
      do i = 1,natoms
        if (.not.geom_cent_get(geom,i,tagi,ci,chgi)) 
     &    call errquit(pname//'geom_cent_get failed:i',911, GEOM_ERR)
        status = geom_tag_to_element(tagi,symi,elemi,atni)
        call inp_lcase(symi)
cinit
        do k = 1,50
           gexpo(i,k) = 1d6
           gcoef(i,k) = 0d0
        end do                  ! k = 1,50
        status=.false.
        do j = 1,mpmaxelem  ! in modelpotential.fh
          status=symi.eq.mpelem(j)
          if(status) then
              gelem(i) = mpelem(j)
              do k = 1,50
                 gexpo(i,k) = mpexpo(j,k)
                 gcoef(i,k) = mpcoef(j,k)
              end do ! k = 1,50
              goto 1984
          end if
        end do ! j = 1,mpmaxelem 
 1984   continue
        if(ga_nodeid().eq.0.and.(.not.status)) write(6,'(a,i3,2x,a)') 
     W       ' modelpotential warning: no entry found for atom ',i,tagi
      end do  ! i = 1,natoms 
c
c     set parameters in rtdb
      status = rtdb_parallel(.false.)
      if (.not. rtdb_cput(rtdb,'modelpotential:gelem',natoms,gelem))
     &  call errquit(pname//'unable to set gelem',0,RTDB_ERR)
      if(.not. rtdb_put(rtdb,'modelpotential:gexpo',mt_dbl,natoms*50,
     & gexpo)) call errquit(pname//'unable to set gexpo',0,RTDB_ERR)
      if(.not. rtdb_put(rtdb,'modelpotential:gcoef',mt_dbl,natoms*50,
     & gcoef)) call errquit(pname//'unable to set gcoef',0,RTDB_ERR)
      status = rtdb_parallel(.true.)
c
c     print model potential information
      write(luout,*)
      write(luout,1000)
      write(luout,*)
      write(luout,*) "Reading Model Potential File:"
      write(luout,*)  mpfile(1:inp_strlen(mpfile))
      write(luout,*)
c
 1000 format(' Using Relativistic ZORA Model Potential Approach'/
     & ' C. van Wullen, J. Chem. Phys. 109, 392 (1998)'/
     & ' C. van Wullen, C. Michauk, J. Chem. Phys. 123, 204113 (2005)')
c
c     test set model potential parameters
c      do i = 1,natoms
c        write(luout,*) "gelem:",gelem(i)
c        do j = 1,50
c           write(luout,*) gexpo(i,j),gcoef(i,j)
c        end do
c      end do
c
      return
      end
