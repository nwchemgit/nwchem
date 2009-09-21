
*     ***************************
*     *			        *
*     *		molecular_init  *
*     *				*
*     ***************************
      subroutine molecular_init(rtdb)
*
* $Id: molecular.F,v 1.8 2009-01-05 17:46:30 bylaska Exp $
*
      implicit none 
      integer   rtdb

#include "rtdb.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical found
      integer       i,j,l
      integer a,b,m,msize
      integer na,nb,cyclic,indxa(20)
      double precision blist(20)

      character*30 filename
      character*255 full_filename,string


      logical value
      common /pspw_molecular/ value

      call Parallel_taskid(taskid)

      !**** SHAKE data file exists ****
      value = rtdb_cget(rtdb,'pspw:molecular_filename',1,filename)
      if (value) then

         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)
          l = index(full_filename,' ') - 1
          if (taskid.eq.MASTER) then
          write(*,*) "Molecular Constraint File:",full_filename(1:l)
          write(*,*)
          end if
          call pspw_molecule_read(full_filename,l)


      !**** input via rtdb ****
      else
        value = rtdb_cget(rtdb,'nwpw:shake_constraint',1,string)
        if (value) then
          full_filename = 'SHAKE'
          call util_file_name_resolve(full_filename, .false.)
c          call util_file_name_noprefix('SHAKE',.false.,
c     >                                .false.,
c     >                                full_filename)
          l = index(full_filename,' ') - 1

          !*** create SHAKE data file ****
          if (taskid.eq.MASTER) then

          inquire(file=full_filename,exist=found)         !*** delete shake file ***
          if (found) call util_file_unlink(full_filename)
         
          !*** create new shake file ***
          open(unit=49,file=full_filename,form='formatted')
          write(49,*) string
          close(49)

          write(*,*) "Molecular Constraint File:",full_filename(1:l)
          write(*,*)
          end if

          call ga_sync() 
          call pspw_molecule_read(full_filename,l)
        end if
      end if

      if (value) then
         call pspw_molecule_msize(msize)
         do m=1,msize
            call pspw_molecule_data(m,na,indxa,blist,cyclic)
            if (taskid.eq.MASTER) then
            if (cyclic.eq.1) then
               write(*,'(A,A,10(A,I3,I3,A,F6.3,A))') 
     >         "Cyclic Constraint:",
     >         "  fixed bonds -", 
     >     (" ( atoms:",indxa(i),indxa(i+1)," dist:",blist(i)," )", 
     >          i=1,(na-1)),
     >      " ( atoms:",indxa(na),indxa(1)," dist:",blist(na)," )" 
            else
               write(*,'(A,A,10(A,I3,I3,A,F6.3,A))') 
     >         "Linear Constraint:",
     >         "  fixed bonds -", 
     >     (" ( atoms:",indxa(i),indxa(i+1)," dist:",blist(i)," )", 
     >          i=1,(na-1))
            end if
            end if
         end do

      end if


      return
      end

*     ***************************
*     *				*
*     *		molecular_end	*
*     *				*
*     ***************************
      subroutine molecular_end()
      implicit none 

      logical value
      common /pspw_molecular/ value

      if (value) then
         call pspw_molecule_end()
      end if

      return
      end 


*     ***************************
*     *				*
*     *	    molecular_fix	*
*     *	 		        *
*     ***************************
      subroutine molecular_fix(r2,r1)
      implicit none
      double precision r2(*),r1(*)

      logical value
      common /pspw_molecular/ value

*     **** local variables ****
      integer maxit
      double precision tol
      parameter (maxit=1000,tol=1.0d-4)

      integer a,b,m,msize
      integer na,nb,cyclic,indxa(20)
      double precision blist(20),mass(20)

*     **** external functions ****
      logical  pspw_qmmm_shake_found
      external pspw_qmmm_shake_found
      double precision ion_amass
      external         ion_amass

      if (pspw_qmmm_shake_found()) call pspw_qmmm_shake(r2,r1)
      if (value) then     

         call pspw_molecule_msize(msize)
         do m=1,msize

            call pspw_molecule_data(m,na,indxa,blist,cyclic)

*           **** set mass variable ****
*           **** should be rewritten since mass is not katm anymore****
            do a=1,na
                mass(a) = ion_amass(indxa(a))
            end do

*           **** set bond length constraints ****
            nb = na-1
            if (cyclic.eq.1) nb = nb+1
            do b=1,nb
               blist(b) = blist(b)**2
            end do
            call shake_chain(na,indxa,nb,
     >                       tol,maxit,
     >                       blist,mass,
     >                       r2,r1)
         end do


      end if

      return
      end 

*     ***********************************
*     *				        *
*     *	    molecular_force_average	*
*     *	 		                *
*     ***********************************
      subroutine molecular_force_average(filename)
      implicit none
      character*(*) filename

#include "mafdecls.fh"
#include "errquit.fh"
      

      logical value
      common /pspw_molecular/ value

*     **** local variables ****
      logical ok
      integer maxit,nion
      double precision tol
      parameter (maxit=1000,tol=1.0d-4)

      integer a,b,m,msize,indx,indx2,indx3,ii,iii,nconst,nframe,l
      integer na,nb,cyclic,indxa(20)
      double precision blist(20),mass(20),E
      double precision a1x,a1y,a1z
      double precision a2x,a2y,a2z
      double precision a3x,a3y,a3z
      integer r1(2),f1(2),af1(2),af2(2),ac1(2)
      character*2 tag2
      character*4 tag4

      if (value) then
         l = index(filename,' ') - 1
         write(*,*)
         write(*,*) "Force Energy File:",filename(1:l)
         open(unit=35,file=filename,form='formatted',status='old')
         read(35,*) nion
         close(unit=35)
         call pspw_molecule_msize(msize)
         nconst = 0
         do m=1,msize
            call pspw_molecule_data(m,na,indxa,blist,cyclic)
            nb = na-1
            if (cyclic.eq.1) nb = nb+1
            nconst = nconst + nb
         end do

         ok =        MA_push_get(mt_dbl,3*nion,'r1',r1(2),r1(1))
         ok = ok.and.MA_push_get(mt_dbl,3*nion,'f1',f1(2),f1(1))
         ok = ok.and.MA_push_get(mt_dbl,nconst,'af1',af1(2),af1(1))
         ok = ok.and.MA_push_get(mt_dbl,nconst,'af2',af2(2),af2(1))
         ok = ok.and.MA_push_get(mt_dbl,nconst,'ac1',ac1(2),ac1(1))
         if (.not.ok)
     >    call errquit('molecular_force_average: out of stack',1,MA_ERR)
         call dcopy(nconst,0.0d0,0,dbl_mb(af1(1)),1)
         call dcopy(nconst,0.0d0,0,dbl_mb(af2(1)),1)
         call dcopy(nconst,0.0d0,0,dbl_mb(ac1(1)),1)

         open(unit=35,file=filename,form='formatted',status='old')
         nframe = 0
         do while (.true.)
            read(35,*,ERR=30,END=30) nion
            read(35,*,ERR=30,END=30) E
            read(35,*,ERR=30,END=30) a1x,a1y,a1z
            read(35,*,ERR=30,END=30) a2x,a2y,a2z
            read(35,*,ERR=30,END=30) a3x,a3y,a3z
            indx = 0
            do ii=1,nion
              read(35,*,ERR=30,END=30) iii,tag2,tag4,
     >                                 dbl_mb(r1(1)+indx),
     >                                 dbl_mb(r1(1)+indx+1),
     >                                 dbl_mb(r1(1)+indx+2),
     >                                 dbl_mb(f1(1)+indx),
     >                                 dbl_mb(f1(1)+indx+1),
     >                                 dbl_mb(f1(1)+indx+2)
              indx = indx + 3
            end do

            nframe = nframe + 1
            call molecular_sub_force_average(dbl_mb(r1(1)),
     >                                       dbl_mb(f1(1)),
     >                                       dbl_mb(af1(1)),
     >                                       dbl_mb(af2(1)),
     >                                       dbl_mb(ac1(1)))
         end do

  30     continue
         close(unit=35)

         if (nframe.eq.0) then
            E = 0.0d0
         else
            E = 1.0d0/dble(nframe)
         end if
         call dscal(nconst,E,dbl_mb(af1(1)),1)
         call dscal(nconst,E,dbl_mb(af2(1)),1)
         call dscal(nconst,E,dbl_mb(ac1(1)),1)

         write(*,110) nframe
         write(*,111) nconst
         write(*,112) nion
  110    format(1x,'frames used           =',I8)
  111    format(1x,'number of constraints =',I8)
  112    format(1x,'number of atoms       =',I8)
         indx  = af1(1)
         indx3 = af2(1)
         indx2 = ac1(1)
         do m=1,msize
            call pspw_molecule_data(m,na,indxa,blist,cyclic)
            nb = na-1
            if (cyclic.eq.1) nb = nb+1
            if (cyclic.eq.1) then
              do b=1,nb-1
                write(*,120) indxa(b),indxa(b+1),blist(b),
     >                       dbl_mb(indx2),dbl_mb(indx),
     >                       dbl_mb(indx3)-dbl_mb(indx)**2
                indx  = indx + 1
                indx2 = indx2 + 1
                indx3 = indx3 + 1
              end do
              write(*,120) indxa(nb),indxa(1),blist(nb),
     >                     dbl_mb(indx2),dbl_mb(indx),
     >                       dbl_mb(indx3)-dbl_mb(indx)**2
                indx  = indx + 1
              indx  = indx + 1
              indx2 = indx2 + 1
              indx3 = indx3 + 1
            else
              do b=1,nb
                write(*,120) indxa(b),indxa(b+1),blist(b),
     >                       dbl_mb(indx2),dbl_mb(indx),
     >                       dbl_mb(indx3)-dbl_mb(indx)**2
                indx  = indx + 1
                indx  = indx + 1
                indx2 = indx2 + 1
                indx3 = indx3 + 1
              end do
            end if
         end do
  120    format(1x,'atoms=',I5,I5,
     >          1x,'dist=',E10.3,
     >          1x,'(shake error=',E10.3,')',
     >          4x,'<F>=',E14.6,
     >          1x,'(<F**2>-<F>**2=',E10.3,')')

         ok =        MA_pop_stack(ac1(2))
         ok = ok.and.MA_pop_stack(af2(2))
         ok = ok.and.MA_pop_stack(af1(2))
         ok = ok.and.MA_pop_stack(f1(2))
         ok = ok.and.MA_pop_stack(r1(2))
         if (.not.ok)
     >    call errquit(
     >         'molecular_force_average: error popping stack',1,MA_ERR)

      end if

      return
      end

      subroutine molecular_sub_force_average(rion,fion,af,af2,ac)
      implicit none
      real*8 rion(3,*)
      real*8 fion(3,*)
      real*8 af(*)
      real*8 af2(*)
      real*8 ac(*)

*     **** local variables
      integer b,m,msize,na,indxa(20),cyclic,indx
      real*8 x,y,z,f1,f2,ff,blist(20),rr

      indx = 1
      call pspw_molecule_msize(msize)
      do m=1,msize
         call pspw_molecule_data(m,na,indxa,blist,cyclic)
         do b=1,na-1
            x = rion(1,indxa(b+1))-rion(1,indxa(b))
            y = rion(2,indxa(b+1))-rion(2,indxa(b))
            z = rion(3,indxa(b+1))-rion(3,indxa(b))
            rr = dsqrt(x*x + y*y + z*z)

            f1 = x*fion(1,indxa(b)) 
     >         + y*fion(2,indxa(b)) 
     >         + z*fion(3,indxa(b))
            f2 = x*fion(1,indxa(b+1)) 
     >         + y*fion(2,indxa(b+1)) 
     >         + z*fion(3,indxa(b+1))
            ff = 0.5*(f2-f1)/rr

            af(indx)  = af(indx)  + ff
            af2(indx) = af2(indx) + ff*ff

            ac(indx) = ac(indx) + dabs(rr-blist(b))

            indx = indx + 1
         end do
         if (cyclic.eq.1) then
            x = rion(1,indxa(na))-rion(1,indxa(1))
            y = rion(2,indxa(na))-rion(2,indxa(1))
            z = rion(3,indxa(na))-rion(3,indxa(1))
            f1 = x*fion(1,indxa(1)) 
     >         + y*fion(2,indxa(1)) 
     >         + z*fion(3,indxa(1))
            f2 = x*fion(1,indxa(na)) 
     >         + y*fion(2,indxa(na)) 
     >         + z*fion(3,indxa(na))
            ff = 0.5*(f2-f1)
            af(indx)  = af(indx)  + ff
            af2(indx) = af2(indx) + ff*ff

            rr = dsqrt(x*x + y*y + z*z)
            ac(indx) = ac(indx) + dabs(rr-blist(na))

            indx = indx + 1
         end if
      end do

      return
      end
