
*     ***************************
*     *			        *
*     *		molecular_init  *
*     *				*
*     ***************************
      subroutine molecular_init(rtdb)
*
* $Id: molecular.F,v 1.2 2001-10-17 23:48:11 bylaska Exp $
*
      implicit none 
      integer   rtdb

#include "rtdb.fh"

      character*30 filename
      character*255 full_filename
      integer       l

      logical value
      common /pspw_molecular/ value


      value = rtdb_cget(rtdb,'pspw:molecular_filename',1,filename)
      if (value) then

         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)
          l = index(full_filename,' ') - 1
          write(*,*) "Molecular Constraint File:",full_filename(1:l)
          write(*,*)
          call pspw_molecule_read(full_filename,l)
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
      double precision r2(3,*),r1(3,*)

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
      double precision ion_amass
      external         ion_amass

      if (value) then     

         call pspw_molecule_msize(msize)
         do m=1,msize

            call pspw_molecule_data(m,na,indxa,blist,cyclic)

*           **** set mass variable ****
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
