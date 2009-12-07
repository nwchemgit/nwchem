*
* $Id: get_psp_type.F,v 1.3 2007-09-24 16:58:10 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           get_psp_type		*
*     *             			*
*     ***********************************

      subroutine get_psp_type(psp_filename,psp_type)
      implicit none 
      character*(*) psp_filename
      integer psp_type

*    *** local variables ***
      logical found
      integer MASTER,taskid
      parameter(MASTER=0)

      integer l,msglen
      character*255 full_filename


      call Parallel_taskid(taskid)

*     **** open psp_filename file ****
      if (taskid.eq.MASTER) then
        call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
        l = index(full_filename,' ') - 1
        open(unit=27,file=full_filename(1:l),form='formatted',
     >             status='old')
        read(27,9090,ERR=9998,END=9998) psp_type
 9090   format(i1)
        go to 9999
 9998   psp_type = 0
 9999   continue

        close(27)
      end if

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,psp_type)

      return
      end

