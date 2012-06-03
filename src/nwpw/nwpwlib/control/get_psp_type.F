*
* $Id$
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

      integer l
      character*255 full_filename
      character*10 test

*     **** external functions ****
      integer  convert_psp_type
      external convert_psp_type

      call Parallel_taskid(taskid)

*     **** open psp_filename file ****
      if (taskid.eq.MASTER) then
        call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
        l = index(full_filename,' ') - 1
        open(unit=27,file=full_filename(1:l),form='formatted',
     >             status='old')
        read(27,'(A)') test
        close(27)

        psp_type = convert_psp_type(test)
      end if

c     **** send header to all nodes ****
      call Parallel_Brdcst_ivalue(MASTER,psp_type)

      return
      end

*     ***********************************
*     *             			*
*     *          convert_psp_type	*
*     *             			*
*     ***********************************
      integer function convert_psp_type(test)
      implicit none
      character*(*) test

      integer psp_type

      psp_type = 0
      if (test(1:1).eq.'0') psp_type = 0
      if (test(1:1).eq.'1') psp_type = 1
      if (test(1:1).eq.'2') psp_type = 2
      if (test(1:1).eq.'3') psp_type = 3
      if (test(1:1).eq.'4') psp_type = 4
      if (test(1:1).eq.'5') psp_type = 5
      if (test(1:1).eq.'6') psp_type = 6
      if (test(1:1).eq.'7') psp_type = 7
      if (test(1:1).eq.'8') psp_type = 8
      if (test(1:1).eq.'9') psp_type = 9

      convert_psp_type = psp_type
      return
      end
