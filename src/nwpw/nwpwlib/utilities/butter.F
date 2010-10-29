*
* $Id$
*

c     ****************************************
c     *                                      *
c     *          Butter_levels               *
c     *                                      *
c     ****************************************
c
c  This routine computes log(dble(np_j))/log(2.0d0)
c
      integer function Butter_levels(np_j)
      implicit none
      integer np_j

*     *** local variables ***
      integer Level,v

      Level = -1
      v     = 1
      do while (v.le.np_j)
        Level = Level + 1
        v     = 2*v
      end do

      Butter_levels = Level
      return
      end


c     ****************************************
c     *                                      *
c     *          Butter_size                 *
c     *                                      *
c     ****************************************

      integer function Butter_size(taskid_j,np_j,na)
      implicit none
      integer taskid_j,np_j,na(*)

*     *** local variables ***
      integer i,l,pr,size,rsize,Level,shift,pfrom

*     *** external functions ****
      integer  Butter_levels
      external Butter_levels

      !Level = (log(dble(np_j))/log(2.0d0))
      Level = Butter_levels(np_j)

      rsize = na(taskid_j+1)
      do l=0,(Level-1)
         shift = 2**l
         pfrom = mod(taskid_j  +shift,np_j)

         size = shift
         if (l.eq.(Level-1)) then 
           size = ((np_j-2**Level)/2) + 1
         end if

         do i=0,(size-1)
            pr = mod(pfrom    + i,np_j)
            rsize = rsize + na(pr+1)
         end do
      end do

      Butter_size = rsize
      return
      end



c     ****************************************
c     *                                      *
c     *          Butter_indexes              *
c     *                                      *
c     ****************************************

      subroutine Butter_indexes(l,taskid_j,np_j,na,
     >                          istart,iend,jstart,jend,
     >                          imodn,imodtask)
      implicit none
      integer l,taskid_j,np_j,na(*)
      integer istart,iend,jstart,jend
      integer imodn,imodtask

*     *** local variables ***
      integer i,shift,size,rsize,rindx,pr,Level,pfrom

*     **** external functions ****
      integer  Butter_levels
      external Butter_levels

      imodn    = 1
      imodtask = 0

      if (l.eq.0) then
         istart =  1
         iend   =  na(taskid_j+1)
         jstart =  1
         jend   =  na(taskid_j+1)
      else
         !Level = (log(dble(np_j))/log(2.0d0))
         Level = Butter_levels(np_j)
         shift = 2**(l-1)
         pfrom = mod(taskid_j + shift,np_j)

         size  = shift
         if (l.eq.Level) then
            size = ((np_j-2**Level)/2) + 1
         endif

         rsize = 0
         do i=0,size-1
           pr = mod(pfrom + i,np_j)
           rsize = rsize + na(pr+1)
         end do

         rindx = 1
         do i=0,shift-1
           pr    = mod(taskid_j + i,np_j)
           rindx = rindx + na(pr+1)
         end do
         istart =  1
         iend   =  na(taskid_j+1)
         jstart = rindx
         jend   = rindx + rsize - 1

*        *** last level and even ***
         if ((l.eq.Level).and.(mod(np_j,2).eq.0)) then

*           *** first np_j/2 procs do even computations, ***
*           *** second np_j/2 proces do odd computations ***
            imodn    = 2
            imodtask = taskid_j/(np_j/2)
            if (imodtask.eq.1) then
              rindx  = istart
              istart = jstart
              jstart = rindx

              rindx = iend
              iend  = jend
              jend  = rindx
            end if
        end if
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *          Butter_indexes_L1           *
c     *                                      *
c     ****************************************

      subroutine Butter_indexes_L1(taskid_j,np_j,na,
     >                             istart,iend,jstart,jend,
     >                             imodn,imodtask)
      implicit none
      integer taskid_j,np_j,na(*)
      integer istart,iend,jstart,jend
      integer imodn,imodtask

*     *** local variables ***
      integer i,shift,size,rsize,rindx,pr,Level,pfrom,rs2

*     **** external functions ****
      integer  Butter_levels
      external Butter_levels

      imodn    = 1
      imodtask = 0

      Level = Butter_levels(np_j)
      shift = 2**(Level-1)
      pfrom = mod(taskid_j + shift,np_j)

      size = ((np_j-2**Level)/2) + 1

      rsize = 0
      do i=0,size-1
        pr = mod(pfrom + i,np_j)
        rsize = rsize + na(pr+1)
      end do

      pr = mod(pfrom + size-1,np_j)
      rs2 = na(pr+1)

      rindx = 1
      do i=0,shift-1
        pr    = mod(taskid_j + i,np_j)
        rindx = rindx + na(pr+1)
      end do
      istart = 1
      iend   = na(taskid_j+1)
      jstart = rindx
      jend   = rindx + rsize - 1
      jend   = jend - rs2 

      return
      end


c     ****************************************
c     *                                      *
c     *          Butter_indexes_L2           *
c     *                                      *
c     ****************************************

      subroutine Butter_indexes_L2(taskid_j,np_j,na,
     >                             istart,iend,jstart,jend,
     >                             imodn,imodtask)
      implicit none
      integer taskid_j,np_j,na(*)
      integer istart,iend,jstart,jend
      integer imodn,imodtask

*     *** local variables ***
      integer i,shift,size,rsize,rindx,pr,Level,pfrom,rs2

*     **** external functions ****
      integer  Butter_levels
      external Butter_levels

      imodn    = 1
      imodtask = 0

      Level = Butter_levels(np_j)
      shift = 2**(Level-1)
      pfrom = mod(taskid_j + shift,np_j)

      size = ((np_j-2**Level)/2) + 1

      rsize = 0
      do i=0,size-1
        pr = mod(pfrom + i,np_j)
        rsize = rsize + na(pr+1)
      end do

      pr = mod(pfrom + size-1,np_j)
      rs2 = na(pr+1)

      rindx = 1
      do i=0,shift-1
        pr    = mod(taskid_j + i,np_j)
        rindx = rindx + na(pr+1)
      end do
      istart = 1
      iend   = na(taskid_j+1)
c      jstart = rindx
      jend   = rindx + rsize - 1
      jstart = jend - rs2 + 1

*     *** last level and even ***
      if (mod(np_j,2).eq.0) then

*        *** first np_j/2 procs do even computations, ***
*        *** second np_j/2 proces do odd computations ***
         imodn    = 2
         imodtask = taskid_j/(np_j/2)
         if (imodtask.eq.1) then
           rindx  = istart
           istart = jstart
           jstart = rindx

           rindx = iend
           iend  = jend
           jend  = rindx
         end if
      end if

      return
      end

