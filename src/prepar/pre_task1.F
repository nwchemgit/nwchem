      logical function pre_task1(irtdb)
c
      implicit none
c
      logical scf,esp
      external scf,esp
c
      integer irtdb
c
      pre_task1=.false.
c
      if(scf(irtdb)) then
      if(esp(irtdb)) then
      pre_task1=.true.
      endif
      endif
c
      return
      end
