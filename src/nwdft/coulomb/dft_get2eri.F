      subroutine dft_get2eri(CD_bas_han, g_2ceri,oskel)
      implicit none
      integer CD_bas_han, g_2ceri(*)
      logical oskel
      external int_2e2c
      call ga_get2eri(cd_bas_han, g_2ceri,oskel,int_2e2c)
      return
      end
