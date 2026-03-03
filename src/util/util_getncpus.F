      integer function util_getncpus()
      USE ISO_C_BINDING
      implicit none
      interface
        subroutine utilc_getncpus(ncpus) bind(C)
          use, intrinsic                  :: ISO_C_BINDING
          integer (C_INT)                 :: ncpus
        end subroutine
      end interface
      integer (C_INT) :: ncpus
      call utilc_getncpus(ncpus)
      util_getncpus=ncpus
      return
      end
