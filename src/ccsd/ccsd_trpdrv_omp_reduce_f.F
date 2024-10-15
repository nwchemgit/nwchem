      subroutine ccsd_trpdrv_omp_fbody_reduce_new (f1n, f1t, f2n, f2t,
     &                                  f3n, f3t, f4n, f4t,
     &                                  eorb,
     &                                  ncor, nocc, nvir,
     &                                  emp4, emp5,
     &                                  i, k,
     &                                  eaijk,
     &                                  dintc1, dintx1, t1v1,
     &                                  dintc2, dintx2, t1v2)

      use omp_lib
      use iso_fortran_env
      implicit none

      double precision, intent(inout) :: emp4, emp5
      double precision, intent(inout) :: f1n(nvir,nvir), f1t(nvir,nvir)
      double precision, intent(inout) :: f2n(nvir,nvir), f2t(nvir,nvir)
      double precision, intent(inout) :: f3n(nvir,nvir), f3t(nvir,nvir)
      double precision, intent(inout) :: f4n(nvir,nvir), f4t(nvir,nvir)
      double precision, intent(in)    :: eorb(*)
      double precision, intent(in)    :: eaijk
      double precision, intent(in)    :: dintc1(nvir), dintx1(nvir)
      double precision, intent(in)    :: dintc2(nvir), dintx2(nvir)
      double precision, intent(in)    :: t1v1(nvir), t1v2(nvir)
      integer, intent(in) :: ncor, nocc, nvir
      integer, intent(in) :: i, k

      double precision  :: emp4i,emp5i,emp4k,emp5k, denom
      double precision  :: f1nbc,f1tbc,f1ncb,f1tcb
      double precision  :: f2nbc,f2tbc,f2ncb,f2tcb
      double precision  :: f3nbc,f3tbc,f3ncb,f3tcb
      double precision  :: f4nbc,f4tbc,f4ncb,f4tcb
      double precision  :: t1v1b,t1v2b,dintx1c,dintx2c,dintc1c,dintc2c
      integer :: b,c

      emp4i = 0.0
      emp5i = 0.0
      emp4k = 0.0
      emp5k = 0.0

      !$omp target teams distribute parallel do collapse(2)
     &      reduction(+:emp5i,emp4i,emp5k,emp4k)
     &      private(f1nbc,f1tbc,f1ncb,f1tcb,f2nbc,f2tbc,f2ncb,f2tcb)
     &      private(f3nbc,f3tbc,f3ncb,f3tcb,f4nbc,f4tbc,f4ncb,f4tcb)
     &      private(t1v1b,t1v2b,dintx1c,dintx2c,dintc1c,dintc2c)
     &      private(denom) firstprivate(eaijk,nvir,ncor,nocc)
      do b=1,nvir
         do c=1,nvir
            denom=-1.0/( eorb(ncor+nocc+b)+eorb(ncor+nocc+c)+eaijk )

            f1nbc = f1n(b,c);
            f1tbc = f1t(b,c);
            f1ncb = f1n(c,b);
            f1tcb = f1t(c,b);

            f2nbc = f2n(b,c);
            f2tbc = f2t(b,c);
            f2ncb = f2n(c,b);
            f2tcb = f2t(c,b);

            f3nbc = f3n(b,c);
            f3tbc = f3t(b,c);
            f3ncb = f3n(c,b);
            f3tcb = f3t(c,b);

            f4nbc = f4n(b,c);
            f4tbc = f4t(b,c);
            f4ncb = f4n(c,b);
            f4tcb = f4t(c,b);

            t1v1b = t1v1(b);
            t1v2b = t1v2(b);

            dintx1c = dintx1(c);
            dintx2c = dintx2(c);
            dintc1c = dintc1(c);
            dintc2c = dintc2(c);

            emp4i = emp4i
     &             + (denom * (f1tbc+f1ncb+f2tcb+f3nbc+f4ncb)
     &             * (f1tbc-f2tbc*2-f3tbc*2+f4tbc)
     &             - denom * (f1nbc+f1tcb+f2ncb+f3ncb)
     &             * (f1tbc*2-f2tbc-f3tbc+f4tbc*2)
     &             + denom * 3 * (f1nbc*(f1nbc+f3ncb+f4tcb*2)
     &             + f2nbc*f2tcb+f3nbc*f4tbc))

            emp4k = emp4k
     &             + (denom * (f1nbc+f1tcb+f2ncb+f3tbc+f4tcb)
     &             * (f1nbc-f2nbc*2-f3nbc*2+f4nbc)
     &             - denom * (f1tbc+f1ncb+f2tcb+f3tcb)
     &             * (f1nbc*2-f2nbc-f3nbc+f4nbc*2)
     &             + denom * 3 * (f1tbc*(f1tbc+f3tcb+f4ncb*2)
     &             + f2tbc*f2ncb+f3tbc*f4nbc))

            emp5i = emp5i
     &             + (denom * t1v1b * dintx1c
     &             * (f1tbc+f2nbc+f4ncb
     &             - (f3tbc+f4nbc+f2ncb+f1nbc+f2tbc+f3ncb)*2
     &             + (f3nbc+f4tbc+f1ncb)*4)
     &             + denom * t1v1b * dintc1c
     &             * (f1nbc+f4nbc+f1tcb -(f2nbc+f3nbc+f2tcb)*2))

           emp5k = emp5k
     &             + (denom * t1v2b * dintx2c
     &             * (f1nbc+f2tbc+f4tcb
     &             - (f3nbc+f4tbc+f2tcb +f1tbc+f2nbc+f3tcb)*2
     &             + (f3tbc+f4nbc+f1tcb)*4)
     &             + denom * t1v2b * dintc2c
     &             * (f1tbc+f4tbc+f1ncb -(f2tbc+f3tbc+f2ncb)*2))
         enddo
      enddo
      !$omp end target teams distribute parallel do

      emp4 = emp4 + emp4i
      emp5 = emp5 + emp5i
      if (i.ne.k) then
         emp4 = emp4 + emp4k
         emp5 = emp5 + emp5k
      end if ! (i.ne.k)

      end
