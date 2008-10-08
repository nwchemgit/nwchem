#include <stdio.h>
#include <math.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>

void ccsd_t_sse2_(int* p_range_p4, int* p_range_p5, int* p_range_p6,
                  int* p_range_h1, int* p_range_h2, int* p_range_h3,
                  double* denom_p4, double* denom_p5, double* denom_p6,
                  double* denom_h1, double* denom_h2, double* denom_h3,
                  double* singles, double* doubles,
                  double* p_factor, double* p_energy){

   int p4,p5,p6,h1,h2,h3;
   int range_p4,range_p5,range_p6;
   int range_h1,range_h2,range_h3;
   double eps_p4,eps_p5,eps_p6;
   double eps_h1,eps_h2,eps_h3;
   double denom_p,denom_h,denom;
   double val_sing,val_doub;
   double factor,energy,energy0;

//   printf("Top of ccsd_t_sse2_\n");

   factor = *(p_factor);

   energy = *(p_energy);

//   printf(" Energy = %18.8f at the top of C function\n",energy);

   range_p4 = *(p_range_p4);
   range_p5 = *(p_range_p5);
   range_p6 = *(p_range_p6);
   range_h1 = *(p_range_h1);
   range_h2 = *(p_range_h2);
   range_h3 = *(p_range_h3);

/*
   printf("range_p4,range_p5,range_p6 = %d %d %d\n",range_p4,range_p5,range_p6);
   printf("range_h1,range_h2,range_h3 = %d %d %d\n",range_h1,range_h2,range_h3);

   printf("*denom_p4,*denom_p5,*denom_p6 = %f %f %f\n",*denom_p4,*denom_p5,*denom_p6);
   printf("*denom_h1,*denom_h2,*denom_h3 = %f %f %f\n",*denom_h1,*denom_h2,*denom_h3);
*/


/*
   int dummy;

   for ( dummy = 1 ; dummy < (range_p4*range_p5*range_p6*range_h1*range_h2*range_h3+1) ; dummy++){
     if ( 
          ( fabs( *(singles+dummy) ) > 1.0e-8 ) | 
          ( fabs( *(doubles+dummy) ) > 1.0e-8 ) 
        ){
       printf("%8d%18.8f%18.8f\n",dummy,*(singles+dummy),*(doubles+dummy));
     }
   }
*/

   energy = 0.0;

   for ( p4 = 0 ; p4 < range_p4 ; p4++ ){
      eps_p4 = *( denom_p4 + p4 + 1 );
      for ( p5 = 0 ; p5 < range_p5 ; p5++ ){
         eps_p5 = *( denom_p5 + p5 + 1 );
         for ( p6 = 0 ; p6 < range_p6 ; p6++ ){
            eps_p6 = *( denom_p6 + p6 + 1 );

            denom_p = ( eps_p4 + eps_p5 + eps_p6 );

            for ( h1 = 0 ; h1 < range_h1 ; h1++ ){
               eps_h1 = *( denom_h1 + h1 + 1 );
               for ( h2 = 0 ; h2 < range_h2 ; h2++ ){
                  eps_h2 = *( denom_h2 + h2 + 1 );
                  for ( h3 = 0 ; h3 < range_h3 ; h3++ ){
                     eps_h3 = *( denom_h3 + h3 + 1 );

                     denom_h = ( eps_h1 + eps_h2 + eps_h3 );

                     denom = 1.0e0 / ( denom_h - denom_p );

                     val_sing = *(singles++);
                     val_doub = *(doubles++);

                     energy0 = factor * denom * val_doub * ( val_doub + val_sing );

                     energy += energy0;

/*
                     if ( fabs(energy0) > 1.0e-8 ){

                     printf(" %3d %3d %3d %3d %3d %3d",p4+1,p5+1,p6+1,h1+1,h2+1,h3+1);
                     printf(" %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f",eps_p4,eps_p5,eps_p6,eps_h1,eps_h2,eps_h3);
                     printf(" %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",factor,denom_p,denom_h,denom,val_sing,val_doub,energy0);

                     }
*/
                     /* Advance pointers on dbl_mb(k_singles) and dbl_mb(k_doubles) */

/*
                     singles++;
                     doubles++;
*/

                  }
               }
            }
         }
      }
   }

//   printf(" Energy = %18.8f at the bottom of C function\n",energy);

   *(p_energy) = energy;

//   printf(" Energy = %18.8f of the pointer to be returned\n",*(p_energy));

//   printf("Bottom of ccsd_t_sse2_\n");

   fflush(stdout);

}
