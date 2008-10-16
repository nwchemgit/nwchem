#ifdef SSE2_VERSION

#ifdef DEBUG_PRINT
  #include <stdio.h>
#endif

#include <math.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#define linesize 64

void ccsd_t_sse2_(int *p_range_p4, int *p_range_p5, int *p_range_p6,
                  int *p_range_h1, int *p_range_h2, int *p_range_h3,
                  double *denom_p4, double *denom_p5, double *denom_p6,
                  double *denom_h1, double *denom_h2, double *denom_h3,
                  double *singles,  double *doubles,
                  double *p_factor, double *p_energy){

   static unsigned int range_p4,range_p5,range_p6,range_h1,range_h2,range_h3;
   void ccsd_t_sse2_unaligned_();

   printf("Top of ccsd_t_sse2_\n");

   range_p4 = *(p_range_p4);
   range_p5 = *(p_range_p5);
   range_p6 = *(p_range_p6);
   range_h1 = *(p_range_h1);
   range_h2 = *(p_range_h2);
   range_h3 = *(p_range_h3);

   // check if aligned
   //if ( range_h3 % (linesize/sizeof(double)) == 0 && (long)singles % 64 == 0 && (long)doubles % 64 == 0 ) {
/*   if ( range_h3 % 2 ) {
      //transpose_aligned, SSE2 used
      ccsd_t_sse2_aligned_(range_p4, range_p5, range_p6, range_h1, range_h2, range_h3,
                           denom_p4, denom_p5, denom_p6, denom_h1, denom_h2, denom_h3,
                           singles,  doubles,  p_factor, p_energy);
   } else { */
      // transpose_misaligned
      ccsd_t_sse2_unaligned_(range_p4, range_p5, range_p6, range_h1, range_h2, range_h3,
                             denom_p4, denom_p5, denom_p6, denom_h1, denom_h2, denom_h3,
                             singles,  doubles,  p_factor, p_energy);
   //}

   printf("Bottom of ccsd_t_sse2_\n");

}

void ccsd_t_sse2_aligned_(int range_p4, int range_p5, int range_p6,
                          int range_h1, int range_h2, int range_h3,
                          double *denom_p4, double *denom_p5, double *denom_p6,
                          double *denom_h1, double *denom_h2, double *denom_h3,
                          double *singles,  double *doubles,
                          double *p_factor, double *p_energy){

   unsigned int p4,p5,p6,h1,h2,h3;
   double eps_p4,eps_p5,eps_p6,eps_h1,eps_h2,eps_h3;
   double energy0;
   double energy1;
   double denom_ph12;

   register __m128d v_factor;
   register __m128d v_energy;
   register __m128d v_denom_p4;
   register __m128d v_denom_p5;
   register __m128d v_denom_p45;
   register __m128d v_denom_p6;
   register __m128d v_denom_p456;
   register __m128d v_denom_h1;
   register __m128d v_denom_ph1;
   register __m128d v_denom_h2;
   register __m128d v_denom_ph12;
   register __m128d v_denom_h3;
   register __m128d v_denom;
   register __m128d v_both;
   register __m128d v_doub;

   printf("Top of ccsd_t_sse2_aligned_\n");

   v_factor = _mm_load1_pd(p_factor);

   _mm_storel_pd(&energy0,v_factor);
   _mm_storeh_pd(&energy1,v_factor);
   printf("v_factor_l = %f\n",energy0);
   printf("v_factor_h = %f\n",energy1);

   v_energy = _mm_setzero_pd();

   _mm_storel_pd(&energy0,v_energy);
   _mm_storeh_pd(&energy1,v_energy);
   printf("v_energy_l = %f\n",energy0);
   printf("v_energy_h = %f\n",energy1);

   printf("Before loops\n");

/*
   for ( p4 = 1 ; p4 <= range_p4 ; p4++ ){
      v_denom_p4 = _mm_load1_pd(denom_p4+p4);
      for ( p5 = 1 ; p5 <= range_p5 ; p5++ ){
         v_denom_p5 = _mm_load1_pd(denom_p5+p5);
         v_denom_p45 = _mm_add_pd(v_denom_p4,v_denom_p5); // eps_p4 + eps_p5
         for ( p6 = 1 ; p6 <= range_p6 ; p6++ ){
            v_denom_p6 = _mm_load1_pd(denom_p6+p6);
            v_denom_p456 = _mm_add_pd(v_denom_p45,v_denom_p6); // eps_p4 + eps_p5 + eps_p6
*/
   for ( p4 = 0 ; p4 < range_p4 ; p4++ ){
      eps_p4 = *( denom_p4 + p4 + 1 );
      for ( p5 = 0 ; p5 < range_p5 ; p5++ ){
         eps_p5 = *( denom_p5 + p5 + 1 );
         for ( p6 = 0 ; p6 < range_p6 ; p6++ ){
            eps_p6 = *( denom_p6 + p6 + 1 );
/*
            for ( h1 = 1 ; h1 <= range_h1 ; h1++ ){
               v_denom_h1 = _mm_load1_pd(denom_h1+h1);
               v_denom_ph1 = _mm_sub_pd(v_denom_h1,v_denom_p456); // eps_h1 - ( eps_p4 + eps_p5 + eps_p6 )
               for ( h2 = 1 ; h2 <= range_h2 ; h2++ ){
                  v_denom_h2 = _mm_load1_pd(denom_h2+h2);
                  v_denom_ph12 = _mm_add_pd(v_denom_ph1,v_denom_h2); // eps_h1 + eps_h2 - ( eps_p4 + eps_p5 + eps_p6 )
*/
            for ( h1 = 0 ; h1 < range_h1 ; h1++ ){
               eps_h1 = *( denom_h1 + h1 + 1 );
               for ( h2 = 0 ; h2 < range_h2 ; h2++ ){
                  eps_h2 = *( denom_h2 + h2 + 1 );
                  denom_ph12 = ( eps_h1 + eps_h2 ) - ( eps_p4 + eps_p5 + eps_p6 );
                  v_denom_ph12 = _mm_load1_pd(&denom_ph12);

/*
   _mm_storel_pd(&energy0,v_denom_ph12);
   _mm_storeh_pd(&energy1,v_denom_ph12);
   printf("v_denom_ph12_l & denom_ph12 = %f %f\n",energy0,denom_ph12);
   printf("v_denom_ph12_h & denom_ph12 = %f %f\n",energy1,denom_ph12);
*/

                  for ( h3 = 1 ; h3 <= range_h3 ; h3+=2 ){
                     v_denom_h3 = _mm_loadu_pd(denom_h3+h3);

                     printf("Indices: %3d %3d %3d %3d %3d %3d\n",p4,p5,p6,h1,h2,h3);
/*
   _mm_storel_pd(&energy0,v_denom_h3);
   _mm_storeh_pd(&energy1,v_denom_h3);
   printf("v_denom_h3_l & *(denom_h3+h3+0) = %f %f\n",energy0,*(denom_h3+h3+0));
   printf("v_denom_h3_h & *(denom_h3+h3+1) = %f %f\n",energy1,*(denom_h3+h3+1));
*/
                     v_denom = _mm_add_pd(v_denom_ph12,v_denom_h3); /* eps_h1 + eps_h2 + eps_h3 - ( eps_p4 + eps_p5 +eps_p6 ) */
/*
   _mm_storel_pd(&energy0,v_denom);
   _mm_storeh_pd(&energy1,v_denom);
   printf("v_denom_l & denom+0) = %f %f\n",energy0,denom_ph12 + *(denom_h3+h3+0));
   printf("v_denom_h & denom+1) = %f %f\n",energy1,denom_ph12 + *(denom_h3+h3+1));
*/
                     v_doub = _mm_loadu_pd(doubles);       /* load doubles into v_doub                             */
                     v_both = _mm_loadu_pd(singles);       /* load singles into v_both (v_sing)                    */

   _mm_storel_pd(&energy0,v_doub);
   _mm_storeh_pd(&energy1,v_doub);
   printf("v_doub_l & *(doubles+0) = %f %f\n",energy0,*(doubles+0));
   printf("v_doub_h & *(doubles+1) = %f %f\n",energy1,*(doubles+1));

   _mm_storel_pd(&energy0,v_both);
   _mm_storeh_pd(&energy1,v_both);
   printf("v_both_l & *(singles+0) = %f %f\n",energy0,*(singles+0));
   printf("v_both_h & *(singles+1) = %f %f\n",energy1,*(singles+1));

                     v_both = _mm_add_pd(v_both,v_doub);   /* v_both = v_sing + v_doub                             */
/*
   _mm_storel_pd(&energy0,v_both);
   _mm_storeh_pd(&energy1,v_both);
   printf("v_both_l & *(singles+0) + *(doubles+0) = %f %f\n",energy0,*(singles+0)+*(doubles+0));
   printf("v_both_h & *(singles+1) + *(doubles+1) = %f %f\n",energy1,*(singles+1)+*(doubles+1));
*/
                     v_doub = _mm_mul_pd(v_doub,v_factor); /* v_doub = v_factor * v_doub                           */
/*
   _mm_storel_pd(&energy0,v_doub);
   _mm_storeh_pd(&energy1,v_doub);
   printf("v_doub_l & *(p_factor) * *(doubles+0) = %f %f\n",energy0,*(p_factor) * *(doubles+0));
   printf("v_doub_h & *(p_factor) * *(doubles+1) = %f %f\n",energy1,*(p_factor) * *(doubles+1));
*/
                     v_both = _mm_mul_pd(v_both,v_doub);   /* v_both = ( v_factor * v_doub ) * ( v_sing + v_doub ) */
/*
   _mm_storel_pd(&energy0,v_both);
   _mm_storeh_pd(&energy1,v_both);
   printf("v_both_l = v_doub * ( v_sing + v_doub ) = %f\n",energy0);
   printf("v_both_h = v_doub * ( v_sing + v_doub ) = %f\n",energy1);
*/
                     v_both = _mm_div_pd(v_both,v_denom);  /* v_both = ( term above ) / v_denom                    */
/*
   _mm_storel_pd(&energy0,v_both);
   _mm_storeh_pd(&energy1,v_both);
   printf("v_both_l = v_doub * v_denom = %f\n",energy0);
   printf("v_both_h = v_doub * v_denom = %f\n",energy1);
*/
                     v_energy = _mm_add_pd(v_energy,v_both);
/*
   _mm_storel_pd(&energy0,v_energy);
   _mm_storeh_pd(&energy1,v_energy);
   printf("v_energy_l = %f\n",energy0);
   printf("v_energy_h = %f\n",energy1);
*/

                     singles += 2;
                     doubles += 2;
                  }
               }
            }
         }
      }
   }

   printf("After loops\n");

   _mm_storel_pd(&energy0,v_energy);
   _mm_storeh_pd(&energy1,v_energy);

   printf("energy0 = %f\n",energy0);
   printf("energy1 = %f\n",energy1);
   *(p_energy) = energy0 + energy1;
   printf("*p_energy = %f\n",*p_energy);

   printf("Bottom of ccsd_t_sse2_aligned_\n");

}

/* This version works
void ccsd_t_sse2_unaligned_(int range_p4, int range_p5, int range_p6,
                            int range_h1, int range_h2, int range_h3,
                            double *denom_p4, double *denom_p5, double *denom_p6,
                            double *denom_h1, double *denom_h2, double *denom_h3,
                            double *singles,  double *doubles,
                            double *p_factor, double *p_energy){

   unsigned int p4,p5,p6,h1,h2,h3;
   double eps_p4,eps_p5,eps_p6,eps_h1,eps_h2,eps_h3;
   double denom_p,denom_h,denom;
   double val_sing,val_doub;
   static double factor;
   double energy;

   printf("Top of ccsd_t_sse2_unaligned_\n");

   factor = *(p_factor);
   energy = *(p_energy);

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

                     energy += factor * denom * val_doub * ( val_doub + val_sing );

                  }
               }
            }
         }
      }
   }

   *(p_energy) = energy;

   fflush(stdout);

}
*/


void ccsd_t_sse2_unaligned_(int range_p4, int range_p5, int range_p6,
                            int range_h1, int range_h2, int range_h3,
                            double *denom_p4, double *denom_p5, double *denom_p6,
                            double *denom_h1, double *denom_h2, double *denom_h3,
                            double *singles,  double *doubles,
                            double *p_factor, double *p_energy){

   unsigned int p4,p5,p6,h1,h2,h3;
   double eps_p4,eps_p5,eps_p6,eps_h1,eps_h2,eps_h3;
   double denom_p,denom_h,denom;
   double val_sing,val_doub;
   static double factor;
   double energy;

//   printf("Top of ccsd_t_sse2_unaligned_\n");

   factor = *(p_factor);
   energy = *(p_energy);

   energy = 0.0;

   #pragma parallel
   #pragma loop count min(4) max(32) avg(16)
   for ( p4 = 0 ; p4 < range_p4 ; p4++ ){
      eps_p4 = *( denom_p4 + p4 + 1 );
      #pragma loop count min(4) max(32) avg(16)
      for ( p5 = 0 ; p5 < range_p5 ; p5++ ){
         eps_p5 = *( denom_p5 + p5 + 1 );
         #pragma loop count min(4) max(32) avg(16)
         for ( p6 = 0 ; p6 < range_p6 ; p6++ ){
            eps_p6 = *( denom_p6 + p6 + 1 );

            denom_p = -1.0*( eps_p4 + eps_p5 + eps_p6 );

            #pragma loop count min(4) max(32) avg(16)
            for ( h1 = 0 ; h1 < range_h1 ; h1++ ){
               eps_h1 = *( denom_h1 + h1 + 1 );
               #pragma loop count min(4) max(32) avg(16)
               for ( h2 = 0 ; h2 < range_h2 ; h2++ ){
                  eps_h2 = *( denom_h2 + h2 + 1 );
                  #pragma loop count min(4) max(32) avg(16)
                  #pragma prefetch
                  #pragma unroll(8)
                  #pragma vector aligned
                  for ( h3 = 0 ; h3 < range_h3 ; h3++ ){
                     eps_h3 = *( denom_h3 + h3 + 1 );

                     denom_h = ( eps_h1 + eps_h2 + eps_h3 );

                     denom = 1.0e0 / ( denom_h + denom_p );

                     val_sing = *(singles++);
                     val_doub = *(doubles++);

                     energy += factor * denom * val_doub * ( val_doub + val_sing );

                  }
               }
            }
         }
      }
   }

   *(p_energy) = energy;

//   fflush(stdout);

}

#endif
