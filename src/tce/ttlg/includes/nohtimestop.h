cudaDeviceSynchronize(); 
tmpTime = omp_get_wtime() - start;
//        printf("\t%6.6lf\t",tmpTime);
        unsigned long int totalsize=8;
        for(int i=0; i < ndim; i++)
        {
                totalsize *= lda[i];
        }
  //      printf("\t%6.2lf\t",2*totalsize/(1000000000*tmpTime));

