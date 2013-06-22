/* $Id$ */
/**
 * NXTVAL Atomic Counter:
 * ----------------------
 * 
 * initalization and termination of atomic counter is collective 
 * operation. The actual get counter operation is atomic RMW.
 *
 * To initialize: If the value of the argument is zero, then nxtval counter 
 * is initialized (In C programs, call as follows)
 *      val=0;
 *      util_gnxtval_(&val); (this is a collective operation)
 *
 * To Terminate: (any negative value passed to the nxtval terminates 
 * this counter)
 *      val=-1;
 *      next_value = util_gnxtval_(&val);
 * 
 * The actual operation (nxtval): If the value is > 0, then nxtval() returns
 * the next value in the counter and increments the counter.
 */
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"

static int g_T;
static long initval=0;
static short int initialized=0;
static int subscript = 0;

Integer util_tcesublock_(Integer *val,Integer *p_handle) {

//    if(*p_handle==0) exit(1);
//ga_error("nxtask: p_handle is zero", 1);

    if(*val > 0) {
//       if(!initialized) exit(1);
//ga_error("nxtask: not yet initialized", 1);
       return (Integer) NGA_Read_inc(g_T, &subscript, 1);
    }
    else if(*val==0) {
       int n = 1;
       initialized=1;
       int p_h = (int)*p_handle;

       /* create task array */
//       g_T = NGA_Create(C_LONG, 1, &n,"Atomic Task", NULL);

       g_T = NGA_Create_config(C_LONG,1,&n,"Atomic Task",NULL,p_h);

       /* Initialize the task array */
       if(GA_Pgroup_nodeid(p_h)==0) {
	  int lo=0, hi=0;
	  NGA_Put (g_T, &lo, &hi, &initval, &hi);
//          printf("PUT %i %i %i\n",sizeof(*p_handle),sizeof(Integer),sizeof(int));
	  initval=0;
       }

       GA_Pgroup_sync(p_h);
//       printf("CREATE %i %i \n",*p_handle,g_T);
       return 0;
    }
    else if (*val < 0) {
        GA_Destroy(g_T);
//        printf("DELETE %i %i \n",*p_handle,g_T);
//        ga_pgroup_sync_(p_handle);
        initialized=0; 
        initval=0; 
        return 0;
   }
    
//    ga_error("nxtval: invalid value passed", 0L);
    return -1;
}

/*void nxtval_initval_(Integer *val) {
    initval=(long) *val;
    }*/
