/*
 $Id$
*/
      Vall_match   = Vall[match];
      dVall_match  = (1.0/(al*r[match]))*Derivative7_4(match,Vall);
      ddVall_match = (-1.0/(r[match]*r[match]*al))   *Derivative7_4(match,Vall)
                   + (+1.0/(r[match]*r[match]*al*al))*Laplacian7_4(match,Vall);

      rc[0]	       = 1.0;
      rc[1]	       = r[match];
      for (i=2; i<13; ++i) 
         rc[i] = rc[1]*rc[i-1];

      /* debug */
      for (i=0; i<Ngrid; ++i)
         upp[i] = (2.0*Vall[i] + l*(l+1.0)/(r[i]*r[i]) - 2.0*el)*ul[i];
      upp_match   = upp[match];
      uppp_match  = (1.0/(al*r[match]))*Derivative7_4(match,upp);
      upppp_match = (-1.0/(r[match]*r[match]*al))   *Derivative7_4(match,upp)
                   + (+1.0/(r[match]*r[match]*al*al))*Laplacian7_4(match,upp);

      printf("dVall_match  = %le\n", dVall_match);
      printf("ddVall_match = %le\n", ddVall_match);

      /******************************************/
      /* Calculate the all-electron core charge */
      /******************************************/
      ae_core_charge = Norm_LogGrid(match,(l+1.0),ul);
      printf("ae_core_charge = %le\n",ae_core_charge);


      /**************************************************************/
      /* define p(rcl), p'(rcl), p''(rcl), p'''(rcl) and p''''(rcl) */
      /**************************************************************/
      poly[0] = log(ul[match]/pow(rc[1],(l+1.0)));
      poly[1] = ul_prime[match]/ul[match] 
              - (l+1.0)/rc[1];
      poly[2] = 2.0*Vall_match
              - 2.0*el 
              - (2.0*(l+1.0)/rc[1])*poly[1] 
              - poly[1]*poly[1];
      poly[3] = 2.0*dVall_match 
              + (2.0*(l+1.0)/rc[2])*poly[1]
              - (2.0*(l+1.0)/rc[1])*poly[2]
              - 2.0*poly[1]*poly[2];
      poly[4] = 2.0*ddVall_match
              - (4.0*(l+1.0)/rc[3])*poly[1]
              + (4.0*(l+1.0)/rc[2])*poly[2]
              - (2.0*(l+1.0)/rc[1])*poly[3]
              - 2.0*poly[2]*poly[2]
              - 2.0*poly[1]*poly[3];
       
