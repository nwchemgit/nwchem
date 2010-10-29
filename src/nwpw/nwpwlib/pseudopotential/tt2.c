/*
 $Id$
*/
      /***************************************************************/
      /* get the coeficients of exp(p(r), i.e. c0,c2,c4,c6,c8,10,c12 */
      /***************************************************************/
      get_cs(c,l,match,rc,ae_core_charge,poly);

      printf("delta     = %le \n",c[0]);

      generate_p(c,match,p);
      generate_dp(c,match,dp);
      generate_ddp(c,match,ddp);
      generate_dddp(c,match,dddp);
      generate_ddddp(c,match,ddddp);

      /***********************************/
      /* generate the pseudowavefunction */
      /***********************************/
      generate_psp_psi(c,l,match,wl);
      for (i=match+1; i<Ngrid; ++i)
         wl[i] = ul[i];
      printf("ul/wl     = %le \n",ul[match]/wl[match]);

     
      /*********************/
      /* generate (d/dr)wl */
      /*********************/
      for (i=0; i<=match; ++i)
      {
        wl_prime[i] = (l+1.0)*pow(r[i],l)*exp(p[i])
                    + pow(r[i],(l+1.0))*dp[i]*exp(p[i]); 
      }
      for (i=match+1; i<Ngrid; ++i)
         wl_prime[i] = ul_prime[i];
      printf("dul/dwl   = %le \n",ul_prime[match]/wl_prime[match]);


      /*************************/
      /* generate (d/dr)**2 wl */
      /*************************/
      for (i=0; i<=match; ++i)
      {
        wlpp[i] = exp(p[i])*( l*(l+1.0)*pow(r[i],l-1.0)
                            + 2.0*(l+1.0)*pow(r[i],l)*dp[i]
                            + pow(r[i],l+1.0)*ddp[i]
                            + pow(r[i],l+1.0)*dp[i]*dp[i]
                            );
      }
      for (i=match+1; i<Ngrid; ++i)
        wlpp[i] = upp[i];
      
      printf("d2ul/dwl2 = %le \n",upp_match/wlpp[match]);


      /*************************/
      /* generate (d/dr)**3 wl */
      /*************************/
      wlppp_match = exp(p[match])*
                 ( (l-1.0)*l*(l+1.0)*pow(r[match],l-2.0)
                 + 3.0*l*(l+1.0)*   pow(r[match],l-1.0)*dp[match]
                 + 3.0*(l+1.0)*     pow(r[match],l    )*dp[match]*dp[match]
                 + 3.0*(l+1.0)*     pow(r[match],l    )*ddp[match]
                 + 3.0*             pow(r[match],l+1.0)*dp[match]*ddp[match]
                 +                  pow(r[match],l+1.0)*dddp[match]
                 +                  pow(r[match],l+1.0)*dp[match]
                                                       *dp[match]*dp[match]
                 );
      printf("d3ul/dwl3 = %le \n",uppp_match/wlppp_match);


      /*************************/
      /* generate (d/dr)**4 wl */
      /*************************/
      wlpppp_match = exp(p[match])*
                 ( (l-2.0)*(l-1.0)*l*(l+1.0)*pow(r[match],l-3.0)
                 + 4.0*(l-1.0)*l*(l+1.0)*    pow(r[match],l-2.0)*dp[match]
                 + 6.0*l*(l+1.0)*   pow(r[match],l-1.0)*dp[match]*dp[match]
                 + 6.0*l*(l+1.0)*   pow(r[match],l-1.0)*ddp[match]
                 + 4.0*(l+1.0)*     pow(r[match],l    )*dp[match]*dp[match]
						       *dp[match]
                 + 12.0*(l+1.0)*    pow(r[match],l    )*dp[match]*ddp[match]
                 + 4.0*(l+1.0)*     pow(r[match],l    )*dp[match]*dp[match]
                 + 4.0*(l+1.0)*     pow(r[match],l    )*dddp[match]
                 + 1.0*             pow(r[match],l+1.0)*dp[match]*dp[match]
						       *dp[match]*dp[match]
                 + 6.0*             pow(r[match],l+1.0)*dp[match]*dp[match]
						       *ddp[match]
                 + 3.0*             pow(r[match],l+1.0)*ddp[match]*ddp[match]
                 + 4.0*             pow(r[match],l+1.0)*dp[match]*dddp[match]
                 +                  pow(r[match],l+1.0)*ddddp[match]
                 );
      printf("d4ul/dwl4 = %le \n",upppp_match/wlpppp_match);



      /*********************************/
      /* invert the pseudowavefunction */
      /*********************************/
      generate_psp(c,l,match,el,Vall,Vl);

      printf("Vall/Vl     = %le \n",Vall[match]/Vl[match]);

      uppp_match = -((l+1.0)/(r[match]*r[match]))*dp[match]
                 +  ((l+1.0)/(r[match]))         *ddp[match]
                 +  0.5*dddp[match]
                 +  1.0*dp[match]*ddp[match];
      printf("dVall = %le  dVl = %le \n",dVall_match,uppp_match);

      uppp_match = 2.0*((l+1.0)/(r[match]*r[match]*r[match]))*dp[match]
                 - 2.0*((l+1.0)/(r[match]*r[match]))         *ddp[match]
                 + ((l+1.0)/(r[match]))*dddp[match]
                 + 0.5*ddddp[match]
                 + ddp[match]*ddp[match]
                 + dp[match]*dddp[match];
      printf("d2Vall = %le  dVl2 = %le \n",ddVall_match,uppp_match);
      
		
      /* debug */
      printf("poly[0] = %le\t p[match]     = %le\n",poly[0],p[match]);
      printf("poly[1] = %le\t dp[match]    = %le\n",poly[1],dp[match]);
      printf("poly[2] = %le\t ddp[match]   = %le\n",poly[2],ddp[match]);
      printf("poly[3] = %le\t dddp[match]  = %le\n",poly[3],dddp[match]);
      printf("poly[4] = %le\t ddddp[match] = %le\n",poly[4],ddddp[match]);
	

