*/gtbce/: defines further paramters for GTBCE calculations
      integer, parameter ::
     &    mxpana = 20
      common /gtbce/ inc_sing(3), inc_doub(5),
     &             igtbcs, iopt_gtb, igtbmod, igtbfusk,
     &             n_eag, n_foo, n_hss, nn_eag,
     &             igtb_test_h1, igtb_prjout,
     &             igtb_closed, i_do_h0,
     &             igtb_disptt, igtb_pad0,
     &             it_eag(mxpana), it_foo(mxpana), it_hss(mxpana),
     &             ng_eag(mxpana), st_eag(mxpana),
     &             en_eag(mxpana), np_eag(mxpana),
     &             expg_thrsh, mxterm_expg
