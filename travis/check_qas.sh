#!/bin/bash
#          echo "The steps context is:"
#          echo "${{ toJson(steps) }}"
          echo "!!!!!!!!!!!!!!!!!!!!!!!!!"
          if [[ -f QA/testoutputs/dft_he2+.out ]]; then \
          grep d= QA/testoutputs/dft_he2+.out || true ;  fi
          if [[ -f QA/testoutputs/h2o_opt.out ]]; then \
          grep @ QA/testoutputs/h2o_opt.out || true ; fi
          if [[ -f QA/testoutputs/h2o2-response.out ]]; then \
          tail -490 QA/testoutputs/h2o2-response.out || true; fi
          if [[ -f QA/testoutputs/prop_mep_gcube.out ]]; then \
          tail -30 QA/testoutputs/prop_mep_gcube.out || true; fi
          if [[ -f QA/testoutputs/pspw.out ]]; then \
          tail -490 QA/testoutputs/pspw.out || true; fi
          if [[ -f QA/testoutputs/dft_smear.out ]]; then \
          grep @ QA/testoutputs/dft_smear.out || true; fi
          if [[ -f QA/testoutputs/dft_smear.out ]]; then \
          egrep d=  QA/testoutputs/dft_smear.out || true; fi
          if [[ -f QA/testoutputs/dft_smear.out ]]; then \
          tail -2000  QA/testoutputs/dft_smear.out || true; fi
          if [[ -f QA/testoutputs/ritddft_h2o.out ]]; then \
          tail -2000  QA/testoutputs/ritddft_h2o.out || true; fi
          if [[ -f QA/testoutputs/ritddft_co.out ]]; then \
          tail -2000  QA/testoutputs/ritddft_co.out || true; fi
          if [[ -f QA/testoutputs/pspw_md.out ]]; then \
          grep 'PSPW energy'  QA/testoutputs/pspw_md.out || true; fi
