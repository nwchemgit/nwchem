# GW

task dft gw

```
gw
  [ ( G0W0 || evGW0 || evGW ) default G0W0 ] 
  [ method ( cdgw || analytic ) default analytic ] 
  [ states ( alpha || beta) occ <integer noqp default 1> vir <integer nvqp default 0> ] 
  [ solver ( graph || newton ) default newton ] 
  [ convergence 
       <real thresh default 0.001> 
       [unit (eV || Ha) default eV]  ]
end
```


