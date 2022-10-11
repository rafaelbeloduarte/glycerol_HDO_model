function [bstrp] = boot(residual, n_boot)
  pkg load statistics;
##  t0 = clock ();
  bstrp =  [];
  bstrp = randsample (residual, n_boot, replacement=true);
##  elapsed_time = etime (clock (), t0)
##  hist(bstrp);
endfunction
