  %% Example for user specified Jacobian.

  %% independents
  x = [1:10:100]';
  %% observed data
  y =[9.2160e-001, 3.3170e-001, 8.9789e-002, 2.8480e-002, 2.6055e-002,...
     8.3641e-003,  4.2362e-003,  3.1693e-003,  1.4739e-004,  2.9406e-004]';
  %% initial values:
  p0=[0.8; 0.05];
  %% bounds
  lb=[0; 0]; ub=[];
  %% Jacobian setting
  opts = optimset ("Jacobian", "on")

  %% model function:
  function [F,J] = myfun (p, x)
    F = p(1) * exp (-p(2) * x);
    if nargout > 1
      J = [exp(- p(2) * x), - p(1) * x .* exp(- p(2) * x)];
    endif
  endfunction

  [c, resnorm, residual, flag, output, lambda, jacob] = ...
      lsqcurvefit (@ (varargin) myfun(varargin{:}), p0, x, y, lb,  ub, opts)
