%{
This program is part of a set of scripts to fit a zeroth order reaction model
to experimental data as well as calculate the parameters confidence intervals.
Copyright (C) 2022  Rafael Belo Duarte

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

contact me at rafaelbeloduarte@pm.me
%}
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
