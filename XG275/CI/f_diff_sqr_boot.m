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
function Y = f_diff_sqr_boot(k, X_calculated, t_range)
  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  X0 = [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, xG0, 0.00001];
  % Resolução das EDOs utilizando o método de Dormand-Prince de ordem 4
  [t X] = ode45(@(t, X) odes_boot(t,X,abs(k)), t_range, X0);
  % Abrindo o vetor solução em suas componentes
  Y = sum(sum((X - X_calculated).^2));
endfunction

% [xmin, fval]=fminsearch(@f_250,[11.3], optimset ("Display", "iter"))
% [x, fval, info, output, grad, hess] = fminunc(@f_300,[10;0.5;0.1], optimset ("TolFun",0.0001))
