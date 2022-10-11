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
function Y = f_XG250(constante)
  intervalo_t = 0:0.01:6;

  % Vetor Condições iniciais
  xG0 = 0.684; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0	0	0.0000705673656836223	0	0.000224334021026903	0.0000263653895395729
0	0.000049194610207348	0.00008038826006635	0.000167626216044625	0.00023443991334025	0.000444875774429263	0.000533540426327624
0	0	0	0	0	0.000230194856854255	0.000233529596418148
0.684442205351682	0.683663031866356	0.681776781563699	0.679726838611213	0.678165959149904	0.674577022915608	0.672741971574695
0	0	0	0	0	0.0000603036664514772	0
0	0.000700825416546395	0.0021718615139465	0.00346159006059108	0.00477848010634989	0.0065692398626041	0.00834782644138291
0	0	0	0	0	0	0];
  x_EG_exp = x(1,:);
  x_12P_exp = x(2,:);
  x_13P_exp = x(3,:);
  x_G_exp = x(4,:);
  x_S_exp = x (5,:);
  x_I_exp = x(6,:);
  x_M_exp = x(7,:);

  % Resolução das EDOs utilizando o método de Dormand-Prince de ordem 4
  [t Y] = ode45(@odes, intervalo_t, Y0);
  % Abrindo o vetor solução em suas componentes
  x_12P = [];
  x_13P = [];
  x_EG = [];
  x_S = [];
  x_I = [];
  x_G = [];
  x_M = [];

  x_12P = Y(:,1);
  x_13P = Y(:,2);
  x_EG = Y(:,3);
  x_S = Y (:,4);
  x_I = Y(:,5);
  x_G = Y(:,6);
  x_M = Y(:,7);

  soma_dif_quad = 0;

  for i = 1:1:length(t_exp)
    for j = 1:1:length(intervalo_t)
      dif = abs(t_exp(i) - intervalo_t(j));
      if dif < 0.001
          soma_dif_quad = soma_dif_quad + real((x_G(j) - x_G_exp(i))^2 + (x_12P(j) - x_12P_exp(i))^2 + (x_13P(j) - x_13P_exp(i))^2 + (x_EG(j) - x_EG_exp(i))^2 + (x_I(j) - x_I_exp(i))^2 + (x_M(j) - x_M_exp(i))^2 + (x_S(j) - x_S_exp(i))^2);
      endif
    endfor
  endfor
  Y = soma_dif_quad
endfunction

% [xmin, fval]=fminsearch(@f_250,[11.3], optimset ("Display", "iter"))
% [x, fval, info, output, grad, hess] = fminunc(@f_300,[10;0.5;0.1], optimset ("TolFun",0.0001))
