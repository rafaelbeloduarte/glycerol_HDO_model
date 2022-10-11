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
function Y = f_XG275(constante)
  intervalo_t = 0:0.01:6;

  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.00194549439845336	0.00087395975525803	0.00149818889892486	0.0024908991976411	0.00504508604594008	0.00387379954231023	0.00362180474518346
0.00216965789727529	0.00137752474610415	0.00352074784393782	0.00443059632909249	0.00664672827892196	0.00707658950004113	0.00881691751137682
0.00120648035047698	0.000741429597178869	0.00258519434767493	0.00311480245840776	0.00424196114852012	0.00333162652721618	0.00356864835517883
0.688822626364104	0.695865769598272	0.670716222663411	0.656100405778356	0.628654017166774	0.624994356029778	0.602407417663549
0	0.0000916192956785067	0.000626363388483587	0.000975518345421958	0.00223768655071472	0.00307931528475682	0.00588481057466947
0.000637319447767249	0.00461144893058251	0.0125505774448675	0.019782580977562	0.0324970220260622	0.0409478188712749	0.0542636268256832
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
