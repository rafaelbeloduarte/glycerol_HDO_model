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
function Y = f_XG265(constante)
  intervalo_t = 0:0.01:6;

  % Vetor Condições iniciais
  xG0 = 0.677; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0.000258008137456242	0.00121734948945218	0.000990533154712822	0.00151925377890684	0.00189741662954417	0.000962807544541986
0	0.000347947079806869	0.00117345315677939	0.00171019875794062	0.00272965691238154	0.00293876797155187	0.00296712067296808
0	0.000349178250508359	0.000884365384512156	0.000918912902739483	0.0016725747998978	0.00118313576345039	0.00151098201015375
0.677031024434055	0.670844940635314	0.663018117833216	0.649846337997237	0.643460785061871	0.642518115451164	0.628074600961928
0	0.000238936100932408	0.000254443701267073	0.000502359904513105	0.000625801411513753	0.000806059076793232	0.00124700350449904
0	0.0033193723371167	0.00846468693777631	0.012777414103993	0.0165508380532851	0.0213796891746793	0.0268823418573825
0	0	0	0	0.0000238830013814539	0.0000291758702947773	0.0000369245562455555];
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
