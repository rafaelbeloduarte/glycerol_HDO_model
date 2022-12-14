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
function Y = f_CXG265(constante)
  intervalo_t = 0:0.01:6;

  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.00319405967940873	0.00164746168966925	0.00435408264677725	0.00589579091206648	0.00992218399149192	0.00937843951014123	0.00891449586718024
0.00115000920151936	0.00380478917780572	0.0124775300589467	0.0190923030315172	0.0433548830767874	0.0423689731542582	0.0394945343094663
0.00100554616685905	0.0014772555551388	0.0046016690585132	0.00506590213688732	0.00808372637745789	0.0078195303291812	0.00879400824960989
0.686869092252988	0.676493129296004	0.638162863496033	0.60625089458323	0.558821105508447	0.552794270849324	0.533776665933563
0.000341420115327105	0.000426603985522671	0.00108735858928218	0.00241640184470745	0.00347481766940102	0.00473875896015277	0.00651928920195604
0.00179131307316157	0.00831071357286639	0.0169422638012227	0.0308377822129207	0.0339635701481284	0.0439271905208685	0.0545461053884707
0	0.000174004664341329	0.000646946161425272	0.00201749537546178	0.00376287562664706	0.0051148183064752	0.00516357939171532];
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
