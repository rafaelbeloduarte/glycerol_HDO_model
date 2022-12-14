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
function Y = f_XG285(constante)
  intervalo_t = 0:0.01:6;

  % Vetor Condições iniciais
  xG0 = 0.687; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.000120736019106005	0.00190542440690286	0.00126244275900532	0.00245650973286063	0.00345920725224455	0.0045534425575047	0.00547419411429379
0.000485247518719119	0.00231941595281414	0.00349664450705343	0.00495123526098734	0.00879527646171398	0.0128643049161132	0.0157702558358661
0.000313476087311038	0.00319273182508872	0.00613915891361497	0.0032510980857634	0.00474166685870984	0.00584445103395859	0.00737038552765253
0.686609344313597	0.674117033507055	0.646049321827957	0.619995278542439	0.578654902103528	0.566372569330732	0.534311242079215
0.000355852100700251	0.000452642255611177	0.00153275471395357	0.00313811173253462	0.007510876639749	0.011367323695359	0.017134530695938
0.00203921003920804	0.00549912899109952	0.0216181875618648	0.0372484791805704	0.054171111606807	0.0561776543209205	0.0675149706050606
0	0	0	0	0	0.0000102198652078027	0.0000176773328838207];
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
