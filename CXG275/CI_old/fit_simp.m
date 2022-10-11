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
%fit_simp.m
function diff = fit_simp(param,t,x)

k1=param(1);
k2=param(2);
k3=param(3);
k4=param(4);
k5=param(5);
k6=param(6);
k7=param(7);

tspan=0:1:6;

xG0 = 0.697; % mol/L
init = [0.00001 0.00001 0.00001 xG0 0.00001 0.00001 0.00001];
[t,Y]=ode45(@odes_nova, tspan, init);
save('test');
Y
x
diff = (Y-x)^2

% Função com o sistema de EDOs
function dY = odes_nova(~, Y)
  V_L = 0.06; % L
  V_G = 0.3 - V_L; % L
  m_T = 1.26*V_L*1000; % g
  MM_G = 92.094;
  MM_12P = 76.095;
  MM_13P = 76.095;
  MM_EG = 62.068;
  MM_I = 158.24;
  MM_M = 16;
  MM_S = 192.16;
  % O vetor coluna das EDOs 'dY'
  % deve ser gerado pois a função 'edo45'aceita apenas
  % duas entradas no campo função (Variável independente, Variável dependente)
  % Y(1) = P_CH3OH, Y(2) = P_H2O, Y(3) = P_H2, Y(4) = P_CO2, Y(5) = P_CO, Y(6) = T
  % P_A = P_CH3OH, P_B = P_H2O, P_C = P_H2, P_D = P_CO2, P_E = P_CO, T=T0
  dY = zeros(7,1);
  n_T = m_T/(Y(2)*MM_12P + Y(3)*MM_13P + Y(1)*MM_EG + Y(5)*MM_S + Y(6)*MM_I + Y(4)*MM_G + Y(7)*MM_M);
  dY(1) = V_L*k1/n_T; % dx_EG/dt
  dY(2) = V_L*k2/n_T; % dx_12P/dt
  dY(3) = V_L*k3/n_T; % dx_13P/dt
  dY(4) = -V_L*k4/n_T; % dx_G/dt
  dY(5) = V_L*k5/n_T; % dx_S/dt
  dY(6) = V_L*k6/n_T; % dx_I/dt
  dY(7) = V_G*k7/n_T; % dx_M/dt
end
end

##  x_EG_exp = x(1,:);
##  x_12P_exp = x(2,:);
##  x_13P_exp = x(3,:);
##  x_G_exp = x(4,:);
##  x_S_exp = x (5,:);
##  x_I_exp = x(6,:);
##  x_M_exp = x(7,:);

