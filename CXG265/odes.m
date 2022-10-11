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
% Função com o sistema de EDOs
function [dY] = odes(t, Y)
  k1 = Y(8);
  k2 = Y(9);
  k3 = Y(10);
  k4 = Y(11);
  k5 = Y(12);
  k6 = Y(13);
  k7 = Y(14);
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
  dY = zeros(14,1);
  n_T = m_T/(Y(1)*MM_12P + Y(2)*MM_13P + Y(3)*MM_EG + Y(4)*MM_S + Y(5)*MM_I + Y(6)*MM_G + Y(7)*MM_M);
  dY(1) = V_L*r12P(Y)/n_T; % dx_12P/dt
  dY(2) = V_L*r13P(Y)/n_T; % dx_13P/dt
  dY(3) = V_L*rEG(Y)/n_T; % dx_EG/dt
  dY(4) = V_L*rS(Y)/n_T; % dx_S/dt
  dY(5) = V_L*rI(Y)/n_T; % dx_I/dt
  dY(6) = V_L*rG(Y)/n_T; % dx_G/dt
  dY(7) = V_G*rM(Y)/n_T; % dx_M/dt
  dY(8) = 0; % k1
  dY(9) = 0; % k2
  dY(10) = 0; % k3
  dY(11) = 0; % k4
  dY(12) = 0; % k5
  dY(13) = 0; % k6
  dY(14) = 0; % k7
end

function resultado_r12P = r12P(Y)
  k1 = Y(8);
  r12P = k1;
  resultado_r12P = r12P;
end

function resultado_r13P = r13P(Y)
  k2 = Y(9);
  r13P = k2;
  resultado_r13P = r13P;
end

function resultado_rEG = rEG(Y)
  k3 = Y(10);
  rEG = k3;
  resultado_rEG = rEG;
end

function resultado_rS = rS(Y)
  k4 = Y(11);
  rS = k4;
  resultado_rS = rS;
end

function resultado_rI = rI(Y)
  k5 = Y(12);
  rI = k5;
  resultado_rI = rI;
end


function resultado_rG = rG(Y)
  k6 = Y(13);
  rG = -k6;
  resultado_rG = rG;
end

function resultado_rM = rM(Y)
  k7 = Y(14);
  rM = k7;
  resultado_rM = rM;
end
