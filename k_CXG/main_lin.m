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
% main_lin##1,2-PDO
##1,3-PDO
##Ethylene glycol
##Scopoletin
##Isobutyl Isovalerate
##Glycerol
##Methane
R = 8.314; % J/mol.K
T = [250 265 275 285];
T = 273.15 + T;
k = [5.98E-04	1.20E-03	2.85E-03	6.81E-03
1.30E-04	4.71E-04	5.61E-04	1.43E-03
2.15E-04	4.20E-04	7.83E-04	1.61E-03
2.66E-05	7.57E-05	2.66E-04	8.05E-04
7.67E-04	1.79E-03	3.70E-03	7.85E-03
2.88E-03	7.02E-03	1.13E-02	2.40E-02
2.78E-05	4.00E-05	9.34E-05	1.89E-04];

k_12PDO = k(1,:);

ln_k_12PDO = log(k_12PDO);
T_inv = 1./T;
P = polyfit(T_inv, ln_k_12PDO, 1);
E_12PDO = -P(1)*R
A_12PDO = exp(P(2))

k_12PDO_predicted = [];
T_predicted = (250+273.15):1:(285+273.15);

for i = 1:1:length(T_predicted)
  k_12PDO_predicted(end + 1) = A_12PDO*exp(-E_12PDO/(R*T_predicted(i)));
  endfor

figure(1);
clf;
plot(T_predicted, k_12PDO_predicted);
hold on;
scatter(T, k_12PDO);
hold on;
title("linearized")
