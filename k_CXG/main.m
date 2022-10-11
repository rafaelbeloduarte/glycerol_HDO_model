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
pkg load statistics
##1,2-PDO
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
k_13PDO = k(2,:);
k_EG = k(3,:);
k_S = k(4,:);
k_I = k(5,:);
k_G = k(6,:);
k_M = k(7,:);

[params_12PDO, fval] = fminsearch(@(params_12PDO) f_objective(params_12PDO,T,k_12PDO),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_13PDO, fval] = fminsearch(@(params_13PDO) f_objective(params_13PDO,T,k_13PDO),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_EG, fval] = fminsearch(@(params_EG) f_objective(params_EG,T,k_EG),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_S, fval] = fminsearch(@(params_S) f_objective(params_S,T,k_S),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_I, fval] = fminsearch(@(params_I) f_objective(params_I,T,k_I),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_G, fval] = fminsearch(@(params_G) f_objective(params_G,T,k_G),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
[params_M, fval] = fminsearch(@(params_M) f_objective(params_M,T,k_M),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));

csvwrite("params_12PDO.csv", params_12PDO);
csvwrite("params_13PDO.csv", params_13PDO);
csvwrite("params_EG.csv", params_EG);
csvwrite("params_S.csv", params_S);
csvwrite("params_I.csv", params_I);
csvwrite("params_G.csv", params_G);
csvwrite("params_M.csv", params_M);


A_12PDO = params_12PDO(1);
E_12PDO = params_12PDO(2);
k_12PDO_predicted = [];
A_13PDO = params_13PDO(1);
E_13PDO = params_13PDO(2);
k_13PDO_predicted = [];
A_EG = params_EG(1);
E_EG = params_EG(2);
k_EG_predicted = [];
A_S = params_S(1);
E_S = params_S(2);
k_S_predicted = [];
A_I = params_I(1);
E_I = params_I(2);
k_I_predicted = [];
A_G = params_G(1);
E_G = params_G(2);
k_G_predicted = [];
A_M = params_M(1);
E_M = params_M(2);
k_M_predicted = [];

T_predicted = (250+273.15):1:(285+273.15);
for i = 1:1:length(T_predicted)
  k_12PDO_predicted(end + 1) = A_12PDO*exp(-E_12PDO/(R*T_predicted(i)));
  k_13PDO_predicted(end + 1) = A_13PDO*exp(-E_13PDO/(R*T_predicted(i)));
  k_EG_predicted(end + 1) = A_EG*exp(-E_EG/(R*T_predicted(i)));
  k_S_predicted(end + 1) = A_S*exp(-E_S/(R*T_predicted(i)));
  k_I_predicted(end + 1) = A_I*exp(-E_I/(R*T_predicted(i)));
  k_G_predicted(end + 1) = A_G*exp(-E_G/(R*T_predicted(i)));
  k_M_predicted(end + 1) = A_M*exp(-E_M/(R*T_predicted(i)));
endfor

res_k_12PDO = [];
res_k_13PDO = [];
res_k_EG = [];
res_k_S = [];
res_k_I = [];
res_k_G = [];
res_k_M = [];

for i = 1:1:length(T)
  for j = 1:1:length(T_predicted)
    if T(i) == T_predicted(j)
      res_k_12PDO(end+1) = k_12PDO(i) - k_12PDO_predicted(j);
      res_k_13PDO(end+1) = k_13PDO(i) - k_13PDO_predicted(j);
      res_k_EG(end+1) = k_EG(i) - k_EG_predicted(j);
      res_k_S(end+1) = k_S(i) - k_S_predicted(j);
      res_k_I(end+1) = k_I(i) - k_I_predicted(j);
      res_k_G(end+1) = k_G(i) - k_G_predicted(j);
      res_k_M(end+1) = k_M(i) - k_M_predicted(j);
    endif
    endfor
  endfor

res_boot_12PDO =  [];
res_boot_13PDO =  [];
res_boot_EG =  [];
res_boot_S =  [];
res_boot_I =  [];
res_boot_G =  [];
res_boot_M =  [];
for i = 1:1:10000
  i
  bstrp = randsample (res_k_12PDO, length(res_k_12PDO), replacement=true);
  res_boot_12PDO(end+1) = mean(bstrp);
  bstrp = randsample (res_k_13PDO, length(res_k_13PDO), replacement=true);
  res_boot_13PDO(end+1) = mean(bstrp);
  bstrp = randsample (res_k_EG, length(res_k_EG), replacement=true);
  res_boot_EG(end+1) = mean(bstrp);
  bstrp = randsample (res_k_S, length(res_k_S), replacement=true);
  res_boot_S(end+1) = mean(bstrp);
  bstrp = randsample (res_k_I, length(res_k_I), replacement=true);
  res_boot_I(end+1) = mean(bstrp);
  bstrp = randsample (res_k_G, length(res_k_G), replacement=true);
  res_boot_G(end+1) = mean(bstrp);
  bstrp = randsample (res_k_M, length(res_k_M), replacement=true);
  res_boot_M(end+1) = mean(bstrp);
endfor

figure(2);
clf;
hist(res_boot_12P);

# now we start bootstraping confidence intervals
n_boot = 4;
n = 2000;
params_12PDO_distribution = [];
params_13PDO_distribution = [];
params_EG_distribution = [];
params_S_distribution = [];
params_I_distribution = [];
params_G_distribution = [];
params_M_distribution = [];
for i = 1:1:n
  i

  e_k_12PDO_boot = boot(res_boot_12PDO, n_boot);
  k_boot_12PDO = k_12PDO + e_k_12PDO_boot;
  [params_12PDO, fval] = fminsearch(@(params_12PDO) f_objective(params_12PDO,T,k_boot_12PDO),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_12PDO_distribution = [params_12PDO_distribution;params_12PDO'];

  e_k_13PDO_boot = boot(res_boot_13PDO, n_boot);
  k_boot_13PDO = k_13PDO + e_k_13PDO_boot;
  [params_13PDO, fval] = fminsearch(@(params_13PDO) f_objective(params_13PDO,T,k_boot_13PDO),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_13PDO_distribution = [params_13PDO_distribution;params_13PDO'];

  e_k_EG_boot = boot(res_boot_EG, n_boot);
  k_boot_EG = k_EG + e_k_EG_boot;
  [params_EG, fval] = fminsearch(@(params_EG) f_objective(params_EG,T,k_boot_EG),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_EG_distribution = [params_EG_distribution;params_EG'];

  e_k_S_boot = boot(res_boot_S, n_boot);
  k_boot_S = k_S + e_k_S_boot;
  [params_S, fval] = fminsearch(@(params_S) f_objective(params_S,T,k_boot_S),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_S_distribution = [params_S_distribution;params_S'];

  e_k_I_boot = boot(res_boot_I, n_boot);
  k_boot_I = k_I + e_k_I_boot;
  [params_I, fval] = fminsearch(@(params_I) f_objective(params_I,T,k_boot_I),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_I_distribution = [params_I_distribution;params_I'];

  e_k_G_boot = boot(res_boot_G, n_boot);
  k_boot_G = k_G + e_k_G_boot;
  [params_G, fval] = fminsearch(@(params_G) f_objective(params_G,T,k_boot_G),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_G_distribution = [params_G_distribution;params_G'];

  e_k_M_boot = boot(res_boot_M, n_boot);
  k_boot_M = k_M + e_k_M_boot;
  [params_M, fval] = fminsearch(@(params_M) f_objective(params_M,T,k_boot_M),[1;1], optimset ("MaxFunEvals", 1e40, "MaxIter", 1e40, "TolFun", 1e-15, "TolX", 1e-15));
  params_M_distribution = [params_M_distribution;params_M'];
  endfor

##figure(2);
##clf;
##hist(E_12PDO_dist);

A_12PDO_dist = params_12PDO_distribution(:,1);
E_12PDO_dist = params_12PDO_distribution(:,2);
A_CI_lw_12PDO = prctile(A_12PDO_dist,2.5);
A_CI_up_12PDO = prctile(A_12PDO_dist,97.5);
E_CI_lw_12PDO = prctile(E_12PDO_dist,2.5);
E_CI_up_12PDO = prctile(E_12PDO_dist,97.5);


A_13PDO_dist = params_13PDO_distribution(:,1);
E_13PDO_dist = params_13PDO_distribution(:,2);
A_CI_lw_13PDO = prctile(A_13PDO_dist,2.5);
A_CI_up_13PDO = prctile(A_13PDO_dist,97.5);
E_CI_lw_13PDO = prctile(E_13PDO_dist,2.5);
E_CI_up_13PDO = prctile(E_13PDO_dist,97.5);

A_EG_dist = params_EG_distribution(:,1);
E_EG_dist = params_EG_distribution(:,2);
A_CI_lw_EG = prctile(A_EG_dist,2.5);
A_CI_up_EG = prctile(A_EG_dist,97.5);
E_CI_lw_EG = prctile(E_EG_dist,2.5);
E_CI_up_EG = prctile(E_EG_dist,97.5);

A_S_dist = params_S_distribution(:,1);
E_S_dist = params_S_distribution(:,2);
A_CI_lw_S = prctile(A_S_dist,2.5);
A_CI_up_S = prctile(A_S_dist,97.5);
E_CI_lw_S = prctile(E_S_dist,2.5);
E_CI_up_S = prctile(E_S_dist,97.5);

A_I_dist = params_I_distribution(:,1);
E_I_dist = params_I_distribution(:,2);
A_CI_lw_I = prctile(A_I_dist,2.5);
A_CI_up_I = prctile(A_I_dist,97.5);
E_CI_lw_I = prctile(E_I_dist,2.5);
E_CI_up_I = prctile(E_I_dist,97.5);

A_G_dist = params_G_distribution(:,1);
E_G_dist = params_G_distribution(:,2);
A_CI_lw_G = prctile(A_G_dist,2.5);
A_CI_up_G = prctile(A_G_dist,97.5);
E_CI_lw_G = prctile(E_G_dist,2.5);
E_CI_up_G = prctile(E_G_dist,97.5);

A_M_dist = params_M_distribution(:,1);
E_M_dist = params_M_distribution(:,2);
A_CI_lw_M = prctile(A_M_dist,2.5);
A_CI_up_M = prctile(A_M_dist,97.5);
E_CI_lw_M = prctile(E_M_dist,2.5);
E_CI_up_M = prctile(E_M_dist,97.5);

CONFIDENCE_INTERVALS_A = [A_CI_lw_12PDO, A_CI_up_12PDO;
                                                       A_CI_lw_13PDO, A_CI_up_13PDO;
                                                       A_CI_lw_EG, A_CI_up_EG;
                                                       A_CI_lw_S, A_CI_up_S;
                                                       A_CI_lw_I, A_CI_up_I;
                                                       A_CI_lw_G, A_CI_up_G;
                                                       A_CI_lw_M, A_CI_up_M];
CONFIDENCE_INTERVALS_E = [E_CI_lw_12PDO, E_CI_up_12PDO;
                                                       E_CI_lw_13PDO, E_CI_up_13PDO;
                                                       E_CI_lw_EG, E_CI_up_EG;
                                                       E_CI_lw_S, E_CI_up_S;
                                                       E_CI_lw_I, E_CI_up_I;
                                                       E_CI_lw_G, E_CI_up_G;
                                                       E_CI_lw_M, E_CI_up_M];

csvwrite("dist_params_12PDO.csv", params_12PDO_distribution);
csvwrite("dist_params_13PDO.csv", params_13PDO_distribution);
csvwrite("dist_params_EG.csv", params_EG_distribution);
csvwrite("dist_params_S.csv", params_S_distribution);
csvwrite("dist_params_I.csv", params_I_distribution);
csvwrite("dist_params_G.csv", params_G_distribution);
csvwrite("dist_params_M.csv", params_M_distribution);
##dlmwrite('dist_k.csv',k_distribution,'delimiter',',','-append');
csvwrite("A_CI.csv", CONFIDENCE_INTERVALS_A);
csvwrite("E_CI.csv", CONFIDENCE_INTERVALS_E);

k_12PDO_lw = [];
k_12PDO_up = [];
k_13PDO_lw = [];
k_13PDO_up = [];
k_EG_lw = [];
k_EG_up = [];
k_S_lw = [];
k_S_up = [];
k_I_lw = [];
k_I_up = [];
k_G_lw = [];
k_G_up = [];
k_M_lw = [];
k_M_up = [];
for i = 1:1:length(T_predicted)
  k_12PDO_lw(end + 1) = A_CI_lw_12PDO*exp(-E_CI_lw_12PDO/(R*T_predicted(i)));
  k_12PDO_up(end + 1) = A_CI_up_12PDO*exp(-E_CI_up_12PDO/(R*T_predicted(i)));
  k_13PDO_lw(end + 1) = A_CI_lw_13PDO*exp(-E_CI_lw_13PDO/(R*T_predicted(i)));
  k_13PDO_up(end + 1) = A_CI_up_13PDO*exp(-E_CI_up_13PDO/(R*T_predicted(i)));
  k_EG_lw(end + 1) = A_CI_lw_EG*exp(-E_CI_lw_EG/(R*T_predicted(i)));
  k_EG_up(end + 1) = A_CI_up_EG*exp(-E_CI_up_EG/(R*T_predicted(i)));
  k_S_lw(end + 1) = A_CI_lw_S*exp(-E_CI_lw_S/(R*T_predicted(i)));
  k_S_up(end + 1) = A_CI_up_S*exp(-E_CI_up_S/(R*T_predicted(i)));
  k_I_lw(end + 1) = A_CI_lw_I*exp(-E_CI_lw_I/(R*T_predicted(i)));
  k_I_up(end + 1) = A_CI_up_I*exp(-E_CI_up_I/(R*T_predicted(i)));
  k_G_lw(end + 1) = A_CI_lw_G*exp(-E_CI_lw_G/(R*T_predicted(i)));
  k_G_up(end + 1) = A_CI_up_G*exp(-E_CI_up_G/(R*T_predicted(i)));
  k_M_lw(end + 1) = A_CI_lw_M*exp(-E_CI_lw_M/(R*T_predicted(i)));
  k_M_up(end + 1) = A_CI_up_M*exp(-E_CI_up_M/(R*T_predicted(i)));
endfor

figure(1);
clf;
% 12pdo
c = [0 0.443 0.737];
plot(T_predicted, k_12PDO_predicted, 'color', c);
hold on;
scatter(T, k_12PDO, 20, c, '+');
hold on;
plot(T_predicted, k_12PDO_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_12PDO_up, 'color', c, "--");
hold on;
% 13pdo
c = [0.847 0.322 0.094];
plot(T_predicted, k_13PDO_predicted, 'color', c);
hold on;
scatter(T, k_13PDO, 20, c, 'o');
hold on;
plot(T_predicted, k_13PDO_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_13PDO_up, 'color', c, "--");
hold on;
% eg
c = [0.925 0.690 0.122];
plot(T_predicted, k_EG_predicted, 'color', c);
hold on;
scatter(T, k_EG, 20, c, '*');
hold on;
plot(T_predicted, k_EG_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_EG_up, 'color', c, "--");
hold on;
% s
c = [0.298 0.741 0.929];
plot(T_predicted, k_S_predicted, 'color', c);
hold on;
scatter(T, k_S, 20, c, 'x');
hold on;
plot(T_predicted, k_S_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_S_up, 'color', c, "--");
hold on;
% i
c = [0.490 0.180 0.553];
plot(T_predicted, k_I_predicted, 'color', c);
hold on;
scatter(T, k_I, 20, c, 's');
hold on;
plot(T_predicted, k_I_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_I_up, 'color', c, "--");
hold on;
% g
c = [1 1 1];
plot(T_predicted, k_G_predicted, 'color', c);
hold on;
scatter(T, k_G, 20, c, 'd');
hold on;
plot(T_predicted, k_G_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_G_up, 'color', c, "--");
hold on;
% m
c = [0.463 0.671 0.184];
plot(T_predicted, k_M_predicted, 'color', c);
hold on;
scatter(T, k_M, 20, c, '^');
hold on;
plot(T_predicted, k_M_lw, 'color', c, ":");
hold on;
plot(T_predicted, k_M_up, 'color', c, "--");
hold on;
title("non linearized")

line_width = 2;
font_size = 16;
xlabel ("T(^oC)");
ylabel ("k (mol g^{-1} h^{-1})");
set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif", "ylim", [-inf inf], 'yscale', 'lin');
