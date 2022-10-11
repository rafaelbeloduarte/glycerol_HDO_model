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
% first simulation to extract initial residuals
t_range = 0:1:6;
% f_val = 0,3367
k0 = [1.5479e-01
   3.2578e-02
   3.6571e-02
   1.8304e-02
   1.7844e-01
   5.4667e-01
   4.3070e-03];

  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(k0)];

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

  x = [x_12P_exp', x_13P_exp',x_EG_exp',x_S_exp',x_I_exp',x_G_exp',x_M_exp'];

  % Resolução das EDOs utilizando o método de Dormand-Prince de ordem 4
  [t Y] = ode45(@(t, Y) odes(t, Y), t_range, Y0);
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

  res_x_12P = [];
  res_x_13P = [];
  res_x_EG = [];
  res_x_S = [];
  res_x_I = [];
  res_x_G = [];
  res_x_M = [];

  for i = 1:1:length(t_exp)
    for j = 1:1:length(t_range)
      if t_exp(i) == t_range(j)
          res_x_12P(end+1) = x_12P_exp(i) - x_12P(j);
          res_x_13P(end+1) = x_13P_exp(i) - x_13P(j);
          res_x_EG(end+1) = x_EG_exp(i) - x_EG(j);
          res_x_S(end+1) = x_S_exp(i) - x_S(j);
          res_x_I(end+1) = x_I_exp(i) - x_I(j);
          res_x_G(end+1) = x_G_exp(i) - x_G(j);
          res_x_M(end+1) = x_M_exp(i) - x_M(j);
      endif
    endfor
endfor

  res_boot_12P =  [];
  res_boot_13P =  [];
  res_boot_EG =  [];
  res_boot_S =  [];
  res_boot_I =  [];
  res_boot_G =  [];
  res_boot_M =  [];
for i = 1:1:10000
  i
  bstrp = randsample (res_x_12P, length(res_x_12P), replacement=true);
  res_boot_12P(end+1) = mean(bstrp);
  bstrp = randsample (res_x_13P, length(res_x_13P), replacement=true);
  res_boot_13P(end+1) = mean(bstrp);
  bstrp = randsample (res_x_EG, length(res_x_EG), replacement=true);
  res_boot_EG(end+1) = mean(bstrp);
  bstrp = randsample (res_x_S, length(res_x_S), replacement=true);
  res_boot_S(end+1) = mean(bstrp);
  bstrp = randsample (res_x_I, length(res_x_I), replacement=true);
  res_boot_I(end+1) = mean(bstrp);
  bstrp = randsample (res_x_G, length(res_x_G), replacement=true);
  res_boot_G(end+1) = mean(bstrp);
  bstrp = randsample (res_x_M, length(res_x_M), replacement=true);
  res_boot_M(end+1) = mean(bstrp);
endfor

figure(4);
clf;
hist(res_boot_12P);
% end of first simulation, residuals extracted

# now we start bootstraping confidence intervals
% load the package for bootstraping
pkg load statistics-bootstrap

n_boot = 7 % use multiples of the total time
##t_boot = (0:6/n_boot:6-6/n_boot)';
t_boot=t_range;
% bootstraping the residuals
% initial conditions for x_j
X0 = [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, xG0, 0.00001];
csvwrite("k0.csv", k0');
k_distribution = [];
n = 2000;
[t X] = ode45(@(t, X) odes_boot(t,X,k0), t_boot, X0);
% recalculating k_j from bootstrapped residuals n times
k = k0;
for i = 1:1:n
  i
  e_x_12P_boot = boot(res_boot_12P, n_boot);
  e_x_13P_boot = boot(res_boot_13P, n_boot);
  e_x_EG_boot = boot(res_boot_EG, n_boot);
  e_x_S_boot = boot(res_boot_S, n_boot);
  e_x_I_boot = boot(res_boot_I, n_boot);
  e_x_G_boot = boot(res_boot_G, n_boot);
  e_x_M_boot = boot(res_boot_M, n_boot);
  E_X = [e_x_12P_boot; e_x_13P_boot; e_x_EG_boot; e_x_S_boot; e_x_I_boot; e_x_G_boot; e_x_M_boot];
  % adding x_j bootstrap residuals (to get x_j*)
  X_boot = X + E_X;
  [k, fval, info, output, grad, hess] = fminunc(@(k0) f_diff_sqr_boot(k0, X_boot, t_boot), k0);
  k
  % salvando os k_j
  k_distribution = [k_distribution;k'];
endfor
k_distribution = k_distribution(all(k_distribution > 0, 2), :);
k_distribution
dist_k_12P = k_distribution(:,1)
dist_k_13P = k_distribution(:,2)
dist_k_EG = k_distribution(:,3)
dist_k_S = k_distribution(:,4)
dist_k_I = k_distribution(:,5)
dist_k_G = k_distribution(:,6)
dist_k_M = k_distribution(:,7)

figure(4);
clf;
hist(dist_k_12P);

CI_lw_12P = prctile(dist_k_12P,2.5);
CI_up_12P = prctile(dist_k_12P,97.5);

CI_lw_13P = prctile(dist_k_13P,2.5);
CI_up_13P = prctile(dist_k_13P,97.5);

CI_lw_EG = prctile(dist_k_EG,2.5);
CI_up_EG = prctile(dist_k_EG,97.5);

CI_lw_S = prctile(dist_k_S,2.5);
CI_up_S = prctile(dist_k_S,97.5);

CI_lw_I = prctile(dist_k_I,2.5);
CI_up_I = prctile(dist_k_I,97.5);

CI_lw_G = prctile(dist_k_G,2.5);
CI_up_G = prctile(dist_k_G,97.5);

CI_lw_M = prctile(dist_k_M,2.5);
CI_up_M = prctile(dist_k_M,97.5);

CONFIDENCE_INTERVALS = [CI_lw_12P, CI_up_12P;
                                                   CI_lw_13P, CI_up_13P;
                                                   CI_lw_EG, CI_up_EG;
                                                   CI_lw_S, CI_up_S;
                                                   CI_lw_I, CI_up_I;
                                                   CI_lw_G, CI_up_G;
                                                   CI_lw_M, CI_up_M]

csvwrite("dist_k.csv", k_distribution);
##dlmwrite('dist_k.csv',k_distribution,'delimiter',',','-append');
csvwrite("CI.csv", CONFIDENCE_INTERVALS);
