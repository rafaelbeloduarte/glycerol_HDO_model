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
k0 = [1.1303e-02
   5.8348e-03
   5.9348e-03
   3.5684e-03
   8.6974e-02
   1.5845e-01
   1.4227e-05];

  % Initial conditions
  xG0 = 0.677; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(k0)];

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

  x = [x_12P_exp', x_13P_exp',x_EG_exp',x_S_exp',x_I_exp',x_G_exp',x_M_exp'];

  % ODEs solution using Dormand-Prince method of order 4
  [t Y] = ode45(@(t, Y) odes(t, Y), t_range, Y0);
  % dividing solution matrix into individual vectors
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
  % salving the k_j
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
