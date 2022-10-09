pkg load statistics
% first simulation to extract initial residuals
t_range = 0:1:6;
% f_val = 0,3367
k0 = [3.2720e-02
   1.2828e-02
   1.1427e-02
   2.0598e-03
   4.8762e-02
   1.9106e-01
   1.0878e-03];

  % Vetor Condições iniciais
  xG0 = 0.663; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(k0)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0.00148966367080551	0.00337922409161915	0.00271950782109699	0.00131112627273551	0.00247708118788349	0.00283805132749683
0	0.00265672349675731	0.00234456734507558	0.00591610976241684	0.00487954217418269	0.00755975198182747	0.0103152823018937
0	0.00259219763974022	0.00113961335119305	0.00438370481534784	0.00157658882294006	0.00291433362518941	0.00308338156307195
0.662842142497291	0.632385167165656	0.637238015782846	0.624655617990531	0.63427215363947	0.616302084101628	0.615219559455091
0	0.000741015866216374	0.000374140361981225	0.000571785436519997	0.000542987895105544	0.000500285647675669	0.000728684765825176
0	0.00383816255262137	0.00546869435959392	0.00602472604137568	0.00950606397390995	0.0127905903010819	0.0134623012467434
0	0.0000928676042083693	0.00022551715907685	0.00036416169332738	0.00125569025026	0.00104583408037889	0.00126523043917613];
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
