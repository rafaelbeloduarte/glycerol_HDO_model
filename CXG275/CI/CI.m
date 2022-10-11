pkg load statistics
% first simulation to extract initial residuals
t_range = 0:1:6;
k0 = [0.0711266099963456
0.0142443243265413
0.0193547518424003
0.00669429418027291
0.0951599505233869
0.283905370037338
0.0023480527582845];

  % Initial conditions
  xG0 = 0.697; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(k0)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.000854591538052104	0.000735741867627219	0.00194497996821644	0.00320424754168804	0.00396335808191175	0.00493929036908437	0.00569302635494354
0.000214868686104798	0.00131685656790968	0.00439722017907256	0.00941648668994965	0.0139112157053546	0.0175640855506391	0.0243574983801011
0	0.000465682969222304	0.00152324554357987	0.00226607982767692	0.002829008888551	0.00380209469987456	0.00409257132061907
0.697410535444331	0.693071469364175	0.674592916711882	0.652972973236454	0.638002912615113	0.62317263780838	0.609259693065981
0.000133923627881811	0.000258008988989206	0.000373434229557013	0.00078300660576858	0.00117820048555334	0.0017452268097832	0.00204860963506344
0.00125001129840931	0.00236740273605756	0.00760322502490318	0.0156877081170168	0.0207096967113278	0.0246597928580262	0.0278815585466593
0.000020913672637073	0.0000921069735930886	0.000374644140212069	0.00155680628740693	0.0015618261590607	0.00232242059337922	0.00334293743232316];
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
n = 5000;
[t X] = ode45(@(t, X) odes_boot(t,X,k0), t_boot, X0);
% recalculating k_j from bootstrapped residuals n times
k = k0;
x_12P_exp_boot = [];
x_13P_exp_boot = [];
x_EG_exp_boot = [];
x_S_exp_boot = [];
x_I_exp_boot = [];
x_G_exp_boot = [];
x_M_exp_boot = [];
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
##  scatter(t_boot, X_boot(:,1))
  x_12P_exp_boot = [x_12P_exp_boot, X_boot(:,1)];
  x_13P_exp_boot = [x_13P_exp_boot, X_boot(:,2)];
  x_EG_exp_boot = [x_EG_exp_boot, X_boot(:,3)];
  x_S_exp_boot = [x_S_exp_boot, X_boot(:,4)];
  x_I_exp_boot = [x_I_exp_boot, X_boot(:,5)];
  x_G_exp_boot = [x_G_exp_boot, X_boot(:,6)];
  x_M_exp_boot = [x_M_exp_boot, X_boot(:,7)];
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

e_x_12P = [];
e_x_13P = [];
e_x_EG = [];
e_x_S = [];
e_x_I = [];
e_x_G = [];
e_x_M = [];
for i = 1:1:length(x_12P_exp_boot(:,1))
  e_x_12P(end+1) = std(x_12P_exp_boot(i,:));
  e_x_13P(end+1) = std(x_13P_exp_boot(i,:));
  e_x_EG(end+1) = std(x_EG_exp_boot(i,:));
  e_x_S(end+1) = std(x_S_exp_boot(i,:));
  e_x_I(end+1) = std(x_I_exp_boot(i,:));
  e_x_G(end+1) = std(x_G_exp_boot(i,:));
  e_x_M(end+1) = std(x_M_exp_boot(i,:));
endfor
hist(e_x_12P);

errors = [e_x_12P; e_x_13P; e_x_EG; e_x_S; e_x_I; e_x_G; e_x_M];

csvwrite("errors.csv", errors);
